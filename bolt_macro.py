"""Bolt Load Calculation Macro
This code extracts axial and shear forces on all beams in the model, performs
bolted joint calculations on them, and exports the results to the project
user_files directory. Each run of script creates a new file with a timestamp.

It requires:
    * ANSYS 23.2 or newer
    * Nodal Forces are exported stored in results
    * Any solver units or user units are allowed. Everything is converted to N mm
    * Do not need to add any command snippets anywhere
    * Beam naming is not required but be set to customize results

Optional Inputs:
    * Beams can be named 'M_<grade>_<length in mm>_<any name>'
        * If name does not follow this convention will assume 8.8 and grip length
        being the length of the beam at solve time

Assumptions:
    * Steel to Steel connections with steel bolts
        * Can change stiffness of bolt and clamp material in `BoltResult.__init__()`
        * Can change joint coefficient of friction in `BoltResult.__init__()`
    * Torsion and moment loads on bolts are insignificant
        * No torsion or moment data is used anywhere
        * Pay extra attention to this if the joint slips!


Release History:
v1: 2023-09-01 Neil Dencklau
    * Initial release
    * Is direct migragation of 'Bolt_Calc_V15.txt' 
"""

import clr
import os
import Ansys.Utilities

clr.AddReferenceToFileAndPath(
    os.path.join(
        Ansys.Utilities.ApplicationConfiguration.DefaultConfiguration.WorkbenchInstallRootDirectoryPath,
        "Addins",
        "ACT",
        "bin",
        "Win64",
        "MathNet.Numerics.dll",
    )
)

import datetime

from System import Array as sys_array
import MathNet
import MathNet.Numerics.LinearAlgebra as la
import mech_dpf
import Ans.DataProcessing as dpf
import units
import wbjn

mech_dpf.setExtAPI(ExtAPI)


def array(*x):
    """Convert arbitrary collection of values to sys_array of Double"""
    return sys_array[float](x)  # float is equivalent to .Net double


def cross_product(a, b):
    """Return Vector of a cross b
    MathNet doesn't have a cross product, so implement one

    Args:
        a (list[float) OR MathNet.Numerics.LinearAlgebra.Double.Vector)
        b (list[float) OR MathNet.Numerics.LinearAlgebra.Double.Vector)
    Returns:
        Cross product of a cross b as a MathNet.Numerics.LinearAlgebra.Double.Vector
    """
    x = a[1] * b[2] - a[2] * b[1]
    y = a[2] * b[0] - a[0] * b[2]
    z = a[0] * b[1] - a[1] * b[0]
    return la.Double.Vector.Build.Dense(array(x, y, z))


def align_force_to_beam(node1, node2, force_global):
    """Converts beam result forces to beam coordinates
    Beam nodes are aligned to global coordinates, but we care
    about the axial and shear forces

    This takes the end nodes of the beam and the force in global CS
    and converts the force to axial tension == +Z
    Note: The X and Y directions alignment is not controlled

    Uses Rodrigues' rotation formula:
        https://en.m.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        Rough idea is rotate force vector about axis normal to the
        global Z vector and beam axis vector.

    Args:
        node1 (list[float]): List of 3 floating point numbers describing the
            first (I) node of the beam in global coordinates
        node2 (list[float]): List of 3 floating point numbers describing the
            second (J) node of the beam in global coordinates
        force_global (list[float]): List of 3 floating point numbers of the beam
            nodal force in global coordinates

    Returns:
        list[float]: 3 Floats indicating force aligned to beam axis
            is +Z is axial tension, X and Y are not aligned to anything in
            particular, just normal to Z
    """
    # Get inputs as vectors
    node1_v = la.Double.Vector.Build.Dense(array(*node1))
    node2_v = la.Double.Vector.Build.Dense(array(*node2))
    force_global_v = la.Double.Vector.Build.Dense(array(*force_global))

    v_source = (node1_v - node2_v).Normalize(2)
    v_target = la.Double.Vector.Build.Dense(array(0, 0, 1))

    # Cosine of angle between source and target
    cos_theta = v_target.DotProduct(v_source)
    theta = MathNet.Numerics.Trig.Acos(cos_theta)
    sin_theta = MathNet.Numerics.Trig.Sin(theta)

    if MathNet.Numerics.Precision.AlmostEqual(cos_theta, 1, 1e-8):
        # Beam axis points +Z, nothing to do
        force_local = force_global_v
    elif MathNet.Numerics.Precision.AlmostEqual(cos_theta, -1, 1e-8):
        # Beam axis points -Z, just flip the vector
        force_local = force_global_v.Negate()
    else:
        # Need to rotate the force
        # Get vector normal to global Z and target Z
        K = cross_product(v_source, v_target).Normalize(2)
        # Apply Rodrigue's
        force_local = (
            force_global_v * cos_theta
            + cross_product(K, force_global_v) * sin_theta
            + K * K.DotProduct(force_global_v) * (1 - cos_theta)
        )

    return [force_local[0], force_local[1], force_local[2]]


class Bolt:
    def __init__(self, ui_beam):
        self.ui_beam = ui_beam
        self.diameter = self._get_bolt_size()
        self._parse_name_grade_length()

    def _get_bolt_size(self):
        """Get the nominal bolt diameter [mm] from UI
        Is the UI value converted to mm rounded to nearest full mm
        """
        d = units.ConvertUnit(
            self.ui_beam.Radius.Value, 
            self.ui_beam.Radius.Unit, 
            "mm"
        )
        d = round(2*d)
        allowed_sizes = [4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 24, 30, 36,]
        if d not in allowed_sizes:
            raise ValueError(
                "Unexpected bolt size {}mm. "
                "Allowed sizes are {} [mm]".format(d, allowed_sizes)
            )
        return d

    def _parse_name_grade_length(self):
        """Parse the UI name of the beam
        Sets values:
            grade (float): If not specified sets to 8.8
            grip_length (float or None): If not specified sets to None
            name (str): If name starts with 'M_' uses name as:
                'M_<grade>_<length>_<name>'. Otherwise uses full name from UI
        """
        # M_088_037_NAME_4
        name_full = self.ui_beam.Name
        if name_full.startswith("M_"):
            name_split = name_full.split("_")
            if len(name_split) <= 3:
                raise ValueError(
                    "Beam named incorrectly. "
                    "Expected format 'M_<grade>_<length>_<name>'"
                )
            
            grade_str = name_split[1]
            try:
                if "." in grade_str: # Entered as '_8.8_', '_08.8_'
                    grade = float(grade_str)
                else: # allother cases assume needs to be divided by 10
                    grade = float(grade_str)/10.0
            except Exception as e:
                raise ValueError(
                    "Could not parse grade, got '{}' "
                    "which could not be parsed into a grade.".format(grade_str)
                )
            allowed_grades = [8.8, 10.9]
            if grade not in allowed_grades:
                raise ValueError(
                    "Bad grade specified, got {} "
                    "but only {} are allowed.".format(grade, allowed_grades)
                )
            self.grade = grade

            try:
                self.grip_length = float(name_split[2])
            except Exception as e:
                raise ValueError(
                    "Could not parse length, got '{}' "
                    "which could not be parsed into a length.".format(name_split[2])
                )

            # Build the beam name back up, rejoining any '_'
            # and removing any commas that will mess with csv output
            name = "_".join(name_split[3:])  # Everything after 3'rd '_' is the name
            name = name.replace(",", "")  
            self.name = name
        else:
            self.grade = 8.8
            self.grip_length = None
            self.name = name_full

    @property
    def area_stress(self):
        """Stress Area in [mm^2]"""
        area_by_diameter = {
            4: 8.78,
            5: 14.2,
            6: 20.1,
            7: 28.8,
            8: 36.6,
            10: 58.0,
            12: 84.0,
            14: 115.0,
            16: 157.0,
            18: 192.0,
            20: 245.0,
            24: 353.0,
            30: 561.0,
            36: 817.0,
        }
        area_stress = area_by_diameter[self.diameter]
        return area_stress

    @property
    def yield_stress(self):
        """Yield Stress in MPa based on grade"""
        if self.grade == 8.8:
            yield_stress = 640  # MPa
        elif self.grade == 10.9:
            yield_stress = 940  # MPa
        else:
            raise ValueError("Invalid grade {}".format(self.grade))
        return yield_stress

    @property
    def proof_load(self):
        """Proof Load in [N]
        Based on 90% of yield stress
        """
        proof_load = 0.9 * self.yield_stress * self.area_stress
        return proof_load

    @property
    def clamp_load(self):
        """Clamp Load in [N]
        Based on 75% of proof load
        """
        clamp_load = 0.75 * self.proof_load
        return clamp_load


class BoltResult:
    def __init__(
        self, bolt, analysis, mesh, elemental_nodal_forces, analysis_time_step
    ):
        # These values are material dependant!
        self.joint_friction_coeff = 0.1
        self.stiffness_bolt = 2.11e5  # Eb, MPa
        self.stiffness_clamp = 2.05e5  # Ec, MPa

        self.bolt = bolt
        self.analysis = analysis
        self.solver_data = analysis.Solution.SolverData
        self.mesh = mesh
        self.elemental_nodal_forces = elemental_nodal_forces
        self.analysis_time_step = analysis_time_step

        self._load_mesh_data()
        self._set_grip_lengths()
        self._load_and_compute_forces()
        self._run_bolt_calc()

    def _load_and_compute_forces(self):
        """Load beam forces, convert to local coordinates, axial and shear

        Sets:
            forces_global: [xi, yi, zi, xj, yj, zj] in solver units global coordinates
            forces: [xi, yi, zi, xj, yj, zj] in solver units with Z axial
            axial_force: Force in N along beam axis, + is tension
            shear_force: Total shear force in N
        """
        self.forces_global = list(
            elemental_nodal_forces.GetEntityDataById(self.element.Id)
        )
        node1_pos = [
            self.nodes[0].X,
            self.nodes[0].Y,
            self.nodes[0].Z,
        ]
        node2_pos = [
            self.nodes[1].X,
            self.nodes[1].Y,
            self.nodes[1].Z,
        ]
        forces_i = align_force_to_beam(node1_pos, node2_pos, self.forces_global[:3])
        forces_j = align_force_to_beam(node1_pos, node2_pos, self.forces_global[3:])

        self.forces = [
            forces_i[0],
            forces_i[1],
            forces_i[2],
            forces_j[0],
            forces_j[1],
            forces_j[2],
        ]
        self.axial_force = self._compute_axial_force()
        self.shear_force = self._compute_shear_force()

    def _load_mesh_data(self):
        element_id = self.solver_data.GetObjectData(self.bolt.ui_beam).ElementId
        self.element = self.mesh.ElementById(element_id)
        self.nodes = self.element.Nodes

    def _set_grip_lengths(self):
        """Set the grip_length and grip_length_nodal in [mm]
        Computes the `grip_length_nodal` from the beam node positions
        Sets `grip_length` to self.bolt.grip_length (if specified)
            otherwise sets it to grip_length_nodal

        """
        grip_solver_units = (
            (self.nodes[0].X - self.nodes[1].X) ** 2
            + (self.nodes[0].Y - self.nodes[1].Y) ** 2
            + (self.nodes[0].Z - self.nodes[1].Z) ** 2
        ) ** 0.5
        self.grip_length_nodal = units.ConvertUnit(
            grip_solver_units, self.mesh.Unit, "mm"
        )  # L, grip length in mm

        if self.bolt.grip_length is not None:
            self.grip_length = self.bolt.grip_length
        else:
            self.grip_length = self.grip_length_nodal

    def _compute_axial_force(self):
        """Compute the axial load on the beam in [N]
        Coordinate system for beams is Z axial
        When looking at the reference end (node I)
            +Z is 'away' from mobile end (node J)
            This means a +Z on I is away from J
            but +Z on J is towards I
            Invert sign on J to make it 'Positive is pulling apart'
        However also want 'load on the bolt by the plate'
            but the reaction force on the node is 'load on the plate by the bolt'
            so invert the whole thing
        This is '-max(zi,-zj)' which is just max(-zi, zj)
        """
        i = self.forces[2]
        j = self.forces[5]
        axial = max(-i, j)
        axial_N = units.ConvertUnit(axial, self.elemental_nodal_forces.Unit, "N")
        return axial_N

    def _compute_shear_force(self):
        """Compute the shear load on the beam in [N]
        Computes the vector sum of shear forces at each end
        and returns the larger of the 2
        """
        xi = self.forces[0]
        yi = self.forces[1]
        xj = self.forces[3]
        yj = self.forces[4]
        shear_i = (xi**2 + yi**2) ** 0.5
        shear_j = (xj**2 + yj**2) ** 0.5
        shear = max(shear_i, shear_j)
        shear_N = units.ConvertUnit(shear, self.elemental_nodal_forces.Unit, "N")
        return shear_N

    def _run_bolt_calc(self):
        d_hole = 1.5 * self.bolt.diameter  # Dh
        d_washer = 2.0 * self.bolt.diameter  # Dw

        d3 = d_washer + 0.5774 * self.grip_length  # d3
        d_clamp = (d3 + d_washer) / 2.0  # Dc
        area_clamp = 0.785 * (d_clamp**2 - d_hole**2)  # Ac
        k_bolt = (self.bolt.area_stress * self.stiffness_bolt) / self.grip_length  # kb
        k_clamp = (area_clamp * self.stiffness_clamp) / self.grip_length  # kc

        k_ratio = k_clamp / k_bolt  # Ratio

        f_clamp = self.bolt.clamp_load - (k_clamp / (k_clamp + k_bolt)) * max(
            self.axial_force, 0
        )  # Fc
        f_bolt = self.bolt.clamp_load + (k_bolt / (k_clamp + k_bolt)) * max(
            self.axial_force, 0
        )  # Fb

        shear_capacity_joint = self.joint_friction_coeff * f_clamp  # Csl, N
        shear_capacity_bolt = (
            self.bolt.yield_stress / 1.732
        ) * self.bolt.area_stress  # Csh, N
        axial_capacity_bolt = self.bolt.yield_stress * self.bolt.area_stress  # Cax, N

        # Compute the safety factors on the joint and bolt
        sf_shear_joint = shear_capacity_joint / self.shear_force
        sf_shear_bolt = shear_capacity_bolt / self.shear_force
        sf_axial_bolt = axial_capacity_bolt / f_bolt

        # If the bolt slips then find the bolt shear stress,
        #   else the bolt shear stress is 0
        if sf_shear_joint < 1.0:
            shear_stress_bolt = self.shear_force / self.bolt.area_stress
        else:
            shear_stress_bolt = 0.0

        axial_stress_bolt = f_bolt / self.bolt.area_stress
        total_stress_bolt = (
            axial_stress_bolt**2 + 3.0 * shear_stress_bolt**2
        ) ** 0.5

        self.f_clamp = f_clamp
        self.f_bolt = f_bolt
        self.k_ratio = k_ratio
        self.shear_capacity_joint = shear_capacity_joint
        self.shear_capacity_bolt = shear_capacity_bolt
        self.axial_capacity_bolt = axial_capacity_bolt
        self.sf_shear_joint = sf_shear_joint
        self.sf_shear_bolt = sf_shear_bolt
        self.sf_axial_bolt = sf_axial_bolt
        self.shear_stress_bolt = shear_stress_bolt
        self.axial_stress_bolt = axial_stress_bolt
        self.total_stress_bolt = total_stress_bolt

    def result_string(self):
        """Return csv formatted string for this result"""
        values = [
            self.analysis.Name.replace(
                ",", " "
            ),  # Analysis Name, clear out any ',' that will mess up csv file
            "{}".format(self.analysis_time_step),  # Time [s]
            self.bolt.name,  # Name
            "{:.1f}".format(self.bolt.diameter),  # Bolt Size [mm]
            "{:.1f}".format(self.bolt.grade),  # Grade
            "{:.2f}".format(self.grip_length),  # Grip Length [mm]
            "{:.2f}".format(self.grip_length_nodal),  # Grip Length (nodal) [mm]
            "{:.1f}".format(self.shear_force),  # Shear Force [N]
            "{:.1f}".format(self.axial_force),  # Axial Force [N]
            "{:.1f}".format(self.shear_capacity_joint),  # Shear Capacity Slip [N]
            "{:.1f}".format(self.shear_capacity_bolt),  # Shear Capacity Bolt [N]
            "{:.1f}".format(self.axial_capacity_bolt),  # Axial Capacity Bolt [N]
            "{:.3f}".format(self.sf_shear_joint),  # SF Slip"
            "{:.3f}".format(self.sf_shear_bolt),  # SF Shear Bolt
            "{:.1f}".format(self.sf_axial_bolt),  # SF Axial Bolt
            "{:.1f}".format(self.shear_stress_bolt),  # Shear Stress Bolt [MPa]
            "{:.1f}".format(self.axial_stress_bolt),  # Axial Stress Bolt [MPa]
            "{:.1f}".format(self.total_stress_bolt),  # Equivalent Stress Bolt [MPa
        ]
        result_csv = ",".join(values)
        return result_csv

    def result_string_header(self):
        """Return csv formatted header for output file"""
        column_names = [
            "Analysis Name",
            "Time [s]",
            "Name",
            "Bolt Size [mm]",
            "Grade",
            "Grip Length [mm]",
            "Grip Length (nodal) [mm]",
            "Shear Force [N]",
            "Axial Force [N]",
            "Shear Capacity Slip [N]",
            "Shear Capacity Bolt [N]",
            "Axial Capacity Bolt [N]",
            "SF Slip",
            "SF Shear Bolt",
            "SF Axial Bolt",
            "Shear Stress Bolt [MPa]",
            "Axial Stress Bolt [MPa]",
            "Equivalent Stress Bolt [MPa]",
        ]
        header_csv = ",".join(column_names)
        return header_csv


bolts = []
ui_beams = ExtAPI.DataModel.Project.Model.Connections.GetChildren[
    Ansys.ACT.Automation.Mechanical.Connections.Beam
](True)
for ui_beam in ui_beams:
    # Dont try and use any suppressed beams
    if ui_beam.Suppressed:
        continue

    try:
        bolts.append(Bolt(ui_beam))
    except Exception as e:
        print("Failed to make a bolt for '{}'. Got error: {}".format(ui_beam.Name, e))
        msg = Ansys.Mechanical.Application.Message(
            "Failed to make a bolt for '{}'. Got error: {}".format(ui_beam.Name, e),
            MessageSeverityType.Error
        )
        ExtAPI.Application.Messages.Add(msg)
        continue


bolt_results = []
for analysis in ExtAPI.DataModel.Project.Model.Analyses:
    # fmt: off
    if analysis.Solution.Status!= Ansys.Mechanical.DataModel.Enums.SolutionStatusType.Done:
        msg = Ansys.Mechanical.Application.Message(
            "Analysis '{}' is not solved, cannot export bolt forces".format(analysis.Name),
            MessageSeverityType.Error
        )
        ExtAPI.Application.Messages.Add(msg)
        continue
    elif analysis.AnalysisSettings.NodalForces != Ansys.Mechanical.DataModel.Enums.OutputControlsNodalForcesType.Yes:
        msg = Ansys.Mechanical.Application.Message(
            "Analysis '{}' does not have 'Analysis Settings -> Output Controls -> Nodal Forces' "
            "turned on, cannot export bolt forces.".format(analysis.Name),
            MessageSeverityType.Error
        )
        ExtAPI.Application.Messages.Add(msg)
        continue
    
    try:
        model = dpf.Model(analysis.ResultFileName)
        mesh = model.Mesh

        if model.ResultInfo.AnalysisType != dpf.enums.AnalysisType.Static:
            msg = Ansys.Mechanical.Application.Message(
                "Analysis '{}' is not 'Static Sturctural, cannot export bolt forces".format(analysis.Name),
                MessageSeverityType.Error
            )
            ExtAPI.Application.Messages.Add(msg)
            continue
        # fmt: on

        beam_elements = dpf.operators.scoping.on_mesh_property()
        beam_elements.inputs.property_name.Connect("beam_elements")
        beam_elements.inputs.mesh.Connect(mesh)

        for step_end_time in analysis.StepsEndTime:
            # Get the data scoped to the expected time
            op_nodal = dpf.operators.result.element_nodal_forces(
                data_sources=model.DataSources,
                mesh_scoping=beam_elements.outputs.mesh_scoping,
                time_scoping=step_end_time,
            )
            # Data is returned as list (by time), want the 1st time
            #   which is the time we just scoped to
            elemental_nodal_forces = op_nodal.outputs.fields_container.GetData()[0]

            for bolt in bolts:
                result = BoltResult(
                    bolt=bolt,
                    analysis=analysis,
                    mesh=mesh,
                    elemental_nodal_forces=elemental_nodal_forces,
                    analysis_time_step=step_end_time,
                )
                bolt_results.append(result)
        
    except Exception as e:
        msg = Ansys.Mechanical.Application.Message(
            "Failed to calculate bolt results for analysis '{}': {}".format(analysis.Name, e),            
            MessageSeverityType.Info,
        )
        ExtAPI.Application.Messages.Add(msg)
    finally:
        # Done with this analysis so release all data streams (ie close the files)
        model.ReleaseStreams()

user_files = wbjn.ExecuteCommand(ExtAPI, "returnValue(GetUserFilesDirectory())")
output_fname = os.path.join(
    user_files,
    "bolt_loads_{}.csv".format(datetime.datetime.now().strftime("%Y%m%dT%H%M%S")),
)
with open(output_fname, "w") as fp:
    fp.write(bolt_results[0].result_string_header())
    fp.write("\n")
    for item in bolt_results:
        fp.write(item.result_string())
        fp.write("\n")

msg = Ansys.Mechanical.Application.Message(
    "Exported {} bolt results to {}".format(len(bolt_results), output_fname),
    MessageSeverityType.Info,
)
ExtAPI.Application.Messages.Add(msg)
