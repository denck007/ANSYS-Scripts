'''
Given a selected object in the tree, add a Total Deformation and Eqv Stress for each time step
Also add a force reaction probe for each support in the tree at each time step
  * Will rename the probes based on the support name
  
 Will look for existing solution objects so it will not add duplicates
'''

fao = Tree.FirstActiveObject

if fao.GetType().ToString().Equals('Ansys.ACT.Automation.Mechanical.Analysis'):
    print("Analysis selected")
    analysis = fao
    solution = analysis.Solution
elif fao.GetType().ToString.Equals('Ansys.ACT.Automation.Mechanical.Solution'):
    print("Solution selected")
    analysis = fao.Parent
    solution = fao
# is result
elif fao.Parent.GetType().ToString().Equals('Ansys.ACT.Automation.Mechanical.Solution'):
    print("Child of solution selected")
    solution = fao.Parent
    analysis = solution.Parent
    
elif fao.Parent.GetType().ToString().Equals('Ansys.ACT.Automation.Mechanical.Analysis'):
    print("Child of analysis selected")
    analysis = fao.Parent
    solution = analysis.Solution
else:
    print("Did not find a match")
    analysis = None
    solution = None
boundary_conditions = {}
for child in analysis.Children:
    ansys_type = child.GetType().ToString()
    if 'Ansys.ACT.Automation.Mechanical.BoundaryConditions.RemoteDisplacement' in ansys_type:
        boundary_conditions[child.ObjectId] = child
        
existing_results = {}
# get all the existing results
for rst in solution.Children:
    ansys_type = rst.GetType().ToString()
    if ansys_type.Equals('Ansys.ACT.Automation.Mechanical.Results.StressResults.EquivalentStress') or ansys_type.Equals('Ansys.ACT.Automation.Mechanical.Results.DeformationResults.TotalDeformation'):
        t = rst.DisplayTime.ToString()
        t = t[:t.find(" [")]
        s = "{} t{}".format(rst.GetType().ToString(),t)
        existing_results[s] = 1
    if 'Ansys.ACT.Automation.Mechanical.Results.ProbeResults' in ansys_type:
        t = rst.DisplayTime.ToString()
        t = t[:t.find(" [")]
        s = "{} t{}".format(rst.BoundaryConditionSelection.ObjectId.ToString(),t)
        existing_results[s] = 1
        
# Add the results as needed
time_steps = analysis.AnalysisSettings.InternalObject.NumberOfSteps
for time_step in range(1,time_steps+1):
    if 'Ansys.ACT.Automation.Mechanical.Results.DeformationResults.TotalDeformation t{}'.format(time_step) not in existing_results:
        rst = solution.AddTotalDeformation()
        rst.DisplayTime = Quantity("{} [sec]".format(time_step))
        rst.Name = rst.Name + "t{}".format(time_step)
    if 'Ansys.ACT.Automation.Mechanical.Results.DeformationResults.TotalDeformation t{}'.format(time_step) not in existing_results:
        rst = solution.AddEquivalentStress()
        rst.DisplayTime = Quantity("{} [sec]".format(time_step))
        rst.Name = rst.Name + "t{}".format(time_step)
    
    for bc in boundary_conditions:
        #print("Adding force reaction {} t{}".format(boundary_conditions[bc].Name,t))
        if "{} t{}".format(bc,time_step) not in existing_results:
            print("Adding force reaction {} t{}".format(boundary_conditions[bc].Name,time_step))
            rst = solution.AddForceReaction()
            rst.Name = "{} t{}".format(boundary_conditions[bc].Name,time_step)
            rst.DisplayTime = Quantity("{} [sec]".format(time_step))
            rst.LocationMethod = LocationDefinitionMethod().BoundaryCondition
            rst.BoundaryConditionSelection = boundary_conditions[bc]
