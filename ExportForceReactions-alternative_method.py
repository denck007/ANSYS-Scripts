"""
Author: Neil Dencklau
Title: Export Force Reactions - Alternative Method
Version: v1.0
Created: 2022-10-26
Revised: 2022-10-26
Ansys version: 2022R2

Description:
Export force reactions for each time step, from every analysis in the current project to a csv file.
Method:
1. Get all the supports in the analysis
2. Identify all the nodes attached to each support
3. Reads the result file, summing the forces on the nodes in each support
4. Formats the result for output to csv

Change log:

Known Issues:
* Only supports that directly attach to nodes works (ie displacement, fixed, nodal). 
    Any supports that attach via MPC (remote displacement) will not work with this code,
    but this code could be updated to output remote displacement as well.
"""


def get_nodes_in_supports(analysis, support_types):
    """
    Given an analysis object and the types of supports to export results for
    Create a dictionary mapping the support name to the node ids in the support
    """
    nodes = {}
    meshObj = analysis.MeshData
    for support_type in support_types:
        for child in analysis.GetChildren(support_type, True):
            nodes[child.Name] = []
            for location_id in child.Location.Ids:
                nodes[child.Name].extend(
                              meshObj.MeshRegionById(location_id).NodeIds
                          )
    return nodes

def read_results_for_analysis(analysis, nodes):
    """
    Given an analysis and a dictionary of {<support name>: <nodes in support>}
    Create a list of results for each support at each time step
    """
    results = []
    for support_name, node_ids in nodes.items():
        with analysis.GetResultsData() as reader:
            times = reader.ListTimeFreq # [1, 2, 3]
            for idx in range(reader.ResultSetCount):
                # ResultSet is indexed starting at 1, not 0
                reader.CurrentResultSet = idx+1 
                force = reader.GetResult("F")
                result = { 
                    "analysis_name": analysis.Name,
                    "support_name": support_name,
                    "time": times[idx],
                    "unit": {},
                    "quantity_type": {},
                    "data": {},
                }

                for component in force.Components:
                    comp_info = force.GetComponentInfo(component)
                    result["unit"][component] = comp_info.Unit
                    result["quantity_type"][component] = comp_info.QuantityName
                    force.SelectComponents([component])
                    result["data"][component] = sum(force.GetNodeValues(node_ids))
                results.append(result)
    return results


def get_unique_keys(data, key):
    """
    Given a list of dictionaries data with format:
        [
            {<key>: {"X": 0, "Y": 323, "Z": 43, "ASDF": -1}},
            {<key>: {"X": 0, "Y": 323, "Z": 43, "QWE}},
        ]
    return: list of unique keys in 'key': ["ASDF", "QWE", "X", "Y", "Z"]
    """
    keys = set()
    for row in data:
        keys.update(row[key].keys())
    return sorted(list(keys))


def get_output_string(results):
    """
    Convert the list of results to a csv style string
    with each component fully listed
    """
    components = get_unique_keys(results, "data")
    output = "Analysis Name,Support Name,Time," 
    output += ",".join(["Result {}".format(c) for c in components]) +  ","
    output += ",".join(["Unit {}".format(c) for c in components]) +  ","
    output += ",".join(["Quantity {}".format(c) for c in components]) +"\n"
    for result in results:
        line = [ 
            result["analysis_name"],
            result["support_name"],
            "{}".format(result["time"])
        ]
        line.extend(["{}".format(result["data"][c]) for c in components])
        line.extend(["{}".format(result["unit"][c]) for c in components])
        line.extend(["{}".format(result["quantity_type"][c]) for c in components])

        output += ",".join(line)
        output += "\n"
    return output


# Now run the functions
support_types = [
    DataModelObjectCategory.FixedSupport, 
    DataModelObjectCategory.Displacement,
    ]
results = []
for analysis in ExtAPI.DataModel.AnalysisList:
    nodes = get_nodes_in_supports(analysis, support_types)
    results.extend(read_results_for_analysis(analysis, nodes))
output = get_output_string(results)

with open("C:\\temp\\forces.csv",'w') as fp:
    fp.write(output)