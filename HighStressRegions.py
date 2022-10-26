"""
Author: Neil Dencklau
Title: Highlight High Stress Elements
Version: v1.1
Created: 2022-10-26
Revised: 2022-10-26
Ansys version: 2022R2

Description:
Create a Named Selection for each distinct region the model that is over a given threshold.
Method:
1. Get the nodal stress values using the Equivalent Stress result object in the tree
2. Filter by nodal stress over the given stress limit
3. Create a named selection containing all the elements over the given limit.
4. Identify unique groups of elements over stress limit 
5. Create named selection for each unique group 

Change Log:
v1.1:
  * Refactored into functions to make code easier to read. No functional or algorithm changes.

"""
from collections import deque

# Set the stress limit, All elements with nodal stresses over this value will be added to the named selection
stressLimit = Quantity("4000. [psi]") # can use units of MPa, Pa, psi
analysis = DataModel.AnalysisList[0] # This is the first analysis, change the index to run on others


def get_stress_result(analysis):
    """
    Return the first Equivalent stress result in the tree for the analysis

    """
    # First find the first vonMises stress in the analysis
    # If it does not exist throw an error
    stress_result = None
    for child in analysis.Solution.Children:
        if child.GetType().Equals(Ansys.ACT.Automation.Mechanical.Results.StressResults.EquivalentStress):
            stress_result = child
    if stress_result is None:
        msg = Ansys.Mechanical.Application.Message("Did not find a Equivalent Stress result in for analysis {}".format(analysis.Name), MessageSeverityType.Error)
        ExtAPI.Application.Messages.Add(msg)
        # For some reason raise does not work with buttons, but does in console.
        # But it still fails, just in a very unclear way.....
        raise Exception("Did not find a Equivalent Stress result in for analysis {}".format(analysis.Name))
    return stress_result


def filter_stress_result(stress_result, stressLimit):
    """
    Given a stress result, return the listing of nodes over the given stressLimit
    
    stress_result: Ansys.ACT.Automation.Mechanical.Results.StressResults.EquivalentStress (a Equivalent Stress result object)
    stressLimit: Quantity (the stress value and unit to filter out values below)

    returns: list[int] (list of node numbers)
    """
    # Next get a listing of the node ids and their corresponding stress value
    # Note that the PlotData property does not exist before ANSYS 2020R1
    nodeIdsInResult = stress_result.PlotData.Values[1]
    stressAtNodes = stress_result.PlotData.Values[2]
    stressUnit = stress_result.PlotData.Dependents.Values[0].Unit

    # Go over every node and see if it exceeds the stressLimit, and add them to a list
    # If there are performance issues, this is probably a good place to start looking
    nodesHighStress = []
    for node,stress in zip(nodeIdsInResult,stressAtNodes):
        if Quantity("{} [{}]".format(stress,stressUnit)) > stressLimit:
            nodesHighStress.append(node)

    return nodesHighStress


def create_named_selection(ids, selection_type, name):
    """
    Given a list of ids (ex element ids), create a named selection with name
    
    ids: list[int] (ex: [1,2,3])
    selection_type: SelectionTypeEnum (ex:SelectionTypeEnum.MeshElements)
    name: str (ex: "my_named_selection")
    """
    sm = ExtAPI.SelectionManager
    # The SelectionManager allows you to do the equivalent of manually selecting
    #   items in the graphics. It can be scoped to any geometry feature type, or mesh feature type
    #   First create a SelectionInfo object with all the data that you want to select, then actually
    #   select the data. This allows for doing many modification to the selection without needing
    #   to update the graphics which will slow things down
    sm.ClearSelection()
    selectionInfo = sm.CreateSelectionInfo(selection_type)
    selectionInfo.Ids = ids
    sm.NewSelection(selectionInfo)

    # Create the named selection. Just like when using Mechanical normally, anything that is selected
    #   when a named selection is created is added to the named selection
    ns = DataModel.Project.Model.AddNamedSelection()
    ns.Name = name


def create_element_groups(elements_high_stress):
    """
    Given a list of elements, create groups where elements that share at least
        one node are in the same group

    elements_high_stress: list[int] (ex: [1, 2, 3, 4])

    returns: list[list[int]] (list containing lists of adjacent elements)
    
    This section using Breadth First Search to find connected components
    The underlying theory is beyond the scope of this example.
    """

    def getNeighborElements(element):
        ''''
        Return the element ids for elements that are adjacent to the element
        
        Note: functions can be defined inside of functions, 
            but they are only available inside of their parent function
        '''
        nodes = meshObj.NodeIdsFromElementIds([element])
        elements = meshObj.ElementIdsFromNodeIds(nodes)
        # list(set()) gets a list of unique values from a list
        elements = list(set(elements))
        return elements

    # Dictionary to allow fast lookup of elements that have high stress
    # This also allows mapping of element id to grouping id
    # Initialize the tracker with -1 for each element, meaning the element has not been explored
    #    this is an implementation detail and there are many different ways to approach this
    elementTracker = {}
    for element in elements_high_stress:
        elementTracker[element] = -1

    groups = [] # The element Ids 
    # Now use breadth first search to find groups of high stress elements
    for element in elements_high_stress:
        if (element in elementTracker) and (elementTracker[element] == -1):
            # Assign the element to the latest group.
            # The new group is not added to the groups list till all items in the group
            #	are found, so the current group is the length of groups, with starting index == 0
            elementTracker[element] = len(groups)
            currentGroup = [element]
    
            q = deque()
            q.append(element)   	 
            while len(q) > 0:
                elem = q.pop()
                # Add all neighboring elements that have not been explored
                for e in getNeighborElements(elem):
                    if (e in elementTracker) and (elementTracker[e] == -1):
                        q.append(e)
                        elementTracker[e] = len(groups)
                        currentGroup.append(e) 
            # all the elements in the group have been added to currentGroup, so add
            #   the current group to the list of groups
            groups.append(currentGroup)
    return groups


# Get the stress values and filter them
stress_result = get_stress_result(analysis)
nodes_high_stress = filter_stress_result(stress_result, stressLimit)

# Convert from nodes to elements, then group adjacent elements
meshObj = analysis.MeshData
elements_high_stress = meshObj.ElementIdsFromNodeIds(nodes_high_stress)
groups = create_element_groups(elements_high_stress)

# Finally create the Named Selections
create_named_selection(elements_high_stress, SelectionTypeEnum.MeshElements, "Elements with nodes over {}".format(stressLimit.ToString()))
for idx, group in enumerate(groups):
    create_named_selection(group, SelectionTypeEnum.MeshElements, "Elements with nodes over {} in group {}".format(stressLimit.ToString(),idx))