"""
Author: Neil Dencklau
Title: Highlight High Stress Elements
Version: v1.0
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

"""


# Set the stress limit, All elements with nodal stresses over this value will be added to the named selection
stressLimit = Quantity("4000. [psi]") # can use units of MPa, Pa, psi

# Get some of the basic objects needed to perform the work.
sm = ExtAPI.SelectionManager
analysis = DataModel.AnalysisList[0] # This is the first analysis, change the index to run on others
meshObj = analysis.MeshData

# First find the first vonMises stress in the analysis
# If it does not exist throw an error
seqv = None
for child in analysis.Solution.Children:
    if child.GetType().Equals(Ansys.ACT.Automation.Mechanical.Results.StressResults.EquivalentStress):
        seqv = child
if seqv is None:
    msg = Ansys.Mechanical.Application.Message("Did not find a Equivalent Stress result in for analysis {}".format(analysis.Name),MessageSeverityType.Error)
    ExtAPI.Application.Messages.Add(msg)
    # For some reason raise does not work with buttons, but does in console.
    # But it still fails, just in a very unclear way.....
    raise Exception("Did not find a Equivalent Stress result in for analysis {}".format(analysis.Name))

# Next get a listing of the node ids and their corresponding stress value
# Note that the PlotData property does not exist before ANSYS 2020R1
nodeIdsInResult = seqv.PlotData.Values[1]
stressAtNodes = seqv.PlotData.Values[2]
stressUnit = seqv.PlotData.Dependents.Values[0].Unit

# Go over every node and see if it exceeds the stressLimit, and add them to a list
nodesHighStress = []
for node,stress in zip(nodeIdsInResult,stressAtNodes):
    if Quantity("{} [{}]".format(stress,stressUnit)) > stressLimit:
        nodesHighStress.append(node)

# Convert the node ids from above to their corresponding element ids
elementsHighStress = meshObj.ElementIdsFromNodeIds(nodesHighStress)

# The SelectionManager allows you to do the equivalent of manually selecting
#   items in the graphics. It can be scoped to any geometry feature type, or mesh feature type
#   First create a SelectionInfo object with all the data that you want to select, then actually
#   select the data. This allows for doing many modification to the selection without needing
#   to update the graphics which will slow things down
sm.ClearSelection()
selectionInfo = sm.CreateSelectionInfo(SelectionTypeEnum.MeshElements)
selectionInfo.Ids = elementsHighStress
sm.NewSelection(selectionInfo)

# Create the named selection. Just like when using Mechanical normally, anything that is selected
#   when a named selection is created is added to the named selection
ns = DataModel.Project.Model.AddNamedSelection()
ns.Name = "Elements with nodes over {}".format(stressLimit.ToString())

# This section using Breadth First Search to find connected components
# The underlying theory is beyond the scope of this example.
from collections import deque

def getNeighborElements(element):
    ''''
    Return the element ids for elements that are adjacent to the element
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
for element in elementsHighStress:
    elementTracker[element] = -1

groups = [] # The element Ids 
# Now use breadth first search to find groups of high stress elements
for element in elementsHighStress:
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

# Creating named selections behaves the same way as in the UI. Whatever is selected when the NS is created is in the NS
#    Of course the values can be edited later if needed
for idx, group in enumerate(groups):
    sm.ClearSelection()
    selectionInfo = sm.CreateSelectionInfo(SelectionTypeEnum.MeshElements)
    selectionInfo.Ids = group
    sm.NewSelection(selectionInfo)
    ns = DataModel.Project.Model.AddNamedSelection()
    ns.Name = "Elements with nodes over {} in group {}".format(stressLimit.ToString(),idx)