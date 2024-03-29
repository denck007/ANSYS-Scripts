"""
Author: Neil Dencklau
Title: Export Force Reactions
Version: v1.1 
Created: 2022-10-26
Revised: 2022-10-26
Ansys version: 2022R2

Description:
Export force reactions for each time step, from every analysis in the current project to a csv file.
Does this by reading Force Reaction probes in the analysis.

Change log:
v1.1: 
  * Outputs to the 'user_files' directory of the current project
  * Use Quantity's Value method instead of manually parsing
"""

import os
working_dir = ExtAPI.DataModel.AnalysisList[0].WorkingDir
user_files = os.path.join(working_dir[:working_dir.rfind("_files")+6],"user_files")
output_fname = os.path.join(user_files,"forces.csv")

# Create the header for the output file
output_data = "Analysis Name, Result Name, Result Time, Unit,X Axis, Y Axis,  Z Axis\n"

for analysis in DataModel.AnalysisList: # go over every analysis in the tree
    for step in analysis.StepsEndTime:# go over every time step in the analysis
        t = Quantity("{} [sec]".format(step))
        for result in analysis.Solution.Children: # Search the results for reactions matching the time step
            # Check to see if the result is a reaction force probe
            # The 'not' keyword is a boolean inversion
            # The 'continue' keyword skips any following code and goes to the next iteration of the loop
            # Done this way to avoid a bunch of nested if statements
            if not result.GetType().Equals(Ansys.ACT.Automation.Mechanical.Results.ProbeResults.ForceReaction):
                continue
    
            # Need to handle case where DisplayTime is last. To get the 'Last' a value of '0' is set as the
            # display time. Check for this and set the correct time.
            result_time = result.DisplayTime
            if result.DisplayTime.Equals(Quantity("0 [sec]")):
                # Get the last time and set the result_time
                last_time = analysis.StepsEndTime[analysis.StepsEndTime.Count-1]
                result_time = Quantity("{} [sec]".format(last_time))

            if result_time.Equals(t):
                msg = Ansys.Mechanical.Application.Message("Exporting reaction {} for time {} in analysis {}".format(result.Name,result_time,analysis.Name), MessageSeverityType.Info)
                ExtAPI.Application.Messages.Add(msg)
                #Results are type Quantity which will print out with the unit which we do not want
                x = result.XAxis.Value
                y = result.YAxis.Value
                z = result.ZAxis.Value
                unit = result.XAxis.Unit
                output_data += "{},{},{},{},{},{},{}\n".format(analysis.Name,result.Name,result_time,unit,x,y,z)

# Open the file in write mode ('w' is write, 'r' is read , 'a'is append)
# keywork with is a context manager. It will automatically close the file at the end of the indented block
#with open(output_fname,'w') as fp:
with open(output_fname,'w') as fp:
    fp.write(output_data)
