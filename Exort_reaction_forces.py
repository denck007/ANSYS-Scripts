'''
Exort all of the reaction forces in the model in the order specified by the 'to_export' variable
Will export to forces.csv in the user_files folder
'''import os

project = DataModel.Project
working_dir = ExtAPI.DataModel.AnalysisList[0].WorkingDir
user_files = os.path.join(working_dir[:working_dir.rfind("_files")+6],"user_files")
fname = os.path.join(user_files,"forces.csv")

to_export = ["front truck","rear truck","front left","front right","rear left","rear right",
            "front buff","front draft","rear buff","rear draft"]

with open(fname,"w") as f:
    f.write("Model Name\tTime step\t")
    for force in to_export:
        for direction in ["x","y","z"]:
            f.write("{} {}\t".format(force,direction))
            pass
    f.write("\n")
for model in project.Model.Analyses:
    time_steps = model.AnalysisSettings.InternalObject.NumberOfSteps
    sol = model.Solution
    for time_step in range(1,time_steps+1):
        forces = [""]*(len(to_export)*3)
        for rst in sol.Children:
            for column_group,reaction_to_find in enumerate(to_export):
                if rst.GetType().ToString().Equals("Ansys.ACT.Automation.Mechanical.Results.ProbeResults.ForceReaction"):
                    if int(rst.InternalObject.Time) == time_step:
                        if reaction_to_find in rst.Name.ToLower():
                            #name = rst.Name
                            #t = rst.Time.ToString()
                            #t = t[:t.find(" ")]
                            x = rst.XAxis.ToString()
                            x= float(x[:x.find(" ")])
                            forces[column_group*3+0] = x
                            y = rst.YAxis.ToString()
                            y = float(y[:y.find(" ")])
                            forces[column_group*3+1] = y
                            z = rst.ZAxis.ToString()
                            z = float(z[:z.find(" ")])
                            forces[column_group*3+2] = z
                            #print("{},{},{},{},{}".format(name,t,x,y,z))
        with open(fname,'a') as f:
            f.write("{}\t{}\t".format(model.Name,time_step))
            for force in forces:
                f.write("{}\t".format(force))
            f.write("\n")
