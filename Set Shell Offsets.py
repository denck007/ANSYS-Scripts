'''
Set the shell offsets for all of the bodies in a model.

Looks through the names of all bodies in the tree and if they end in 
  "BOTTOM" or "MID" set them to the corresponding offset. 
    Otherwise set them  to TOP offset
'''

def set_offset(obj):
    '''
    Set the offset of the part
    '''
    if not obj.GetType().Equals(Ansys.ACT.Automation.Mechanical.Body):
        return
    if obj.Name[-6:] == "BOTTOM":
        obj.InternalObject.ShellOffsetType =  2# 2=Bottom
    elif obj.Name[-3:] == "MID":
        obj.InternalObject.ShellOffsetType = 1 # 1=Midsurface
    else:
        obj.InternalObject.ShellOffsetType = 0 # 0=Top 

def set_offset_in_children(obj):
    '''
    Set the offset of all of the children of obj in the tree
    '''
    for c in obj.Children:
        #print("Name: {} Type: {}".format(c.Name,c.InternalObject.GeometryType))
        if len(c.Children)>0:
            set_offset_in_children(c)
        elif  c.InternalObject.GeometryType == 1: # type 0 is solid, type 1 is surface
            set_offset(c)
set_offset_in_children(Model.Geometry)
