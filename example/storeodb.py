############### Important ################
#Add the odbhelper module to Abaqus python packages 
#(Abaqus --> version --> tools --> SMApy --> Lib --> site-packages)
############### Important ################

import odbhelper

# ------------------ GLOBAL VARIABLES ----------------------------
# GENERAL VARIABLE NAMES:
# -----------------------
file_name = 'odbexample.odb'
instance_name = 'GRANULATION-1'
step_name = 'high_load1' # steps can be defined with name or number
frame_num = 'all' # -1 stands for "last frame"


odbh = odbhelper.OdbHelper(file_name,step_name,frame_num,instance_name)

var = 'NT11'
all_nodes_tuple = odb.set_node_list_instance[0]  # all nodes
odb.csvVariableWriter(variable = var, regiontuple = all_nodes_tuple)

var = 'SDV1'
all_elements_tuple = odb.set_element_list_instance[0]  # all elements
odb.csvVariableWriter(variable = var, regiontuple = all_elements_tuple)

var = 'SDV1'
odbh.export2vtk(var)
odbh.export2vtk(var, instances = 'all')

step_name = 0 # steps can be defined with name or number
odb.defineStep(step_name, frame_num)
var = 'S'
odb.export2vtk(var)