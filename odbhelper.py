import numpy as np
import odbAccess
import csv
import os
from operator import add
from abaqusConstants import *
import warnings


class OdbHelper:
    """ Class containing methods to help reading .odb files and extracting the desired data.
    After calling the constructor, user should call defineStepAndFrame and defineInstance methods
    before using any other methods.

    IMPORTANT: for a complete use of this module, the first elementSet/nodeSet of all the instances
    must include all elements/nodes (sets are listed in alphabetical order by Abaqus).
    Creating "ALL_ELEMENTS" and "ALL_NODES" sets in all instances is usually enough.
    """
    # ------------------ GLOBAL VARIABLES ----------------------------
    # OUTPUT VARIABLES AT NODES:
    # -----------------------
    ov_nodes = {
        "CF": "Point loads at nodes",
        "CM3": "Point moments at nodes",
        "RF": "Reaction force at nodes",
        "RM3": "Reaction moment at nodes",
        "U": "Spatial displacement at nodes",
        "UR": "Rotational displacement at nodes"
    }
    # OUTPUT VARIABLES AT INTEGRATION POINTS:
    # -----------------------
    ov_int_points = {
        "CENER": "Creep dissipation energy density at integration points",
        "DMENER": "Damage dissipation energy density at integration points",
        "E": "Strain components force at integration points",
        "EENER": "Electrostatic energy density at integration points",
        "JENER": "Electrical current dissipation per unit volume at integration points",
        "PENER": "Plastic dissipation energy density at integration points",
        "S": "Stress components at integration points",
        "SENER": "Strain energy density at integration points",
        "VENER": "Viscous dissipation energy at integration points"
    }

    element_vtk_cell_type = {
        # trusss
        "T2D2": 3,              # vtk_line
        "T2D2H": 3,
        "T2D3": 4,              # vtk_poly_line
        "T2D3H": 4,
        "T3D2": 3,
        "T3D2H": 3,
        "T3D3": 4,
        "T3D3H": 4,
        # beams
        "B21": 3,
        "B21H": 3,
        "B22": 4,
        "B22H": 4,
        "B31": 3,
        "B31H": 3,
        "B32": 4,
        "B32H": 4,
        "B33": 4,
        "B33H": 4,
        # surfaces
        "S4": 9,                # vtk_quad
        "S4R": 9,
        "S4RS": 9,
        "S4RSW": 9,
        "S4R5": 9,
        "S8R": 23,              # vtk_quadratic_quad
        "S8R5": 23,
        "STRI3": 5,             # vtk_triangle
        "S3": 5,
        "S3R": 5,
        "S3RS": 5,
        "CPE3": 5,
        "CPE3T": 5,
        # volumes
        "C3D8": 12,             # vtk_hexahedron
        "C3D8H": 12,
        "C3D8I": 12,
        "C3D8IH": 12,
        "C3D8R": 12,
        "C3D8RH": 12,
        "C3D20": 25,
        "C3D20H": 25,           # vtk_quadratic_hexahedron
        "C3D20R": 25,
        "C3D20RH": 25,
        "C3D4": 10,             # vtk_tetrahedron
        "C3D4T": 10,
        "C3D4H": 10,
        "C3D10": 24,            # vtk_quadratic_tetrahedron
        "C3D10H": 24,
        "C3D10I": 24,
        "C3D10M": 24,
        "C3D10MH": 24,
        "C3D6": 13,             # vtk_wedge
    }

    def __init__(self, filename, step = [], frame_num = [], instance_name = []):
        self.file_name = filename
        self.odb = odbAccess.openOdb(filename)
        self.assembly = self.odb.rootAssembly
        self.step_name = []
        self.frame = []
        self.instance_name = []
        if step != []:
            self.defineStep(step, frame_num)
        if (step == [] and frame_num != []):
            warnings.warn('Frame cannot be defined if no step is selected')
        if instance_name != []:
            self.defineInstance(instance_name)

    def defineStep(self, step, frame_num = []): # by default get the last available frame
        self.frame = [] # delete frames of previous steps defined
        if isinstance(step, str):
            self.step_name = step
            self.step = self.odb.steps[step]
            print("Step name defined: " + step)
        elif isinstance(step, int):
            self.step_name = self.odb.steps.items()[step][0]
            self.step = self.odb.steps.items()[step][1]
            print("Step number defined: " + self.step_name)

        if frame_num != []:
            self.defineFrame(frame_num)
        else:
            warnings.warn('Step defined but no frame selected')

    def defineFrame(self, frame_num=-1):  # by default get the last available frame
        if self.step_name == []:
            raise AttributeError('No step defined before defining frame')

        if frame_num == 'all':
            self.frame = self.step.frames
        else:
            # to avoid errors if the user inputs a single frame (integer) or a list of frames
            if type(frame_num) is not list:
                frame_num = [frame_num]

            try:
                for i in frame_num:
                    self.frame.append(self.step.frames[i])
            except:
                raise ValueError('Wrong input frames')

        print("Frame defined: " + str(frame_num))

    def defineInstance(self, instance_name):
        try:
            self.instance_name = instance_name
            self.instance = self.assembly.instances[self.instance_name]
        except:
            self.instance_name = instance_name.upper()  # sometimes in the Abaqus odb file the instance name gets automatically uppercased
            self.instance = self.assembly.instances[self.instance_name]
            print("Instance name has been UPPERCASED")
        print("Instance name defined: " + self.instance_name)
        self.extractNodeSetNames()
        self.extractElementSetNames()

    def extractNodeSetNames(self):
        # nset = self.set_node_list_assembly[0]  # all elements
        # nset[0] set name, nset[1] actual set object

        self.set_node_list_instance = self.instance.nodeSets.items()
        self.set_node_names_instance = [nt[0] for nt in
                                        self.set_node_list_instance]  # get node set name from each tuple
        self.set_node_list_assembly = self.assembly.nodeSets.items()
        self.set_node_names_assembly = [nt[0] for nt in self.set_node_list_assembly]

    def extractElementSetNames(self):
        # elset = self.set_element_list_assembly[0]  # all elements
        # elset[0] set name, elset[1] actual set object

        self.set_element_list_instance = self.instance.elementSets.items()
        self.set_element_names_instance = [nt[0] for nt in
                                           self.set_element_list_instance]  # get node set name from each tuple
        self.set_element_list_assembly = self.assembly.elementSets.items()
        self.set_element_names_assembly = [nt[0] for nt in self.set_element_list_assembly]

    def csvVariableStore(self, fo, frame_info):
        # store output variable as .csv file
        # input: fielOutputs file (ex. fo = odbh.frame[0].fieldOutputs['E'])
        #        info of the current frame to store - list

        # store the name of the variables in fieldOutputs.values.data
        complabel = fo.componentLabels
        if complabel == ():
            complabel = [fo.name]
        else:
            complabel = [i.lower() for i in complabel]

        # store the name of the variables in fieldOutputs.values such as mises, tresca, etc...
        validinv = []
        if fo.validInvariants != ():
            switcher = {
                'min_inplane_principal': 'minInPlanePrincipal',
                'max_inplane_principal': 'maxInPlanePrincipal',
                'max_principal': 'maxPrincipal',
                'min_principal': 'minPrincipal',
                'mid_principal': 'midPrincipal',
                'outofplane_principal': 'outOfPlanePrincipal',
                }
            # need to convert "min_inplane_principal" to "minInplanePrincipal" attribute
            for i in fo.validInvariants:
                attribute = str(i).lower()
                if attribute in switcher.keys():
                    attribute = switcher[attribute]
                validinv.append(attribute)

        # variable position
        pos = str(fo.values[0].position)

        # header for the .csv file identifying the node/integration point the value refers to.
        if pos == 'INTEGRATION_POINT':
            id = ['elementLabel']
            fo = fo.getSubset(position= CENTROID)
        elif pos == 'NODAL':
            id = ['nodeLabel']

        # second line of the .csv file
        csvhead = id + complabel + validinv

        # ex. csvhead = ['elementLabel', 'integrationPoint', 's11', 's22', 's33, 's12', ' 'mises',
        # 'max_inplane_principal', 'min_inplane_principal','outofplane_principal', 'max_principal',
        # 'mid_principal', 'min_principal', 'tresca', 'press', 'inv3']

        output = []
        for val in fo.values:
            if isinstance(val.data, float):
                data_complabel = [val.data]
            else:
                data_complabel = list(val.data)

            output.append([val.__getattribute__(param) for param in id] +
                          data_complabel +
                          [val.__getattribute__(param) for param in validinv])

        #store variable properties: first line in csv file
        description = [fo.name, fo.description, fo.type, 'position: ' + str(fo.values[0].position)] \
                      + frame_info

        finalcsv = [description] + [csvhead] + output
        return finalcsv

    def csvVariableWriter(self, variable, regiontuple = []):
        """ write output variable in multiple .csv files (one per frame)
         Input:
                - variable: Abaqus variable to write in .csv file (see list above)
                - regiontuple: tuple from set_element_list or set_node_list methods
         """

        if regiontuple == []:
            pos = str(self.frame[0].fieldOutputs[variable].values[0].position)
            if pos == 'INTEGRATION_POINT':
                region_name = 'ALL_ELEMENTS'
            elif pos == 'NODAL':
                region_name = 'ALL_NODES'
            region_set = self.instance
        else:
            # [0] set name, [1] actual set object
            region_name = regiontuple[0]
            region_set = regiontuple[1]

        directory_result = 'exported_results'
        if not os.path.exists(directory_result):
            os.makedirs(directory_result)

        directory_csv = directory_result + '\\csv_results'
        if not os.path.exists(directory_csv):
            os.makedirs(directory_csv)

        instname = self.instance_name
        directory_instance = directory_csv + '\\' + instname
        if not os.path.exists(directory_instance):
            os.makedirs(directory_instance)

        directory_region = directory_instance + '\\' + region_name.lower()
        if not os.path.exists(directory_region):
            os.makedirs(directory_region)

        directory_variable = directory_region + '\\' + variable.upper()
        if not os.path.exists(directory_variable):
            os.makedirs(directory_variable)

        directory_step = directory_variable + '\\' + self.step_name.lower()
        if not os.path.exists(directory_step):
            os.makedirs(directory_step)


        for fram in self.frame:

            print('Exporting variable ' + variable + ' at frame ' + '%d' % fram.incrementNumber + ' in .csv format  ...')

            info = [fram.description, 'Frame Value: ' + '%.2f' % fram.frameValue,
                    'Region Name: ' + region_name]

            fo = fram.fieldOutputs[variable].getSubset(region=region_set)

            csvoutput = self.csvVariableStore(fo, info)

            csvfilename = instname + '_' + region_name + '_' + variable + '_' + str(int(fram.frameValue)) + '.csv'
            csvfilepath = directory_step + '\\' + csvfilename
            with open(csvfilepath, 'wb') as writeFile:
                writer = csv.writer(writeFile)
                writer.writerows(csvoutput)

    def vtkVariableStore(self, fo):
        # store output variable as .vtk file
        # input: fielOutputs file (ex. fo = odbh.frame[0].fieldOutputs['S'])
        #        info of the current frame to store - list

        varname = fo.name

        # store the name of the variables in fieldOutputs.values.data
        complabel = fo.componentLabels
        if complabel == ():
            complabel = [fo.name]
        else:
            complabel = [i.lower() for i in complabel]

        # store the name of the variables in fieldOutputs.values such as mises, tresca, etc...
        validinv = []
        if fo.validInvariants != ():
            switcher = {
                'min_inplane_principal': 'minInPlanePrincipal',
                'max_inplane_principal': 'maxInPlanePrincipal',
                'max_principal': 'maxPrincipal',
                'min_principal': 'minPrincipal',
                'mid_principal': 'midPrincipal',
                'outofplane_principal': 'outOfPlanePrincipal',
                }
            # need to convert "min_inplane_principal" to "minInplanePrincipal" attribute
            for i in fo.validInvariants:
                attribute = str(i).lower()
                if attribute in switcher.keys():
                    attribute = switcher[attribute]
                validinv.append(attribute)

        # array names of the field data in .vtk file
        fieldarrays = complabel + validinv

        # ex. fieldarrays = ['s11', 's22', 's33', 's12', 'mises',
        # 'max_inplane_principal', 'min_inplane_principal','outofplane_principal', 'max_principal',
        # 'mid_principal', 'min_principal', 'tresca', 'press', 'inv3']

        n_tuples = len(fo.values)
        out_arrays = ['<DataArray type="Float32" Name="' + array + '" format="ascii">\n'
                      for array in fieldarrays]

        for value in fo.values:
            if isinstance(value.data, float):
                data_complabel = [value.data]
            else:
                data_complabel = list(value.data)

            val_complabel = [str(val) + '\n' for val in data_complabel]
            val_indiv = [str(value.__getattribute__(param)) + '\n' for param in validinv]

            val_out= val_complabel + val_indiv
            # element-wise addition between lists
            out_arrays = list(map(add, out_arrays, val_out))

        end_data_arrays = ['</DataArray>\n'] * len(fieldarrays)
        out_arrays = list(map(add, out_arrays, end_data_arrays))

        return ''.join(out_arrays)

    def vtkVariableStore_DataFile(self, fo):
        # store output variable as vtk DataFile
        # input: fielOutputs file (ex. fo = odbh.frame[0].fieldOutputs['S'])
        #        info of the current frame to store - list

        varname = fo.name

        # store the name of the variables in fieldOutputs.values.data
        complabel = fo.componentLabels
        if complabel == ():
            complabel = [fo.name]
        else:
            complabel = [i.lower() for i in complabel]

        # store the name of the variables in fieldOutputs.values such as mises, tresca, etc...
        validinv = []
        if fo.validInvariants != ():
            switcher = {
                'min_inplane_principal': 'minInPlanePrincipal',
                'max_inplane_principal': 'maxInPlanePrincipal',
                'max_principal': 'maxPrincipal',
                'min_principal': 'minPrincipal',
                'mid_principal': 'midPrincipal',
                'outofplane_principal': 'outOfPlanePrincipal',
                }
            # need to convert "min_inplane_principal" to "minInplanePrincipal" attribute
            for i in fo.validInvariants:
                attribute = str(i).lower()
                if attribute in switcher.keys():
                    attribute = switcher[attribute]
                validinv.append(attribute)

        # array names of the field data in .vtk file
        fieldarrays = complabel + validinv

        # ex. fieldarrays = ['s11', 's22', 's33', 's12', 'mises',
        # 'max_inplane_principal', 'min_inplane_principal','outofplane_principal', 'max_principal',
        # 'mid_principal', 'min_principal', 'tresca', 'press', 'inv3']

        n_tuples = len(fo.values)
        out_arrays = [array + ' 1 ' + str(n_tuples) + ' float\n' for array in fieldarrays]

        for value in fo.values:
            if isinstance(value.data, float):
                data_complabel = [value.data]
            else:
                data_complabel = list(value.data)

            val_complabel = [str(val) + '\n' for val in data_complabel]
            val_indiv = [str(value.__getattribute__(param)) + '\n' for param in validinv]

            val_out= val_complabel + val_indiv
            # element-wise addition between lists

            out_arrays = list(map(add, out_arrays, val_out))

        out_vtk = 'FIELD ' + varname + ' ' + str(len(fieldarrays)) + '\n' + ''.join(out_arrays)

        return out_vtk

    def export2vtk(self, variable, instances = []):
        """         This method exports output data variable to vtk file (XML format).
         WARNINGS: -only exports the instance defined by the OdbHelper object
                   -it assumes all elements have the same connectivity
        Inputs:
                   - variable keyword
                   - instances to store in .vtk file. Might be only the instance
                     selected when created the OdbHelper ojbect (instances property empty)
                     or all the instances in the model (instaces = 'all').
                     If only one istance is selected, the first elementSet/nodeSet
                     of the instance must include all elements/nodes
                     (sets are listed in alphabetical order by Abaqus)
                   """

        directory_result = 'exported_results'
        if not os.path.exists(directory_result):
            os.makedirs(directory_result)

        directory_vtk = directory_result + '\\vtk_results'
        if not os.path.exists(directory_vtk):
            os.makedirs(directory_vtk)

        inst = []
        if instances == 'all':
            instname = 'all_instances'
            inst = [i[1] for i in self.assembly.instances.items()]

        else:
            instname = self.instance_name
            inst = [self.instance]

        directory_instance = directory_vtk + '\\' + instname
        if not os.path.exists(directory_instance):
            os.makedirs(directory_instance)

        directory_step = directory_instance + '\\' + self.step_name.lower()
        if not os.path.exists(directory_step):
            os.makedirs(directory_step)

        directory_variable = directory_step + '\\' + variable.upper()
        if not os.path.exists(directory_variable):
            os.makedirs(directory_variable)

        for num, fram in enumerate(self.frame):

            print('Exporting variable ' + variable + ' at frame ' + '%d' % fram.incrementNumber + ' in .vtu format  ...')

            fo = fram.fieldOutputs[variable]

            # variable position
            pos = str(fo.values[0].position)

            out_pointdata = []
            out_celldata = []
            # header for the .csv file identifying the node/integration point the value refers to.
            if pos == 'INTEGRATION_POINT':
                # need to check how to properly assign integration point data in a vtk file.
                # At the moment, integration point data cannot be assigned.
                # In this code, data at integration points is interpolated and assigned to the centroid of the element.
                fo = fo.getSubset(position = CENTROID)
                for i in inst:
                    fo = fo.getSubset(region=i)
                    out_celldata.append(self.vtkVariableStore(fo))
            elif pos == 'NODAL':
                for i in inst:
                    fo = fo.getSubset(region=i)
                    out_pointdata.append(self.vtkVariableStore(fo))

            vtkfilename = instname + '_' + variable + '_' + '%04d' % num + '.vtu'
            vtkfilepath = directory_variable + '\\' + vtkfilename

            out = '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">\n'
            out += '\n<UnstructuredGrid>\n'
            out += '\n<FieldData>\n'
            out += '<DataArray type="Float64" Name="TIME" NumberOfTuples="1" format="ascii" ' + \
                   'RangeMin="'+str(fram.frameValue)+'" RangeMax="'+str(fram.frameValue)+'">\n'
            out += str(fram.frameValue) + '\n'
            out += '</DataArray>\n'
            out += '<DataArray type="Int32" Name="CYCLE" NumberOfTuples="1" format="ascii" ' + \
                   'RangeMin="'+str(fram.frameValue)+'" RangeMax="'+str(fram.frameValue)+'">\n'
            out += str(fram.incrementNumber) + '\n'
            out += '</DataArray>\n'
            out += '<DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii" ' + \
                   'RangeMin="'+str(fram.frameValue)+'" RangeMax="'+str(fram.frameValue)+'">\n'
            out += str(fram.frameValue) + '\n'
            out += '</DataArray>\n'
            out += '</FieldData>\n'

            # Store instance geometries in .vtk format

            nnodes = [] # number of nodes (of all instances selected)
            nelems = [] # number of elements (of all instances selected)
            out_nodes = ''  # node list in .vtk format (of all instances selected)
            node_label = [] # node label list
            out_elements = '' # elements list in .vtk format (of all instances selected)
            out_connect = '' # connectivity list in .vtk format (of all instances selected)
            out_offset = '' # vtk files in xml formats needs the offset to define cells
            offset = 0
            node_count = 0
            for i in inst:
                nnodes.append(len(i.nodes))
                nelems.append(len(i.elements))

                for node in i.nodes:
                    node_label.append(node.label)
                    out_nodes += ' '.join('%10.5E' % n for n in node.coordinates)
                    out_nodes += '\n'
                    
                connectivity = []
                for element in i.elements:
                    # IMPORTANT!!!: node index in paraview start in 0 (that's why i-1)
                    vtk_con = [int(e) - 1 + node_count for e in element.connectivity]
                    connectivity.append(vtk_con)
                    out_connect += str(self.element_vtk_cell_type[element.type]) + '\n'
                    offset += len(element.connectivity)
                    out_offset += str(offset) + '\n'

                node_count += int(nnodes[-1])
                new_label = list(range(len(node_label)))
                node_vtk = [n - 1 for n in node_label]
                # renumber nodes in connectivity matrix if Abaqus node labels and node order do not correspond
                if node_vtk != new_label:
                    nel = len(connectivity)
                    import numpy as np
                    cmatrix = np.asarray(connectivity)

                    def matrix_replace(a, val_old, val_new):
                        d = dict(zip(val_old, val_new))
                        return [d.get(e, e) for e in a]

                    new_connectivity_flat = matrix_replace(cmatrix.flatten(), node_vtk, new_label)
                    connectivity = np.reshape(new_connectivity_flat, (nel, -1)).tolist()

                out_connectivity = [''.join(str(x)) for x in connectivity]
                out_elements += '\n'.join(out_connectivity).replace('[', '').replace(']', '').replace(',', '')
                out_elements += '\n'

            # Report Abaqus nodel labels as string
            out_node_label_str = '<DataArray type="Float32" Name="Abaqus_node_label" format="ascii">\n'
            out_node_label_str += '\n'.join(map(str, node_label))
            out_node_label_str += '</DataArray>\n'

            # Writing instance geometries in final .vtk file

            out += '\n<Piece NumberOfPoints = "' + str(sum(nnodes)) + '" NumberOfCells = "' + str(sum(nelems)) + '">\n'

            out += '\n<Points>\n'
            out += '<DataArray type="Float64" NumberOfComponents=" ' + str(len(inst[0].nodes[0].coordinates)) + '" format="ascii">\n'
            out += out_nodes
            out += '</DataArray>\n'
            out += '</Points>\n'

            out += '\n<Cells>\n'

            out += '<DataArray type = "Int64" Name = "connectivity" format = "ascii">\n'
            out += out_elements
            out += '</DataArray>\n'

            out += '<DataArray type="Int64" Name="offsets" format="ascii">'
            out += out_offset
            out += '</DataArray>\n'

            out += '<DataArray type="UInt8" Name="types" format="ascii">'
            out += out_connect
            out += '</DataArray>\n'

            out += '\n</Cells>\n'

            out += '<CellData>\n'
            out += ''.join(out_celldata)
            out += '</CellData>\n'

            out += '\n<PointData>\n'
            out += out_node_label_str
            out += ''.join(out_pointdata)
            out += '</PointData>\n'

            out += '\n</Piece>\n</UnstructuredGrid>\n</VTKFile>'

            with open(vtkfilepath, 'w') as f:
                f.write(out)

    def export2vtk_DataFile(self, variable, instances = []):
        """         This method exports output data variable to vtk file (DataFile format).
         WARNINGS: -only exports the instance defined by the OdbHelper object
                   -it assumes all elements have the same connectivity
        Inputs:
                   - variable keyword
                   - instances to store in .vtk file. Might be only the instance
                     selected when created the OdbHelper ojbect (instances property empty)
                     or all the instances in the model (instaces = 'all').
                     If only one istance is selected, the first elementSet/nodeSet
                     of the instance must include all elements/nodes
                     (sets are listed in alphabetical order by Abaqus)
                   """

        directory_result = 'exported_results'
        if not os.path.exists(directory_result):
            os.makedirs(directory_result)

        directory_vtk = directory_result + '\\vtk_results'
        if not os.path.exists(directory_vtk):
            os.makedirs(directory_vtk)

        inst = []
        if instances == 'all':
            instname = 'all_instances'
            inst = [i[1] for i in self.assembly.instances.items()]

        else:
            instname = self.instance_name
            inst = [self.instance]

        directory_instance = directory_vtk + '\\' + instname
        if not os.path.exists(directory_instance):
            os.makedirs(directory_instance)

        directory_step = directory_instance + '\\' + self.step_name.lower()
        if not os.path.exists(directory_step):
            os.makedirs(directory_step)

        directory_variable = directory_step + '\\' + variable.upper()
        if not os.path.exists(directory_variable):
            os.makedirs(directory_variable)

        for num, fram in enumerate(self.frame):

            print('Exporting variable ' + variable + ' at frame ' + '%d' % num + 'in DataFile.vtu format  ...')

            fo = fram.fieldOutputs[variable]

            # variable position
            pos = str(fo.values[0].position)

            out_pointdata = []
            out_celldata = []
            # header for the .csv file identifying the node/integration point the value refers to.
            if pos == 'INTEGRATION_POINT':
                # need to check how to properly assign integration point data in a vtk file.
                # At the moment, integration point data cannot be assigned.
                # In this code, data at integration points is interpolated and assigned to the centroid of the element.
                fo = fo.getSubset(position=CENTROID)

            for i in inst:
                fo = fo.getSubset(region=i)
                out_pointdata.append(self.vtkVariableStore(fo))

            vtkfilename = instname + '_' + variable + '_' + '%04d' % num + '.vtk'
            vtkfilepath = directory_variable + '\\' + vtkfilename

            out = '# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n'

            # Store instance geometries in .vtk format

            nnodes = [] # number of nodes (of all instances selected)
            nelems = [] # number of elements (of all instances selected)
            conn = len(self.instance.elements[0].connectivity) # assumed the same for all instances
            out_nodes = '' # node list in .vtk format (of all instances selected)
            node_label = [] # node label list
            out_elements = '' # elements list in .vtk format (of all instances selected)
            out_connect = '' # connectivity list in .vtk format (of all instances selected)
            node_count = 0
            for i in inst:
                nnodes.append(len(i.nodes))
                nelems.append(len(i.elements))

                for node in i.nodes:
                    node_label.append(node.label)
                    out_nodes += ' '.join('%10.5E' % n for n in node.coordinates)
                    out_nodes += '\n'

                for element in i.elements:
                    # IMPORTANT!!!: node index in paraview start in 0 (that's why i-1)
                    vtk_con = [int(e) - 1 + node_count for e in element.connectivity]
                    out_elements += str(conn) + ' ' + ' '.join('%8d' % c for c in vtk_con)
                    out_elements += '\n'

                    out_connect += str(self.element_vtk_cell_type[element.type]) + '\n'

                node_count += int(nnodes[-1])

            # Writing instance geometries in final .vtk file

            out += 'POINTS {0} float\n'.format(sum(nnodes))
            out += out_nodes

            out += 'CELLS {0} {1}\n'.format(sum(nelems), sum(nelems) * (conn + 1))
            out += out_elements

            out += 'CELL_TYPES {0}\n'.format(sum(nelems))
            out += out_connect

            out += 'CELL_DATA {0}\n'.format(sum(nelems))
            out += ''.join(out_celldata)

            out += 'POINT_DATA {0}\n'.format(sum(nnodes))
            out += ''.join(out_pointdata)

            with open(vtkfilepath, 'w') as f:
                f.write(out)
