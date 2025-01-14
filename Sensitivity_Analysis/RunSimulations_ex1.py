from mpi4py import MPI 
from uq_physicell.uq_physicell import PhysiCell_Model, get_xml_element_value, get_rule_index_in_csv
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from shutil import rmtree
import pcdl
from SALib import ProblemSpec
import os


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def custom_summary_func(OutputFolder,SummaryFile, dic_params, SampleID, ReplicateID):
    mcds_ts = pcdl.TimeSeries(OutputFolder, microenv=True, graph=False, settingxml=None, physiboss=False, verbose=False)
    for mcds in mcds_ts.get_mcds_list():
        df_cell = mcds.get_cell_df() 
        df_tumor_live = df_cell[ (df_cell['cell_type'] == 'tumor') & (df_cell['dead'] == False) ]
        df_tumor_dead = df_cell[ (df_cell['cell_type'] == 'tumor') & (df_cell['dead'] == True) ]
        df_motile_live = df_cell[ (df_cell['cell_type'] == 'motile_tumor') & (df_cell['dead'] == False) ]
        df_motile_dead = df_cell[ (df_cell['cell_type'] == 'motile_tumor') & (df_cell['dead'] == True) ]
        # distance from the center (0,0,0)
        distances_tumor_live = np.sqrt( (df_tumor_live['position_x'])**2 + (df_tumor_live['position_y'])**2 + (df_tumor_live['position_z'])**2 )
        distances_tumor_dead = np.sqrt( (df_tumor_dead['position_x'])**2 + (df_tumor_dead['position_y'])**2 + (df_tumor_dead['position_z'])**2 )
        distances_motile_live = np.sqrt( (df_motile_live['position_x'])**2 + (df_motile_live['position_y'])**2 + (df_motile_live['position_z'])**2 )
        distances_motile_dead = np.sqrt( (df_motile_dead['position_x'])**2 + (df_motile_dead['position_y'])**2 + (df_motile_dead['position_z'])**2 )
        data = {'time': mcds.get_time(), 'replicate': ReplicateID, 'sample': SampleID, 'runtime': mcds.get_runtime(),
                'tumor_live': len(df_tumor_live), 'tumor_dead': len(df_tumor_dead), 'motile_live': len(df_motile_live), 'motile_dead': len(df_motile_dead),
                'dist_tumor_live': distances_tumor_live.to_numpy(), 'dist_tumor_dead': distances_tumor_dead.to_numpy(),
                'dist_motile_live': distances_motile_live.to_numpy(), 'dist_motile_dead': distances_motile_dead.to_numpy()}      
        data_conc = {**data,**dic_params} 
        if ( mcds.get_time() == 0 ): 
            # create the dataframe
            df = pd.DataFrame([data_conc]) 
        else:
            # append the dictionary to the dataframe
            df = pd.concat([df, pd.DataFrame([data_conc])])
    # remove replicate output folder
    rmtree( OutputFolder )
    df.to_feather(SummaryFile.replace('.csv','.feather'))
    return

if __name__ == '__main__':
    PhysiCellModel = PhysiCell_Model("Sensitivity_Analysis/ConfigFile.ini", 'model_hypoxia')
    # Define reference values for the parameters
    ref_parameters = {'xml': [], 'rules': []}
    # Parameters from xml
    tree = ET.parse(PhysiCellModel.configFile_ref)
    xml_root = tree.getroot()
    for key_xml in PhysiCellModel.keys_variable_params:
        text_elem = get_xml_element_value(xml_root, key_xml)
        ref_parameters['xml'].append(float(text_elem))
    # Parameters from rules
    for key_rule, list_rule in PhysiCellModel.parameters_rules.items():
        id_rule = get_rule_index_in_csv(PhysiCellModel.rules, key_rule)
        parameter_rule = key_rule.split(',')[-1]
        ref_parameters['rules'].append( float(PhysiCellModel.rules[id_rule][parameter_rule]) )
    ref_parameters_xml = np.array(ref_parameters['xml'])
    ref_parameters_rules = np.array(ref_parameters['rules'])
    print('Reference parameters XML: ', ref_parameters_xml, '\nReference parameters Rules: ', ref_parameters_rules)

    # Define the samples for the local sensitivity analysis +- 1%, 5%, 10%, 20%
    parameterSamplesXML = [ref_parameters_xml] # reference value first sample
    parameterSamplesRules = [ref_parameters_rules] # reference value first sample
    # For XML parameters
    for i in range(len(ref_parameters_xml)):
        # change only one parameter at a time
        for var in [-0.01,0.01,-0.05,0.05,-0.1,0.1,-0.2,0.2]:
            sample_xml = ref_parameters_xml.copy()
            sample_xml[i] = ref_parameters_xml[i] * (1 + var)
            parameterSamplesXML.append(sample_xml)
            parameterSamplesRules.append(ref_parameters_rules.copy())
    # For Rules parameters
    for i in range(len(ref_parameters_rules)):
        # change only one parameter at a time
        for var in [-0.01,0.01,-0.05,0.05,-0.1,0.1,-0.2,0.2]:
            sample_rules = ref_parameters_rules.copy()
            sample_rules[i] = ref_parameters_rules[i] * (1 + var)
            parameterSamplesRules.append(sample_rules)
            parameterSamplesXML.append(ref_parameters_xml.copy())
    
    # Size of parameterSamples and parameterSamplesRules are the same
    parameterSamplesXML = np.array(parameterSamplesXML)
    print('Number of samples: ', len(parameterSamplesXML))

    # Generate a three list with size NumSimulations = len(Samples) or len(Replicates or len(Parameters)
    ParametersXML = []; ParametersRules = []; Samples = []; Replicates = []
    for sampleID in range(len(parameterSamplesXML)):
        for replicateID in np.arange(PhysiCellModel.numReplicates):
            # check if the file already exists
            if ( os.path.isfile(PhysiCellModel.outputs_folder+'SummaryFile_%06d_%02d.feather'%(sampleID,replicateID)) ) : continue
            ParametersXML.append(parameterSamplesXML[sampleID])
            ParametersRules.append(parameterSamplesRules[sampleID])
            Samples.append(sampleID)
            Replicates.append(replicateID)
    
    SplitIndexes = np.array_split(np.arange(len(Samples)),size, axis=0) # split [0,1,...,NumSimulations-1] in size arrays equally (or +1 in some ranks)
    
    print(f"Total number of simulations: {len(Samples)} Simulations per rank: {int(len(Samples)/size)}")

    for ind_sim in SplitIndexes[rank]:
        print('Rank: ',rank, ', Simulation: ', ind_sim, ', Sample: ', Samples[ind_sim],', Replicate: ', Replicates[ind_sim])
        PhysiCellModel.RunModel(Samples[ind_sim], Replicates[ind_sim], Parameters=ParametersXML[ind_sim], ParametersRules = ParametersRules[ind_sim],SummaryFunction=custom_summary_func)

    MPI.Finalize()
