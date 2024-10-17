from mpi4py import MPI 
from uq_physicell.uq_physicell import PhysiCell_Model, get_rule_index_in_csv
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
    mcds_ts = pcdl.TimeSeries(OutputFolder, microenv=True, graph=False, settingxml=None, verbose=False)
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
    ref_parameters = []
    for key_rule, list_rule in PhysiCellModel.parameters_rules.items():
        id_rule = get_rule_index_in_csv(PhysiCellModel.rules, key_rule)
        parameter_rule = key_rule.split(',')[-1]
        ref_parameters.append(float(PhysiCellModel.rules[id_rule][parameter_rule])*0.8)
    ref_parameters = np.array(ref_parameters)
    print('Reference parameters: ', ref_parameters)

    # Define the samples for the local sensitivity analysis +- 1%, 5%, 10%, 20%
    parameterSamples = [ref_parameters] # reference value first sample
    for i in range(len(ref_parameters)):
        # change only one parameter at a time
        for var in [-0.01,0.01,-0.05,0.05,-0.1,0.1,-0.2,0.2]:
            sample = ref_parameters.copy()
            sample[i] = ref_parameters[i]*(1+var)
            parameterSamples.append(sample)
    parameterSamples = np.array(parameterSamples)
    print('Number of samples: ', len(parameterSamples))

    # Generate a three list with size NumSimulations = len(Samples) or len(Replicates or len(Parameters)
    Parameters = []; Samples = []; Replicates = []
    for sampleID in range(len(parameterSamples)):
        for replicateID in np.arange(PhysiCellModel.numReplicates):
            # check if the file already exists
            if ( os.path.isfile(PhysiCellModel.outputs_folder+'SummaryFile_%06d_%02d.feather'%(sampleID,replicateID)) ) : continue
            Parameters.append(parameterSamples[sampleID])
            Samples.append(sampleID)
            Replicates.append(replicateID)
    
    SplitIndexes = np.array_split(np.arange(len(Samples)),size, axis=0) # split [0,1,...,NumSimulations-1] in size arrays equally (or +1 in some ranks)
    
    print(f"Total number of simulations: {len(Samples)} Simulations per rank: {int(len(Samples)/size)}")

    for ind_sim in SplitIndexes[rank]:
        print('Rank: ',rank, ', Simulation: ', ind_sim, ', Sample: ', Samples[ind_sim],', Replicate: ', Replicates[ind_sim])
        PhysiCellModel.RunModel(Samples[ind_sim], Replicates[ind_sim], Parameters=np.array([]), ParametersRules = Parameters[ind_sim],SummaryFunction=custom_summary_func)

    MPI.Finalize()
