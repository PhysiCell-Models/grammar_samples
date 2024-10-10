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
        df_m2_live = df_cell[ (df_cell['cell_type'] == 'M2_macrophage') & (df_cell['dead'] == False) ]
        df_m2_dead = df_cell[ (df_cell['cell_type'] == 'M2_macrophage') & (df_cell['dead'] == True) ]
        # distance from the center (0,0,0)
        distances_tumor_live = np.sqrt( (df_tumor_live['position_x'])**2 + (df_tumor_live['position_y'])**2 + (df_tumor_live['position_z'])**2 )
        distances_tumor_dead = np.sqrt( (df_tumor_dead['position_x'])**2 + (df_tumor_dead['position_y'])**2 + (df_tumor_dead['position_z'])**2 )
        distances_m2_live = np.sqrt( (df_m2_live['position_x'])**2 + (df_m2_live['position_y'])**2 + (df_m2_live['position_z'])**2 )
        distances_m2_dead = np.sqrt( (df_m2_dead['position_x'])**2 + (df_m2_dead['position_y'])**2 + (df_m2_dead['position_z'])**2 )
        data = {'time': mcds.get_time(), 'replicate': ReplicateID, 'sample': SampleID, 'runtime': mcds.get_runtime(),
                'tumor_live': len(df_tumor_live), 'tumor_dead': len(df_tumor_dead), 'm2_live': len(df_m2_live), 'm2_dead': len(df_m2_dead),
                'dist_mean_tumor_live': distances_tumor_live.mean(), 'dist_mean_tumor_dead': distances_tumor_dead.mean(),
                'dist_mean_m2_live': distances_m2_live.mean(), 'dist_mean_m2_dead': distances_m2_dead.mean(),
                'dist_std_tumor_live': distances_tumor_live.std(), 'dist_std_tumor_dead': distances_tumor_dead.std(),
                'dist_std_m2_live': distances_m2_live.std(), 'dist_std_m2_dead': distances_m2_dead.std()}      
        data_conc = {**data,**dic_params} 
        if ( mcds.get_time() == 0 ): 
            # create the dataframe
            df = pd.DataFrame([data_conc]) 
        else:
            # append the dictionary to the dataframe
            df = pd.concat([df, pd.DataFrame([data_conc])])
    # remove replicate output folder
    rmtree( OutputFolder )
    df.to_csv(SummaryFile, sep='\t', encoding='utf-8')
    return

if __name__ == '__main__':
    PhysiCellModel = PhysiCell_Model("Sensitivity_Analysis/ConfigFile.ini", 'model_tam_egf')
    # Define parameters of rules in range +/- 20% of the reference value
    names_parameters = []
    bounds_parameters = []
    for key_rule, list_rule in PhysiCellModel.parameters_rules.items():
        id_rule = get_rule_index_in_csv(PhysiCellModel.rules, key_rule)
        parameter_rule = key_rule.split(',')[-1]
        names_parameters.append(list_rule[1])
        bounds_parameters.append([float(PhysiCellModel.rules[id_rule][parameter_rule])*0.8, float(PhysiCellModel.rules[id_rule][parameter_rule])*1.2])
        # id_rule, name, parametername (saturation, half max, or hill power), value of reference
        # print(id_rule, list_rule[1], parameter_rule, PhysiCellModel.rules[id_rule][parameter_rule])
    
    # Define SA problem
    sa_sobol = ProblemSpec({'names': names_parameters, 'bounds': bounds_parameters})
    
    # Sample parameters 
    sa_sobol.sample_sobol(2**6, calc_second_order=True, seed=42) # False: N*(D+2) True: N*(2D+2)
    
    # Generate a three list with size NumSimulations = len(Samples) or len(Replicates or len(Parameters)
    Parameters = []; Samples = []; Replicates = []
    for sampleID in range(sa_sobol.samples.shape[0]):
        for replicateID in np.arange(PhysiCellModel.numReplicates):
            # check if the file already exists
            if ( os.path.isfile(PhysiCellModel.outputs_folder+'SummaryFile_%06d_%02d.csv'%(sampleID,replicateID)) ) : continue
            Parameters.append(sa_sobol.samples[sampleID])
            Samples.append(sampleID)
            Replicates.append(replicateID)
    
    SplitIndexes = np.array_split(np.arange(len(Samples)),size, axis=0) # split [0,1,...,NumSimulations-1] in size arrays equally (or +1 in some ranks)
    
    print(f"Total number of simulations: {len(Samples)} Simulations per rank: {int(len(Samples)/size)}")

    for ind_sim in SplitIndexes[rank]:
        print('Rank: ',rank, ', Simulation: ', ind_sim, ', Sample: ', Samples[ind_sim],', Replicate: ', Replicates[ind_sim])
        PhysiCellModel.RunModel(Samples[ind_sim], Replicates[ind_sim], Parameters=np.array([]), ParametersRules = Parameters[ind_sim],SummaryFunction=custom_summary_func)

    MPI.Finalize()
