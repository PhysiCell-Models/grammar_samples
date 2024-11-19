# DRB python script to do Nelder-Mead for pc-cortex

import numpy as np
import pandas as pd
import scipy.io as sio
import xml.etree.ElementTree as ET

from scipy.optimize import minimize

import os
import subprocess
import time
import copy

home_dir = os.path.expanduser("~")
path_to_physicell = "." # e.g. path_to_physicell = home_dir + "/grammar_samples"
path_to_sbatch = f"{path_to_physicell}/user_projects/neuro_dev/Calibration/cortex_calibration.sbat"

layer_celldef_id = {6: 4, 5: 5, 4: 6, 3: 7, 2: 9}
layer_list = [2, 3, 4, 5, 6, 'total']


user_name = input("Please enter your username: ")
region = "AUD"

# initial values for the parameters based on the region being parameterized
if region=="AUD":
    original_start_time = {6: 1440.0, 5: 3100.0, 4: 6700.0, 3: 7800.0}
    initial_rgc_ec50 = 4700.0
    layer_counts_data = {2: 338, 4: 121, 5: 362, 6: 220} # AUD in 1.086
elif region=="SOM":
    original_start_time = {6: 1440.0, 5: 5200.0, 4: 7200.0, 3: 10800.0}
    initial_rgc_ec50 = 6200.0
    layer_counts_data = {2: 331, 4: 275, 5: 190, 6: 441} # SOM in ????
else:
    raise ValueError(f"Error: region {region} not recognized")

layer_counts_data['total'] = sum(layer_counts_data.values())

def main():
    result, order = optimizeParameters(layer_counts_data, min_replicates=1, maxfev = 100)
    print(result)
    print(order)

    # print results and order to out.txt
    with open(f"{region}Results.txt", "w") as f:
        f.write("Result:\n")
        f.write(str(result))
        f.write("\n")
        f.write("Parameter Order:\n")
        f.write(str(order))

def layerCountsAtEnd( path_to_output ):
    if not os.path.exists(path_to_output):
        raise FileNotFoundError(f"Error: output folder {path_to_output} does not exist")
    # read the layer counts from the output file
    path_to_filename = path_to_output + '/final_cells.mat'
    cells_df = pd.DataFrame(sio.loadmat(path_to_filename)['cells'].T, index=None, columns=None)

    # count 
    layer_counts_simulated = {}
    for layer, celldef_id in layer_celldef_id.items():
        layer_counts_simulated[layer] = (cells_df[5] == celldef_id).sum()

    return layer_counts_simulated

def initializeParameters():
    parameters = {}
    # rgc ec50 of cycle rate decrease due to time
    initialzeRGCEC50(parameters)
    for layer_start in [6, 5, 4]:
        initializeGap(parameters, layer_start)
    for p in parameters.values():
        p["current_value"] = p["initial_value"] # ensure that the current value is set to the initial value when initializing

    return parameters

def initialzeRGCEC50(parameters):
    p = {}
    p["name"] = "rgc_ec50"
    p["initial_value"] = initial_rgc_ec50
    p["min_value"] = 0.0
    p["max_value"] = 10000.0
    p["min_step_size"] = 1.0
    p["max_step_size"] = 1000.0
    sub = {}
    sub["location"] = "rules"
    sub["line_number_0_based"] = 1
    sub["col_number_0_based"] = 5
    p["subs"] = [sub]
    parameters[p["name"]] = p

def initializeGap(parameters, layer_start):
    # these "gaps" span the time layer i is forming before i-1 begins to form. at the end of the gap, the asymmetric division probabilities shift towards the next layer
    # layer_start is i
    p = {}
    p["name"] = f"gap_{layer_start}"
    p["initial_value"] = original_start_time[layer_start-1] - original_start_time[layer_start]
    p["min_value"] = 0.0
    p["max_value"] = 10000.0
    p["min_step_size"] = 1.0
    p["max_step_size"] = 1000.0
    p["subs"] = setUpTimeSubs(layer_start)
    parameters[p["name"]] = p

def getPreviousEndTime(config_tree, layer_start):
    previous_start_time = config_tree.find(".//user_parameters").find(f".//layer_{layer_start+1}_end_time").text
    return float(previous_start_time)
    
def setUpTimeSubs(layer_start):
    layer_6_decrease_line_0_based = 3 # first line to have asymmetric division rule. they will go layer_6, layer_5, etc. and within each layer go increase then decrease
    prev_layer_decrease_line_0_based = layer_6_decrease_line_0_based - 2*layer_start + 2*6

    sub_previous_layer_decrease = {}
    sub_previous_layer_decrease["location"] = "rules"
    sub_previous_layer_decrease["line_number_0_based"] = prev_layer_decrease_line_0_based
    sub_previous_layer_decrease["col_number_0_based"] = 5

    sub_next_layer_increase = {}
    sub_next_layer_increase["location"] = "rules"
    sub_next_layer_increase["line_number_0_based"] = prev_layer_decrease_line_0_based + 1
    sub_next_layer_increase["col_number_0_based"] = 5

    if layer_start == 6:
        sub_previous_layer_decrease["fn"] = lambda _, delta : 1440.0 + delta
        sub_next_layer_increase["fn"] = lambda _, delta : 1440.0 + delta
    elif layer_start == 4:
        sub_previous_layer_decrease["fn"] = lambda df, delta : float(df.iloc[prev_layer_decrease_line_0_based - 1, 5]) + delta
        sub_next_layer_increase["fn"] = lambda df, delta : float(df.iloc[prev_layer_decrease_line_0_based - 1, 5]) + delta
        # set the transition from 3->2 at the halfway point since we don't have data separating 2 and 3
        sub_layer_3_decrease = {}
        sub_layer_3_decrease["location"] = "rules"
        sub_layer_3_decrease["line_number_0_based"] = prev_layer_decrease_line_0_based + 2
        sub_layer_3_decrease["col_number_0_based"] = 5
        sub_layer_3_decrease["fn"] = lambda df, delta : 0.5 * ((float(df.iloc[prev_layer_decrease_line_0_based - 1, 5]) + delta) + 9*1440) # layer 3 stops halfway from the start of it to the end of the layer formation (at t=9*1440)

        sub_layer_2_increase = {}
        sub_layer_2_increase["location"] = "rules"
        sub_layer_2_increase["line_number_0_based"] = prev_layer_decrease_line_0_based + 3
        sub_layer_2_increase["col_number_0_based"] = 5
        sub_layer_2_increase["fn"] = lambda df, delta : 0.5 * ((float(df.iloc[prev_layer_decrease_line_0_based - 1, 5]) + delta) + 9*1440) # layer 2 starts halfway from the start of layer 3 to the end of layer 3 (at t=9*1440)
        return [sub_previous_layer_decrease, sub_next_layer_increase, sub_layer_3_decrease, sub_layer_2_increase]
    else:
        sub_previous_layer_decrease["fn"] = lambda df, delta : float(df.iloc[prev_layer_decrease_line_0_based - 1, 5]) + delta
        sub_next_layer_increase["fn"] = lambda df, delta : float(df.iloc[prev_layer_decrease_line_0_based - 1, 5]) + delta

    return [sub_previous_layer_decrease, sub_next_layer_increase]

def writeNewParameters( parameters, idx ):
    path_to_config = f"{path_to_physicell}/config/PhysiCell_settings.xml"
    if os.path.exists(path_to_config) == False:
        raise FileNotFoundError(f"Error: config file {path_to_config} does not exist")
    path_to_rules = f"{path_to_physicell}/config/cell_rules.csv"
    if os.path.exists(path_to_rules) == False:
        raise FileNotFoundError(f"Error: rules file {path_to_rules} does not exist")
    config_tree = ET.parse(path_to_config)
    rules_df = pd.read_csv( path_to_rules, dtype = str, header = None)

    p_values_list = copy.deepcopy(list(parameters.keys())) # deepcopy just to make sure I'm not changing the original list
    p_values_list.sort(reverse=True) # make sure to start with gap 6 then gap 5 etc. rgc ec50 can be anywhere
    for par_name in p_values_list:
        p = parameters[par_name]
        for sub in p["subs"]:
            if sub["location"] == "rules":
                if "fn" in sub:
                    new_value = sub["fn"](rules_df, p["current_value"])
                else:
                    new_value = p["current_value"]
                rules_df.iloc[sub["line_number_0_based"], sub["col_number_0_based"]] = str(new_value)
            elif sub["location"] == "config":
                if "fn" in sub:
                    new_value = sub["fn"](config_tree, p["current_value"])
                else:
                    new_value = p["current_value"]
                final_node = config_tree
                xml_path = sub["xml_path"].split("//")
                for xp in xml_path:
                    tokens = xp.split(":")
                    if len(tokens)==1:
                        final_node = final_node.find(f".//{tokens[0]}")
                        continue
                    for var in final_node.findall(tokens[0]):
                        if var.attrib[tokens[1]] == tokens[2]:
                            final_node = var
                            break
                final_node.text = str(new_value)
            else:
                print(f"Error: unknown location {sub['location']} for parameter {sub['name']}")

    path_to_new_config = f"{path_to_config[:-4]}_{idx}.xml"
    path_to_new_rules = f"{path_to_rules[:-4]}_{idx}.csv"
    path_to_new_rules_folder = os.path.dirname(path_to_new_rules)
    filename_new_rules = os.path.basename(path_to_new_rules)

    ruleset_node = config_tree.find(".//cell_rules").find("rulesets").find("ruleset")
    ruleset_node.find("folder").text = path_to_new_rules_folder
    ruleset_node.find("filename").text = filename_new_rules

    save_folder_node = config_tree.find("save").find("folder")
    save_folder_node.text = f"output_{idx}"

    config_tree.write(path_to_new_config)
    rules_df.to_csv(path_to_new_rules, header = False, index = False)

def jobIDInQueue( jobid ):
    # check if jobid is in the queue for the user to decide if it is time to process or just keep waiting
    queue_line_split = subprocess.check_output(["squeue", "-u", user_name, "-j", jobid]).decode("utf-8").split(" ")
    for i, line in enumerate(queue_line_split):
        if i==0:
            continue # first line is headings
        line_words = line.split(" ")
        for word in line_words:
            # check if word starts with jobid
            if word.startswith(jobid):
                return True
    return False

def cleanUpSimulations( min_replicates ):
    for i in range(min_replicates):
        path_to_output = f"{path_to_physicell}/output_{i}"
        os.system(f"rm -rf {path_to_output}")

def printCurrentParameters( parameters ):
    print("Current parameters:")
    for p in parameters.values():
        print(f"{p['name']}: {p['current_value']}")
    print("\n", flush=True)

def runSimulationsAndError( x, parameters, layer_counts_data, parameter_order, min_replicates ):
    for name, i in parameter_order.items():
        parameters[name]["current_value"] = x[i]

    printCurrentParameters(parameters)
    
    for i in range(min_replicates):
        writeNewParameters(parameters, i)
        path_to_output = f"{path_to_physicell}/output_{i}"
        os.system(f"rm -rf {path_to_output}")
        os.makedirs(path_to_output, exist_ok = True)
        
    time.sleep(1) # make sure the files are written before running the simulations

    runReplicates(range(min_replicates))

    # check if layer_counts_data is a dictionary
    if isinstance(layer_counts_data, dict):
        layer_counts_simulated_all = {i: [] for i in layer_celldef_id.keys()}
        layer_counts_simulated_all['total'] = []
        for i in range(min_replicates):
            path_to_output = f"{path_to_physicell}/output_{i}"
            layer_counts_simulated = layerCountsAtEnd(path_to_output)
            layer_counts_simulated['total'] = sum(layer_counts_simulated.values())
            for layer, count in layer_counts_simulated.items():
                layer_counts_simulated_all[layer].append(count)
        
        # calculate mean for each layer
        simulated_layer_means = {i: np.mean(counts) for i, counts in layer_counts_simulated_all.items()}
        simulated_layer_stds = {i: np.std(counts) for i, counts in layer_counts_simulated_all.items()}
        temp_dict = {i: (simulated_layer_means[i], simulated_layer_stds[i]) for i in layer_list}
        for layer in layer_list:
            print(f"- Layer {layer}:")
            print(f"\tcounts: {layer_counts_simulated_all[layer]}")
            print(f"\t(mean, std): {temp_dict[layer]}\n")
        if 3 not in layer_counts_data.keys():
            # in this case, layers 2 and 3 are combined in the data
            simulated_layer_means[2] = simulated_layer_means[2] + simulated_layer_means[3]
        # calculate error
        temp_dict = {i: (simulated_layer_means[i], layer_counts_data[i], simulated_layer_means[i] - layer_counts_data[i]) for i in layer_counts_data.keys()}
        print("\n(Model, Data, Model - Data):")
        for layer, v in temp_dict.items():
            print(f"\t- Layer {layer}: {v}")
        # important that we loop over the data keys, since the simulated keys will include 2 and 3
        # since we our now adding the total to the dict, we can ignore one of the layer errors; choose 2 since it is affected by the most parameters
        error = sum([np.square(simulated_layer_means[i] - layer_counts_data[i]) for i in layer_counts_data.keys() if i!=2])
    else: # assume that it is the total cell count
        layer_counts_simulated_all = []
        for i in range(min_replicates):
            path_to_output = f"{path_to_physicell}/output_{i}"
            layer_counts_simulated = layerCountsAtEnd(path_to_output)
            layer_counts_simulated_all.append(sum(layer_counts_simulated.values()))
        
        layer_means_simulated = np.mean(layer_counts_simulated_all)
        # calculate error
        error = np.square(layer_means_simulated - layer_counts_data)
    
    print(f"Error: {error}\n", flush=True)
    return error

def runReplicates(replicate_ids):
    setUpSlurmScript(replicate_ids)
    
    # run sbatch script
    job = subprocess.check_output(["sbatch", path_to_sbatch])
    # get jobid
    jobid = job.decode("utf-8").split(" ")[-1].strip()

    while jobIDInQueue(jobid) == True:
        time.sleep(5)

    ids_left = []
    for id in replicate_ids:
        path_to_output = f"{path_to_physicell}/output_{id}/final_cells.mat"
        if not os.path.exists(path_to_output):
            ids_left.append(id)

    if len(ids_left) > 0:
        runReplicates(ids_left)
        
def setUpSlurmScript( replicate_ids ):
    with open(path_to_sbatch, "r") as f:
        lines = f.readlines()
    
    # find line beginning with #SBATCH --array
    for i, line in enumerate(lines):
        if line.startswith("#SBATCH --array"):
            break

    # replace the number of jobs to go 0 thru min_replicates
    id_str = ",".join(map(str, replicate_ids))
    lines[i] = f"#SBATCH --array={id_str}\n"

    with open(path_to_sbatch, "w") as f:
        f.writelines(lines)

def optimizeParameters( layer_counts_data, initial_parameters = None, min_replicates = 6, maxfev = None ):
    setUpSlurmScript(range(min_replicates))
    if initial_parameters is None:
        initial_parameters = initializeParameters()
    parameter_order = {}
    x0 = np.zeros(len(initial_parameters))
    mins = np.zeros(len(initial_parameters))
    maxs = np.zeros(len(initial_parameters))
    for i, p in enumerate(initial_parameters.values()):
        parameter_order[p["name"]] = i
        x0[i] = p["current_value"]
        mins[i] = p["min_value"]
        maxs[i] = p["max_value"]
    f = lambda x : runSimulationsAndError(x, initial_parameters, layer_counts_data, parameter_order, min_replicates)
    bounds = [(low, high) for low, high in zip(mins, maxs)]
    options = {}
    if maxfev is not None:
        options["maxfev"] = maxfev
    result = minimize(f, x0, method='Nelder-Mead', bounds=bounds, options=options)
    return result, parameter_order
    
if __name__ == "__main__":
    main()