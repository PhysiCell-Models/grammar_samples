import numpy as np
from physicool.config import ConfigFileParser
from LF_MCMC import ABC_MCMC
import subprocess
import sys, os
import scipy.io as sio

# cd ./documents/github/PC_inverse/cell_rules2/physicell/user_projects/BCH_slate/custom_modules; python BCH_Calibration.py
# parameters: [1] non-hypoxic, [2] hypoxic, [3] post-hypoxic, [4] boundary cell uptake, [5] oxygen diffusion
def HypoxiaModel(parameters):
    xml_data = ConfigFileParser("/Users/morgana/Documents/GitHub/BCH_inverse/PhysiCell/user_projects/BCH_slate/config/PhysiCell_settings.xml")
    cell_data = []
    # Add cell definitions here as necessary
    cell_data.append(xml_data.read_cell_data("non-hypoxic"))
    cell_data.append(xml_data.read_cell_data("hypoxic"))
    cell_data.append(xml_data.read_cell_data("post-hypoxic"))
    cell_data.append(xml_data.read_cell_data("boundary"))
    for i in range(0, len(cell_data)):
        cell_data[i].secretion[0].uptake_rate = parameters[i]
        xml_data.write_secretion_params(cell_data[i].name, cell_data[i].secretion, update_file = True)
    microenvironment_data = []
    # Add microenvironment parameters here as necessary (^^^)
    microenvironment_data.append(xml_data.read_me_params()[0])
    for i in range(0, len(microenvironment_data)):
        microenvironment_data[i].diffusion_coefficient = parameters[len(cell_data) + i]
        xml_data.write_substance_params(microenvironment_data[i])
    os.system("cd ~/documents/github/BCH_inverse/physicell; make load PROJ=clean_rules; make load PROJ=BCH_slate; ./rules_modeling")
    new_oxygen_data = sio.loadmat("/Users/morgana/Documents/GitHub/BCH_inverse/PhysiCell/output_BCH_slate/final_microenvironment0.mat")['multiscale_microenvironment'][4]
    return new_oxygen_data
    
# custom cell data
base_oxygen_data = sio.loadmat("/Users/morgana/Documents/GitHub/BCH_inverse/PhysiCell/user_projects/BCH_slate/custom_modules/inverse/LF_MCMC_oxygen_prior_switchexport")['multiscale_microenvironment'][4]
    
# Gaussian Prior
mu = np.array([7.5, 15.0, 7.5, 0, 100000.0])
sigma = np.array([0.6, 0.8, 0.6, 0, 5000.0])
ABC_MCMC(HypoxiaModel, base_oxygen_data, mu, sigma, "/Users/morgana/Documents/GitHub/BCH_inverse/PhysiCell/user_projects/BCH_slate/custom_modules/LF_MCMC_output/uptake_diffusion_calibration.csv")
