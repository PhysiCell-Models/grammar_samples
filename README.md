# Sample models for the PhysiCell grammar / rules paper

## Latest Release

We have made our latest release available on Zenodo. You can find it at the following link:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14900735.svg)](https://doi.org/10.5281/zenodo.14900735)

## Manuscript Release

The version synced with the published manuscript is here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16285252.svg)](https://doi.org/10.5281/zenodo.16285252)

This release includes all the necessary files and documentation for _Human interpretable grammar encodes multicellular systems biology models to democratize virtual cell laboratories_ (with preprint named [_Digitize your Biology! Modeling multicellular systems through interpretable cell behavior_](https://doi.org/10.1101/2023.09.17.557982)). This release is dated July 21, 2025. Please refer to the [release notes](https://github.com/PhysiCell-Models/grammar_samples/releases/tag/v2.0.1) for more details.

## Getting started
### Windows users
Note: these instructions are written for Linux/MacOS systems.
For Windows, the commands may need to be adjusted.
In particular, the path separators need to be changed from `/` to `\` and the `./project` executable updated to `project.exe`.

### How to run the models
To  use any of these models:
1. clone or download this repository
2. from a terminal, navigate to the repository directory
3. run the model-specific commands to load the model and run the simulation

Generally, the commands are given by:
``` sh
make load PROJ=<project_name> && make && ./project
```
where `<project_name>` is the name of the project directory in `user_projects`.
Some of the models have multiple configuration files, requiring a second argument to specify the desired configuration file:
```sh
make load PROJ=<project_name> && make && ./project <path/to/config>
```
where `<path/to/config>` is the path to the desired configuration file, typically `./config/<filename.xml>`.

### Viewing the results in PhysiCell-Studio
Finally, [PhysiCell-Studio](https://github.com/PhysiCell-Tools/PhysiCell-Studio) can be used to interactively explore the model results.
First, follow the [instructions](https://github.com/PhysiCell-Tools/Studio-Guide/blob/main/README.md#installing-and-running) to download the studio.
After starting the simulation (the call to `./project`), you may immediately open a new terminal and run the following command:
```sh
python <path/to/studio/bin/studio.py> -c <path/to/config> -e ./project
``` 
where `<path/to/studio/bin/studio.py>` is the path to the studio executable (located in the studio directory in `bin/studio.py`).
The config and project arguments are optional, with studio defaulting to `config/PhysiCell_settings.xml` and `./project`, respectively.
Navigate to the `Plots` tab and use the GUI to explore the model results.

## Model 1: Modeling the progression of hypoxia in a metastatic tumor ([link](user_projects/hypoxia))

### Description: 
The model represents oxygen-dependent tumor growth, incorporating mechano-feedback mechanisms and hypoxia-induced cell migration.

### Load and run simulation:
``` bash
make load PROJ=hypoxia && make && ./project
```
### Cell Hypothesis Rules:
<details open>
    <summary>In tumor cells:</summary>
    <ul>
    <li>oxygen increases cycle entry from 0.00072 towards 0.0072 with a Hill response, with half-max 21.5 and Hill power 4.</li>
    <li>pressure decreases cycle entry from 0.00072 towards 0 with a Hill response, with half-max 0.25 and Hill power 3.</li>
    <li>oxygen decreases necrosis from 0.0028 towards 0 with a Hill response, with half-max 3.75 and Hill power 8.</li>
    <li>oxygen decreases transform to motile tumor from 0.001 towards 0 with a Hill response, with half-max 6.75 and Hill power 8.</li>
    </ul>
</details>
<details open>
    <summary>In motile tumor cells:</summary>
    <ul>
    <li>oxygen increases cycle entry from 0.00072 towards 0.0072 with a Hill response, with half-max 21.5 and Hill power 4.</li>
    <li>pressure decreases cycle entry from 0.00072 towards 0 with a Hill response, with half-max 0.25 and Hill power 3.</li>
    <li>oxygen decreases necrosis from 0.0028 towards 0 with a Hill response, with half-max 3.75 and Hill power 8.</li>
    <li>oxygen increases transform to tumor from 0 towards 0.005 with a Hill response, with half-max 6.75 and Hill power 8.</li>
    </ul>
</details>

### Sensitivity Analysis
A local sensitivity analysis was performed on 24 parameters associated with the cell hypothesis rules to evaluate their impact on the model's behavior. [View the analysis script](Sensitivity_Analysis/SA_script_ex1.ipynb).

## Model 2: Fibroblast-mediated invasion of neoplastic cells ([link](user_projects/epi_caf_invasion))
### Description:
The model represents the interactions between neoplastic cancer cells and fibroblasts in the context of cancer invasion.
Two scenarios are considered: a coculture experiment and one based on patient samples.

For the coculture experiment, use the `PhysiCell_settings_coculture.xml` configuration file.
In that file, set the `<cell_positions>` element to set the desired initial cancer:fibroblast ratio.
The options are located in `config/ics/coculture`.

For the PDAC condition, two samples are provided in `config/ics/`: `PDAC01` and `PDAC02`.
Each of these has an associated `cells.csv` and `substrates.csv` file setting initial conditions for the cells and substrates, respectively.
To run the model with one of these samples, set the `<cell_positions>` element in `PhysiCell_settings_PDAC.xml` to the desired `cells.csv` file and also set the `microenvironment_setup//options//initial_condition` to the desired `substrates.csv` file.
Note: these XML locations can both be found in the file by searching for `PDAC`.

Finally, these simulations take longer to run than many others in this repository.
Each of the two configuration files saves output to a specified folder.
Thus, if you run the coculture, for example, on two different ratios, the second run will overwrite the first.
Plan accordingly!

### Load and run simulation:
```sh
make load PROJ=epi_caf_invasion && make && ./project ./config/PhysiCell_settings_coculture.xml
```
```sh
make load PROJ=epi_caf_invasion && make && ./project ./config/PhysiCell_settings_PDAC.xml
```

### Cell Hypothesis Rules:
<details open>
    <summary>In epithelial_normal cells:</summary>
    <ul>
    <li>pressure decreases cycle entry from 0 towards 0 with a Hill response, with half-max 0.5 and Hill power 4.</li>
    <li>ecm increases transform to mesenchymal_normal from 0 towards 0.01 with a Hill response, with half-max 0.01 and Hill power 4.</li>
    </ul>
</details>
<details open>
    <summary>In mesenchymal_normal cells:</summary>
    <ul>
    <li>ecm decreases migration speed from 0.249909 towards 0.249678 with a Hill response, with half-max 0.0595916 and Hill power 1.22413.</li>
    <li>ecm increases migration speed from 0.249909 towards 3.78565 with a Hill response, with half-max 3.11941 and Hill power 9.99815.</li>
    <li>inflammatory_signal decreases transform to epithelial_normal from 0.01 towards 0 with a Hill response, with half-max 0.2 and Hill power 4.</li>
    </ul>
</details>
<details open>
    <summary>In fibroblast cells:</summary>
    <ul>
    <li>ecm decreases migration speed from 8.03971e-06 towards 6.90516e-06 with a Hill response, with half-max 2.07833 and Hill power 1.1551.</li>
    <li>ecm increases migration speed from 8.03971e-06 towards 3.47055 with a Hill response, with half-max 9.99968 and Hill power 1.15302.</li>
    </ul>
</details>
<details open>
    <summary>In epithelial_tumor cells:</summary>
    <ul>
    <li>pressure decreases cycle entry from 0.001 towards 0 with a Hill response, with half-max 1 and Hill power 4.</li>
    <li>ecm increases transform to mesenchymal_tumor from 0 towards 0.01 with a Hill response, with half-max 0.01 and Hill power 4.</li>
    </ul>
</details>
<details open>
    <summary>In mesenchymal_tumor cells:</summary>
    <ul>
    <li>ecm decreases migration speed from 0.249909 towards 0.249678 with a Hill response, with half-max 0.0595916 and Hill power 1.22413.</li>
    <li>ecm increases migration speed from 0.249909 towards 3.78565 with a Hill response, with half-max 3.11941 and Hill power 9.99815.</li>
    <li>inflammatory_signal decreases transform to epithelial_tumor from 0.01 towards 0 with a Hill response, with half-max 0.2 and Hill power 4.</li>
    </ul>
</details>

## Model 3: Base tumor-immune example ([link](user_projects/tumor_immune_base))
### Description: 
The model represents the basic interactions between the tumor and the immune system, incorporating immune responses to tumor groth, such as macrophage polarization and CD8 T cell cytotoxicity.

### Load and run simulation:
``` bash
make load PROJ=tumor_immune_base && make && ./project
```
### Cell Hypothesis Rules:
<details open>
    <summary>In tumor cells:</summary>
    <ul>
    <li>oxygen increases cycle entry from 0.00072 towards 0.0072 with a Hill response, with half-max 21.5 and Hill power 4.</li>
    <li>pressure decreases cycle entry from 0.00072 towards 0 with a Hill response, with half-max 0.25 and Hill power 3.</li>
    <li>oxygen decreases necrosis from 0.0028 towards 0 with a Hill response, with half-max 3.75 and Hill power 8.</li>
    <li>damage increases apoptosis from 5.31667e-05 towards 0.023 with a Hill response, with half-max 45 and Hill power 16.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
    </ul>
</details>
<details open>
    <summary>In macrophage cells:</summary>
    <ul>
    <li>oxygen increases pro-inflammatory factor secretion from 0 towards 10 with a Hill response, with half-max 12 and Hill power 16.</li>
    <li>oxygen decreases anti-inflammatory factor secretion from 10 towards 0 with a Hill response, with half-max 12 and Hill power 16.</li>
    </ul>
</details>
<details open>
    <summary>In CD8 T cell cells:</summary>
    <ul>
    <li>anti-inflammatory factor decreases attack tumor from 0.1 towards 0 with a Hill response, with half-max 0.5 and Hill power 8.</li>
    <li>pro-inflammatory factor increases attack tumor from 0.1 towards 1 with a Hill response, with half-max 0.5 and Hill power 8.</li>
    <li>anti-inflammatory factor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.5 and Hill power 8.</li>
    <li>contact with tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.5 and Hill power 8.</li>
    </ul>
</details>

### Sensitivity Analysis
A local sensitivity analysis was conducted on 26 parameters associated with the cell hypothesis rules to evaluate their impact on tumor-immune interactions. [View the analysis script](Sensitivity_Analysis/SA_script_ex3.ipynb).

## Model 4: Extended tumor-immune example ([link](user_projects/tumor_immune_extended))
### Description:
The model extends the base tumor-immune model by including distinct macrophage cell states and naive CD8 T cells that can becom activated.

### Load and run simulation:
```sh
make load PROJ=tumor_immune_extended && make && ./project
```

### Cell Hypothesis Rules:
<details open>
    <summary>In tumor cells:</summary>
    <ul>
    <li>pressure decreases cycle entry from 0.00072 towards 0 with a Hill response, with half-max 1 and Hill power 4.</li>
    <li>oxygen increases cycle entry from 0.00072 towards 0.00072 with a Hill response, with half-max 21.5 and Hill power 4.</li>
    <li>oxygen decreases necrosis from 0.0028 towards 0 with a Hill response, with half-max 3.75 and Hill power 8.</li>
    <li>damage increases apoptosis from 5.31667e-05 towards 0.072 with a Hill response, with half-max 180 and Hill power 2.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    </ul>
</details>
<details open>
    <summary>In M0 macrophage cells:</summary>
    <ul>
    <li>contact with dead cell decreases migration speed from 1 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 4.</li>
    <li>contact with dead cell increases transform to M1 macrophage from 0 towards 0.05 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    </ul>
</details>
<details open>
    <summary>In M1 macrophage cells:</summary>
    <ul>
    <li>contact with dead cell decreases migration speed from 1 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 4.</li>
    <li>oxygen decreases transform to M2 macrophage from 0.01 towards 0 with a Hill response, with half-max 5 and Hill power 4.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    </ul>
</details>
<details open>
    <summary>In M2 macrophage cells:</summary>
    <ul>
    <li>contact with dead cell decreases migration speed from 1 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 4.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    </ul>
</details>
<details open>
    <summary>In naive T cell cells:</summary>
    <ul>
    <li>IL-10 decreases transform to CD8 T cell from 0.001 towards 0 with a Hill response, with half-max 0.25 and Hill power 2.</li>
    <li>IFN-gamma increases transform to CD8 T cell from 0.001 towards 0.01 with a Hill response, with half-max 0.25 and Hill power 2.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    </ul>
</details>
<details open>
    <summary>In CD8 T cell cells:</summary>
    <ul>
    <li>IL-10 decreases attack tumor from 0.1 towards 0 with a Hill response, with half-max 0.25 and Hill power 2.</li>
    <li>IL-10 decreases migration speed from 1 towards 0.1 with a Hill response, with half-max 0.25 and Hill power 2.</li>
    <li>contact with tumor decreases migration speed from 1 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>IFN-gamma increases cycle entry from 7.2e-05 towards 0.00041 with a Hill response, with half-max 0.25 and Hill power 2.</li>
    <li>IL-10 increases transform to exhausted T cell from 0 towards 0.005 with a Hill response, with half-max 0.25 and Hill power 4.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    </ul>
</details>
<details open>
    <summary>In exhausted T cell cells:</summary>
    <ul>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10.</li>
    </ul>
</details>

## Model 5: Tumor-associated macrophage and EGFR signaling example ([link](user_projects/tam_egf))
### Description:
The model represents the interactions between tumor cells and tumor-associated macrophages (TAMs), incorporating the effect of EGFR signaling in tumor growth and motility.

Three hypotheses are tested: EGFR signaling increasing (1) tumor growth, (2) tumor motility, and (3) both.
These correspond to the three config files and the three rules files with filenames ending with `_grow`, `_go`, and `_go_and_grow`, respectively.

Furthermore, the cell initial conditions is provided in `config/cells_with_M2.csv` (used to create panels 5C-E in the paper).
An alternative initial condition with M2-like macrophages replaced with naive M0 macrophages is given in `config/cells_with_M0.csv` (used to create panel 3H with the `systemic_therapy` parameter set accordingly).
A third alternative initial condition with CD4+ T cells and M0 macrophages is provided in `config/cells_with_M0_and_CD4.csv` (used to create panel 3I).
For the desired hypothesis, change the config file (e.g., `PhysiCell_settings_grow.xml`) to use the desired cells file in the `<cell_positions>` element.

Finally, systemic therapy is controlled via the `<user_parameter>` `systemic_therapy`in each of the config files.
Set this to one of `IFN-gamma`, `IL-4`, or `none` to simulate the effect of the respective therapy.

### Load and run simulation:
```sh
make load PROJ=tam_egf && make && ./project ./config/PhysiCell_settings_grow.xml
```
```sh
make load PROJ=tam_egf && make && ./project ./config/PhysiCell_settings_go.xml
```
```sh
make load PROJ=tam_egf && make && ./project ./config/PhysiCell_settings_go_and_grow.xml
```

### Cell Hypothesis Rules:
<details open>
    <summary>Rules in common for all hypotheses:</summary>
    <details open>
        <summary>In tumor cells:</summary>
        <ul>
        <li>pressure decreases cycle entry from 0.0001 towards 0 with a Hill response, with half-max 1 and Hill power 4.</li>
        <li>oxygen decreases necrosis from 0.0028 towards 0 with a Hill response, with half-max 3.75 and Hill power 8.</li>
        <li>damage increases apoptosis from 7.2e-05 towards 0.072 with a Hill response, with half-max 180 and Hill power 2.</li>
        <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
        <li>IFN-gamma decreases migration speed from 0.17 towards 0 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        </ul>
    </details>
    <details open>
    <summary>In M0 macrophage cells:</summary>
        <ul>
        <li>contact with dead cell increases transform to M1 macrophage from 0 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 10.</li>
        <li>debris increases transform to M1 macrophage from 0 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 10.</li>
        <li>contact with dead cell decreases migration speed from 0.5 towards 0 with a Hill response, with half-max 0.01 and Hill power 4.</li>
        <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
        </ul>
    </details>
    <details open>
        <summary>In M1 macrophage cells:</summary>
        <ul>
        <li>contact with dead cell decreases migration speed from 0.5 towards 0 with a Hill response, with half-max 0.01 and Hill power 4.</li>
        <li>oxygen decreases transform to M2 macrophage from 0.01 towards 0 with a Hill response, with half-max 5 and Hill power 4.</li>
        <li>IL-4 increases transform to M2 macrophage from 0.01 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 4.</li>
        <li>IFN-gamma increases cycle entry from 7.2e-05 towards 0.00036 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>IFN-gamma increases phagocytose apoptotic cell from 0.001 towards 0.05 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>IFN-gamma increases phagocytose necrotic cell from 0.001 towards 0.05 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
        </ul>
    </details>
    <details open>
        <summary>In M2 macrophage cells:</summary>
        <ul>
        <li>IFN-gamma increases transform to M1 macrophage from 0 towards 0.001 with a Hill response, with half-max 1 and Hill power 4.</li>
        <li>contact with dead cell decreases migration speed from 0.5 towards 0 with a Hill response, with half-max 0.01 and Hill power 4.</li>
        <li>IFN-gamma decreases cycle entry from 7.2e-05 towards 0 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>IL-4 increases cycle entry from 7.2e-05 towards 0.001 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>IFN-gamma increases phagocytose apoptotic cell from 0.001 towards 0.05 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>IFN-gamma increases phagocytose necrotic cell from 0.001 towards 0.05 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
        </ul>
    </details>
    <details open>
        <summary>In Th2 CD4 T cell cells:</summary>
        <ul>
        <li>IL-4 increases cycle entry from 0 towards 0.00036 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>IFN-gamma decreases cycle entry from 0 towards 0 with a Hill response, with half-max 0.25 and Hill power 2.</li>
        <li>contact with M0 macrophage decreases migration speed from 0.5 towards 0.25 with a Hill response, with half-max 0.1 and Hill power 2.</li>
        <li>contact with M1 macrophage decreases migration speed from 0.5 towards 0.25 with a Hill response, with half-max 0.1 and Hill power 2.</li>
        <li>contact with M2 macrophage decreases migration speed from 0.5 towards 0.25 with a Hill response, with half-max 0.1 and Hill power 2.</li>
        <li>contact with tumor decreases migration speed from 0.5 towards 0.25 with a Hill response, with half-max 0.1 and Hill power 2.</li>
        <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
        </ul>
    </details>
</details>
<details open>
    <summary>Grow hypothesis:</summary>
    <details open>
        <summary>In tumor cells:</summary>
        <ul>
        <li>EGF increases cycle entry from 0.0001 towards 0.001 with a Hill response, with half-max 2 and Hill power 10.
        </li>
        </ul>
    </details>
 </details>
<details open>
    <summary>Go hypothesis:</summary>
    <details open>
        <summary>In tumor cells:</summary>
        <ul>
        <li>EGF increases migration speed from 0.17 towards 0.5 with a Hill response, with half-max 2 and Hill power 10.</li>
        <li>EGF increases migration persistence time from 5 towards 10 with a Hill response, with half-max 2 and Hill power 10.</li>
        <li>EGF decreases cell-cell adhesion from 0.4 towards 0.2 with a Hill response, with half-max 2 and Hill power 10.</li>
        </ul>
    </details>
 </details>

## Model 6: PDAC immunotherapy example ([link](user_projects/pdac_therapy))
### Description:
The model represents the interactions between the tumor and the immune system in pancreatic ductal adenocarcinoma (PDAC), incorporating the effect of immunotherapies.

The model includes initial conditions representing 16 different samples, each treated with one of eight therapies representing any combination of GVAX (`GVAX`), nivolumab (`ICI`), and urelumab (`URU`).
These 128 files are located in `config/ic_cells/` and follow the naming convention `PDAC_TISSUE_<number>_<therapy>.csv` where `<number>` is a sample index and `<therapy>` is one of `GVAX`, `GVAX_ICI`, etc.
To run the model with different samples and therapies, change the `config/PhysiCell_settings.xml` file inside the `<cell_positions>` element to the desired file.

### Load and run simulation:
```sh
make load PROJ=pdac_therapy && make && ./project
```

### Cell Hypothesis Rules:
<details open>
    <summary>In PD-L1lo_tumor cells:</summary>
    <ul>
    <li>pressure decreases cycle entry from 0.00072 towards 0 with a Hill response, with half-max 0.25 and Hill power 3.</li>
    <li>damage increases apoptosis from 7.2e-05 towards 0.023 with a Hill response, with half-max 45 and Hill power 16.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
    </ul>
</details>
<details open>
    <summary>In PD-L1hi_tumor cells:</summary>
    <ul>
    <li>pressure decreases cycle entry from 0.00072 towards 0 with a Hill response, with half-max 0.25 and Hill power 3.</li>
    <li>damage increases apoptosis from 7.2e-05 towards 0.023 with a Hill response, with half-max 45 and Hill power 16.</li>
    <li>dead increases debris secretion from 0 towards 0.017 with a Hill response, with half-max 0.1 and Hill power 10. Rule applies to dead cells.</li>
    </ul>
</details>
<details open>
    <summary>In macrophage cells:</summary>
    <ul>
    <li>oxygen decreases anti-inflammatory factor secretion from 1 towards 0 with a Hill response, with half-max 5 and Hill power 4.</li>
    <li>oxygen increases pro-inflammatory factor secretion from 0 towards 0.1 with a Hill response, with half-max 5 and Hill power 4.</li>
    </ul>
</details>
<details open>
    <summary>In PD-1hi_CD137lo_CD8_Tcell cells:</summary>
    <ul>
    <li>contact with PD-L1hi_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>contact with PD-L1lo_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases migration speed from 1 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1hi_tumor from 0 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1lo_tumor from 3.33e-06 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    </ul>
</details>
<details open>
    <summary>In PD-1lo_CD137lo_CD8_Tcell cells:</summary>
    <ul>
    <li>contact with PD-L1hi_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>contact with PD-L1lo_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases migration speed from 1 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1hi_tumor from 3.33e-05 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>pro-inflammatory factor increases attack PD-L1hi_tumor from 3.33e-05 towards 0.01 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1lo_tumor from 3.33e-05 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>pro-inflammatory factor increases attack PD-L1lo_tumor from 3.33e-05 towards 0.01 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    </ul>
</details>
<details open>
    <summary>In PD-1hi_CD137hi_CD8_Tcell cells:</summary>
    <ul>
    <li>contact with PD-L1hi_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>contact with PD-L1lo_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases migration speed from 1 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1hi_tumor from 0 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1lo_tumor from 3.33e-05 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    </ul>
</details>
<details open>
    <summary>In PD-1lo_CD137hi_CD8_Tcell cells:</summary>
    <ul>
    <li>contact with PD-L1hi_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>contact with PD-L1lo_tumor decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases migration speed from 1 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1hi_tumor from 0.000333 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>pro-inflammatory factor increases attack PD-L1hi_tumor from 0.000333 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    <li>anti-inflammatory factor decreases attack PD-L1lo_tumor from 0.000333 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>pro-inflammatory factor increases attack PD-L1lo_tumor from 0.000333 towards 0.1 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    </ul>
</details>
<details open>
    <summary>In PD-1hi_CD4_Tcell cells:</summary>
    <ul>
    <li>anti-inflammatory factor decreases migration speed from 1 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>contact with macrophage decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    </ul>
</details>
<details open>
    <summary>In PD-1lo_CD4_Tcell cells:</summary>
    <ul>
    <li>anti-inflammatory factor decreases migration speed from 1 towards 0 with a Hill response, with half-max 2.5 and Hill power 2.</li>
    <li>contact with macrophage decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.1 and Hill power 2.</li>
    </ul>
</details>

## Model 7: Neuro development example ([link](user_projects/neuro_dev))
### Description:
The model represents the formation of the laminar structure of the cortex, incorporating asymmetric division of a progenitor layer and migration of neurons to form the cortical layers.

Two configuration files are included to represent two regions of the cortex: the somatosensory cortex (SOM) and the auditory cortex (AUD).

### Load and run simulation:
```sh
make load PROJ=neuro_dev && make && ./project ./config/PhysiCell_settings_SOM.xml
```
```sh
make load PROJ=neuro_dev && make && ./project ./config/PhysiCell_settings_AUD.xml
```

### Cell Hypothesis Rules:
<details open>
    <summary>SOM:</summary>
    <details open>
        <summary>In rgc cells:</summary>
        <ul>
        <li>contact with apical decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.5 and Hill power 8.</li>
        <li>time decreases cycle entry from 0.0013 towards 0.00085 with a Hill response, with half-max 9206 and Hill power 10.</li>
        <li>time increases asymmetric division to layer_6 from 0 towards 1 with a Hill response, with half-max 1440 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_6 from 0 towards 0 with a Hill response, with half-max 5007 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_5 from 0 towards 1 with a Hill response, with half-max 5007 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_5 from 0 towards 0 with a Hill response, with half-max 6554 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_4 from 0 towards 1 with a Hill response, with half-max 6554 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_4 from 0 towards 0 with a Hill response, with half-max 9010 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_3 from 0 towards 1 with a Hill response, with half-max 9010 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_3 from 0 towards 0 with a Hill response, with half-max 10985 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_2 from 0 towards 1 with a Hill response, with half-max 10985 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_2 from 0 towards 0 with a Hill response, with half-max 12960 and Hill power 256.</li>
        <li>time increases apoptosis from 0 towards 1 with a Hill response, with half-max 12960 and Hill power 256.</li>
        </ul>
    </details>
</details>
<details open>
    <summary>AUD:</summary>
    <details open>
        <summary>In rgc cells:</summary>
        <ul>
        <li>contact with apical decreases migration speed from 1 towards 0 with a Hill response, with half-max 0.5 and Hill power 8.</li>
        <li>time decreases cycle entry from 0.0013 towards 0.00085 with a Hill response, with half-max 4808 and Hill power 10.</li>
        <li>time increases asymmetric division to layer_6 from 0 towards 1 with a Hill response, with half-max 1440 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_6 from 0 towards 0 with a Hill response, with half-max 3309 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_5 from 0 towards 1 with a Hill response, with half-max 3309 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_5 from 0 towards 0 with a Hill response, with half-max 7068 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_4 from 0 towards 1 with a Hill response, with half-max 7068 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_4 from 0 towards 0 with a Hill response, with half-max 8638 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_3 from 0 towards 1 with a Hill response, with half-max 8638 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_3 from 0 towards 0 with a Hill response, with half-max 10799 and Hill power 256.</li>
        <li>time increases asymmetric division to layer_2 from 0 towards 1 with a Hill response, with half-max 10799 and Hill power 256.</li>
        <li>time decreases asymmetric division to layer_2 from 0 towards 0 with a Hill response, with half-max 12960 and Hill power 256.</li>
        <li>time increases apoptosis from 0 towards 1 with a Hill response, with half-max 12960 and Hill power 256.</li>
        </ul>
    </details>
</details>
