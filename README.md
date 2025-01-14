# Sample models for the PhysiCell grammar / rules paper

## Model 1: Modeling the progression of hypoxia in a metastatic tumor ([link](user_projects/example1_hypoxia))

### Description: 
The model represents oxygen-dependent tumor growth, incorporating mechano-feedback mechanisms and hypoxia-induced cell migration.

### Load and run simulation:
``` bash
make load PROJ=example1_hypoxia && make && ./project
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
