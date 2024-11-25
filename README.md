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


## Model 3: Base tumor-immune example ([link](user_projects/example3_immune))
### Description: 
The model represents the basic interactions between the tumor and the immune system, incorporating immune responses to tumor groth, such as macrophage polarization and CD8 T cell cytotoxicity.

### Load and run simulation:
``` bash
make load PROJ=example3_immune && make && ./project
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