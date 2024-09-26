# Sample models for the PhysiCell grammar / rules paper

## Model 1: Modeling the progression of hypoxia in a metastatic tumor ([link](user_projects/example1_hypoxia))

### Description: 
The model represents oxygen-dependent tumor growth, incorporating mechano-feedback mechanisms and hypoxia-induced cell migration.

### Cell Hypothesis Rules:

| Cell        | Signal           | Direction  | Behavior                  | Saturation value  | Half max value | Hill power | Apply to dead |
|-------------|------------------|------------|---------------------------|-------------------|----------------|------------|---------------|
| tumor       | oxygen           | increases  | cycle entry               | 7.00E-04          | 21.5           | 4          | 0             |
| tumor       | pressure         | decreases  | cycle entry               | 0                 | 0.25           | 3          | 0             |
| tumor       | oxygen           | decreases  | necrosis                  | 0                 | 3.75           | 8          | 0             |
| tumor       | oxygen           | decreases  | transform to motile tumor | 0.0               | 6.75           | 8          | 0             |
| motile tumor| oxygen           | increases  | cycle entry               | 7.00E-04          | 21.5           | 4          | 0             |
| motile tumor| pressure         | decreases  | cycle entry               | 0                 | 0.25           | 3          | 0             |
| motile tumor| oxygen           | decreases  | necrosis                  | 0                 | 3.75           | 8          | 0             |
| motile tumor| oxygen           | increases  | transform to tumor        | 0.005             | 6.75           | 8          | 0             |



## Model 2: 