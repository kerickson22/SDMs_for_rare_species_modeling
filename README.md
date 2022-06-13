# Species distribution modeling for very rare species

Code associated with K. D. Erickson & A. B. Smith. **Modeling the rarest of the rare : A comparison between joint species distribution models, ensembles of small models, and single-species models at extremely low sample sizes**. 


1. Simulate Species 

2. Run models: 

| File              | Outputs |
| ----------------- | ----------- |
| `HMSC_Simple.R`   |  720: "../data/models/HMSC_Simple/SpeciesType/NumPresences/model_ReplicateNumber.RData" |    
| `HMSC_Joint_Simple.R` |  720: "../data/models/HMSC_Joint/SpeciesType/NumPresences/model_ReplicateNumber.RData"  |
| `ESM.R`               |  720: "../data/models/ESM_bivariate/SpeciesType/NumPresences/model_ReplicateNumber.RData" |
|                          720: "../data/models/ESM/SpeciesType/NumPresences/model_ReplicateNumber.RData            |
|                          'results.RData'                                                                          |


3a. Process Models 

| File              | Outputs |
| ----------------- | ----------- |
| `Process_Models_HMSC_Single_Simple.R`   | `results.RData`    |    
| `Process_Models_HMSC_Joint_Simple.R`    | `results.RData`    |

 *(ESM model processing takes place within the ESM script*


3b. PC2 Calibration
| File              | Outputs |
| ----------------- | ----------- |
| `PC2 Calibration.R`   | `results2.RData`    |    



4. Figures 




