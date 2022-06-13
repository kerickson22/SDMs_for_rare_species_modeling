# Species distribution modeling for very rare species

Code associated with K. D. Erickson & A. B. Smith. **Modeling the rarest of the rare : A comparison between joint species distribution models, ensembles of small models, and single-species models at extremely low sample sizes**. 


1. Simulate Species: Creates 4 species types x 6 Numbers of Presences x 30 replicates = 720 simulated species

2. Run models: Runs the models for each of the simulated species. 
*`HMSC_Simple.R`:    
* `HMSC_Joint_Simple.R` 
* `ESM.R`: Runs both simple and complex ESMs in the same script                                                                          |

3a. Process Models: Loops through each of the models that were run in step 2 and gathers summary statistics. 
* `Process_Models_HMSC_Single_Simple.R`   
* `Process_Models_HMSC_Joint_Simple.R`   

 *(Note: ESM model processing takes place within the ESM script*

3b. PC2 Calibration: Calculates PC2 calibration for each model type

 *`PC2 Calibration.R`   

4. Figures: Creates the figures that appear in the manuscript




