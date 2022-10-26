# Species distribution modeling for very rare species

Code associated with K. D. Erickson & A. B. Smith. **Modeling the rarest of the rare : A comparison between joint species distribution models, ensembles of small models, and single-species models at extremely low sample sizes**. 


1. Assemble data and simulate species 
* `00a_Run_Once.R`: Loads data from @Norberg_etal_2019, subsets to SE states, calculates statistics associated with species niches, determine how many archetypes to use for the Species Archetype Model, and pick "marginal" niche location as well as specialist niche breadth. 

* `00b_Constants.R`: Script that is sourced in all subsequent files. Loads packages, creates directories for storing modeling results, loads constants. 

* `01_Simulate_Species.R`: Creates 4 species types {Central Generalist, Central Specialist, Marginal Generalist, Marginal Specialist} x 6 Numbers of Presences {2, 4, 8, 32, 64} x 100 replicates = 2400 Simulated Species 


2. Run models: Runs the models for each of the simulated species.
 
* `02_Model_01_glm.R`:   Single-species, deterministic glm {for n > 8}
* `02_Model_02_ESM_linear.R`: Single-species, deterministic Ensembles of Small Models, using only linear terms {for n >2}
* `02_Model_03_ESM_polynomial.R`: Single-species, deterministic Ensembles of Small Models, using both linear and quadratic terms {for n > 8}
* `02_Model_04_Hmsc:   JSDM, Bayesian Hierarchical Models of Species Communities {for n > 2}
* `02_Model_05_SAM: JSDM, deterministic Species Archetype Model 

3a. Process Models: Loops through each of the models that were run in step 2 and gathers summary statistics. 

* `03_Collate_Results.R` 

4. Display Items: Creates the figures and tables that appear in the manuscript

* `04_Figures_and_Tables.R`




