# LogicModel
This is a small introduction on how to use the logic model.



(1) MAIN FUNCTION

main.m is the main function of the logic model. Here the glucose and the nitrogen levels, the active crosstalks and the knockouts are set. Then the model runs and saves all the plots and data in the folder specified by the foldername. 



(2) THE MODEL FUNCTION

runLogicModel.m intializes the model and loops through all boolean operations, specified in runUntilSteadyState.m and crosstalk.m, until the logial steady state is reached. Knockouts are proteins whose presence is set to 0 after the initialization and stays 0 throughout the whole procedure. 
If the input is more than one step (e.g. glucoseLevels = [1 0 1 0]) the model also successively loops over this sequence, with the output of each step being the input to the next step. Thus, the initalization is only done in the very beginning to provide a starting point.

Note that the functions runUntilSteadyState.m and crosstalks.m include all boolean operations that are contained in the model.

Involved functions: initializeModel.m, knockouts.m, runUntilSteadyState.m, crosstalks.m



(3) SAVING AND PLOTTING

In the last step, saveAndPlot.m saves the data and schematic pictures are created and saved in the Data folder. The functions are mainly for plotting and have nothing to do with the model itself.

Involved functions: arrowCoordinates,.m, ellipseMaker.m, getProteinPosition.m, getNameTag.m, includeProtein.m, makeModelPicture.m, plotEmptyCell.m, rectangleMaker.m, saveAndPlot.m



(4) EXAMPLE: CROSSTALK TESTING

In the folder Crosstalk Simulation the files main.m and mainTOR.m are provided which reproduce comparisons between the model with and without crosstalk and gapfilling, as well as the effect of all possible combinations of crosstalks on the gene expression in the model. Output files are created in the Data folder. To save computation time the pictures of the individual pathways are not created. 



(5) HOW TO ADAPT THE CODE

In case modules or single components want to be added, they have to be initalized in the initalizeModel.m function. All boolean operations can be added in the runUntilSteadyState.m and the crosstalk.m function. If a new pathway should also be plotted they have to be added in the saveAndPlot.m function, both writing the txt file and in the names that are listed for creating the figures. If the table for the added pathway is in the same format as the others, the plots are generated automatically.










