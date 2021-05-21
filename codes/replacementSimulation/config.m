%% config
clear; home;
load('../BehavNeuralTbl.mat')

%%
res_dmpfcinrdm_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_RDM, 'dmpfc', 1e4);
res_dmpfcinoc_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_OC, 'dmpfc', 1e4);
res_dmpfcinsc_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_SC, 'dmpfc', 1e4);

%%
res_daccinrdm_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_RDM, 'dacc', 1e4);
res_daccinoc_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_OC, 'dacc', 1e4);
res_daccinsc_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_SC, 'dacc', 1e4);

%%
res_fpcinrdm_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_RDM, 'fpc', 1e4);
res_fpcinoc_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_OC, 'fpc', 1e4);
res_fpcinsc_simulation_replacement_10k = DidDoSimulation(BehavNeuraltbl_SC, 'fpc', 1e4);