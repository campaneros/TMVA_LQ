# TMVA_LQ

This repository allow to create reduced tree with all the weight and run a BDT on them.

## Create reduced tree
To create a reduced tree run 
```
pytho3 Classificator_create_dataset_reduced.py -c makePlots_datasetlists.txt -o ntuples_training
```


## Run BDT on reduced Tree
To run the BDT from the reduced trees produced before, make a input file list
```
python3 Classificator_1muon.py -c makePlots_datasetlists_reduced.txt -f TMVA_BDT_1muon
```

where the option -f must be followed by the name of the desired root output file for the BDT and the weight file. The weight file should be in the directory dataset/weights

## Apply the BDT results from xml file
**To be optimised still**\
once you have the output from the BDT you should create new root file with a branch with the output of the BDT, to do that run 
```
python3 App_Classificator_1muons.py -c makePlots_datasetlists_reduced.txt -f TMVA_BDT_1muon -o BDT
```

where the file to be used should be the weight file created from the BDT at the previous step.
