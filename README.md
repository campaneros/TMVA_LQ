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
python3 Classificator_1muon.py -c makePlots_datasetlists -f TMVA_BDT_1muon
```

where the option -f must be followed by the name of the desired root output file for the BDT and the weight file
