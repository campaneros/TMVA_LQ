# TMVA_LQ

This repository allow to create reduced tree with all the weight and run a BDT on them.
## Setup the environment
Download the following CMSSW arch
```
scram p -n CMSSW_12_1_0_TMVA CMSSW CMSSW_12_1_0_pre4_ROOT624
cd CMSSW_12_1_0_pre4_ROOT624/src
cmsenv
```

Download the code
```
git clone git@github.com:campaneros/TMVA_LQ.git
cd TMVA_LQ
git checkout Snakemake
```

## Run With Snakemake
To download snakmake
```
pip3 install snakemake --user
```
[Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/index.html)

Edit the file config_new.yaml
```
vim config_new.yaml 
```
and modify all the directory according to your need

To view the rules that will be excecuted run
```
nohup -s Snakefile_new -n -c30
```

To Create a reduced dataset, Run a BDT and Apply the BDT Results on MC and Signal Sample
```
cmsenv
nohup -s Snakefile_new -c30
```


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
