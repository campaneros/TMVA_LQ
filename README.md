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
snakemake -s Snakefile_new -n -c30
```

To Create a reduced dataset, Run a BDT and Apply the BDT Results on MC and Signal Sample
```
cmsenv
nohup snakemake -s Snakefile_new -c30&
```
The nohup option will detach the process from the terminal so you can close the terminal without killing it. The output of the process wil be in the file nohup.out

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

once you have the output from the BDT you should create new root file with a branch with the output of the BDT, to do that run 
```
python3 App_Classificator_1muons.py -c makePlots_datasetlists_reduced.txt -f TMVA_BDT_1muon -o BDT
```

where the file to be used should be the weight file created from the BDT at the previous step.

## Do Significance Scan on the variable
To run a Significance Scan on the desired variable run (Example: Met_Significance)
```
for M in 1000 2000 3000;do echo ${M}; for p in 0p1 1p0; do python3 Scan_Met_2muons.py -c makePlots_datasetlists_2muons_BDT_noPFisp_real.txt -o Signal_BDT_2Muons_Met_dE0p5 --step 0.5 --mass ${M} --coupling ${p}; done;done
```
where the -o option is the output directory of the plots produced by the scripts and the -c is the file containing the list of all root files. The --step option is the step of the variable scan. **For the two variable scan the second variable step has to be changed inside the code**
To unite all the plots on one image run
```
montage -tile 3x2x -geometry -2-2 -resize -0.5-0.5 Signal_BDT_2Muons_Met_dE0p5/LQ*.png Signal_BDT_2Muons_Met_dE0p5/monatge.jpg
```

### Variable Significance Scan for 1 Muon
1. $\Delta\phi_{\mu jet}$                                     [Scan_dPhi.py]
2. Met Significance                                           [Scan_Met.py]
3. $P_{T}^{\mu}/P_{T}^{jet}$ && $P_{T}^{jet}/M_{\mu jet}$     [Scan_MuonPT.py]

### Variable Significance Scan for 2 Muons category
1. $P_{T}^{\mu}/P_{T}^{jet}$ && $P_{T}^{\mu}/M_{\mu jet}$     [Scan_PtOvM.py]
2. $P_{T}^{\mu ,2}/M_{\mu jet}$ && $\Delta\phi_{\mu jet}$     [Scan_MuonPT_dPhi_2muons.py]
3. $|\delta \eta_{\mu\mu}|$                                   [Scan_Eta1_m_Eta2.py]
4. Met Significance                                           [Scan_Met_2muons.py]




