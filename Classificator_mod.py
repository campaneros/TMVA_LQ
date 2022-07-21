#!/usr/bin/env python                                                                                                                                             

from ROOT import TMVA, TFile, TTree, TCut, gROOT
from subprocess import call
from os.path import isfile
 
import torch
from torch import nn


import os
import sys
import optparse
import datetime
import subprocess
from glob import glob
from collections import defaultdict
from collections import OrderedDict
import tdrstyle
import CMS_lumi

from ROOT import *

usage = "usage: To be run from RootTreeAnalyzer:  python analysis/makePlots.py -c analysis/makePlots.txt -o output/TEST_cms_lq_prod2/plots -f TestOutput --fracMC 1.0 --fracData 1.0"

parser = optparse.OptionParser(usage)

parser.add_option("-c", "--config", dest="configFile",
                  help="config file")

#parser.add_option("-o", "--output", dest="outputdir",
#                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

#parser.add_option("-f", "--outFileName", dest="outFileName",
#                  help="the output filename")

#parser.add_option("--fracMC", dest="fracMC", default=1,
#                  help="fraction of data events to draw")

#parser.add_option("--fracData", dest="fracData", default=1,
#                  help="fraction of data events to draw")

(opt, args) = parser.parse_args()

if not opt.configFile:
    parser.error('config file not provided')

#if not opt.outputdir:
#    parser.error('output dir not provided')

#if not opt.outFileName:
#    parser.error('output filename not provided')


#os.system("mkdir -p "+opt.outputdir)

## Read input files
cfg = open(opt.configFile, "r")

treeName = ""
histoCounterName = ""
histoCounterBin = -1
cutSequence = ""
lumi = -1
pileupvar = ""
pileupdata = ""
pileupmc = ""
h_pileupdata = ""
h_pileupmc = ""
applyPUweight = 0

listSamples = []
listVars = []
listDraw = []

inputFile = {}
background = {}
signal = {}


# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
 
output = TFile.Open('TMVA_BDT_single_event_2.root', 'RECREATE')
factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification')
			#!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification')

dataloader = TMVA.DataLoader('dataset')


dataloader.AddVariable('dPhi_muj_ak4')
dataloader.AddVariable('dEta_muj_ak4')
dataloader.AddVariable('nAK4jets')
dataloader.AddVariable('metSig')
dataloader.AddVariable('muonpt_over_jet_pt')
dataloader.AddVariable('jet_pt_over_m_muj')
dataloader.AddVariable('muon_pt_over_m_muj')

for line in cfg:

    #skip commented out lines or empty lines
    if (line.startswith("#")):
        continue
    if line in ['\n','\r\n']:
        continue

    line = line.rstrip('\n')
    splitline = line.split(" ")
    linetype = splitline[0]
    linesize = len(line.split(" "))

    if (linetype == "TREENAME=" and linesize==2):
        treeName = splitline[1]

    if (linetype == "HISTONAME=" and linesize==3):
        histoCounterName = splitline[1]
        histoCounterBin = int(splitline[2])

    if (linetype == "CUTS=" and linesize==2):
        cuts = splitline[1]
        cutSequence = cuts

    if (linetype == "LUMI=" and linesize==2):
        lumi = splitline[1]

    if (linetype == "SAMPLE=" and linesize==6):
        listSamples.append(line)

    if (linetype == "VAR=" and linesize==8):
        listVars.append(line)

    if (linetype == "DRAW=" and linesize==4):
        listDraw.append(line)

    if (linetype == "APPLYPUWEIGHT=" and linesize==2):
        applyPUweight = int(splitline[1])

    if (linetype == "PILEUPVAR=" and linesize==2):
        pileupvar = splitline[1]

    if (linetype == "PILEUPDATA=" and linesize==3):
        pileupdata = splitline[1]
        h_pileupdata = splitline[2]

    if (linetype == "PILEUPMC=" and linesize==3):
        pileupmc = splitline[1]
        h_pileupmc = splitline[2]


#lumieff = float(opt.fracData) * float(lumi)
#change the CMS_lumi variables (see CMS_lumi.py)
#CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(lumieff/1000.)
#CMS_lumi.writeExtraText = 1
#CMS_lumi.extraText = "Preliminary"
#CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4


dictDraw = OrderedDict() #list of all histograms to be plotted for this variable
for lineD in listDraw:
    lineD = lineD.rstrip('')
    splitlineD = lineD.split(" ")
    dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[2])
    dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[3])

print(dictDraw)


print(listSamples)
for lineS in listSamples:
        lineS = lineS.rstrip('\n')
        splitlineS = lineS.split(" ")

        sample = splitlineS[1]
        files = splitlineS[2]
        xsec = splitlineS[3]
        kfac = splitlineS[4]
        name = splitlineS[5]
        ngen = 0

        #print(kfac)
        #print("Processing "+sample)

        chain = TChain(treeName)
        listFiles = files.split(";")
        for subFiles in listFiles:
          #  print("adding "+subFiles)
            chain.Add('%s' % subFiles)
            inputFile['%s'%sample] = TFile.Open(subFiles)

            if(float(xsec)<0):
                thisNgen = inputFile['%s'%sample].Get(histoCounterName).GetBinContent(1)
            else:
                thisNgen = inputFile['%s'%sample].Get(histoCounterName).GetBinContent(int(histoCounterBin))
        #    print("Current Ngen: "+str(thisNgen))
            ngen=ngen+thisNgen
            #inputFile.Close()
            entries = chain.GetEntries()
         #   print("Total number of generated events: "+ str(ngen))
         #   print("Total number of selected events: "+str(entries))

            maxentries = 0

            lumieq=float(ngen)/float(xsec)
            weight=(float(lumi)*float(xsec))/(float(ngen))

            if 'LQ' in sample:
                signal['%s'%sample] = inputFile['%s'%sample].Get(treeName)
                print(sample)             
                for i,event in enumerate(signal['%s'%sample]):
                      weight_sample=event.weight_PU*event.weight_Generator*event.Muon1_weight_all*event.Muon1_weight_Trigger*float(kfac)
                      tot_weight=float(weight_sample)*float(weight)
                      v = std.vector('double')() 
                      #for hh in ['dPhi_muj_ak4']:
                      v.push_back(event.dPhi_muj_ak4) 	
                      v.push_back(event.dEta_muj_ak4)
                      v.push_back(event.nAK4jets)
                      metSig=event.met/sqrt(event.sumEt)
                      ratio_pt=event.Muon1_Pt/event.AK4Jet1_Pt
                      muon_over_muj=event.Muon1_Pt/event.m_muj_ak4
                      jet_over_muj=event.AK4Jet1_Pt/event.m_muj_ak4
                      v.push_back(metSig)
                      v.push_back(ratio_pt)
                      v.push_back(jet_over_muj)
                      v.push_back(muon_over_muj)
                      if (i<signal['%s'%sample].GetEntries()/2.0):
                          dataloader.AddSignalTrainingEvent(v, tot_weight)
                      else:
                           dataloader.AddSignalTestEvent(v, tot_weight)

            else:
                background['%s'%sample] = inputFile['%s'%sample].Get(treeName)
                print(sample)             
                for i,event in enumerate(background['%s'%sample]):
                      if (i>background['%s'%sample].GetEntries()*(weight/0.2)):
                          print(i,background['%s'%sample].GetEntries(),background['%s'%sample].GetEntries()*(weight/0.2)) 
                          break
                      if 'WJets' in sample and (event.LHE_Vpt<100):
                          break
                      elif 'NLO' in sample and (event.LHE_Vpt>100):
                          break
                      weight_sample=event.weight_PU*event.weight_Generator*event.Muon1_weight_all*event.Muon1_weight_Trigger*float(kfac)
                      tot_weight=float(weight_sample)*float(weight)
                      v = std.vector('double')() 
                      #for hh in ['dPhi_muj_ak4']:
                      v.push_back(event.dPhi_muj_ak4) 	
                      v.push_back(event.dEta_muj_ak4)
                      v.push_back(event.nAK4jets)
                      metSig=event.met/sqrt(event.sumEt)
                      ratio_pt=event.Muon1_Pt/event.AK4Jet1_Pt
                      muon_over_muj=event.Muon1_Pt/event.m_muj_ak4
                      jet_over_muj=event.AK4Jet1_Pt/event.m_muj_ak4
                      v.push_back(metSig)
                      v.push_back(ratio_pt)
                      v.push_back(jet_over_muj)
                      v.push_back(muon_over_muj)
                      if (weight/0.2>=1) and (i<background['%s'%sample].GetEntries()/2.0): 
                          dataloader.AddBackgroundTrainingEvent(v, tot_weight)
                      elif (weight/0.2<=1) and (i<(background['%s'%sample].GetEntries()*(weight/0.2))/2.0):
                          dataloader.AddBackgroundTrainingEvent(v, tot_weight) 
                      else:
                          dataloader.AddBackgroundTestEvent(v, tot_weight)   
                #dataloader.AddBackgroundTree(background['%s'%sample], weight)	


#dataloader.AddVariable('met/sqrt(sumEt)')
#dataloader.AddVariable('Muon1_Pt/AK4Jet1_Pt')
#dataloader.AddVariable('nAK4jets')
#dataloader.AddVariable('AK4Jet1_Pt/m_muj_ak4')
#dataloader.AddVariable('Muon1_Pt/m_muj_ak4')


# Load data
#if not isfile('tmva_class_example.root'):
#    call(['curl', '-L', '-O', 'http://root.cern.ch/files/tmva_class_example.root'])
 
#if not isfile('tmva_example_multiple_background.root'):
#    createDataMacro = str(gROOT.GetTutorialDir()) + '/tmva/createData.C'
#    print(createDataMacro)
#    gROOT.ProcessLine('.L {}'.format(createDataMacro))
#    gROOT.ProcessLine('create_MultipleBackground(4000)')

#data = TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/SingleLQ_umuLQumu_M3000_Lambda1p0_reduced_skim.root')
#signal = data.Get('rootTupleTree/tree')
#dataB= TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/WJetsToLNu_HT-100To200_md_2018_reduced_skim.root')
#background = dataB.Get('rootTupleTree/tree')
#dataB1= TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/WJetsToLNu_HT-200To400_md_2018_reduced_skim.root')
#background1 = dataB1.Get('rootTupleTree/tree')
#dataB2= TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/WJetsToLNu_HT-400To600_md_2018_reduced_skim.root')
#background2 = dataB2.Get('rootTupleTree/tree')
 

#for branch in signal.GetListOfBranches():
#    dataloader.AddVariable(branch.GetName())

 



#dataloader.AddSignalTree(signal, 1.0)
#dataloader.AddBackgroundTree(background, 0.5)
#dataloader.AddBackgroundTree(background1, 0.5)
#dataloader.AddBackgroundTree(background2, 0.5)

#mycutb ='passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05'

#mycuts='PassJSON==1&&passHLTMuon==1&&fabs(Muon1_Eta)<2.4&&Muon1_Pt>55&&AK4Jet1_Pt>90&&Muon1_Pt/AK4Jet1_Pt>0.6&&Muon1_Pt/AK4Jet1_Pt<1.6&&dPhi_muj_ak4>2.5&&met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05'
#mycutb='PassJSON==1&&passHLTMuon==1&&fabs(Muon1_Eta)<2.4&&Muon1_Pt>55&&AK4Jet1_Pt>90&&Muon1_Pt/AK4Jet1_Pt>0.6&&Muon1_Pt/AK4Jet1_Pt<1.6&&dPhi_muj_ak4>2.5&&met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05'
mycuts=''
mycutb=''

dataloader.PrepareTrainingAndTestTree(mycuts,mycutb,
                                      'nTrain_Signal=:nTrain_Background=:SplitMode=Random:NormMode=EqualNumEvents:!V')
 
 
# Generate model
 
# Define model
 
 
# Book methods
#factory.BookMethod(dataloader, TMVA.Types.kCuts, 'Cuts',
 #                  '!H:!V:FitMethod=GA:EffSel:VarProp=FSmart')
#factory.BookMethod(dataloader, TMVA.Types.kFisher, 'Fisher',
#                   '!H:!V:Fisher:VarTransform=D,G')
#factory.BookMethod(dataloader, TMVA.Types.kPyTorch, 'PyTorch',
#          	   'H:!V:VarTransform=D,G:FilenameModel=model.pt:NumEpochs=20:BatchSize=32')
factory.BookMethod(dataloader, TMVA.Types.kBDT, 'BDT',
		    '!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=25')
                   #'!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20')
#		   '!H:!V:VarTransform=D,G:NTrees=1000:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=1:nCuts=20:MaxDepth=4') 
 
# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
 
 
# Plot ROC Curves
roc = factory.GetROCCurve(dataloader)
roc.SaveAs('ROC_ClassificationPyTorch_BDT_2.png')
