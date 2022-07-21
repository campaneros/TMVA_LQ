#!/usr/bin/env python                                                                                                                                             

from ROOT import TMVA, TFile, TTree, TCut, gROOT, RDataFrame, Numba
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

parser.add_option("-o", "--output", dest="outputdir",
                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

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
os.system("mkdir -p "+opt.outputdir)

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

dfs={}
dfs_1muon={}
dfs_2muon={}
dfb={}
dfb_1muon={}
dfb_2muon={}

# Setup TMVA
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
            weight=(float(lumi)*float(xsec)*float(kfac))/(float(ngen))

            if 'LQ' in sample:
                print(sample)
                #signal['%s'%sample] = inputFile['%s'%sample].Get(treeName)
                dfs['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_all","weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{}".format(weight))
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_sample","{}".format(weight))
                dfs_1muon['%s'%sample]=dfs['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05").Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                dfs_2muon['%s'%sample]=dfs['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1").Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)

            else:
                print(sample)
                #background['%s'%sample] = inputFile['%s'%sample].Get(treeName)
                dfb['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
                dfb['%s'%sample]=dfb['%s'%sample].Define("weight_all","weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{}".format(weight))
                dfb['%s'%sample]=dfb['%s'%sample].Define("weight_test","{}".format(weight))
                if 'WJets' in sample:
                	dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("LHE_Vpt>=100 && passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05").Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                	dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1").Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)
                elif 'NLO' in sample:
                	dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("LHE_Vpt<100 && passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05").Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                	dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1").Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)
                else:
                	dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05").Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                	dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.1 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1").Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)




