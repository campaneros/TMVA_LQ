#!/usr/bin/env python                                                                                                                                             

from ROOT import TMVA, TFile, TTree, TCut, gROOT, TString
from subprocess import call
from os.path import isfile
 
import torch
from torch import nn

from array import array


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

parser.add_option("-f", "--FileName", dest="FileName",
                  help="the output filename")

#parser.add_option("--fracMC", dest="fracMC", default=1,
#                  help="fraction of data events to draw")

#parser.add_option("--fracData", dest="fracData", default=1,
#                  help="fraction of data events to draw")

(opt, args) = parser.parse_args()

if not opt.configFile:
    parser.error('config file not provided')

#if not opt.outputdir:
#    parser.error('output dir not provided')

if not opt.FileName:
    parser.error('output filename not provided')


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
reader = TMVA.Reader("Color:!Silent")

 

#dataloader = TMVA.DataLoader('dataset')





branches={}


for key in ['dPhi_muj_ak4','dEta_muj_ak4','met/sqrt(sumEt)','Muon1_Pt/AK4Jet1_Pt','nAK4jets','AK4Jet1_Pt/m_muj_ak4','Muon1_Pt/m_muj_ak4']:
	branches['%s' %key] = array('f', [-999])
	reader.AddVariable(key, branches['%s' %key])

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


            if 'LQ' in sample:
                signal['%s'%sample] = inputFile['%s'%sample].Get(treeName)
            else:
                background['%s'%sample] = inputFile['%s'%sample].Get(treeName)



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

 



#book methods
reader.BookMVA('BDT', TString('dataset/weights/'+opt.FileName+'_BDT.weights.xml'))
 
 
# Print some example classifications
print('Some signal example classifications:')
for i in range(200):
	for key in signal:
		signal['%s'%key].GetEntry(i)
		print(reader.EvaluateMVA('BDT'))
print('')
 
print('Some background example classifications:')
for i in range(2000):
	for key in background:
		background['%s'%key].GetEntry(i)
		print(reader.EvaluateMVA('BDT'))
 
 
