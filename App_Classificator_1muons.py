#!/usr/bin/env python                                                                                                                                             

from ROOT import TMVA, TFile, TTree, TCut, gROOT, TString, TH1F, RDataFrame, VecOps
from subprocess import call
from os.path import isfile


import ROOT.VecOps as RV 

import numpy as np

import torch
from torch import nn

from array import array

import math as m 

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

parser.add_option("-f", "--FileName", dest="FileName",
                  help="the output filename")

#parser.add_option("--fracMC", dest="fracMC", default=1,
#                  help="fraction of data events to draw")

#parser.add_option("--fracData", dest="fracData", default=1,
#                  help="fraction of data events to draw")

(opt, args) = parser.parse_args()

if not opt.configFile:
    parser.error('config file not provided')

if not opt.outputdir:
    parser.error('output dir not provided')

if not opt.FileName:
    parser.error('output filename not provided')


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

gROOT.SetBatch(True)

# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
#reader = TMVA.Reader("Color:!Silent")

 

#dataloader = TMVA.DataLoader('dataset')




df={}


model=TMVA.Experimental.RReader('dataset/weights/'+opt.FileName+'_BDT.weights.xml')
variables = model.GetVariableNames()
lenght=len(variables)
print(variables)
print(lenght)



gInterpreter.ProcessLine('''TMVA::Experimental::RReader model("dataset/weights/%s_BDT.weights.xml");
computeModel = TMVA::Experimental::Compute<%d, float>(model);
'''%(opt.FileName,lenght))



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





dictDraw = OrderedDict() #list of all histograms to be plotted for this variable
for lineD in listDraw:
    lineD = lineD.rstrip('')
    splitlineD = lineD.split(" ")
    dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[2])
    dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[3])


c = TCanvas("c_","c_",800,700)

histo={}
histobd={}

histos = TH1D('Signal', 'Signal', 100, -1, 1)
histob = TH1D('background', 'background', 100, -1, 1)


#ROOT.gInterpreter.Declare("""


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
                df['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
                #histo['%s'%sample]=TH1F('Signal_%s'%sample, 'Signal_%s'%sample, 100, -1, 1)
                df['%s'%sample]=df['%s'%sample].Define("BDT",computeModel, model.GetVariableNames()).Snapshot(treeName,opt.outputdir+'/%s_1_muon_dataset.root'%sample)
                histo['%s'%sample]=df['%s'%sample].Histo1D(ROOT.RDF.TH1DModel("Signal","Signal", 100, -1., 1.),"BDT")
                histo['%s'%sample].Draw()
                for key in ['png','pdf']:
                	c.SaveAs(opt.outputdir+'/Signal_%s.%s'%(sample,key))

            else:
                background['%s'%sample] = inputFile['%s'%sample].Get(treeName)
                if (background['%s'%sample].GetEntries()>0):
                	df['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
                	df['%s'%sample]=df['%s'%sample].Define("BDT",computeModel, model.GetVariableNames()).Snapshot(treeName,opt.outputdir+'/%s_1_muon_dataset.root'%sample)
                	histobd['%s'%sample]=df['%s'%sample].Histo1D(ROOT.RDF.TH1DModel("Signal","Signal", 100, -1., 1.),"BDT","weight_all_1muon")
                	histobd['%s'%sample].Draw()
                	for key in ['png','pdf']:
                		c.SaveAs(opt.outputdir+'/Signal_%s.%s'%(sample,key))
#                	x=TMVA.Experimental.AsTensor['float'](df['%s'%sample], variables)

test={}

for i,hist in enumerate(histo):
	#histo['%s'%hist].Draw()
	if i == 0:
		histos=histo['%s'%hist].Clone()
	else:
		test['%s'%hist]=TH1F('Signal_%s'%sample, 'Signal_%s'%sample, 100, -1, 1)
		test['%s'%hist]=histo['%s'%hist].Clone()
		histos.Add(test['%s'%hist])

for i,hist in enumerate(histobd):
	#histo['%s'%hist].Draw()
	if i == 0:
		histob=histobd['%s'%hist].Clone()
	else:
		test['%s'%hist]=TH1F('Signal_%s'%sample, 'Signal_%s'%sample, 100, -1, 1)
		test['%s'%hist]=histobd['%s'%hist].Clone()
		histob.Add(test['%s'%hist])



histos.Draw('HIST')
for key in ['png','pdf']:
	c.SaveAs(opt.outputdir+'/Signal_df.%s'%key)

histob.Draw('HIST')
for key in ['png','pdf']:
	c.SaveAs(opt.outputdir+'/background.%s'%key)

histob.SetLineColor(kRed)
histob.Draw('HIST')
histos.Draw('SAME')
for key in ['png','pdf']:
	c.SaveAs(opt.outputdir+'/S_and_B.%s'%key) 
