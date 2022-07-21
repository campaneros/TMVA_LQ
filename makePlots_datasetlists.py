#! /usr/bin/env python

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

parser.add_option("-f", "--outFileName", dest="outFileName",
                  help="the output filename")

parser.add_option("--fracMC", dest="fracMC", default=1,
                  help="fraction of data events to draw")

parser.add_option("--fracData", dest="fracData", default=1,
                  help="fraction of data events to draw")

(opt, args) = parser.parse_args()

if not opt.configFile:   
    parser.error('config file not provided')

if not opt.outputdir:   
    parser.error('output dir not provided')

if not opt.outFileName:   
    parser.error('output filename not provided')

################################################

## ROOT settings
#gROOT.Reset()
gROOT.SetBatch(True)

## Set style
tdrstyle.setTDRStyle()

## ???????????????????
dump = TH1
# necessary,  otherwise it crashes with 
#Traceback (most recent call last):
#  File "scripts/makePlots.py", line 271, in <module>
#    pad1.SetTopMargin(0.1)
#NameError: name 'pad1' is not defined
## ???????????????????

## Create output dir
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


lumieff = float(opt.fracData) * float(lumi) 
#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = "%.1f fb^{-1}"%(lumieff/1000.)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

## Pileup weights
if(applyPUweight):
    gROOT.ProcessLine('TFile *inputFileData=TFile::Open("%s");'%pileupdata)
    gROOT.ProcessLine('TFile *inputFileMC=TFile::Open("%s");'%pileupmc)
    gROOT.ProcessLine('TH1F* h_data = (TH1F*) inputFileData->Get("%s");'%h_pileupdata)
    gROOT.ProcessLine('TH1F* h_mc = (TH1F*) inputFileMC->Get("%s");'%h_pileupmc)
    gROOT.ProcessLine('h_data->Scale(1/h_data->Integral());')
    gROOT.ProcessLine('h_mc->Scale(1/h_mc->Integral());')
    gROOT.ProcessLine('TH1F* h_data_over_mc = (TH1F*) h_data->Clone();')
    gROOT.ProcessLine('h_data_over_mc->Divide(h_data,h_mc);')
    gROOT.ProcessLine('double getPUweight(double x){return h_data_over_mc->GetBinContent(x+1);};')

## Sum of two 2D vectors
gROOT.ProcessLine('double pTsum(double pt1, double phi1, double pt2, double phi2){TVector2 v1;TVector2 v2;TVector2 vsum;v1.SetMagPhi(pt1,phi1);v2.SetMagPhi(pt2,phi2);vsum=v1+v2;return vsum.Mod();};')
#gROOT.ProcessLine('double mTsum(double Et1, double phi1, double Et2, double phi2){double Mt;Mt=sqrt(2*Et1*Et2(-1*cos(phi1-phi2)));return Mt;};')
gROOT.ProcessLine('double mTsum(double Et1, double phi1, double Et2, double phi2){double Mt;Mt=sqrt(2*Et1*Et2*(1-cos(phi1-phi2)));return Mt;};')

## Make plots
outFile = open (opt.outFileName+".txt","w")

outFile.write('SAMPLE							NAME			cross_section			Ngen		Lumi		Lumieq			weights\n')

## Loop over draw options
dictDraw = OrderedDict() #list of all histograms to be plotted for this variable
for lineD in listDraw:
    lineD = lineD.rstrip('')
    splitlineD = lineD.split(" ")    
    dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[2])
    dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[3])

print dictDraw

## Loop over variables
print listSamples
for lineS in listSamples:    
        lineS = lineS.rstrip('\n')
        splitlineS = lineS.split(" ")

        sample = splitlineS[1]
        files = splitlineS[2]
        xsec = splitlineS[3]
        kfac = splitlineS[4]
	name = splitlineS[5]
        ngen = 0

	print kfac
        print "Processing "+sample

	
        chain = TChain(treeName)        
        listFiles = files.split(";")
        print listFiles
	for subFiles in listFiles:
            print "adding "+subFiles
            chain.Add('%s' % subFiles)
            inputFile = TFile.Open(subFiles)
            
	    if(float(xsec)<0):
	    	thisNgen = inputFile.Get(histoCounterName).GetBinContent(1)
	    else:
	    	thisNgen = inputFile.Get(histoCounterName).GetBinContent(int(histoCounterBin)) 
            print "Current Ngen: "+str(thisNgen)
            ngen=ngen+thisNgen
            inputFile.Close()

        entries = chain.GetEntries()
        print "Total number of generated events: "+ str(ngen)
        print "Total number of selected events: "+str(entries)

        maxentries = 0
		
	lumieq=float(ngen)/float(xsec)
	weight=(float(lumi)*float(xsec))/(float(ngen))
	outFile.write(sample+'					'+name+'		'+str(xsec)+'		'+str(ngen)+'		'+str(lumi)+'		'+str(lumieq)+'		'+ str(weight))
	outFile.write('\n')

outFile.close()

#cfg.close()

## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]

