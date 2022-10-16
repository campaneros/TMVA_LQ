
from ROOT import TMVA, TFile, TTree, TCut, gROOT, RDataFrame, Numba
from subprocess import call
from os.path import isfile


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




usage = "usage: to run -f filename -o outputdir -c configfile"

parser = optparse.OptionParser(usage)

parser.add_option("-c", "--config", dest="configFile",
                  help="config file")

parser.add_option("-o", "--output", dest="outputdir",
                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

parser.add_option("-f", "--filename", dest="outfile",
                  help="the filename with all the correct path")

parser.add_option("-d","--directory", dest="rootdir")
(opt, args) = parser.parse_args()

if opt.outputdir:
	os.system("mkdir -p "+opt.outputdir)

## Read input files
cfg = open(opt.configFile, "r")


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


	with open(opt.outfile, 'w') as f:
		f.write('TREENAME= '+ treeName+'\n')	
		#f.write('HISTONAME= '+ histoCounterName +' '+ histoCounterBin +'\n')	
		for lineS in listSamples:
			lineS = lineS.rstrip('\n')
			splitlineS = lineS.split(" ")

			sample = splitlineS[1]
			files = splitlineS[2]
			xsec = splitlineS[3]
			kfac = splitlineS[4]
			name = splitlineS[5]

			root_before_f=files.split("/")[-1]
			root_before_f=root_before_f.split(".")
				
			if '1muon' in opt.outfile:	
				root_f=sample+'_1_muon_dataset.'+root_before_f[1]
			else:
				root_f=sample+'_2_muon_dataset.'+root_before_f[1]

			path= os.path.abspath(os.getcwd())
			join_path=os.path.join(path,opt.rootdir,root_f)
			chain = TChain(treeName)
			listFiles = files.split(";")
			for subFiles in listFiles:
				chain.Add('%s' % subFiles)
				inputFile['%s'%sample] = TFile.Open(subFiles)
				if (inputFile['%s'%sample].Get(treeName).GetEntries()>0):
					f.write('SAMPLE= '+sample+' '+join_path+' '+xsec+' '+kfac+' '+name+'\n')


cfg.close()
	


