#!/usr/bin/env python                                                                                                                                             

from ROOT import TMVA, TFile, TTree, TCut, gROOT, TString, TH1F, RDataFrame, VecOps
from subprocess import call
from os.path import isfile



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
from multiprocessing import Queue, Process, current_process, Pool

usage = "usage: To be run from RootTreeAnalyzer:  python analysis/makePlots.py -c analysis/makePlots.txt -o output/TEST_cms_lq_prod2/plots -f TestOutput --fracMC 1.0 --fracData 1.0"

parser = optparse.OptionParser(usage)

parser.add_option("-c", "--config", dest="configFile",
                  help="config file")

parser.add_option("-o", "--output", dest="outputdir",
                  help="the directory contains the output of the program. Can be AFS or EOS directory.")

#parser.add_option("-f", "--FileName", dest="FileName",
#                  help="the output filename")

parser.add_option("--step", dest="step", default=0.1,
                  help="BDT scale step")
parser.add_option("--mass", dest="mass", default=3000,
                  help="signal mass")
parser.add_option("--coupling", dest="coupling", default="1p0",
                  help="signal coupling")




#parser.add_option("--fracData", dest="fracData", default=1,
#                  help="fraction of data events to draw")

(opt, args) = parser.parse_args()

if not opt.configFile:
    parser.error('config file not provided')

if not opt.outputdir:
    parser.error('output dir not provided')

#if 2/float(opt.step) %2 != 0:
#    parser.error('wrong step choise, 2/step should be an integer, you obtain '+str(2/float(opt.step)))


#if not opt.FileName:
#    parser.error('output filename not provided')

def variable(bdt):
	#global n_sig
	#global n_back
	n_signal=0
	n_background=0
	
	#n_sig['%s'%bdt]=0
	#n_back['%s'%bdt]=0


	#print(listSamples)
	for i,lineV in enumerate(listSamples):
			lineV = lineV.rstrip('\n')
			splitlineV = lineV.split(" ")
	
			sample = splitlineV[1]
			files = splitlineV[2]
			xsec = splitlineV[3]
			kfac = splitlineV[4]
			name = splitlineV[5]
			ngen = 0
	#		print(bdt,sample)

			chain = TChain(treeName)
			listFiles = files.split(";")
			for subFiles in listFiles:
	#			print(subFiles)
				chain.Add('%s' % subFiles)
				inputFile['%s'%sample] = TFile.Open(subFiles)
				if 'LQ' in sample:	
					if (str(opt.mass) in sample) and (str(opt.coupling) in sample):
						signal['%s'%sample] = inputFile['%s'%sample].Get(treeName)
						df['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
						df['%s_%s'%(sample,bdt)]=df['%s'%sample].Filter("(m_muj_ak4>({}-{}*2*0.05))&&(m_muj_ak4<({}+{}*2*0.05))&&DPhi_muj_ak4>3.08&&met_over_sq_sumEt<3&&Muon1_Pt_over_AK4Jet1_Pt>0.9&&AK4Jet1_Pt_over_m_muj_ak4>0.3&&Muon1_PFIso<{}".format(int(opt.mass),int(opt.mass),int(opt.mass),int(opt.mass),bdt))
						histo['%s_%s'%(sample,bdt)]=df['%s_%s'%(sample,bdt)].Histo1D(ROOT.RDF.TH1DModel("Signal","Signal", 100, int(opt.mass)-int(opt.mass)*0.10, int(opt.mass)+int(opt.mass)*0.10),"m_muj_ak4","weight_all")
						histo['clone_%s_%s'%(sample,bdt)]= histo['%s_%s'%(sample,bdt)].Clone()
						n_signal+=float(histo['clone_%s_%s'%(sample,bdt)].GetSumOfWeights())
						#if i==0:
							#print(n_signal)

				else:
					background['%s'%sample] = inputFile['%s'%sample].Get(treeName)
					if (background['%s'%sample].GetEntries()>0):
						df['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
						df['%s_%s'%(sample,bdt)]=df['%s'%sample].Filter("(m_muj_ak4>({}-{}*2*0.05))&&(m_muj_ak4<({}+{}*2*0.05))&&DPhi_muj_ak4>3.08&&met_over_sq_sumEt<3&&Muon1_Pt_over_AK4Jet1_Pt>0.9&&AK4Jet1_Pt_over_m_muj_ak4>0.3&&Muon1_PFIso<{}".format(int(opt.mass),int(opt.mass),int(opt.mass),int(opt.mass),bdt))
						histo['%s_%s'%(sample,bdt)]=df['%s_%s'%(sample,bdt)].Histo1D(ROOT.RDF.TH1DModel("Background","Background", 100, int(opt.mass)-int(opt.mass)*0.10, int(opt.mass)+int(opt.mass)*0.10),"m_muj_ak4","weight_all")
						histo['clone_%s_%s'%(sample,bdt)]= histo['%s_%s'%(sample,bdt)].Clone()
						n_background+= float(histo['clone_%s_%s'%(sample,bdt)].GetSumOfWeights())

			#print (n_signal,n_background)
			#return (n_signal,n_background)
	#print(n_signal,n_background)
	#n_sig['%s'%bdt]=n_signal
	#n_back['%s'%bdt]=n_background
	print(n_signal,n_background,bdt)
	return (n_signal,n_background)					

#def main():
	#process = Pool(10)
	#BDT=[]
	#for bdt in np.arange(-1,1,0.1):
	#	BDT.append(bdt)	



if __name__ == '__main__':
	output=opt.outputdir
	os.system("mkdir -p "+output)

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
	df={}


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
        
        
 		



	c = TCanvas("c_","c_",800,700)

	histo={}
	histobd={}

	histos = TH1D('Signal', 'Signal', 100, -1, 1)
	histob = TH1D('background', 'background', 100, -1, 1)

	n_back={}
	n_sig={}
	nbin=int(0.1/float(opt.step))
	histo['Significance']=TH1D('Significance','Significance',nbin,0,0.1)
		#print(i)
		#n_back=0
		#n_sig=0
		#n_signal=0
		#n_background=0
	print('###################')
	pool=Pool(processes=25)
	n = pool.map(variable,np.arange(0,0.1,float(opt.step)))
	pool.close()
	pool.join()
	print('###################')
	#n.reverse()
	for i,k in enumerate(n):
	#print(i,n_sig,n_back)
		#processes = [Process(target=variable, args=(lineV,)) for lineV in listSamples]
		#for p in processes:
        	#        p.start()
	
		#for p in processes:
                #	p.join()
		#print(n_signal,n_background)
		if k[1]<0:
			back=0
		else:
			back=k[1]
			print(k[0],back, k[0]/(3/2+sqrt(back)), i)
		histo['Significance'].SetBinContent(i+1,k[0]/(3/2+sqrt(back)))
#		for mass in [1000,2000,3000]:
#			n_signal+=n_sig['%s_%s'%(mass,bdt)]
#			n_background+=n_back['%s_%s'%(mass,bdt)]
		

#       print dictDraw
        #outFile.Close()
#       print "Output file created: "+opt.outputdir+"/"+opt.outFileName+".root"
	#for i,bdt in enumerate(np.arange(-1,1,0.1)):
	#	print(i,n_sig['%s'%bdt],n_back['%s'%bdt], n_sig['%s'%bdt]/(0.4990/2+sqrt(n_back['%s'%bdt])))
	#	histo['Significance'].SetBinContent(i,n_sig['%s'%bdt]/(0.4990/2+sqrt(n_back['%s'%bdt])))
	histo['Significance'].Draw()
	c.SaveAs(output+"/LQ_"+str(opt.mass)+"_"+str(opt.coupling)+".png")


	cfg.close()

	with open('log_job.txt', 'w') as f:
		f.write("OK")
