#!/usr/bin/env python                                                                                                                                             

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
from multiprocessing import Queue, Process, current_process

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


# Setup TMVA



def variable(lineV,histoCounterBini,histoCounterName):
        lineV = lineV.rstrip('\n')
        splitlineV = lineV.split(" ")

        sample = splitlineV[1]
        files = splitlineV[2]
        xsec = splitlineV[3]
        kfac = splitlineV[4]
        name = splitlineV[5]
        ngen = 0
        realngen=0
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
                realNgen = inputFile['%s'%sample].Get(histoCounterName).GetBinContent(int(1))
            else:
                #thisNgen = inputFile['%s'%sample].Get(histoCounterName).GetBinContent(int(3))
                thisNgen = inputFile['%s'%sample].Get(histoCounterName).GetBinContent(int(histoCounterBin))
                realNgen = inputFile['%s'%sample].Get(histoCounterName).GetBinContent(int(1))
        #    print("Current Ngen: "+str(thisNgen))
            ngen=ngen+thisNgen
            realngen=realngen+realNgen
            #inputFile.Close()
            entries = chain.GetEntries()
         #   print("Total number of generated events: "+ str(ngen))
         #   print("Total number of selected events: "+str(entries))

            maxentries = 0

            lumieq=float(ngen)/float(xsec)
            weight=(float(lumi)*float(xsec)*float(kfac))/(float(ngen))
            weight_real=(float(lumi)*float(xsec)*float(kfac))/(float(realngen))
  



            if 'LQ' in sample:
                #signal['%s'%sample] = inputFile['%s'%sample].Get(treeName)
                dfs['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_all","float(weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{})".format(weight))
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_sample","float({})".format(weight))
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_sample_nogen","float({})".format(weight_real))
                dfs['%s'%sample]=dfs['%s'%sample].Define("met_over_sq_sumEt","float(met/sqrt(sumEt))")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon1_Pt_over_AK4Jet1_Pt","float(Muon1_Pt/AK4Jet1_Pt)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon1_Pt_over_m_muj_ak4","float(Muon1_Pt/m_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("AK4Jet1_Pt_over_m_muj_ak4","float(AK4Jet1_Pt/m_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon2_Pt_over_m_muj_ak4","float(Muon2_Pt/m_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon1_Eta_m_Muon2_Eta","float(Muon1_Eta-Muon2_Eta)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("DPhi_muj_ak4","float(dPhi_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("DEta_muj_ak4","float(dEta_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("NAK4jets","float(nAK4jets)")
                dfs_1muon['%s'%sample]=dfs['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05").Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                dfs_2muon['%s'%sample]=dfs['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110").Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)

            elif "DAT" in name:
                print(sample)
                dfs['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
                weight=-1 
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_all","float({})".format(weight))
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_sample","float({})".format(weight))
                dfs['%s'%sample]=dfs['%s'%sample].Define("weight_sample_nogen","float({})".format(weight))
                dfs['%s'%sample]=dfs['%s'%sample].Define("met_over_sq_sumEt","float(met/sqrt(sumEt))")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon1_Pt_over_AK4Jet1_Pt","float(Muon1_Pt/AK4Jet1_Pt)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon1_Pt_over_m_muj_ak4","float(Muon1_Pt/m_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("AK4Jet1_Pt_over_m_muj_ak4","float(AK4Jet1_Pt/m_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon2_Pt_over_m_muj_ak4","float(Muon2_Pt/m_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("Muon1_Eta_m_Muon2_Eta","float(Muon1_Eta-Muon2_Eta)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("DPhi_muj_ak4","float(dPhi_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("DEta_muj_ak4","float(dEta_muj_ak4)")
                dfs['%s'%sample]=dfs['%s'%sample].Define("NAK4jets","float(nAK4jets)")
                dfs_1muon['%s'%sample]=dfs['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05").Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                dfs_2muon['%s'%sample]=dfs['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110").Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)




            else:
                test=int(0.2/weight_real)
                #background['%s'%sample] = inputFile['%s'%sample].Get(treeName)
                dfb['%s'%sample]=RDataFrame(treeName,inputFile['%s'%sample])
                if (int(0.2/weight_real)<=1):
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_all","float(weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{})".format(weight))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_sample","float({})".format(weight))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_sample_real","float({})".format(weight_real))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_all_1muon","float(weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{})".format(weight))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_sample_1muon","float({})".format(weight))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("met_over_sq_sumEt","float(met/sqrt(sumEt))")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("Muon1_Pt_over_AK4Jet1_Pt","float(Muon1_Pt/AK4Jet1_Pt)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("Muon1_Pt_over_m_muj_ak4","float(Muon1_Pt/m_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("AK4Jet1_Pt_over_m_muj_ak4","float(AK4Jet1_Pt/m_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("Muon2_Pt_over_m_muj_ak4","float(Muon2_Pt/m_muj_ak4)")
               		dfb['%s'%sample]=dfb['%s'%sample].Define("Muon1_Eta_m_Muon2_Eta","float(Muon1_Eta-Muon2_Eta)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("DPhi_muj_ak4","float(dPhi_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("DEta_muj_ak4","float(dEta_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("NAK4jets","float(nAK4jets)")
                	if 'WJets' in sample:
                		dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("LHE_Vpt>=100 && passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05")
                		dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110")
                	elif 'NLO' in sample:
                		dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("LHE_Vpt<100 && passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05")
                		dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110")
                	else:
                		dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05")
                		dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110")
                	dfb_1muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                	dfb_1muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset_reduced.root'%sample)
                	dfb_2muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)
                else:
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_all","float(weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{})".format(weight))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_sample","float({})".format(weight))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("weight_sample_real","float({})".format(weight_real))
                	dfb['%s'%sample]=dfb['%s'%sample].Define("met_over_sq_sumEt","float(met/sqrt(sumEt))")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("Muon1_Pt_over_AK4Jet1_Pt","float(Muon1_Pt/AK4Jet1_Pt)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("Muon1_Pt_over_m_muj_ak4","float(Muon1_Pt/m_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("AK4Jet1_Pt_over_m_muj_ak4","float(AK4Jet1_Pt/m_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("Muon2_Pt_over_m_muj_ak4","float(Muon2_Pt/m_muj_ak4)")
               		dfb['%s'%sample]=dfb['%s'%sample].Define("Muon1_Eta_m_Muon2_Eta","float(Muon1_Eta-Muon2_Eta)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("DPhi_muj_ak4","float(dPhi_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("DEta_muj_ak4","float(dEta_muj_ak4)")
                	dfb['%s'%sample]=dfb['%s'%sample].Define("NAK4jets","float(nAK4jets)")
                	if 'WJets' in sample:
                		dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("LHE_Vpt>=100 && passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05")
                		dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110")
                	elif 'NLO' in sample:
                		dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("LHE_Vpt<100 && passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05")
                		dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110")
                	else:
                		dfb_1muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&nMuonId==1&&Muon1_PFIso<0.05")
                		dfb_2muon['%s'%sample]=dfb['%s'%sample].Filter("passNoiseFilters==1 && PassJSON==1 && passHLTMuon==1 && fabs(Muon1_Eta)<2.4 && Muon1_PFIso<0.05 && Muon1_Pt>55 && AK4Jet1_Pt>90 && Muon1_Pt/AK4Jet1_Pt>0.6 && Muon1_Pt/AK4Jet1_Pt<1.6 && dPhi_muj_ak4>2.5 && met/sqrt(sumEt)<10&&nEleIdIso==0&&Muon1_PFIso<0.05&&nMuonId==2&&Muon2_TrkIso<0.1&&sqrt((Muon1_Phi-Muon2_Phi)*(Muon1_Phi-Muon2_Phi)+(Muon1_Eta-Muon2_Eta)*(Muon1_Eta-Muon2_Eta))>0.1&&m_mumu>110")
                	if ('WJets' in sample) or ('DY' in sample) or ('TT'in sample):
                		dfb_1muon['%s'%sample]=dfb_1muon['%s'%sample].Define("weight_all_1muon","float(weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{})".format(weight*int(0.2/weight_real)))
                		dfb_1muon['%s'%sample]=dfb_1muon['%s'%sample].Define("weight_sample_1muon","float({})".format(weight*int(0.2/weight_real)))
                		dfb_1muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                		int_we=int(0.2/weight_real)	
                		dfb_1muon['%s'%sample].Range(0,0,int_we).Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset_reduced.root'%sample)
                		dfb_2muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)
                	else:
                		dfb_1muon['%s'%sample]=dfb_1muon['%s'%sample].Define("weight_all_1muon","float(weight_PU*weight_Generator*Muon1_weight_all*Muon1_weight_Trigger*{})".format(weight))
                		dfb_1muon['%s'%sample]=dfb_1muon['%s'%sample].Define("weight_sample_1muon","float({})".format(weight))
                		dfb_1muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset.root'%sample)
                		dfb_1muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_1_muon_dataset_reduced.root'%sample)
                		dfb_2muon['%s'%sample].Snapshot(treeName,opt.outputdir+'/'+'%s_2_muon_dataset.root'%sample)
                print(sample)
                print(int(0.2/weight_real))


def main():
        processes = [Process(target=variable, args=(lineV,histoCounterBin,histoCounterName)) for lineV in listSamples]

        for p in processes:
                p.start()

        for p in processes:
                p.join()





if __name__ == '__main__':

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
        inputFile = {}
        background = {}
        signal = {}

        dfs={}
        dfs_1muon={}
        dfs_2muon={}
        dfb={}
        dfb_1muon={}
        dfb_2muon={}
#	dfb_reduced={}

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



        iPeriod = 4


        dictDraw = OrderedDict() #list of all histograms to be plotted for this variable
        for lineD in listDraw:
                lineD = lineD.rstrip('')
                splitlineD = lineD.split(" ")
                dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[2])
                dictDraw.setdefault(splitlineD[1],[]).append(splitlineD[3])


#       print dictDraw
        main()
        #outFile.Close()
#       print "Output file created: "+opt.outputdir+"/"+opt.outFileName+".root"


        cfg.close()
                    
