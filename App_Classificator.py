
from ROOT import TMVA, TFile, TString
from array import array
from subprocess import call
from os.path import isfile


from ROOT import * 
 
# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
reader = TMVA.Reader("Color:!Silent")
 
 
# Load data
if not isfile('tmva_class_example.root'):
    call(['curl', '-L', '-O', 'http://root.cern.ch/files/tmva_class_example.root'])
 
data = TFile.Open('tmva_class_example.root')
signal = data.Get('TreeS')
background = data.Get('TreeB')
 
branches = {}
for branch in signal.GetListOfBranches():
    branchName = branch.GetName()
    branches[branchName] = array('f', [-999])
    reader.AddVariable(branchName, branches[branchName])
    signal.SetBranchAddress(branchName, branches[branchName])
    background.SetBranchAddress(branchName, branches[branchName])
 
 
 
 
# Book methods
reader.BookMVA('BDT', TString('dataset/weights/TMVAClassification_BDT.weights.xml'))
 
 
# Print some example classifications
print('Some signal example classifications:')
for i in range(20):
    signal.GetEntry(i)
    print(reader.EvaluateMVA('BDT'))
print('')
 
print('Some background example classifications:')
for i in range(20):
    background.GetEntry(i)
    print(reader.EvaluateMVA('BDT'))



histoss = TH1F('Signal', 'Signal', 100, -1, 1)
histob = TH1F('background', 'background', 100, -1, 1)
for i,event in enumerate(signal):
	histoss.Fill(reader.EvaluateMVA('BDT'))	
	
for i,event in enumerate(background):
	histob.Fill(reader.EvaluateMVA('BDT'))	

c = TCanvas("c_","c_",800,700)


histob.SetLineColor(kRed)
histob.Draw()
histoss.Draw('SAME')

c.SaveAs('Example.png')
