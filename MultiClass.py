from ROOT import TMVA, TFile, TTree, TCut, gROOT
from os.path import isfile
 
import torch
from torch import nn
 
 
# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()
 
output = TFile.Open('TMVA_noBDT.root', 'RECREATE')
factory = TMVA.Factory('TMVAClassification', output,
    '!V:!Silent:Color:DrawProgressBar:Transformations=D,G:AnalysisType=multiclass')
 
 
# Load data
if not isfile('tmva_example_multiple_background.root'):
    createDataMacro = str(gROOT.GetTutorialDir()) + '/tmva/createData.C'
    print(createDataMacro)
    gROOT.ProcessLine('.L {}'.format(createDataMacro))
    gROOT.ProcessLine('create_MultipleBackground(4000)')
 
data = TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/SingleLQ_umuLQumu_M3000_Lambda1p0_reduced_skim.root')
signal = data.Get('rootTupleTree/tree')
dataB= TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/WJetsToLNu_HT-100To200_md_2018_reduced_skim.root')
background = dataB.Get('rootTupleTree/tree')
dataB1= TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/WJetsToLNu_HT-200To400_md_2018_reduced_skim.root')
background1 = dataB1.Get('rootTupleTree/tree')
dataB2= TFile.Open('/data/mcampana/CMS/ntuples/cms_lq_prod1_mu_5/merged/WJetsToLNu_HT-400To600_md_2018_reduced_skim.root')
background2 = dataB2.Get('rootTupleTree/tree')
 

dataloader = TMVA.DataLoader('dataset')
#for branch in signal.GetListOfBranches():
#    dataloader.AddVariable(branch.GetName())

 
dataloader.AddVariable('dPhi_muj_ak4')
dataloader.AddVariable('AK4Jet1_Pt')
dataloader.AddVariable('Muon1_Pt')
dataloader.AddVariable('AK4Jet1_Pt/Muon1_Pt')



dataloader.AddTree(signal,'Signal')
dataloader.AddTree(background, 'bck_0')
dataloader.AddTree(background1, 'bck_1')
dataloader.AddTree(background2, 'bck_2')

dataloader.PrepareTrainingAndTestTree(TCut(''),
                                      'SplitMode=Random:NormMode=EqualNumEvents:!V')
 
#for branch in signal.GetListOfBranches():
#    dataloader.AddVariable(branch.GetName())
 
 
# Generate model
# Define model
model = nn.Sequential()
model.add_module('linear_1', nn.Linear(in_features=4, out_features=256))
model.add_module('relu', nn.ReLU())
model.add_module('linear_2', nn.Linear(in_features=256, out_features=4))
model.add_module('softmax', nn.Softmax(dim=1))
 
 
# Set loss and optimizer
loss = nn.CrossEntropyLoss()
optimizer = torch.optim.SGD
 
 
# Define train function
def train(model, train_loader, val_loader, num_epochs, batch_size, optimizer, criterion, save_best, scheduler):
    trainer = optimizer(model.parameters(), lr=0.01)
    schedule, schedulerSteps = scheduler
    best_val = None
 
    for epoch in range(num_epochs):
        # Training Loop
        # Set to train mode
        model.train()
        running_train_loss = 0.0
        running_val_loss = 0.0
        for i, (X, y) in enumerate(train_loader):
            trainer.zero_grad()
            output = model(X)
            target = torch.max(y, 1)[1]
            train_loss = criterion(output, target)
            train_loss.backward()
            trainer.step()
 
            # print train statistics
            running_train_loss += train_loss.item()
            if i % 256 == 255:    # print every 32 mini-batches
                print("[{}, {}] train loss: {:.3f}".format(epoch+1, i+1, running_train_loss / 256))
                running_train_loss = 0.0
 
        if schedule:
            schedule(optimizer, epoch, schedulerSteps)
 
        # Validation Loop
        # Set to eval mode
        model.eval()
        with torch.no_grad():
            for i, (X, y) in enumerate(val_loader):
                output = model(X)
                target = torch.max(y, 1)[1]
                val_loss = criterion(output, target)
                running_val_loss += val_loss.item()
 
            curr_val = running_val_loss / len(val_loader)
            if save_best:
               if best_val==None:
                   best_val = curr_val
               best_val = save_best(model, curr_val, best_val)
 
            # print val statistics per epoch
            print("[{}] val loss: {:.3f}".format(epoch+1, curr_val))
            running_val_loss = 0.0
 
    print("Finished Training on {} Epochs!".format(epoch+1))
 
    return model
 
 
# Define predict function
def predict(model, test_X, batch_size=256):
    # Set to eval mode
    model.eval()
   
    test_dataset = torch.utils.data.TensorDataset(torch.Tensor(test_X))
    test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
 
    predictions = []
    with torch.no_grad():
        for i, data in enumerate(test_loader):
            X = data[0]
            outputs = model(X)
            predictions.append(outputs)
        preds = torch.cat(predictions)
   
    return preds.numpy()
 
 
load_model_custom_objects = {"optimizer": optimizer, "criterion": loss, "train_func": train, "predict_func": predict}
 
 
# Store model to file
# Convert the model to torchscript before saving
m = torch.jit.script(model)
torch.jit.save(m, "model.pt")
print(m)
 
 
# Book methods
factory.BookMethod(dataloader, TMVA.Types.kFisher, 'Fisher',
        '!H:!V:Fisher:VarTransform=D,G')
factory.BookMethod(dataloader, TMVA.Types.kPyTorch, "PyTorch",
        'H:!V:VarTransform=D,G:FilenameModel=model.pt:NumEpochs=20:BatchSize=256')
factory.BookMethod(dataloader, TMVA.Types.kBDT, 'BDT',
	'!H:!V:VarTransform=D,G:NTrees=1000:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=4') 
 
 
# Run TMVA
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
 
# Plot ROC Curves
roc = factory.GetROCCurve(dataloader)
roc.SaveAs('ROC_MulticlassPyTorch.png')
