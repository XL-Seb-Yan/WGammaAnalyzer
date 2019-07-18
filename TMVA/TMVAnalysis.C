void TMVAnalysis( )
{
  TFile* outFile = TFile::Open("TMVA.root", "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis", outFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

  TFile *infile_s = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/DEV/CMSSW_9_4_9/src/WGammaAnalyzer/TMVA/SignalMC800_WGamma_select.root");
  TFile *infile_b = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/DEV/CMSSW_9_4_9/src/WGammaAnalyzer/TMVA/SinglePhoton2017C_WGamma_select.root");

  dataloader->AddSignalTree((TTree*)infile_s->Get("Events"));
  dataloader->AddBackgroundTree((TTree*)infile_b->Get("Events"));
  dataloader->AddVariable("ak8puppijet_tau21", 'F');
  dataloader->AddVariable("sys_costhetastar", 'F');
  dataloader->AddVariable("sys_ptoverm", 'F');

  dataloader->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  //factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=-1" );

  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outFile->Close();
  delete factory;
}

// Use TMVA::TMVAGui(¡°TMVA.root¡±) in ROOT prompt to launch GUI
