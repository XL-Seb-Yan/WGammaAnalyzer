void TMVAnalysis( )
{
  TFile* outFile = TFile::Open("MC1600.root", "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis", outFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("MC1400");

  TFile *infile_s = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC1600N_WGamma_full_full.root");
  TFile *infile_b = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017BCD.root");

  dataloader->AddSignalTree((TTree*)infile_s->Get("Events"));
  dataloader->AddBackgroundTree((TTree*)infile_b->Get("Events"));
  //dataloader->AddVariable("ak8puppijet_tau21", 'F');
  dataloader->AddVariable("sys_costhetastar", 'F');
  dataloader->AddVariable("sys_ptoverm", 'F');
  dataloader->AddVariable("sys_seperation", 'F');
  dataloader->AddVariable("photon_pt", 'F');
  dataloader->AddVariable("photon_eta", 'F');
  dataloader->AddVariable("ak8puppijet_pt", 'F');
  dataloader->AddVariable("ak8puppijet_eta", 'F');
  dataloader->AddVariable("ak8puppijet_e", 'F');

  TCut cut_s = "(ak8puppijet_masssoftdropcorr > 65 && ak8puppijet_masssoftdropcorr < 105) && (ak8puppijet_eta > -2 && ak8puppijet_eta < 2) && (photon_eta > -1.44 && photon_eta < 1.44)";
  
  TCut cut_b = "(sys_invmass > 1500 && sys_invmass < 1600) && (ak8puppijet_masssoftdropcorr > 40 && ak8puppijet_masssoftdropcorr < 65) && (ak8puppijet_eta > -2 && ak8puppijet_eta < 2) && (photon_eta > -1.44&& photon_eta < 1.44)";

  dataloader->PrepareTrainingAndTestTree(cut_s,cut_b,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  //factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=-1" );

  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );



  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outFile->Close();
  delete factory;
}

// Use TMVA::TMVAGui(¡°TMVA.root¡±) in ROOT prompt to launch GUI
