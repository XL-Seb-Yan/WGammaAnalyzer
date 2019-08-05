void TMVAnalysisCombineMC( )
{
  TFile* outFile = TFile::Open("MC8001600.root", "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis", outFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("MC8001600");

  TFile *infile_s = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC8001600_WGamma_6595_conbine.root");
  TFile *infile_b = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017_WGamma_50105_full.root");

  dataloader->AddSignalTree((TTree*)infile_s->Get("Events"));
  dataloader->AddBackgroundTree((TTree*)infile_b->Get("Events"));
  //dataloader->AddVariable("ak8puppijet_tau21", 'F');
  dataloader->AddVariable("sys_costhetastar", 'F');
  dataloader->AddVariable("sys_ptoverm", 'F');
  dataloader->AddVariable("sys_seperation", 'F');
  //dataloader->AddVariable("photon_pt", 'F');
  dataloader->AddVariable("photon_eta", 'F');
  //dataloader->AddVariable("ak8puppijet_pt", 'F');
  dataloader->AddVariable("ak8puppijet_eta", 'F');
  dataloader->AddVariable("ak8puppijet_pt / sys_invmass", 'F');
  dataloader->AddVariable("ak8puppijet_e / sys_invmass", 'F');

  TCut cut_s = "(ak8puppijet_masssoftdropcorr > 65 && ak8puppijet_masssoftdropcorr < 95) && (ak8puppijet_eta > -2 && ak8puppijet_eta < 2) && (photon_eta > -1.4 && photon_eta < 1.4) && ak8puppijet_tau21 < 0.4";
  TCut cut_b = "(ak8puppijet_masssoftdropcorr > 50 && ak8puppijet_masssoftdropcorr < 65) && (ak8puppijet_eta > -2 && ak8puppijet_eta < 2) && (photon_eta > -1.4 && photon_eta < 1.4) && ak8puppijet_tau21 < 0.4";

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
