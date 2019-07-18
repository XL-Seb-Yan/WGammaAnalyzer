void TMVAnalysis( )
{
  TFile* outFile = TFile::Open("TMVA.root", "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis", outFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TFile *infile_s = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SignalMC800_WGamma_select.root");
  TFile *infile_b = TFile::Open("/home/xyan13/WGProj/CMSSW_9_4_9/src/WGammaAnalyzer/Selection/SelOutPut/ntuples/SinglePhoton2017C_WGamma_select.root");

  factory->AddSignalTree((TTree*)infile_s->Get("Events"));
  factory->AddBackgroundTree((TTree*)infile_b->Get("Events"));
  factory->AddVariable("ak8puppijet_tau21", 'F');
  factory->AddVariable("ak8puppijet_massdiff", 'F');
  factory->AddVariable("sys_costhetastar", 'F');
  factory->AddVariable("sys_ptoverm", 'F');

  factory->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  //factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=-1" );

  factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outFile->Close();
  delete factory;
}
