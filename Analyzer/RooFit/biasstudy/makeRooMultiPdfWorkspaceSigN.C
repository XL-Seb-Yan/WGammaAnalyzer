
void makeRooMultiPdfWorkspaceSigN(){

  gROOT->SetBatch(1);
  // Load the combine Library 
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  RooRealVar *x = new RooRealVar("m","m",600,5000,"");
  RooPlot *frame = x->frame();

  // Open anchor workspace narrow
  int sigmass_N[16] = {700,800,900,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3500,4000,5000};
  TFile *f_1 = NULL;
  TFile *f_2 = NULL;
  for(int i = 0; i<15; i++){
    TString sig_type_1 = std::to_string(sigmass_N[i])+"N";
    TString sig_type_2 = std::to_string(sigmass_N[i+1])+"N";
	f_1 = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/Analyzer/CMSSW_9_4_13/src/WGammaAnalyzer/Analyzer/RooFit/biasstudy/"+sig_type_1+"-shapes-Unbinned-CBGaus.root");
    f_2 = TFile::Open("/afs/cern.ch/work/x/xuyan/work5/PROD17/Analyzer/CMSSW_9_4_13/src/WGammaAnalyzer/Analyzer/RooFit/biasstudy/"+sig_type_2+"-shapes-Unbinned-CBGaus.root");
    if(f_1 == NULL) continue;
    if(f_2 == NULL) continue;
    cout<<"Processing..."<<endl;
    RooWorkspace *w_1 = (RooWorkspace*)f_1->Get("w");
    RooWorkspace *w_2 = (RooWorkspace*)f_2->Get("w");
    RooDataSet *data_1 = (RooDataSet*)w_1->data("signal_MC");
    RooDataSet *data_2 = (RooDataSet*)w_2->data("signal_MC");

    // Open output workspace for anchor pdfs
    RooRealVar* m = NULL;
    RooRealVar* CB_mean = NULL;
    RooRealVar* CB_sigma = NULL;
    RooRealVar* CB_alpha = NULL;
    RooRealVar* CB_n = NULL;
    RooCBShape* CB_model = NULL;
    RooRealVar* Gaus_mean = NULL;
    RooRealVar* Gaus_sigma = NULL;
    RooGaussian* Gaus_model = NULL;
    RooRealVar* frac = NULL;
    RooAddPdf* com_model = NULL;
    // TFile *fout_1 = new TFile("signal_pdfs_"+sig_type_1+".root","RECREATE");
    // RooWorkspace wout_1("signals","signals");
    // m = new RooRealVar("m","m",00,5000,"");
    // CB_mean = new RooRealVar("CB_mean","CB_mean",w_1->var("CB_mean")->getValV());
    // CB_sigma = new RooRealVar("CB_sigma","CB_sigma",w_1->var("CB_sigma")->getValV());
    // CB_alpha = new RooRealVar("CB_alpha","CB_alpha",w_1->var("CB_alpha")->getValV());
    // CB_n = new RooRealVar("CB_n","CB_n",w_1->var("CB_n")->getValV());
    // CB_model = new RooCBShape("CB","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
    // Gaus_mean = new RooRealVar("Gaus_mean","Gaus_mean",w_1->var("Gaus_mean")->getValV());
    // Gaus_sigma = new RooRealVar("Gaus_sigma","Gaus_sigma",w_1->var("Gaus_sigma")->getValV());
    // Gaus_model = new RooGaussian("Gaussian","Gaussian Function",*m,*Gaus_mean,*Gaus_sigma);
    // frac = new RooRealVar("frac","frac",w_1->var("frac")->getValV());
    // com_model = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model),RooArgList(*frac));
    // wout_1.import(*data_1);
	// wout_1.import(*com_model);
    // wout_1.Write();
    // fout_1->Close();
    // delete m;
    // delete CB_mean;
    // delete CB_sigma;
    // delete CB_alpha;
    // delete CB_n;
    // delete CB_model;
    // delete Gaus_mean;
    // delete Gaus_sigma;
    // delete Gaus_model;
    // delete frac;
    // delete com_model;
    
    // TFile *fout_2 = new TFile("signal_pdfs_"+sig_type_2+".root","RECREATE");
    // RooWorkspace wout_2("signals","signals");
    // m = new RooRealVar("m","m",00,5000,"");
    // CB_mean = new RooRealVar("CB_mean","CB_mean",w_2->var("CB_mean")->getValV());
    // CB_sigma = new RooRealVar("CB_sigma","CB_sigma",w_2->var("CB_sigma")->getValV());
    // CB_alpha = new RooRealVar("CB_alpha","CB_alpha",w_2->var("CB_alpha")->getValV());
    // CB_n = new RooRealVar("CB_n","CB_n",w_2->var("CB_n")->getValV());
    // CB_model = new RooCBShape("CB","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
    // Gaus_mean = new RooRealVar("Gaus_mean","Gaus_mean",w_2->var("Gaus_mean")->getValV());
    // Gaus_sigma = new RooRealVar("Gaus_sigma","Gaus_sigma",w_2->var("Gaus_sigma")->getValV());
    // Gaus_model = new RooGaussian("Gaussian","Gaussian Function",*m,*Gaus_mean,*Gaus_sigma);
    // frac = new RooRealVar("frac","frac",w_2->var("frac")->getValV());
    // com_model = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model),RooArgList(*frac));
    // wout_2.import(*data_2);
	// wout_2.import(*com_model);
    // wout_2.Write();
    // fout_2->Close();
    // delete m;
    // delete CB_mean;
    // delete CB_sigma;
    // delete CB_alpha;
    // delete CB_n;
    // delete CB_model;
    // delete Gaus_mean;
    // delete Gaus_sigma;
    // delete Gaus_model;
    // delete frac;
    // delete com_model;
    
    //interpolated signal shapes
    int step = 50;
    int npoints = (sigmass_N[i+1] - sigmass_N[i]) / step;
    float CB_mean_low = w_1->var("CB_mean")->getValV();
    float CB_sigma_low = w_1->var("CB_sigma")->getValV();
    float CB_alpha_low = w_1->var("CB_alpha")->getValV();
    float CB_n_low = w_1->var("CB_n")->getValV();
    float Gaus_mean_low = w_1->var("Gaus_mean")->getValV();
    float Gaus_sigma_low = w_1->var("Gaus_sigma")->getValV();
    float frac_low = w_1->var("frac")->getValV();
    float CB_mean_high = w_2->var("CB_mean")->getValV();
    float CB_sigma_high = w_2->var("CB_sigma")->getValV();
    float CB_alpha_high = w_2->var("CB_alpha")->getValV();
    float CB_n_high = w_2->var("CB_n")->getValV();
    float Gaus_mean_high = w_2->var("Gaus_mean")->getValV();
    float Gaus_sigma_high = w_2->var("Gaus_sigma")->getValV();
    float frac_high = w_2->var("frac")->getValV();
    cout<<"========================================================================================="<<endl;
    cout<<"Anchor variables low:  "<<CB_mean_low<<" "<<CB_sigma_low<<" "<<CB_alpha_low<<" "<<CB_n_low<<" "<<Gaus_mean_low<<" "<<Gaus_sigma_low<<" "<<frac_low<<endl;
    cout<<"Anchor variables high:  "<<CB_mean_high<<" "<<CB_sigma_high<<" "<<CB_alpha_high<<" "<<CB_n_high<<" "<<Gaus_mean_high<<" "<<Gaus_sigma_high<<" "<<frac_high<<endl;
    cout<<"========================================================================================="<<endl;
    for(int ip=0; ip<npoints+1; ip++){
        m = new RooRealVar("m","m",600,6000,"");
        float CB_mean_p = CB_mean_low + ip * (CB_mean_high - CB_mean_low) / npoints;
        float CB_sigma_p = CB_sigma_low + ip * (CB_sigma_high - CB_sigma_low) / npoints;
        float CB_alpha_p = CB_alpha_low + ip * (CB_alpha_high - CB_alpha_low) / npoints;
        float CB_n_p = CB_n_low + ip * (CB_n_high - CB_n_low) / npoints;
        float Gaus_mean_p = Gaus_mean_low + ip * (Gaus_mean_high - Gaus_mean_low) / npoints;
        float Gaus_sigma_p = Gaus_sigma_low + ip * (Gaus_sigma_high - Gaus_sigma_low) / npoints;
        float frac_p = frac_low + ip * (frac_high - frac_low) / npoints;
        cout<<"-------------------------------------------------------------------------------------"<<endl;
        cout<<"Interpolated variables:  "<<CB_mean_p<<" "<<CB_sigma_p<<" "<<CB_alpha_p<<" "<<CB_n_p<<" "<<Gaus_mean_p<<" "<<Gaus_sigma_p<<" "<<frac_p<<endl;
        cout<<"-------------------------------------------------------------------------------------"<<endl;
        CB_mean = new RooRealVar("CB_mean","CB_mean",CB_mean_p);
        CB_sigma = new RooRealVar("CB_sigma","CB_sigma",CB_sigma_p);
        CB_alpha = new RooRealVar("CB_alpha","CB_alpha",CB_alpha_p);
        CB_n = new RooRealVar("CB_n","CB_n",CB_n_p);
        CB_model = new RooCBShape("CB","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
        Gaus_mean = new RooRealVar("Gaus_mean","Gaus_mean",Gaus_mean_p);
        Gaus_sigma = new RooRealVar("Gaus_sigma","Gaus_sigma",Gaus_sigma_p);
        Gaus_model = new RooGaussian("Gaussian","Gaussian Function",*m,*Gaus_mean,*Gaus_sigma);
        frac = new RooRealVar("frac","frac",frac_p);
        com_model = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model),RooArgList(*frac));
        
        TString sig_type_p = std::to_string(sigmass_N[i] + ip * step)+"N";
        TFile *fout_p = new TFile("signal_pdfs_"+sig_type_p+".root","RECREATE");
        RooWorkspace wout_p("signals","signals");
        wout_p.import(*com_model);
        wout_p.Write();
        fout_p->Close();
        delete m;
        delete CB_mean;
        delete CB_sigma;
        delete CB_alpha;
        delete CB_n;
        delete CB_model;
        delete Gaus_mean;
        delete Gaus_sigma;
        delete Gaus_model;
        delete frac;
        delete com_model;
    }
  }
}