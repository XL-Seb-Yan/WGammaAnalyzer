
void makeRooMultiPdfWorkspaceSigW(){

  gROOT->SetBatch(1);
  // Load the combine Library 
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  RooRealVar *x = new RooRealVar("m","m",600,7500,"");
  RooPlot *frame = x->frame();
  
  // Open anchor workspace wide
  int sigmass_W[15]={700,800,900,1000,1200,1400,1600,1800,2000,2400,2600,2800,3500,5000,6000};
  TFile *f_1 = NULL;
  TFile *f_2 = NULL;
  for(int i = 0; i<14; i++){
    TString sig_type_1 = std::to_string(sigmass_W[i])+"W";
    TString sig_type_2 = std::to_string(sigmass_W[i+1])+"W";
	f_1 = TFile::Open(sig_type_1+"-shapes-Unbinned-CB2Gaus.root");
    f_2 = TFile::Open(sig_type_2+"-shapes-Unbinned-CB2Gaus.root");
    cout<<"Processing..."<<endl;
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
    RooRealVar* Gaus_sigma_1 = NULL;
    RooRealVar* Gaus_sigma_2 = NULL;
    RooGaussian* Gaus_model_1 = NULL;
    RooGaussian* Gaus_model_2 = NULL;
    RooRealVar* frac_1 = NULL;
    RooRealVar* frac_2 = NULL;
    RooAddPdf* com_model_1 = NULL;
    RooAddPdf* com_model = NULL;

    //interpolated signal shapes
    int step = 10;
    int npoints = (sigmass_W[i+1] - sigmass_W[i]) / step;
    float CB_mean_low = w_1->var("CB_mean")->getValV();
    float CB_sigma_low = w_1->var("CB_sigma")->getValV();
    float CB_alpha_low = w_1->var("CB_alpha")->getValV();
    float CB_n_low = w_1->var("CB_n")->getValV();
    float Gaus_sigma_1_low = w_1->var("Gaus_sigma_1")->getValV();
    float Gaus_sigma_2_low = w_1->var("Gaus_sigma_2")->getValV();
    float frac_1_low = w_1->var("frac1")->getValV();
    float frac_2_low = w_1->var("frac2")->getValV();
    float CB_mean_high = w_2->var("CB_mean")->getValV();
    float CB_sigma_high = w_2->var("CB_sigma")->getValV();
    float CB_alpha_high = w_2->var("CB_alpha")->getValV();
    float CB_n_high = w_2->var("CB_n")->getValV();
    float Gaus_sigma_1_high = w_2->var("Gaus_sigma_1")->getValV();
    float Gaus_sigma_2_high = w_2->var("Gaus_sigma_2")->getValV();
    float frac_1_high = w_2->var("frac1")->getValV();
    float frac_2_high = w_2->var("frac2")->getValV();
    cout<<"========================================================================================="<<endl;
    cout<<"Anchor variables low:  "<<CB_mean_low<<" "<<CB_sigma_low<<" "<<CB_alpha_low<<" "<<CB_n_low<<" "<<Gaus_sigma_1_low<<" "<<Gaus_sigma_2_low<<" "<<frac_1_low<<" "<<frac_2_low<<endl;
    cout<<"Anchor variables high:  "<<CB_mean_high<<" "<<CB_sigma_high<<" "<<CB_alpha_high<<" "<<CB_n_high<<" "<<Gaus_sigma_1_high<<" "<<Gaus_sigma_2_high<<" "<<frac_1_high<<" "<<frac_2_high<<endl;
    cout<<"========================================================================================="<<endl;
    for(int ip=0; ip<npoints+1; ip++){
        m = new RooRealVar("m","m",(sigmass_W[i]+step*ip)*0.75,(sigmass_W[i]+step*ip)*1.25,"");
        float CB_mean_p = CB_mean_low + ip * (CB_mean_high - CB_mean_low) / npoints;
        float CB_sigma_p = CB_sigma_low + ip * (CB_sigma_high - CB_sigma_low) / npoints;
        float CB_alpha_p = CB_alpha_low + ip * (CB_alpha_high - CB_alpha_low) / npoints;
        float CB_n_p = CB_n_low + ip * (CB_n_high - CB_n_low) / npoints;
        float Gaus_sigma_1_p = Gaus_sigma_1_low + ip * (Gaus_sigma_1_high - Gaus_sigma_1_low) / npoints;
        float Gaus_sigma_2_p = Gaus_sigma_2_low + ip * (Gaus_sigma_2_high - Gaus_sigma_2_low) / npoints;
        float frac_1_p = frac_1_low + ip * (frac_1_high - frac_1_low) / npoints;
        float frac_2_p = frac_2_low + ip * (frac_2_high - frac_2_low) / npoints;
        cout<<"-------------------------------------------------------------------------------------"<<endl;
        cout<<"Interpolated variables:  "<<CB_mean_p<<" "<<CB_sigma_p<<" "<<CB_alpha_p<<" "<<CB_n_p<<" "<<Gaus_sigma_1_p<<" "<<Gaus_sigma_2_p<<" "<<frac_1_p<<" "<<frac_2_p<<endl;
        cout<<"-------------------------------------------------------------------------------------"<<endl;
        CB_mean = new RooRealVar("CB_mean","CB_mean",CB_mean_p);
        CB_sigma = new RooRealVar("CB_sigma","CB_sigma",CB_sigma_p);
        CB_alpha = new RooRealVar("CB_alpha","CB_alpha",CB_alpha_p);
        CB_n = new RooRealVar("CB_n","CB_n",CB_n_p);
        CB_model = new RooCBShape("CB","Cystal Ball Function",*m,*CB_mean,*CB_sigma,*CB_alpha,*CB_n);
        Gaus_sigma_1 = new RooRealVar("Gaus_sigma_1","Gaus_sigma_1",Gaus_sigma_1_p);
        Gaus_model_1 = new RooGaussian("Gaussian_1","Gaussian Function",*m,*CB_mean,*Gaus_sigma_1);
        Gaus_sigma_2 = new RooRealVar("Gaus_sigma_2","Gaus_sigma_2",Gaus_sigma_2_p);
        Gaus_model_2 = new RooGaussian("Gaussian_2","Gaussian Function",*m,*CB_mean,*Gaus_sigma_2);
        frac_1 = new RooRealVar("frac1","frac1",frac_1_p);
        com_model_1 = new RooAddPdf("CBGaus","CBGaus",RooArgList(*CB_model,*Gaus_model_1),RooArgList(*frac_1));
        frac_2 = new RooRealVar("frac2","frac2",frac_2_p);
        com_model = new RooAddPdf("CB2Gaus","CB2Gaus",RooArgList(*com_model_1,*Gaus_model_2),RooArgList(*frac_2));
        
        TString sig_type_p = std::to_string(sigmass_W[i] + ip * step)+"W";
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
        delete Gaus_sigma_1;
        delete Gaus_sigma_2;
        delete Gaus_model_1;
        delete Gaus_model_2;
        delete frac_1;
        delete frac_2;
        delete com_model_1;
        delete com_model;
    }
  }
}