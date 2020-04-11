#include <TMath.h>
void plot(){
	TF1* f = new TF1("experf","TMath::Exp([0]*x)*(1.+TMath::Erf((x-[1])/[2]))/2.",0,4000);
	f->SetParameters(-1,500,50);
	f->Draw();
}