
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"

class angularCorr{

public:
	angularCorr();
	TH3D* GetData();
	vector<TH1D*> AngCut(TH3D*,Int_t,Double_t, Double_t);	
	vector<Double_t> FitHisto(vector<TH1D*>, Int_t, Double_t,Double_t, Double_t);
	vector<Double_t> WeightAdjust(vector<Double_t>,Int_t);
	vector<Double_t> ErrorAdjust(vector<Double_t>,Int_t);
	vector<Double_t> AltError(vector<Double_t>,Int_t);
	void AngularCorrHisto(vector<Double_t>, TString,vector<Double_t>,Int_t gndstate);
	void CountFold(vector<Double_t>,vector<Double_t>,Int_t gndstate);


private:
	
};
