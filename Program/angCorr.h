
#include <iostream>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TPeak.h"
#include "GCanvas.h"
#include "TStyle.h"


using namespace std;

class angCorr{

public:
	angCorr();
	~angCorr();
	void GetData(TH3D*);
    void IndexProject(TH3D*, Double_t, TH1D* [52]);
    void FitPeak(TH1D* [52], Double_t, Double_t [52]);
    void Normalize(Double_t [52], Double_t [52]);
    void Error(Double_t [52], Double_t [52]);
    void CorrelationPlot(Double_t [52], Double_t [52]); 

private:
};
