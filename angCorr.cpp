/**
    @file angCorr.cpp

    @author Connor Bray
    @author Connor Natzke

    @date 25 Apr. 2018

    @brief A flexible GRSISort macro for extracting angular correlations from (GRIFFIN) data.

    ### Dependencies:

    * GRSISort (Commit 6423d91 or newer).

    * Root 6 (With minuit2 enabled).

    ### Usage instructions:

    After you have converted the g4out.root file using NTuple, launch GRSISort. Load the commands using `.L angCorr.cpp+` (or the full path if it is in another directory).

    Then you can use the angCorr() command to calculate angular correlations. See documentation below for specific syntax, but in general you can either provide a set of energies to find all angular correlations between, or tell the command to find the n largest peaks and calculate angular correlations between those.

    ### TODO:

    * Further testing and QA

    * Determine if multithreading is worthwhile, and, if so, if there is a way to do so without massive memory usage
*/

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TPeak.h"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

using namespace std;

/// Gate energy tolerance (keV)
const Double_t gateTolerance = 6;

/// Second gamma-ray energy tolerance (keV)
const Double_t energyTolerance = 6;

/// Save extra root files flag
int save=0;

/// Base address for file IO
string baseAddr="";

/// Output file stream
ofstream out;

/// Extraction status indicator
int status=0;

/**
@brief Cuts the 3D Histogram at each angular index and projects/plots a histogram for each index

   ## Returns a vector of histograms representing the histogram cut for each index

   ### Arguments:
   * `TH3D * xyzProj` - 3D histogram projection from converted.root
   * `Double_t firstPeak` - Energy value in keV to cut at
*/
vector<TH1D*> AngCut(TH3D * xyzProj, Double_t firstPeak){
    vector<TH1D*> histVec(52);
    for(Int_t i=1;i<52;i++){
        delete histVec[i]; // Clear memory space
        xyzProj->GetZaxis()->SetRangeUser(i,i+1); // First cut
        TH2D * ggProj=(TH2D*) xyzProj->Project3D("xy");
        ggProj->GetXaxis()->SetRangeUser(firstPeak-gateTolerance,firstPeak+gateTolerance); //Second Cut
        histVec[i]=ggProj->ProjectionY((TString)("index"+to_string(i)),firstPeak-gateTolerance,firstPeak+gateTolerance);
        histVec[i]->SetTitle((TString)("Projection Histogram for Index "+to_string(i)));
        histVec[i]->SetMinimum(-1.0);
        histVec[i]->SetMaximum(1000000000);
    }
    return histVec;
}

/**
@brief Extracts counts for each index

   ## Takes in vector of histograms,fits a skewed Gaussian to the photo-peak of interest, integrates fit to find counts, returns counts with background subtracted

   ### Arguments:
   * `vector<Double_t> &counts` - Return value, vector of 52 length representing counts at each index
   * `vector<TH1D*> histVec` - Histograms from AngCut to fit to
   * `Double_t firstPeak, Double_t secPeak` - First and second peak energies. firstPeak energy is only used for printing, the peak is fit around secPeak
*/
void FitHisto(vector<Double_t> &counts, vector<TH1D*> histVec, Double_t firstPeak, Double_t secPeak){
    // Using the TPEAK method for fitting
    TPeak* peak=new TPeak(secPeak, secPeak-energyTolerance, secPeak+energyTolerance);
    peak->SetParameter("centroid",secPeak);
    peak->SetParameter("sigma",1);
    peak->SetParameter("Height",1e+04);

    for(Int_t i=1;i<52;i++){
        histVec[i]->GetXaxis()->SetRangeUser(secPeak-energyTolerance-5, secPeak+energyTolerance+5); //Setting Range of Histogram to allow data fitting
        if(!peak->Fit(histVec[i],"ERMQS")){
            cout << "First Peak: " << firstPeak << ", Second Peak: " << secPeak << ", No data error." << endl;
            status=-1;
            if(out.is_open()) out << "First Peak: " << firstPeak << ", Second Peak: " << secPeak << ", No data error." << endl;
            return;
        }
        counts[i]=peak->GetArea(); // Integrating peak to find counts
        histVec[i]->Draw("SAME");
    }

    if(save){
        TString hpath(baseAddr+"ggac"+to_string(firstPeak)+"-"+to_string(secPeak)+".root");
        cout<<hpath.Data() << endl;
        TFile histroot(hpath,"RECREATE");
        for(int i=0;i<52;i++) histroot.WriteTObject(histVec[i]);
    }
    delete peak;
}


/**
@brief Normalizes the counts for the photopeak

   ### Arguments:
   * `vector<Double_t> &counts` - Return value and argument representing the counts at each index
   * `vector<Double_t> &counts` - Return value representing error in counts at each index
*/
void WeightAdjust(vector<Double_t> &counts, vector<Double_t> &error){
    for(Int_t i=0;i<52;i++){
        Int_t factor;
        if(i==1||i==6||i==8||i==25||i==26||i==43||i==45||i==50) factor = 128;
        else if(i==5||i==10||i==12||i==24||i==27||i==39||i==41||i==46) factor = 48;
        else if(i==7||i==9||i==11||i==14||i==18||i==20||i==31||i==33||i==37||i==40||i==42||i==44) factor = 96;
        else factor = 64;
        error[i]=TMath::Sqrt(counts[i])/factor;
        counts[i]/=factor;
    }
}

/**
@brief Plots and Fits Angular Correlations

   ### Arguments:
   * `vector<Double_t> wcounts` - Vector of 52 length representing counts at each index
   * `vector<Double_t> countErrors` - Vector of 52 length representing error in counts at each index
   * `Double_t firstPeak, Double_t secPeak` - First and second peak energies
*/
void AngularCorrHisto(vector<Double_t> wcounts, vector<Double_t> countErrors, Double_t firstPeak, Double_t secPeak){
    //Angular Indices
    Double_t index[]={0.0,18.79097,25.60153, 26.69036,31.94623,33.65414, 44.36426, 46.79372, 48.57554, 49.79788, 53.83362, 60.15106, 62.70487, 63.08604, 65.01569, 66.46082, 67.45617, 69.86404, 70.86009, 73.08384, 76.38138, 78.66898, 83.04252, 86.22840, 86.23761, 88.47356, 91.52644, 93.76239, 93.77160, 96.95749, 101.33102, 103.61822, 106.91616, 109.13991, 110.13596, 112.54383, 113.53918, 114.98431, 116.91396, 117.29513, 119.84894, 126.16638, 130.20212, 131.42446, 133.20628, 135.63574, 146.34586, 148.05377, 153.30964, 154.39847, 161.21315,180.0};
    vector<Double_t> indexVec(52);
    vector<Double_t> indexErrors(52); // Initialize vector for counts

    for(Int_t i=0;i<52;i++){
        indexVec[i]=TMath::Cos(index[i]*TMath::Pi()/(180));
        indexErrors[i]=0;
    }

    TGraphErrors * g=new TGraphErrors(52,&indexVec[0],&wcounts[0],&indexErrors[0],&countErrors[0]);
    TString c_name(to_string(firstPeak)+" to "+to_string(secPeak)+" keV decay;Cos(#theta);Normalized Counts");
    TF1* efit=new TF1("efit","(1+[0]*(3*x*x-1)/2+[1]*(35*x*x*x*x-30*x*x+3)/8)*[2]",-1,1); // Fit Legendre polynomial
    efit->SetParNames("a2","a4","scale");
    efit->SetParameters(1,1,1000);

    gStyle->SetOptFit(2);
    g->Fit("efit","ERQM");

    Double_t chi2,a22,a44,scale,NDF;
    Double_t a22e,a44e,scalee;
    chi2=efit->GetChisquare();
    a22=efit->GetParameter("a2");
    a44=efit->GetParameter("a4");
    scale=efit->GetParameter("scale");
    a22e=efit->GetParError(0);
    a44e=efit->GetParError(1);
    scalee=efit->GetParError(2);
    NDF=efit->GetNDF();
    cout << "First Peak: " << firstPeak << ", Second Peak: " << secPeak << ", χ^2: " << chi2 << "/" << NDF << ", a22: " << a22 << "±" << a22e << ", A44: " << a44 << "±" << a44e << ", scale: " << scale << "±" << scalee << endl;
    if(out.is_open()) out << "First Peak: " << firstPeak << ", Second Peak: " << secPeak << ", χ^2: " << chi2 << "/" << NDF << ", a22: " << a22 << "±" << a22e << ", A44: " << a44 << "±" << a44e << ", scale: " << scale << "±" << scalee << endl;
    // Finds min and max counts for graph scaling
    Double_t max=TMath::MaxElement(52,g->GetY());
    Double_t min=efit->GetMinimum();

    g->SetTitle(c_name);
    g->SetMaximum(max+70);
    g->SetMinimum(min-70);
    g->SetMarkerStyle(7);
    g->GetYaxis()->SetTitleOffset(1.4);
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->Draw("APE");

    TString rtpath(baseAddr+"ggac_histo"+to_string(firstPeak)+"-"+to_string(secPeak)+".root");//Root File  name
    cout<<rtpath.Data() << endl;
    TFile ggacroot(rtpath,"RECREATE"); // Create .root file
    ggacroot.WriteTObject(g);

    delete g;
    delete efit;
}

/**
@brief Extracts all angular correlations from list of peaks in GRIFFIN converted data

   ### Arguments:
   * `vector<Double_t> peaks` - Vector of energy levels (in keV) to fit peaks between. All sets of peaks are analyzed for angular correlations.
   * `bool addback=0` - Defaults to false, which uses the crystal specific data for the angular correlations (better angular resolution, worse count resolution). If set to true, then it uses the addback spectra (worse angular resolution, better count resolution)
   * `string outputDir=baseAddr` - Base address to read data from and write output files to
   * `bool saveExtra=save` - Defaults to false, but if set to true then it saves extra histograms from intermediary steps
*/
void angCorr(vector<Double_t> peaks, bool addback=0, string outputDir=baseAddr, bool saveExtra=save){
    save=saveExtra; // Flag to save extra files (such as intermediary histograms or pdfs)
    baseAddr=outputDir;
    gROOT->SetBatch(kTRUE);

    TFile* isoData = new TFile((baseAddr+"converted.root").c_str());
    isoData->cd("GriffinND");
    THnSparse * dataFull;
    if(!addback) dataFull=(THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
    else dataFull=(THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_addback_sparse");
    out.open("fit.txt",std::ios_base::app);
    TH3D * xyzProj=((TH3D*) dataFull->Projection(0,1,2)); //3D Histogram;
    gStyle->SetOptStat(0);

    for(unsigned int i=0;i<peaks.size()-1;i++){ // Because of the way root and GRSISort work, multithreading this has been ineffective and lead to huge ram usage.
        vector<TH1D*> histVec=AngCut(xyzProj,peaks[i]);
        for(unsigned int j=i;++j<peaks.size();){
            status=0;
            vector<Double_t> error(52);
            vector<Double_t> counts(52);
            if(!status) FitHisto(counts,histVec,peaks[i],peaks[j]);
            if(!status) WeightAdjust(counts,error); // Only use if you are using the full GRIFFIN Detector array, weights are specific
            if(!status) AngularCorrHisto(counts,error,peaks[i],peaks[j]);
        }
    }

    delete isoData;
    out.close();
}

/**
@brief Finds n highest peaks in GRIFFIN converted data, then extracts all angular correlations between them

   ### Arguments:
   * `int maxPeaks=10` - Max number of peaks to find. Defaults to 10.
   * `bool addback=0` - Defaults to false, which uses the crystal specific data for the angular correlations (better angular resolution, worse count resolution). If set to true, then it uses the addback spectra (worse angular resolution, better count resolution)
   * `string outputDir=baseAddr` - Base address to read data from and write output files to
   * `bool saveExtra=save` - Defaults to false, but if set to true then it saves extra histograms from intermediary steps
*/
void angCorr(int maxPeaks=10, bool addback=0, string outputDir=baseAddr, bool saveExtra=save){
    TFile* converted = new TFile((baseAddr+"converted.root").c_str());
    converted->cd("Griffin1D");
    TH1F * peakPlot=(TH1F*) gDirectory->Get("griffin_crystal_unsup_edep_sum");
    TSpectrum *s = new TSpectrum(maxPeaks); // Max number of peaks (takes the most significant)
    s->SetAverageWindow(1);
    Int_t nfound = s->Search(peakPlot,1,"nodraw",0.000000001);
    Double_t *p = s->GetPositionX();
    vector<Double_t> peaks(p,p+nfound);
    delete converted;
    delete s;
    for(Double_t e : peaks) cout << e << ",";
    cout << endl << peaks.size() << endl;
    angCorr(peaks, addback, outputDir, saveExtra);
}
