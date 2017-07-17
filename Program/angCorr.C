//
// July 17, 2017
// Connor Natzke
// ROOT v6.10.02
// GRSISort v3.1.3.2
//

#include "angCorr.h"

angCorr::angCorr() {
    Double_t firstPeak = 1172;      // Gate Peak
    Double_t secondPeak = 1332;     // Coincident Peak
    TString isotope = "^{60}Co";    // Parent Isotope
    Int_t gndState = 0;             // Parent Groundstate

    TH3D* xyzProj = new TH3D();
    GetData(xyzProj);

    vector<TH1D*> histVec(52);
    IndexProject(xyzProj, firstPeak, histVec);

}

///////////////////////////////////////////////////////////////
// Deconstructor
///////////////////////////////////////////////////////////////

angCorr::~angCorr() {

}

///////////////////////////////////////////////////////////////
// Reads input ROOT file and builds 3D Coincidence Matrix
///////////////////////////////////////////////////////////////

void angCorr::GetData(TH3D* xyzProj){

    TFile* isoData = new TFile("/home/data/cnatzke/SimulationResults/Converted50mm.root");
    isoData->cd("GriffinND");
    THnSparse* dataFull = (THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
    cout << "File Loaded" << endl;
    xyzProj = (TH3D*) dataFull->Projection(0,1,2);
    cout << "3D Histogram Generated" << endl;

    delete isoData;
    delete dataFull;

    return;
}

///////////////////////////////////////////////////////////////
// Isolates indivudual indices and projects coincidence spectra
///////////////////////////////////////////////////////////////

void angCorr::IndexProject(TH3D* xyzProj, Double_t firstPeak, vector<TH1D*> histVec) {

    Double_t gateLow;
    Double_t gateHigh;

    gateLow = firstPeak-2.;         // Lower bound of gate peak
    gateHigh = firstPeak+2.;        // Upper bound of gate peak

    GCanvas c1("c1","c1",600,600);
    gStyle->SetOptStat(0);          // Turn off histogram statistics
    c1.SetLogy();                   // Set Log scale on y-axis

    for(Int_t i=0;i<52;i++){                                        // Loop through indices

        if(i!=0){
            Double_t indLow=i;
            Double_t indHigh=indLow+1;

            TString hTitle("Index ");                               // Sets histogram titles
            hTitle += indLow;
            hTitle += " Projection Histogram;Energy (keV);Index ";
            cout << hTitle.Data() << endl;

            TString hName("index");
            hName += i;
            cout << hName.Data() << endl;

            delete histVec[i];                                      // Clear memory
            xyzProj->GetZaxis()->SetRangeUser(indLow, indHigh);
            TH2D* ggProj = (TH2D*) xyzProj->Project3D("xy");
            histVec[i] = ggProj->ProjectionY(hName, gateLow, gateHigh);
            histVec[i]->SetTitle(hTitle);
            histVec[i]->SetMinimum(1.0);
            histVec[i]->SetMaximum(1.0);
            histVec[i]->Draw();

            c1.Update();
            TString path("/home/data/cnatzke/SimulationResults/Data/RadiusVariation/AngularHistograms/NonFitted/");
            path += hName;
            path += ".png";
            cout << path.Data() << endl;
            c1.SaveAs(path);
            }
        else{}
        }
        return;
}







