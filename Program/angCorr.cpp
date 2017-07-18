//
// July 17, 2017
// Connor Natzke
// ROOT v6.10.02
// GRSISort v3.1.3.2
//

#include "angCorr.h"

angCorr::angCorr() {
    Double_t firstPeak = 1172;                  // Gate Peak
    Double_t secondPeak = 1332;                 // Coincident Peak
    TString isotope = "^{60}Co";                // Parent Isotope
    Int_t gndState = 0;                         // Parent Groundstate

    TH3D* xyzProj = new TH3D();
    GetData(xyzProj);                           // Loads data into ROOT

    TH1D* histVec[52];
    IndexProject(xyzProj, firstPeak, histVec);  // Projects coincidence matrix

    Double_t counts[52];
    FitPeak(histVec, secondPeak, counts);       // Fits peak 

    Double_t errors[52];
    Error(counts, errors);                      // Calculates Errors

    Double_t normCounts[52];
    Normalize(counts, normCounts);              // Normalizes Counts

    Double_t normErrors[52];
    Normalize(errors, normErrors);


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

    TFile* isoData = new TFile("/home/data/cnatzke/SimulationResults/Converted20mm.root");
    isoData->cd("GriffinND");
    THnSparse* dataFull = (THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
    cout << "File Loaded" << endl;
    *xyzProj = *dataFull->Projection(0,1,2);
    cout << "3D Histogram Generated" << endl;

    delete dataFull;
//    delete isoData;       //This seems to crash everything, no idea why

    return;
}

///////////////////////////////////////////////////////////////
// Isolates indivudual indices and projects coincidence spectra
///////////////////////////////////////////////////////////////

void angCorr::IndexProject(TH3D* xyzProj, Double_t firstPeak, TH1D* histVec[52]) {

    cout << "Creating projection histograms" << endl;

    Double_t gateLow;
    Double_t gateHigh;

    gateLow = firstPeak-2.;         // Lower bound of gate peak
    gateHigh = firstPeak+2.;        // Upper bound of gate peak

    GCanvas* c1 = new GCanvas();
    gStyle->SetOptStat(0);          // Turn off histogram statistics
    c1->SetLogy();                   // Set Log scale on y-axis

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

            xyzProj->GetZaxis()->SetRangeUser(indLow, indHigh);
            TH2D* ggProj = (TH2D*) xyzProj->Project3D("xy");
            histVec[i] = ggProj->ProjectionY(hName, gateLow, gateHigh);
            histVec[i]->SetTitle(hTitle);
            histVec[i]->SetMinimum(1.0);
            histVec[i]->SetMaximum(20000.);
            histVec[i]->Draw();

            c1->Update();
            TString path("/home/data/cnatzke/SimulationResults/Data/RadiusVariation/AngularHistograms/NonFitted/");
            path += hName;
            path += ".png";
            cout << path.Data() << endl;
            c1->SaveAs(path);

            }
        else{}
        }

    return;
}


///////////////////////////////////////////////////////////////
// Fits Coincedence Spectra and Returns Counts
///////////////////////////////////////////////////////////////

void angCorr::FitPeak(TH1D* histVec[52], Double_t secondPeak, Double_t counts[52]) {

    Double_t enLow = secondPeak-10.;
    Double_t enHigh = secondPeak+10.;

    GCanvas* c2 = new GCanvas();
    gStyle->SetOptStat(0);                                  // Turn off histogram statistics
    c2->SetLogy();                                          // Set Log scale on y-axis

    TPeak* peak = new TPeak(secondPeak, enLow, enHigh);     // Defining new TPeak for Fitting
    peak->SetParameter("centroid", secondPeak);             // Setting intial guesses for fitting
    peak->SetParameter("sigma", 1.);
    peak->SetParameter("Height", 6.885e+03);

    TString hPath("/home/data/cnatzke/SimulationResults/Data/RadiusVariation/RootFiles/ggac");
    hPath += 50;
    hPath += "mmHistos.root";
    cout << hPath.Data() << endl;
    TFile histRoot(hPath, "RECREATE");                      // ROOT file with fitted histograms

    for(Int_t i=0;i<52;i++){

        if(i!=0){
            histVec[i]->GetXaxis()->SetRangeUser(enLow-5., enHigh+5.);
            histVec[i]->SetMinimum(1.0);
            histVec[i]->SetMaximum(20000);

            peak->Fit(histVec[i]);                          // Fit coincidence spectra
            counts[i] = peak->GetArea();

            TString hName("index");
            hName += i;
            cout << hName.Data() << endl;

            c2->Update();
            TString savePath("/home/data/cnatzke/SimulationResults/Data/RadiusVariation/AngularHistograms/Fitted/");
            savePath += hName;
            savePath += ".png";
            cout << savePath.Data() << endl;
            c2->SaveAs(savePath);                           // Saves fitted histos for review

            histRoot.WriteTObject(histVec[i]);              // Writes ROOT file
            }
        else{}
    }
    return;
}

///////////////////////////////////////////////////////////////
// Normalizes input data 
///////////////////////////////////////////////////////////////

void angCorr::Normalize(Double_t array[52], Double_t norm[52]) {

    for(Int_t i=0;i<52;i++){

        if(i==1||i==6||i==8||i==25||i==26||i==43||i==45||i==50){
            norm[i]=(array[i])/(128);
        }
        else if(i==5||i==10||i==12||i==24||i==27||i==39||i==41||i==46){
            norm[i]=(array[i])/(48);
        }
        else if(i==7||i==9||i==11||i==14||i==18||i==20||i==31||i==33||i==37||i==40||i==42||i==44){
            norm[i]=(array[i])/(96);
        }
        else{
            norm[i]=(array[i])/(64);
        }
    }
    return;

}

///////////////////////////////////////////////////////////////
// Calculates Error Bars
///////////////////////////////////////////////////////////////

void angCorr::Error(Double_t counts[52], Double_t errors[52]) {

    for(Int_t i=0;i<52;i++){

        errors[i]=TMath::Sqrt(counts[i]);
    }
    return;
}



