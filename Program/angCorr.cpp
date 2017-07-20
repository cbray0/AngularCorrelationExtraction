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

/*    Double_t errors[52];
    Error(counts, errors);                      // Calculates Errors

    Double_t normCounts[52];
    Normalize(counts, normCounts);              // Normalizes Counts

    Double_t normErrors[52];
    Normalize(errors, normErrors);
*/
    for(Int_t i=0;i<52;i++){
        cout << *counts << endl;
//        cout << normErrors << endl;
    }
//    CorrelationPlot(normCounts, normErrors);

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
//            counts[i] = peak->GetArea();
            counts[i] = 23.;

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


///////////////////////////////////////////////////////////////
// Plots counts vs angle
///////////////////////////////////////////////////////////////

void angCorr::CorrelationPlot(Double_t normCounts[52], Double_t normError[52]){

    Double_t indexDeg[52] = {0.0, 18.79097, 25.60153, 26.69036, 31.94623, 33.65414, 44.36426, 46.79372, 48.57554, 49.79788, 53.83362, 60.15106, 62.70487, 63.08604, 65.01569, 66.46082, 67.45617, 69.86404, 70.86009, 73.08384, 76.38138, 78.66898, 83.04252, 86.22840, 86.23761, 88.47356, 91.52644, 93.76239, 93.77160, 96.95749, 101.33102, 103.61822, 106.91616, 109.13991, 110.13596, 112.54383, 113.53918, 114.98431, 116.91396, 117.29513, 119.84894, 126.16638, 130.20212, 131.42446, 133.20628, 135.63574, 146.34586, 148.05377, 153.30964, 154.39847, 161.21315, 180.0};

    Double_t indexCos[52];
    Double_t indexError[52];


    for(Int_t i=0;i<52;i++){
        indexCos[i] = TMath::Cos(indexDeg[i]*TMath::Pi()/180);
        indexError[i] = 0.;
    }

    GCanvas* c3 = new GCanvas();
    gStyle->SetOptStat(0);                                  // Turn off histogram statistics

    TGraphErrors* g = new TGraphErrors(52, indexCos, normCounts, indexError, normError);

    TString gName("1172 to 1332 keV Decay (4^{+} #rightarrow 2^{+} #rightarrow 0^{+});Cos(#theta);Normalized Counts");
    
//    TF1* fit = new TF1("fit", "(1+0.5*[0]*(3*x*x-1)+0.125*[1]*(35*x*x*x*x-30*x*x+3))*[2]",-1,1);
//    fit->SetParameter(0,0.01);
//    fit->SetParameter(1,1.0);
//    fit->SetParameter(2,150);
//    fit->SetParName(0,"a2");
//    fit->SetParName(1,"a4");
//    fit->SetParName(2,"scale");

//    g->Fit(fit);
    gStyle->SetOptFit(1);               // Displays fit information`

    Double_t max = TMath::MaxElement(52, g->GetY());
//    Double_t min = fit->GetMinimum();

    // Graph Styling
    g->SetTitle(gName);
    g->SetMaximum(max+20);
//    g->SetMinimum(min-20);
    g->SetMarkerStyle(7);
    g->GetYaxis()->SetTitleOffset(1.4);
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();

    g->Draw("ap");

    // Saving Graph
    c3->Update();
    TString gPath("/home/data/cnatzke/SimulationResults/Data/RadiusVariation/Graphics/corr");
    gPath += 50;
    gPath += "mm.pdf";
    cout << gPath.Data() << endl;
    c3->SaveAs(gPath);

    // ROOT file
    TString rtPath("/home/data/cnatzke/SimulationResults/Data/RadiusVariation/RootFiles/ggac");
    rtPath += 50;
    rtPath += "mm.root";
    cout << rtPath.Data() << endl;
    
    TFile ggacRoot(rtPath, "RECREATE");
    ggacRoot.WriteTObject(g);

}

