//Rewrite of Angular Analysis Code
//Nov 13, 2016
//Should be compatible with addback spectra now


#include <iostream>
#include <vector>
#include "angCorr.H"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"

using namespace std;

angCorr::angCorr(){
    Double_t firstPeak=1172;
    Double_t secPeak=1332;
	TString isotope="^{60}Co";      //Isotope
	Int_t gndstate=0;	        	//Ground State

	Double_t gateLow=firstPeak-2;	//First gate lower bound
	Double_t gateHigh=firstPeak+2;	//First gate upper bound
	Double_t ELow=secPeak-10;       //Second gate lower bound
	Double_t EHigh=secPeak+10;      //Second gate upper bound
	Int_t angSplits=52;	        	//Number of Angular Indices

	TH3D *xyzProj=GetData();
	vector<TH1D*> histVec=AngCut(xyzProj,gateLow, gateHigh);	
	vector<Double_t> rawcounts=FitHisto(histVec, ELow, EHigh, secPeak);
//Only use if you are using the full GRIFFIN Detector array, weights are specific
    vector<Double_t> wcounts=WeightAdjust(rawcounts);
	vector<Double_t> error=Error(rawcounts);
	AngularCorrHisto(wcounts,isotope,error,gndstate);
//	CountFold(wcounts,errorAlt,gndstate);

}


//**********************************************************************************
//Generates a 3D histogram from the Converted.root file and the THnSparse
//**********************************************************************************
TH3D* angCorr::GetData(){
  
    TFile* isoData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Converted.root");
    isoData->cd("GriffinND");
    THnSparse * dataFull=(THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_addback_sparse");
    TH3D * xyzProj=(TH3D*) dataFull->Projection(0,1,2); //3D Histogram

    return xyzProj;
}

//**********************************************************************************
//Cuts the 3D Histogram at each angular index and projects/plots a histogram for each index
//**********************************************************************************
vector<TH1D*> angCorr::AngCut(TH3D* xyzProj, Double_t gateLow, Double_t gateHigh){

    vector<TH1D*> histVec(52); //Initialize Vector of Histograms

    c1=TCanvas("c1","c1",600,600);
    gStyle->SetOptStat(0);
    c1.SetLogy();

    for(Int_t i=0;i<52;i++){ //loop through angular indices

        if(i!=0&&i!=1&&i!=3){                   //Removes Index 0
        Double_t indLow=i;           //Starts at bin0
        Double_t indHigh=indLow+1;   //cut at 1 bin higher

       //histogram titles
        TString h_title("Angular Histogram for Index ");
        h_title+=indLow;
        h_title+=";Energy (keV);Index ";
        h_title+=indLow;
        cout << h_title.Data() << endl;

       //Projection titles
        TString p_title("Projection Histogram for Index ");
        p_title+=indLow;
        p_title+=";Energy (keV);Counts";
        cout << p_title.Data() << endl;
        
       //Histogram Names
        TString h_name("h_angIndex_");
        h_name+=i;
        cout << h_name.Data() << endl;
        

        delete histVec[i]; //clear memory space
        xyzProj->GetZaxis()->SetRangeUser(indLow,indHigh);      //First Cut
        TH2D * ggProj=(TH2D*) xyzProj->Project3D("xy");
        ggProj->GetXaxis()->SetRangeUser(gateLow,gateHigh);         //Second Cut
        //histVec[i]=ggProj->ProjectionY(h_name,gateLow,gateHigh,"e");
        histVec[i]=ggProj->ProjectionY(h_name,gateLow,gateHigh);
        histVec[i]->SetTitle(p_title);
        histVec[i]->SetMinimum(1.0);                            //Set y axis minimum
        histVec[i]->SetMaximum(20000);                          //Set y axis maximum
        histVec[i]->Draw();

        //Save Histogram
        c1.Update();
        TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Angular_Histograms/NonFitted/");
        path+=h_name;
        path+=".png";
        cout << path.Data() << endl;
        c1.SaveAs(path);
        }
    else{}
    }
    return histVec;
}

//**********************************************************************************
//Takes in vector of histograms,fits a skewed gaussian to the photopeak of interest, integrates fit to find counts, returns counts with background subtracted
//**********************************************************************************
vector<Double_t> angCorr::FitHisto(vector<TH1D*> histVec,Double_t ELow,Double_t EHigh, Double_t secPeak){
	
	c1=TCanvas("c1","c1",600,600);
	gStyle->SetOptStat(0);
	c1.SetLogy();

    //Initialize Vector for fitted histograms
	vector<TF1*> fitVec(52); 
	vector<Double_t> counts(52); 

    //Using the TPEAK Method for fitting
    TPeak* peak=new TPeak(secPeak, ELow, EHigh);
	


	for(Int_t i=0;i<52;i++){ 

        if(i!=0&&i!=1&&i!=3){                   //Removes Index 0
		histVec[i]->GetXaxis()->SetRangeUser(ELow-5,EHigh+5);   //Setting Range of Histogram to allow data fitting
		histVec[i]->SetMinimum(1.0);							//Set y-axis minimum
		histVec[i]->SetMaximum(20000);							//Set y-axis maximum




        delete fitVec[i];
//      fitVec[i]=peak->Fit(histVec[i]);                        //NEED TO CHANGE STORAGE OR VECTOR TYPE
        peak->Fit(histVec[i]);
        counts[i]=peak->GetArea();	    				//Integrating Peak to Find Counts
        histVec[i]->Draw("SAME"); 
	   //Histogram Names
		TString h_name("h_angIndex_");
		h_name+=i;
		h_name+="_fitted";
		cout << h_name.Data() << endl;

		//Save Histogram
		c1.Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Angular_Histograms/Fitted/");
		path+=h_name;
		path+=".png";
		cout << path.Data() << endl;
		c1.SaveAs(path);
        }
        else{}
	}
	return counts;
}

//**********************************************************************************
//Normalizes the counts for the photopeak
//**********************************************************************************
vector<Double_t> angCorr::WeightAdjust(vector<Double_t> counts){

	vector<Double_t> wcounts(52);

	for(Int_t i=0;i<52;i++){
		if(i==1||i==6||i==8||i==25||i==26||i==43||i==45||i==50){
				wcounts[i]=counts[i]/(128);
		}
		else if(i==5||i==10||i==12||i==24||i==27||i==39||i==41||i==46){
				wcounts[i]=counts[i]/(48);
		}
		else if(i==7||i==9||i==11||i==14||i==18||i==20||i==31||i==33||i==37||i==40||i==42||i==44){
				wcounts[i]=counts[i]/(96);
		}
		else{
				wcounts[i]=counts[i]/(64);
		}
	}
	return wcounts;
}

//**********************************************************************************
//Statistical Error Calculation
//**********************************************************************************
vector<Double_t> angCorr::Error(vector<Double_t> counts){

	vector<Double_t> errorAlt(52);

	for(Int_t i=0;i<52;i++){
		if(i==1||i==6||i==8||i==25||i==26||i==43||i==45||i==50){
				errorAlt[i]=TMath::Sqrt(counts[i])/(128);
		}
		else if(i==5||i==10||i==12||i==24||i==27||i==39||i==41||i==46){
				errorAlt[i]=TMath::Sqrt(counts[i])/(48);
		}
		else if(i==7||i==9||i==11||i==14||i==18||i==20||i==31||i==33||i==37||i==40||i==42||i==44){
				errorAlt[i]=TMath::Sqrt(counts[i])/(96);
		}
		else{
				errorAlt[i]=TMath::Sqrt(counts[i])/(64);
		}
	}
	return errorAlt;
}

//**********************************************************************************
//Plots Counts vs Angles
//**********************************************************************************
void angCorr::AngularCorrHisto(vector<Double_t> wcounts, TString isotope,vector<Double_t> countErrors,Int_t gndstate){


	//Angular Indices
	Double_t index[]={0.0,18.79097,25.60153, 26.69036,31.94623,33.65414, 44.36426, 46.79372, 48.57554, 49.79788, 53.83362, 60.15106, 62.70487, 63.08604, 65.01569, 66.46082, 67.45617, 69.86404, 70.86009, 73.08384, 76.38138, 78.66898, 83.04252, 86.22840, 86.23761, 88.47356, 91.52644, 93.76239, 93.77160, 96.95749, 101.33102, 103.61822, 106.91616, 109.13991, 110.13596, 112.54383, 113.53918, 114.98431, 116.91396, 117.29513, 119.84894, 126.16638, 130.20212, 131.42446, 133.20628, 135.63574, 146.34586, 148.05377, 153.30964, 154.39847, 161.21315,180.0};


	vector<Double_t> indexVec(52);
	vector<Double_t> errorCounts(52); //Initalize Vector for Counts
	vector<Double_t> indexErrors(52); //Initalize Vector for Counts

	for(Int_t i=0;i<52;i++){
		indexVec[i]=TMath::Cos(index[i]*TMath::Pi()/(180));
		errorCounts[i]=TMath::Sqrt(wcounts[i]);
		indexErrors[i]=0;
	}

	TCanvas * c1=new TCanvas("c1","c1",600,600);

	TGraphErrors * g=new TGraphErrors(52,&indexVec[0],&wcounts[0],&indexErrors[0],&countErrors[0]);

	//Plot name
		//TString c_name("Angular Correlation for ");
//		c_name+=isotope;
//		c_name+=";Cos(#theta);Normalized Counts";
//		cout << c_name.Data() << endl;
		TString c_name("1173 to 1332 keV decay (4^{+} #rightarrow 2^{+} #rightarrow 0^{+})");
//		c_name+=isotope;
		c_name+=";Cos(#theta);Normalized Counts";
		cout << c_name.Data() << endl;
	
	//Defining Fit Function
//		TF1* th=new TF1("th","(1+1/2*0.357143*(3*x*x-1)+1/8*1.142857*(35*x*x*x*x-30*x*x+3))*13.2",-1,1);
		TF1* efit=new TF1("efit","(1+1/2*[0]*(3*x*x-1)+1/8*[1]*(35*x*x*x*x-30*x*x+3))*[2]",-1,1);
		efit->SetParameter(0,0.5);
		efit->SetParameter(1,1.5);
		efit->SetParameter(2,150);
		efit->SetParName(0,"a2");
		efit->SetParName(1,"a4");
		efit->SetParName(2,"scale");

//		TF1* efitalt=new TF1("efitalt","(1+1/2*[0]*(3*x*x-1))*[1]",0,1);
//		efitalt->SetParameter(0,0.5);
//		efitalt->SetParameter(1,30);
//		efitalt->SetParName(0,"a2");
//		efitalt->SetParName(1,"scale");
//		efitalt->SetLineColor(8);

		g->Fit("efit");

		Double_t chi2;
		chi2=efit->GetChisquare();
        cout << chi2 << endl;

    //Finds min and max counts for graph scaling
    Double_t max=TMath::MaxElement(52,g->GetY());
    Double_t min=efit->GetMinimum();

	g->SetTitle(c_name);
	g->SetMaximum(max+20);
	g->SetMinimum(min-20);
//	g->SetMinimum(20);
	g->SetMarkerStyle(7);
	g->GetYaxis()->SetTitleOffset(1.4);
	g->GetXaxis()->CenterTitle();
	g->GetYaxis()->CenterTitle();
	g->Draw("AP");
//	efit->Draw("SAME");
//	th->Draw("SAME");


	 //Save Histogram
		c1->Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/");
		path+="corr";
		path+=".pdf";
		cout << path.Data() << endl;
		c1->SaveAs(path);
	
//	//Root File  name
//		TString rtpath("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac_state_");
//		rtpath+=gndstate;
//		rtpath+="_full.root";
//		cout<<rtpath.Data() << endl;
	//Root File  name
		TString rtpath("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac");
		rtpath+=50;
		rtpath+="mm-grsi.root";
		cout<<rtpath.Data() << endl;

	//Initialize .root file
		TFile ggacroot(rtpath,"RECREATE");
		ggacroot.WriteTObject(g);
	
}


