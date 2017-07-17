//Connor Natzke
//July 28 2016
//
//This is an alternate macro to generate gg corr plots. It splits by angular index then by energy. 
//

#include <iostream>
#include <vector>
#include "AngularCorrAlt.H"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"

using namespace std;

AngularCorrAlt::AngularCorrAlt(){
	Double_t gateLow=1171;	//First Gate
	Double_t gateHigh=1173;	//First Gate
	Double_t ELow=1330;		//Low Energy Gate
	Double_t EHigh=1334;	//High Energy Gate
	Int_t angSplits=52;		//Number of Angular Indices
	TString isotope="^{60}Co"; //Isotope
	Int_t gndstate=0;		//Ground State

	TH3D *xyzProj=GetData();
	vector<TH1D*> histVec=AngCut(xyzProj,angSplits,gateLow, gateHigh);	
	vector<Double_t> rawcounts=FitHisto(histVec, angSplits, ELow, EHigh);
	vector<Double_t> errorFit=ErrorFits(histVec, angSplits, ELow, EHigh);
	//Only use if you are using the full GRIFFIN Detector array, weights are specific
	vector<Double_t> wcounts=WeightAdjust(rawcounts,angSplits);
	vector<Double_t> countsError=ErrorAdjust(errorFit,angSplits);
	vector<Double_t> errorAlt=AltError(rawcounts,angSplits);
	AngularCorrHisto(wcounts,isotope,errorAlt,gndstate);
	CountFold(wcounts,errorAlt,gndstate);
}

//**********************************************************************************
//Generates a 3D histogram from the Converted.root file and the THnSparse
//**********************************************************************************
TH3D* AngularCorrAlt::GetData(){
  
	TFile* isoData = new TFile("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Converted.root");
	isoData->cd("GriffinND");
	THnSparse * dataFull=(THnSparse*) gDirectory->Get("griffin_crystal_unsup_gamma_gamma_corr_edep_cry_sparse");
	TH3D * xyzProj=(TH3D*) dataFull->Projection(0,1,2);	//3D Histogram

	return xyzProj;
}

//**********************************************************************************
//Cuts the 3D Histogram at each angular index and projects/plots a histogram for each index
//**********************************************************************************
vector<TH1D*> AngularCorrAlt::AngCut(TH3D* xyzProj, Int_t angSplits,Double_t gateLow, Double_t gateHigh){

	vector<TH1D*> histVec(angSplits); //Initialize Vector of Histograms

	c1=TCanvas("c1","c1",600,600);
	gStyle->SetOptStat(10);
	c1.SetLogy();

	for(Int_t i=0;i<angSplits;i++){ //loop through angular indices

		Double_t indLow=i; //Starts at bin0
		Double_t indHigh=indLow+1; //cut at 1 bin higher

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
		xyzProj->GetZaxis()->SetRangeUser(indLow,indHigh);		//First Cut
		TH2D * ggProj=(TH2D*) xyzProj->Project3D("xy");
		ggProj->GetXaxis()->SetRangeUser(gateLow,gateHigh);			//Second Cut
		//histVec[i]=ggProj->ProjectionY(h_name,gateLow,gateHigh,"e");
		histVec[i]=ggProj->ProjectionY(h_name,gateLow,gateHigh);
		histVec[i]->SetTitle(p_title);
		histVec[i]->SetMinimum(1.0);							//Set y axis minimum
		histVec[i]->SetMaximum(20000);							//Set y axis maximum
		histVec[i]->Draw();

		//Save Histogram
		c1.Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/Angular_Histograms/NonFitted/");
		path+=h_name;
		path+=".png";
		cout << path.Data() << endl;
		c1.SaveAs(path);
	}
	return histVec;
}

//**********************************************************************************
//Takes in vector of histograms,fits a gaussian to the photopeak of interest, integrates fit to find counts, returns counts
//**********************************************************************************
vector<Double_t> AngularCorrAlt::FitHisto(vector<TH1D*> histVec,Int_t angSplits, Double_t ELow,Double_t EHigh){
	
	c1=TCanvas("c1","c1",600,600);
	gStyle->SetOptStat(10);
	c1.SetLogy();
	
	//Describing the Fit
	TF1 * gFit=new TF1("gFit","gaus",ELow, EHigh);
	vector<Double_t> counts(angSplits); //Initalize Vector for Counts


	for(Int_t i=0;i<angSplits;i++){ 
		histVec[i]->GetXaxis()->SetRangeUser(ELow-5,EHigh+5);	//Setting Range of Histogram to allow data fitting
		histVec[i]->SetMinimum(1.0);							//Set y-axis minimum
		histVec[i]->SetMaximum(20000);							//Set y-axis maximum
		histVec[i]->Fit("gFit","Q");							//Fitting
		counts[i]=gFit->Integral(ELow,EHigh);					//Integrating Peak to Find Counts
//		histVec[i]->GetXaxis()->SetRangeUser(0,2500);			//Resetting Range for Cosmetic Purposes
		histVec[i]->Draw();
		gFit->Draw("SAME");

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
	return counts;
}
//**********************************************************************************
//Calculates errors from histogram fit and integral for error bars on the final correlation plot
//**********************************************************************************
vector<Double_t> AngularCorrAlt::ErrorFits(vector<TH1D*> histVec,Int_t angSplits, Double_t ELow,Double_t EHigh){
	
	c1=TCanvas("c1","c1",600,600);
	gStyle->SetOptStat(10);
	c1.SetLogy();
	
	//Describing the Fit
	TF1 * gFit=new TF1("gFit","gaus",ELow, EHigh);
	vector<Double_t> countsError(angSplits); //Initalize Vector for Counts


	for(Int_t i=0;i<angSplits;i++){ 
		histVec[i]->GetXaxis()->SetRangeUser(ELow-5,EHigh+5);	//Setting Range of Histogram to allow data fitting
		histVec[i]->SetMinimum(1.0);							//Set y-axis minimum
		histVec[i]->SetMaximum(20000);							//Set y-axis maximum
		histVec[i]->Fit("gFit","Q");							//Fitting
		countsError[i]=gFit->IntegralError(ELow,EHigh);			//Finds Error in Counts
//		histVec[i]->GetXaxis()->SetRangeUser(0,2500);			//Resetting Range for Cosmetic Purposes
		histVec[i]->Draw();
		gFit->Draw("SAME");

		
//		cout << countsError[i] << " ";
	}
//	cout << endl;
	return countsError;
}
//**********************************************************************************
//Normalizes the counts for the photopeak
//**********************************************************************************
vector<Double_t> AngularCorrAlt::WeightAdjust(vector<Double_t> counts, Int_t angSplits){

	vector<Double_t> wcounts(angSplits);

	for(Int_t i=0;i<angSplits;i++){
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
		cout << wcounts[i] << " ";
	}
	cout << endl;
	return wcounts;
}
//**********************************************************************************
//Alternate Error Calculation
//**********************************************************************************
vector<Double_t> AngularCorrAlt::AltError(vector<Double_t> counts, Int_t angSplits){

	vector<Double_t> errorAlt(angSplits);

	for(Int_t i=0;i<angSplits;i++){
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
//Normalizes Errors
//**********************************************************************************
vector<Double_t> AngularCorrAlt::ErrorAdjust(vector<Double_t> errorfit, Int_t angSplits){

	vector<Double_t> wcountsError(angSplits);

	for(Int_t i=0;i<angSplits;i++){
		if(i==1||i==6||i==8||i==25||i==26||i==43||i==45||i==50){
				wcountsError[i]=errorfit[i]/(128);
		}
		else if(i==5||i==10||i==12||i==24||i==27||i==39||i==41||i==46){
				wcountsError[i]=errorfit[i]/(48);
		}
		else if(i==7||i==9||i==11||i==14||i==18||i==20||i==31||i==33||i==37||i==40||i==42||i==44){
				wcountsError[i]=errorfit[i]/(96);
		}
		else{
				wcountsError[i]=errorfit[i]/(64);
		}
//		cout << wcountsError[i] << " ";
	}
//	cout << endl;
	return wcountsError;
}
//**********************************************************************************
//Plots Counts vs Angles
//**********************************************************************************
void AngularCorrAlt::AngularCorrHisto(vector<Double_t> wcounts, TString isotope,vector<Double_t> countErrors,Int_t gndstate){


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
		TF1* th=new TF1("th","(1+1/2*0.357143*(3*x*x-1)+1/8*1.142857*(35*x*x*x*x-30*x*x+3))*13.2",-1,1);
		TF1* efit=new TF1("efit","(1+1/2*[0]*(3*x*x-1)+1/8*[1]*(35*x*x*x*x-30*x*x+3))*[2]",-1,1);
		efit->SetParameter(0,0.5);
		efit->SetParameter(1,1.5);
		efit->SetParameter(2,30);
		efit->SetParName(0,"a2");
		efit->SetParName(1,"a4");
		efit->SetParName(2,"scale");

		TF1* efitalt=new TF1("efitalt","(1+1/2*[0]*(3*x*x-1))*[1]",0,1);
		efitalt->SetParameter(0,0.5);
		efitalt->SetParameter(1,30);
		efitalt->SetParName(0,"a2");
		efitalt->SetParName(1,"scale");
		efitalt->SetLineColor(8);

//		g->Fit("efit");

		Double_t chi2;
		chi2=efitalt->GetChisquare();
		cout << chi2 << endl;

	g->SetTitle(c_name);
	g->SetMaximum(160);
	g->SetMinimum(100);
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
		rtpath+="mm.root";
		cout<<rtpath.Data() << endl;

	//Initialize .root file
		TFile ggacroot(rtpath,"RECREATE");
		ggacroot.WriteTObject(g);
	
}
//**********************************************************************************
//Folding The counts
//**********************************************************************************

void AngularCorrAlt::CountFold(vector<Double_t> counts,vector<Double_t> errors,Int_t gndstate){

	vector<Double_t> countsFoldedFull(26);
	vector<Double_t> countsFolded(13);
	vector<Double_t> errorVec(26);
	vector<Double_t> indexError(26);

  for(Int_t i=0;i<26;i++){
  	countsFoldedFull[i]=counts[i]+counts[52-i];
	errorVec[i]=TMath::Sqrt(TMath::Power(errors[i],2)+TMath::Power(errors[52-i],2));
	indexError[i]=0;
// 	cout << errorVec[i] << " ";
  }
//  cout << endl;
  
  //grouping similar angles
	countsFolded[0]=(countsFoldedFull[1]);
	countsFolded[1]=(countsFoldedFull[2]+countsFoldedFull[3])/2;
	countsFolded[2]=(countsFoldedFull[4]+countsFoldedFull[5])/2;
	countsFolded[3]=(countsFoldedFull[6]+countsFoldedFull[7]i+countsFoldedFull[8]+countsFoldedFull[9])/4;
	countsFolded[4]=countsFoldedFull[10];
	countsFolded[5]=(countsFoldedFull[11]+countsFoldedFull[12]+countsFoldedFull[13])/3;
	countsFolded[6]=(countsFoldedFull[14]+countsFoldedFull[15]+countsFoldedFull[16])/3;
	countsFolded[7]=(countsFoldedFull[17]+countsFoldedFull[18]+countsFoldedFull[19])/3;
	countsFolded[8]=(countsFoldedFull[20]+countsFoldedFull[21]+countsFoldedFull[22])/3;
	countsFolded[9]=(countsFoldedFull[18]+countsFoldedFull[19])/2;
	countsFolded[10]=(countsFoldedFull[20]+countsFoldedFull[21])/2;
	countsFolded[11]=(countsFoldedFull[22]+countsFoldedFull[23])/2;
	countsFolded[12]=(countsFoldedFull[24]+countsFoldedFull[25])/2;

	//Defining the reduced angles 
	Double_t indexRed[]={0.0,18.0,25.5,32.0,46.5,53.0,61.5,67.0,73.0,80.5,87.0};
	Double_t index[]={0.0,18.79097,25.60153, 26.69036,31.94623,33.65414, 44.36426, 46.79372, 48.57554, 49.79788, 53.83362, 60.15106, 62.70487, 63.08604, 65.01569, 66.46082, 67.45617, 69.86404, 70.86009, 73.08384, 76.38138, 78.66898, 83.04252, 86.22840, 86.23761, 88.47356};
	
	vector<Double_t> indexVecRed(26);
	for(Int_t i=0;i<26;i++){
		indexVecRed[i]=TMath::Cos(index[i]*TMath::Pi()/(180));
	}
//	for(Int_t i=0;i<13;i++){
//		indexVecRed[i]=TMath::Cos(index[i]*TMath::Pi()/(180));
//	}


	TCanvas * c2=new TCanvas("c2","c2",600,600);
	TGraphErrors* f=new TGraphErrors(26,&indexVecRed[0],&countsFoldedFull[0],&indexError[0],&errorVec[0]);
		TF1* efit=new TF1("efit","(1+1/2*[0]*(3*x*x-1)+1/8*[1]*(35*x*x*x*x-30*x*x+3))*[2]",0,1);
		efit->SetParameter(0,0.5);
		efit->SetParameter(1,1.5);
		efit->SetParameter(2,30);
		efit->SetParName(0,"a2");
		efit->SetParName(1,"a4");
		efit->SetParName(2,"scale");

		TF1* efitalt=new TF1("efitalt","(1+1/2*[0]*(3*x*x-1))*[1]",0,1);
		efitalt->SetParameter(0,0.5);
		efitalt->SetParameter(1,30);
		efitalt->SetParName(0,"a2");
		efitalt->SetParName(1,"scale");
		efitalt->SetLineColor(8);

		f->Fit("efit");

		Double_t chi2;
		chi2=efitalt->GetChisquare();
		cout << chi2 << endl;
//	TGraph* f=new TGraph(13,&indexVecRed[0],&countsFolded[0]);

		TString c_name("1173 to 1332 keV decay (4^{+} #rightarrow 2^{+} #rightarrow 0^{+})");
//		c_name+=isotope;
		c_name+=";Cos(#theta);Normalized Counts";
		cout << c_name.Data() << endl;

	f->SetMarkerStyle(7);
	f->Draw("AP");
	f->SetTitle(c_name);
	f->SetMaximum(300);
	f->SetMinimum(200);
	f->SetMarkerStyle(7);
	f->GetYaxis()->SetTitleOffset(1.4);
	f->GetXaxis()->CenterTitle();
	f->GetYaxis()->CenterTitle();
	f->Draw("AP");

	//Root File  name
		TString rtpath("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/ggac_state_");
		rtpath+=gndstate;
		rtpath+=".root";
		cout<<rtpath.Data() << endl;

	//Initialize .root file
		TFile ggacroot(rtpath,"RECREATE");
		ggacroot.WriteTObject(f);

	efitalt->Draw("SAME");



     //Save Histogram
		c2->Update();
		TString path("/home/cnatzke/Griffin_Sim_v10/detectorSimulations_v10-build/Simulation_Results/");
		path+="corrfolded";
		path+=".pdf";
		cout << path.Data() << endl;
		c2->SaveAs(path);
	

}


