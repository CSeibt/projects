#include <iostream>
using namespace std;
#include <fstream>
#include "stdio.h"
#include <string>
#include <cstdio>

#include <iomanip>
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "TLine.h"
#include "TSystem.h"
#include "TString.h"
#include <TUUID.h>
#include "TLegend.h"
#include "TPaveStats.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "THStack.h"

#define TWOTHIRTYONE 2147483648 // 2^31
#define ROIS 4
#define MAXRONLINES 10000000
#define nentriesMAX 10000000
/*/
Newest Version of MakeHist, with every analysis step split up into own respective functions 
/*/


//Open the root file containing the respective histograms
TFile* File(
    TString run = "run021",
    TString prefixPath = "TU5/BeaBA/"
){
	TString folder = "root/";//"";
	TString list = "DataR_";
    TString Suffix = ".bin";
    
    // Variable definitions.
    TString listModeSuffix = list+run+Suffix;
    TString prefix = prefixPath + run + "/" + folder;
    TString input = prefix;

    
    // Open ROOT file
    TFile* rootFile = new TFile(prefix+listModeSuffix+".root","read");
    //if (rootFile->IsZombie()) return;
    
	return rootFile;
	}


//Draw the time difference and energy histograms directly from the root file, it can be choosen between all histograms or only histograms from one detector
void rawHistograms(
	TFile *rootFile,
	TString run,
	bool all = true,    //true = draw every histogram, false = draw histogram chosen with det
	int det = 0         //Choose detector: TU5 = 0, Sz1 = 1, Sz2 = 2; only necessary if all = false
){
	const int numchannels = 3;
	
	TString detectors[numchannels] = {"TU5","Sz1","Sz2"};    			//Detector names
	
	
	TH1I* histoTiming[numchannels][numchannels];
    TH1D* histo_raw[numchannels];
        
    //Open the raw histograms
    for(int i=0;i<numchannels;i++) {
    histo_raw[i] = (TH1D*)rootFile->Get(Form("histo_raw[%d]",i));
       for(int j=0;j<numchannels;j++) {
           histoTiming[i][j] = (TH1I*)rootFile->Get(Form("histoTiming[%d][%d]",i,j));
       }
    }
    if (all){                                                           //true                                         
		TCanvas* ctime[numchannels];									//create time-canvas array
		TCanvas* craw[numchannels];										//create raw-canvas array
		
		for (int j=0;j<numchannels;j++) {								//create time canvas with three time diff histograms
			ctime[j] = new TCanvas(run+"_"+detectors[j]+"_timing",run+"_"+detectors[j]+"_timing",100,10,1200,700);  //900
			ctime[j]->Divide(1,3);
			for (int i=0;i<3;i++) {
				ctime[j]->cd(i+1);
				gPad->SetLogy();
				histoTiming[i][j]->Draw();
				}
			//create raw canvas with raw histogram
			craw[j] = new TCanvas(run+"_"+detectors[j]+"_raw",run+"_"+detectors[j]+"_raw",200,100,1200,700);
			//c0->SetTitle(detectors[j]+" "+run+ " raw");
			gPad->SetLogy();
			histo_raw[j]->SetTitle(detectors[j]+" raw");
			histo_raw[j]->SetXTitle("Channels");
			histo_raw[j]->SetYTitle("Counts");
			//histo_raw[j]->SetBinContent(1, 0);
			histo_raw[j]->Draw();
			}
		}
	else{
		TString detector = "";
		TCanvas* c0 = new TCanvas(run+"_"+detectors[det]+"_timing",run+"_"+detectors[det]+"_timing",100,10,1200,700);  //900
		c0->Divide(1,3);

		for (int i=0;i<3;i++) {
			c0->cd(i+1);
			gPad->SetLogy();
			histoTiming[i][det]->Draw();
			}
		
		TCanvas* c1 = new TCanvas(run+"_"+detectors[det]+"_raw",run+"_"+detectors[det]+"_raw",200,100,1200,700);
		//c1->SetTitle(run + " raw");
		gPad->SetLogy();
		histo_raw[det]->SetTitle(detectors[det] + " raw");
		histo_raw[det]->SetXTitle("Channels");
		histo_raw[det]->SetYTitle("Counts");
		histo_raw[det]->SetBinContent(1, 0);
		histo_raw[det]->Draw();
		
		}
    //cout << detectors[2] << "and" << detectors[0] << endl;
	}

//Draw the Energy histogram of one detector, regular (blue), with pile-up events (red) and saturation events (green) separated
void Flags(
	TTree* data,
	Int_t det = 1,
	Int_t nch = 16384,
	bool pileup = true,
	bool saturation = true
){
	const int numchannels = 3;
	
	TString detectors[numchannels] = {"TU5","Sz1","Sz2"};    			//Detector names
	
	
	TH1I* histo_adc = new TH1I("histo_adc", detectors[det]+" ADC spectrum; Channel; Counts", nch, 0, nch);
	//The previous line sets up: 			Title,						   XTitle 	YTitle
    data->Draw("adc>>histo_adc", Form("det==%i", det), "goff");
    
    TCanvas* c1_pileup_saturation = new TCanvas(detectors[det]+"_adc_spectra_with_flags", detectors[det]+"_adc_spectra_with_flags", 100,10,1200,700);
    gPad->SetLogy(true);
    histo_adc->SetStats(0);
    histo_adc->Draw();

    TLegend* leg = new TLegend(0.6,0.8,0.9,0.9); //new legend with x_min,y_min,x_max,y_max
    leg->SetFillColor(0); //fill color white
    leg->SetTextSize(0.03); 
  
    leg->AddEntry(histo_adc, "raw","l");
    if (pileup){
		TH1I* histo_pileup = new TH1I("histo_pileup", "spectrum; Channel; Counts", nch, 0, nch);
		data->Draw("adc>>histo_pileup", Form("det==%i && pileup==1", det), "goff");
		histo_pileup->SetLineColor(kRed);
		histo_pileup->Draw("same");	
		leg->AddEntry(histo_pileup, "pile-up","l");
	}
	if (saturation){
		TH1I* histo_saturation = new TH1I("histo_saturation", "spectrum; Channel; Counts", nch, 0, nch);
		data->Draw("adc>>histo_saturation", Form("det==%i && saturation!=0", det), "goff");
		histo_saturation->SetLineColor(kGreen);
		histo_saturation->Draw("same");
		leg->AddEntry(histo_saturation, "saturation","l");
	}
    leg->Draw("same");
    gPad->Update();
	}

//Draw the Energy histogram of one detector, regular (blue), with time and energy cuts (green), and reverse time cuts (red)
void ADCcuts(
	TTree* data,
	Int_t detA,															//Detector for spectrum
	Int_t detB,															//Detector for coincidence cuts
	Double_t thrTimeDiff,												//Time difference threshold for coincidence, in ps
	Double_t thrEnergy,													//Energy threshold for pile-up events, in channels
	Int_t nch = 16384												
){
	
	const int numchannels = 3;
	
	TString detectors[numchannels] = {"TU5","Sz1","Sz2"};    			//Detector names
	
	
    TH1I* histoTcut0 = new TH1I("histoTcut0", detectors[detA]+" ADC spectrum; Channel; Counts", nch, 0, nch);
	//The previous line sets up: 			  Title,						 XTitle   YTitle
    data->Draw("adc>>histoTcut0", Form("det==%i", detA), "goff");
    TH1I* histoTcut1 = new TH1I("histoTcut1", "spectrum; Channel; Counts", nch, 0, nch);
    data->Draw("adc>>histoTcut1", Form("det==%i && adc>%f &&((TimeDiff_before%i<%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<%f && EnergyDep_after%i>%f))", detA, thrEnergy, detB, thrTimeDiff, detB, thrEnergy, detB, thrTimeDiff, detB, thrEnergy), "goff");
    //prev line: Draw energy histogram of Detector A only with: Events with energy >1 and an event in Detector B within the time threshold, excluding pile-up 
    TH1I* histoTcut2 = new TH1I("histoTcut2", "spectrum; Channel; Counts", nch, 0, nch);
    data->Draw("adc>>histoTcut2", Form("det==%i && !(TimeDiff_before%i<%f || TimeDiff_after%i<%f)", detA, detB, thrTimeDiff, detB, thrTimeDiff), "goff");
    //prev line: Draw energy histogram of Detector A without coincidence (Pile-Up is ignored here, must be fixed)
    // draw energy spectra
    TCanvas* c1_Tcut = new TCanvas(detectors[detA]+"_Energy_spectra_with_TimeDiff_cuts", detectors[detA]+"_Energy_spectra_with_TimeDiff_cuts", 100,10,1200,700);
    gPad->SetLogy(true);
    histoTcut0->SetStats(0);
    histoTcut0->Draw();
    histoTcut1->SetLineColor(kGreen);
    histoTcut1->Draw("same");
    histoTcut2->SetLineColor(kRed);
    histoTcut2->Draw("same");
    
    TLegend* leg1 = new TLegend(0.6,0.8,0.9,0.9); //new legend with x_min,y_min,x_max,y_max
    leg1->SetFillColor(0); //fill color white
    leg1->SetTextSize(0.03); 
 
    leg1->AddEntry(histoTcut0, "energy dep w/out cuts","l");
    leg1->AddEntry(histoTcut1, "energy dep with cuts","l");
    leg1->AddEntry(histoTcut2, "energy dep with cuts reverse","l");
    leg1->Draw("same");
    gPad->Update();
	}

//Draw the time difference histogram of one detector in respect to a second detector, regular (blue), with time diff and energy cuts (green) and with time diff reverse (red)
void TimeDiffcuts(
	TTree* data,
	Int_t detA,															//Detector for spectrum
	Int_t detB,															//Detector for coincidence cuts
	Double_t thrTimeDiff,												//Time difference threshold for coincidence, in ps
	Double_t thrEnergy,													//Energy threshold for pile-up events, in channels
	Int_t nch = 16384
){
	const int numchannels = 3;
	
	TString detectors[numchannels] = {"TU5","Sz1","Sz2"};    			//Detector names
	
	TH1I* histoTtime0 = new TH1I("histoTtime0", "t("+detectors[detA]+")-t("+detectors[detB]+"); dt / ps; Counts", 2001, -4E6-2000, 4E6+2000);
	//prev line:								Title											XTitle	 YTitle
    data->Draw(Form("TimeDiff_before%i>>histoTtime0", detB), Form("det==%i", detA), "goff");
    data->Draw(Form("-TimeDiff_after%i>>+histoTtime0", detB), Form("det==%i", detA), "goff");
    //prev lines (2): Time differences for each event in detA: draw time differences labeled with detB
    TH1I* histoTtime1 = new TH1I("histoTtime1", "spectrum; Timediff; Counts", 2001, -4E6-2000, 4E6+2000);
    data->Draw(Form("TimeDiff_before%i>>histoTtime1", detB), Form("det==%i && adc>%f && ((TimeDiff_before%i<=%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<=%f && EnergyDep_after%i>%f))", detA, thrEnergy, detB, thrTimeDiff, detB, thrEnergy, detB, thrTimeDiff, detB, thrEnergy), "goff");		//Time diff before
    data->Draw(Form("-TimeDiff_after%i>>+histoTtime1", detB), Form("det==%i && adc>%f && ((TimeDiff_before%i<=%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<=%f && EnergyDep_after%i>%f))", detA, thrEnergy, detB, thrTimeDiff, detB, thrEnergy, detB, thrTimeDiff, detB, thrEnergy), "goff");		//Time diff after
    //prev lines (2): Time differences for each event in det A labeled with detB 
    //only with E(detA)>energy threshold, time difference detB-detA<=time threshold, E(detB)>energy threshold 
    TH1I* histoTtime2 = new TH1I("histoTtime2", "spectrum; Timediff; Counts", 2001, -4E6-2000, 4E6+2000);
    data->Draw(Form("TimeDiff_before%i>>histoTtime2", detB), Form("det==%i && !(TimeDiff_before%i<%f || TimeDiff_after%i<%f)", detA, detB, thrTimeDiff, detB, thrTimeDiff), "goff");		//Time diff before
    data->Draw(Form("-TimeDiff_after%i>>+histoTtime2", detB), Form("det==%i && !(TimeDiff_before%i<%f || TimeDiff_after%i<%f)", detA, detB, thrTimeDiff, detB, thrTimeDiff), "goff");		//Time diff after
    //prev lines (2): Time differences for each event in detA labeled with detB, excluding coincidence events (dT<time threshold)
    // draw energy spectra
    TCanvas* c1_Ttime = new TCanvas(detectors[detA]+"-"+detectors[detB]+"_coincidence_spectrum", detectors[detA]+"-"+detectors[detB]+"_coincidence_spectrum", 100,10,1200,700);
    gPad->SetLogy(true);
    histoTtime0->SetStats(0);
    histoTtime0->Draw();
    histoTtime1->SetLineColor(kGreen);
    histoTtime1->Draw("same");
    histoTtime2->SetLineColor(kRed);
    histoTtime2->Draw("same");
    
    TLegend* leg2 = new TLegend(0.6,0.8,0.9,0.9); //new legend with x_min,y_min,x_max,y_max
    leg2->SetFillColor(0); //fill color white
    leg2->SetTextSize(0.03); 
  
    leg2->AddEntry(histoTtime0, "time diff w/out cuts","l");
    leg2->AddEntry(histoTtime1, "time diff with cuts","l");
    leg2->AddEntry(histoTtime2, "time diff with cuts reverse","l");
    leg2->Draw("same");
    gPad->Update();
	}

//Draw the time differences of detector A in respect to the detector B, but only considering coincidence events of detector B, which are above energy threshold
TH1I* TimeDiffMuons(
	TTree* data,
	Int_t detA,
	Int_t detB,
	Int_t detC,
	Double_t thrTimeDiff,												//Time difference threshold for coincidence, in ps
	Double_t thrEnergy,													//Energy threshold for pile-up events, in channels
	Int_t nch = 16384
){
	const int numchannels = 3;
	
	TString detectors[numchannels] = {"TU5","Sz1","Sz2"};    			//Detector names
	TH1I* histoCoincidence = new TH1I("histoCoincidence", "t("+detectors[detA]+")-t("+detectors[detB]+"); dt / ps; Counts", 2001, -4E6-2000, 4E6+2000);
    data->Draw(Form("-TimeDiff_before%i>>histoCoincidence", detA), Form("det==%i && adc>%f && EnergyDep_before%i>%f && ((TimeDiff_before%i<%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<%f && EnergyDep_after%i>%f))", detB, thrEnergy, detA, thrEnergy, detC, thrTimeDiff, detC, thrEnergy, detC, thrTimeDiff, detC, thrEnergy), "goff");
    data->Draw(Form("TimeDiff_after%i>>+histoCoincidence", detA), Form("det==%i && adc>%f && EnergyDep_after%i>%f && ((TimeDiff_before%i<%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<%f && EnergyDep_after%i>%f))", detB, thrEnergy, detA, thrEnergy, detC, thrTimeDiff, detC, thrEnergy, detC, thrTimeDiff, detC, thrEnergy), "goff");
    //prev lines (2): Draw time difference of detector B labeled with detector A, 
    //only including entries with E(detB)>energy threshold, E(detA)>energy threshold and coincidence between detector B and detector C [E(detC)>energy threshold]
    TH1I* histoTtime = new TH1I("histoTtime", "t("+detectors[detA]+")-t("+detectors[detB]+"); dt / ps; Counts", 2001, -4E6-2000, 4E6+2000);
    data->Draw(Form("-TimeDiff_before%i>>histoTtime", detB), Form("det==%i && adc>%f && EnergyDep_before%i>%f", detA, thrEnergy, detB, thrEnergy), "goff");
    data->Draw(Form("TimeDiff_after%i>>+histoTtime", detB), Form("det==%i && adc>%f && EnergyDep_after%i>%f", detA, thrEnergy, detB, thrEnergy), "goff");
    //prev lines (2): Draw time differences of detector B labeled with detector a, only including entries with E(detB)>energy threshold and E(detA)>energy threshold
    TCanvas* c1_Coin = new TCanvas("Coincidence "+detectors[detA]+" - "+detectors[detB]+" with muons only", "Coincidence "+detectors[detA]+" - "+detectors[detB]+" with muons only",100,10,1200,700);
    gPad->SetLogy(true);
    histoTtime->SetStats(0);
    histoTtime->Draw();    
    histoCoincidence->SetStats(1);
    histoCoincidence->SetLineColor(kGreen);
    histoCoincidence->Draw("sames");
    
    TLegend* leg2 = new TLegend(0.5,0.8,0.7,0.9); 						//new legend with x_min,y_min,x_max,y_max
    leg2->SetFillColor(0); //fill color white
    leg2->SetTextSize(0.03); 
  
    leg2->AddEntry(histoTtime, "time diff w/out cuts","l");
    leg2->AddEntry(histoCoincidence, "time diff with cuts","l");
    leg2->Draw("same");
    gPad->Update();
    
	return histoCoincidence;
	
	} 

TH1I* MuonVeto(
	TTree* data,
	Int_t detA,
	Int_t detB,
	Int_t detC,
	Double_t thrTimeDiff,
	Double_t thrEnergy,
	Int_t nch = 16384
){
	const int numchannels = 3;
	
	TString detectors[numchannels] = {"TU5","Sz1","Sz2"};    			//Detector names
	
	TH1I* histo =  new TH1I("histo", detectors[detA]+"; Channel; Counts", nch, 0 ,nch);
	data->Draw("adc>>histo",
	Form("det==%i && adc>%f", 
	detA, thrEnergy), "goff"); //&& adc>%f
	
	TH1I* histoCuts = new TH1I("histoCuts", detectors[detA]+" Veto Cuts; Channel; Counts", nch, 0, nch);
	data->Draw("adc>>histoCuts", 
	Form("det==%i && adc>%f && ((TimeDiff_before%i<=%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<=%f && EnergyDep_after%i>%f) || (TimeDiff_before%i<=%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<=%f && EnergyDep_after%i>%f))", 
	detA, thrEnergy, detB, thrTimeDiff, detB, thrEnergy, detB, thrTimeDiff, detB, thrEnergy, detC, thrTimeDiff, detC, thrEnergy, detC, thrTimeDiff, detC, thrEnergy), "goff"); //&& adc>%f
	
	TH1I* histoVeto =  new TH1I("histoVeto", detectors[detA]+" - Muon Veto Cuts; Channel; Counts", nch, 0 ,nch);
	data->Draw("adc>>histoVeto",
	Form("det==%i && adc>%f && !((TimeDiff_before%i<=%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<=%f && EnergyDep_after%i>%f) || (TimeDiff_before%i<=%f && EnergyDep_before%i>%f) || (TimeDiff_after%i<=%f && EnergyDep_after%i>%f))", 
	detA, thrEnergy, detB, thrTimeDiff, detB, thrEnergy, detB, thrTimeDiff, detB, thrEnergy, detC, thrTimeDiff, detC, thrEnergy, detC, thrTimeDiff, detC, thrEnergy), "goff"); //&& adc>%f
	
	TCanvas* c1_Veto = new TCanvas(detectors[detA]+"_Energy_spectra_with_veto", detectors[detA]+"_Energy_spectra_with_veto", 100, 10, 1200, 700);
	gPad->SetLogy(true);
	histo->SetStats(0);
	histo->Draw();
	histoVeto->SetStats(0);
	histoVeto->SetLineColor(kCyan);
	histoVeto->Draw("same");
	histoCuts->SetStats(0);
	histoCuts->SetLineColor(kRed);
	histoCuts->Draw("same");
	cout << "Integral ohne Veto: "<< histo->Integral(0, nch) << endl;
	cout << "Integral mit Veto: "<< histoVeto->Integral(0, nch) << endl;
	cout << "Integral d. Veto cuts: "<< histoCuts->Integral(0, nch) <<endl;
	return histoVeto;
	}
	
	
//Rate-function: Timing/Rate of events with pile-up (red) and saturation (green) flags, 
//!!!CURRENTLY NOT WORKING!!!
void Rate(
	TTree* data,
	Int_t det,
	bool pileup = true,
	bool saturation = true
){
	const int numchannels = 3;
	
	TString detectors[numchannels] = {"TU5","Sz1","Sz2"};    			//Detector names
	
	TBranch* t = data->GetBranch("Time");
	if (!t) cout << "Could not get branch Time" << endl;
	
	Int_t tmax = 300000;												//Measuring time (there may be a possibility f determining this)
	
	TH1I* histoRate = new TH1I("histoRate", detectors[det]+" signal rate; Time / s; Counts", tmax, 0, tmax);
	data->Draw("Time>>histoRate", Form("det==%i", det), "goff");		//Time - Events of det
	
	
	TH1I* histoRPileup = new TH1I("histoRPileup", detectors[det]+" signal rate of pile-up; Time / s; Counts", tmax, 0, tmax);
	data->Draw("Time>>histoRPileup", Form("det==%i && pileup==1", det), "goff");
	TH1I* histoRSaturation = new TH1I("histoRSaturation", detectors[det]+" signal rate of saturation; Time / s; Counts", tmax, 0, tmax);
	data->Draw("Time>>histoRSaturation", Form("det==%i && saturation!=0", det), "goff");
	
	TCanvas* c1_Time = new TCanvas(detectors[det]+"_signal_rate", detectors[det]+"_signal_rate", 100, 10, 1200, 700);
	gPad->SetLogy(true);
//	histoRate->SetStats(0); 
	histoRate->Draw();
	histoRPileup->SetLineColor(kRed);
	histoRPileup->Draw("same");
	histoRSaturation->SetLineColor(kGreen);
	histoRSaturation->Draw("same");
	
	TLegend* leg = new TLegend(0.6,0.8,0.9,0.9); //new legend with x_min,y_min,x_max,y_max
    leg->SetFillColor(0); //fill color white
    leg->SetTextSize(0.03); 
  
    leg->AddEntry(histoRate, "raw","l");
    leg->AddEntry(histoRPileup, "pile-up","l");
    leg->AddEntry(histoRSaturation, "saturation","l");
    leg->Draw("same");
    gPad->Update(); 
	}
	
void MakeHist_220704(){
	//General cosmetics
	gStyle->SetOptStat(1001111);
    gStyle->SetOptFit(0);
    gStyle->SetStripDecimals(kFALSE);

    //Variables
    TString run = "run021";			//run
    Int_t TU5 = 0;					//TU5 (X-ray detector)
    Int_t Sz1 = 1;					//Muon panel 1
    Int_t Sz2 = 2;					//Muon Panel 2
    Double_t thrTimediff = 1.E4;	//time difference threshold
    Double_t thrEnergy = 1;			//energy threshold (pile-up)
    Int_t nch = 16384;
    
    
	TFile* file = File();			//Open File 
	TTree* data = (TTree*)file->Get("Data");		//Get data from file
	
	rawHistograms(file, run, true, TU5);									//Open raw histograms
	
	Flags(data, TU5);														//Energy historgam with pile-up and saturation
	
	ADCcuts(data, TU5, Sz1, thrTimediff, thrEnergy, nch);				//Energy histogram with coincidence cuts

	TimeDiffcuts(data, TU5, Sz1, thrTimediff, thrEnergy, nch);			//Time difference histogram with coincidence cuts
		
	TH1I* coincidence = TimeDiffMuons(data, TU5, Sz2, Sz1, thrTimediff, thrEnergy, nch);	//Time difference histogram detA - detB with coincidence cuts detB - detC
	
	TH1I* veto = MuonVeto(data, TU5, Sz1, Sz2, thrTimediff, thrEnergy, nch);

	Rate(data, TU5);
	}
