#include <iostream>
using namespace std;
#include <fstream>
#include "stdio.h"
#include <string>
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

#include <dirent.h>

#define ROIS 4
#define MAXRONLINES 10000000
#define nentriesMAX 2000000000


int* arrayay(int* arr[]){
	return arr;
	}



void MakeTree_220707(){
	TString run = "run015";
	int* y[10];
	int* a = arrayay(y);
	cout << a[0] << endl;
	}


