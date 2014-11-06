/*
To compile: g++ -o lookForGainDrift lookForGainDrift.cc `root-config --cflags --libs`
I have this set up so I have a text file that lists all the runs I want to look over. 
I do it this way because I've already done some data selection on the data set I'm looking 
at when I run this code, so I want to omit some runs from a given range. If you just want to use
a set range this code can be easily changed to do that.
So to run I would do: ./lookForGainDrift -c [detector channel number] -n [number of runs in my txt file]
On PDSF you can just do: qsub process.csh -c [detector channel number] -n [number of runs in my txt file]
and in process.csh have the command: ./lookForGainDrift $1 $2 $3 $4
*/

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TAxis.h"
#include "TF1.h"
#include "TKey.h"
#include <algorithm>
#include <fstream>
#include <iostream>

int lookForGainDrift(int channel, int ct){  

  int start_run_no = 4494;
  int end_run_no =   5331;
  int base = 40000000;
  char partID[7] = "P3CK3";
  
  double x_val[ct], pulserAvg[ct], pulserArea[ct], pulserOverArea[ct];
  
  //I have this text file in the directory I'm running in and it has a list of all 
  //the runs I want to use for the data set.
  char txtName[100];
  sprintf(txtName, "accepted_Runs_Ch%d_%d_%d.txt", channel, start_run_no, end_run_no);

  std::ifstream inputTxt(txtName, std::ifstream::in);
  if(!(inputTxt.is_open())){
    printf("Input file is not open\n");
    return 1;
    }
  for(int i=0; i<ct; i++) inputTxt>>x_val[i];
  inputTxt.close();
  
  TFile *f;
  TTree *t;
  TH1D *htemp;
  double jmin = 0.;
  double jmax = 3000.*pow(10.,3.);
  int jn = 30000;

  const int chCt = 6;
  int wCh[chCt] =    { 112, 			 114,			    118, 
  				       144, 			 146, 			    148 };
  //approximate channel where the pulser is located for each detector
  double pX[chCt] =  { 280.*pow(10.,3.), 1344.*pow(10.,3.), 1393.*pow(10.,3.),
  				       37.5*pow(10.,3.), 1350.*pow(10.,3.), 1361.*pow(10.,3.) };
  //approximate channel where the "pulser overshoot" is located for each detector
  double poX[chCt] = {  22.*pow(10.,3.), 108.*pow(10.,3.),  116.*pow(10.,3.),
  					   3.15*pow(10.,3.), 110.*pow(10.,3.),  113.*pow(10.,3.) };

  int whichCh=999;
  for(int i=0; i<chCt; i++) if(wCh[i]==channel) whichCh=i;
  if(whichCh==999){
    printf("Invalid Channel. Only can use:");
    for(int i=0; i<(chCt-1); i++) printf(" %d,",wCh[i]);
    printf(" %d\n",wCh[chCt-1]);
    return 1;
    }
  
  double pL = pX[whichCh] - 3.*pow(10.,3.);
  double pU =  pX[whichCh] + 3.*pow(10.,3.);
  
  double poL = poX[whichCh] - 3.*pow(10.,3.);
  double poU = poX[whichCh] + 3.*pow(10.,3.);
 
  int binL, binU;
  double runningAvg;

  TChain chain("mjdTree");
  for(int i=0; i<ct; i++){
    chain.Add( Form("${MJDDATADIR}/surfprot/data/gatified/%s/mjd_run%.0f.root",partID,x_val[i]+base) );
    runningAvg = 0.;
    f = TFile::Open( Form("${MJDDATADIR}/surfprot/data/gatified/%s/mjd_run%.0f.root",partID,x_val[i]+base) );
    t = (TTree*)f->Get("mjdTree");
    t->Draw( Form("energy>>htemp(%d, %f, %f)",jn,jmin,jmax), Form("((channel==%d)&&(timestamp<(300.*pow(10.,8.))))",channel), "goff");
    htemp = (TH1D*)gDirectory->Get("htemp");
    pulserOverArea[i] = htemp->Integral( htemp->FindBin(poL),htemp->FindBin(poU) );
    binL = htemp->FindBin(pL);
    binU = htemp->FindBin(pU);
    pulserArea[i] = htemp->Integral( binL, binU );
    for(int j=binL; j<=binU; j++) runningAvg += (htemp->GetBinContent(j)*htemp->GetXaxis()->GetBinCenter(j));
    pulserAvg[i] = runningAvg/pulserArea[i];
    f->Close();
  	}
      
  TFile *newFile = new TFile( Form("PulserStability_Runs_%d_%d.root", start_run_no, end_run_no),"update" );
  TGraph **g = new TGraph*[3];
  char gDesc[3][50];
  sprintf(gDesc[0],"pulserAvg_Ch%d",channel); 
  sprintf(gDesc[1],"pulserArea_Ch%d",channel);
  sprintf(gDesc[2],"pulserOvershootArea_Ch%d",channel);
  TKey **gKey = new TKey*[3];
  
  for(int i=0; i<3; i++){
    gKey[i] = (TKey*)newFile->FindKey(gDesc[i]);
    if(gKey[i]!=0){
      gKey[i]->Delete();
      printf("Removed TGraph %s and will re-write\n",gDesc[i]);
      }
    if(i==0) g[0] = new TGraph(ct,x_val,pulserAvg);
    if(i==1) g[1] = new TGraph(ct,x_val,pulserArea);
    if(i==2) g[2] = new TGraph(ct,x_val,pulserOverArea);
    g[i]->SetMarkerStyle(23);
    g[i]->Write(gDesc[i]);
	}
	
  chain.Draw( Form("energy>>htemp(%d, %f, %f)",jn,jmin,jmax), Form("((channel==%d)&&(timestamp<(300.*pow(10.,8.))))",channel), "goff");
  htemp = (TH1D*)gDirectory->Get("htemp");
  htemp->Write( Form("spectrum_Ch%d",channel) );
  
  newFile->Close();
  
  return 0;
}

int main(int argc, char *argv[]){

  int oVal, index;
  int channel, acceptedCount;
  while( (oVal=getopt(argc, argv,"c:n:")) != -1 )
    switch (oVal){
      case 'c':
        channel = atoi(optarg);
        break;
      case 'n':
        acceptedCount = atoi(optarg);
        break;
      case '?':
        if((optopt == 'c') || (optopt == 'n')){
          printf("Option -%c requires an argument. Usage:\n", optopt);
          printf(" -c [channel]\n -n [number of accepted runs]\n");
          }
        else if(isprint(optopt))
          printf("Unknown option '-%c'\n",optopt);
        else
          printf("Unknown option character '\\x%x'\n",optopt);
        return 1;
      default:
        abort ();
      }
  for(index=optind; index<argc; index++)
	printf("Non-option argument %s\n",argv[index]);

  printf("channel = %d\n",channel);
  printf("acceptedCount = %d\n",acceptedCount);
  
  int returnVal = lookForGainDrift(channel,acceptedCount);
  
  return returnVal;
  }
