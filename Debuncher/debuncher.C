
// std includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include "math.h"
#include <fstream>

//  root includes
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2.h"


//  fasterac includes
#include "fasterac/fasterac.h"
#include "fasterac/fast_data.h"
#include "fasterac/qdc_caras.h"
#include "fasterac/qtdc.h"

using namespace std;
char filefast1[100];
float ToF = 0.015; // [ms] Time of flight from a debuncher to the MCP
bool firstEvent;
bool EOF1, EOF2;
int counterwhile, eventCounter, trigCounter, IonSum;
Float_t mean, rms, timestamp, timestampRef, charge;


TString fName1;
TH1I *bunchHist, *bunchSumHist;
TH2F *bunchStatHist;

  //  faster file reader
  faster_file_reader_p   reader;
  //  faster data
  faster_data_p          data;
  unsigned char          alias;
  unsigned short         label;
Double_t		clockt;
//  double               clock_a;
//  group data
//  unsigned short         lsize;

  //  qdc tdc data  (from faster group)
//  qdc_t_x1               qdc1;
//  qdc_t_x1               qdc3;
//  qdc_x1				 qdcPM4;

  //  qtdc
  //qtdc                   qt1;
  //qtdc_counter           qtdc_cnt1;
  	qdc_t_x1    q_tdc;
	qdc_counter q_count;



/////////////////////////////////////////////////////
////////////////// MAIN ////////////////////////////
///////////////////////////////////////////////////

int main (int argc, char** argv) {

      reader = faster_file_reader_open (argv[1]);
	float histLeft=0.001 ; // in ms
	float histRight=10. ;// in ms
// take these two from the argv[]
if (argc<4) { cout<<" The use of debuncher: debuncher <file path/file> <start debunch time [ms]> <stop debunch time [ms]>"<<endl; 
			return EXIT_FAILURE;
}

	histLeft = atof(argv[2]);
	histRight= atof(argv[3]);

bunchHist = new TH1I("Bunch", "Bunch", 1000, histLeft , histRight);
bunchSumHist = new TH1I("BunchTotal", "BunchTotal", 1000, histLeft , histRight);
bunchStatHist = new TH2F("mean-RMS","mean-RMS", 1000, histLeft , histRight, 1000, histLeft, histRight);

cout<<"b:"<<flush<<endl;
cout<<"reader:"<<reader<<flush<<endl;

      if (reader == NULL) {
        reader = faster_file_reader_open ("");
cout<<"reader==NULL:"<<flush<<endl;
            if (reader == NULL) {
                printf ("error opening file %s\n",argv[1]);
                printf ("error opening file %s\n", "");
                return EXIT_FAILURE;
            }
      }
      else strcpy(filefast1, argv[1]);

cout<<"c:"<<flush<<endl;

      fName1 =  filefast1;
      fName1.ReplaceAll (".fast", ".root");

		TFile *f = new TFile (fName1, "recreate"); // root file

cout<<"d:"<<flush<<endl;
		ofstream myfile; // textual analysis file
		string str1(filefast1);
		str1.erase(str1.end()-5, str1.end());
		str1.append("_analis.txt");
		const char * c1 = str1.c_str();
		myfile.open (c1); // fstream doesn't accept string as argument

		myfile<<"Ev.No \t Number of ev \t mean \t RMS "<<endl;





	eventCounter=trigCounter=IonSum=0;

cout<<"1:"<<flush<<endl;
      EOF1=((data = faster_file_reader_next (reader)) == NULL) ;

      if (!EOF1){
        alias = faster_data_type_alias (data);  
        label = faster_data_label      (data);
        clockt = faster_data_clock_sec  (data);
      }

//cout<<" before while "<<endl<<flush;
  // main loop
while (EOF1 == false) 
{                    //  read each data
//cout<<"while:"<<flush<<endl;
    counterwhile++;
	if (counterwhile==int(float(counterwhile)/1000)*1000) cout<<"counter "<<counterwhile<<endl;

    if ((alias == QDC_TDC_X1_TYPE_ALIAS)&&(label==1))
	{
		trigCounter++;
		if (eventCounter>0) {
			mean=bunchHist->GetMean();
			rms =bunchHist->GetRMS();
			myfile<<trigCounter<<"\t"<<eventCounter<<"\t"<<mean<<"\t"<<rms<<endl;
			
			bunchStatHist->Fill(mean,rms);
			
			bunchHist->Reset();
			IonSum+=eventCounter;
			
		}
      	faster_data_load (data, &q_tdc);
      	timestampRef = clockt * 1e9 + qdc_conv_dt_ns (q_tdc.tdc);   // long double
		charge = qdc_conv_q_mVns (q_tdc.q1); // in mV.ns
//cout<<"alias: "<<"QDC_TDC_X1"<<", label:"<<label<<", clockt:"<<clockt<<", timestamp:"<<timestampRef<<endl;

		eventCounter=0; // reset the number of events - counting started
		firstEvent=true; // after this we have a first event

//    if (q_tdc.q1_saturated) printf ("        saturated : q1");
//      printf ("\n");
	}





    if ((alias == QDC_TDC_X1_TYPE_ALIAS)&&(label==2))
	{
      	faster_data_load (data, &q_tdc);
      	timestamp = clockt * 1e9 + qdc_conv_dt_ns (q_tdc.tdc);   // long double
//cout<<"alias: "<<"QDC_TDC_X1"<<", label:"<<label<<", clockt:"<<clockt<<", timestamp:"<<timestamp<<endl;

if (((timestamp/1e6-timestampRef/1e6)>(histLeft+histRight+ToF))||((timestamp/1e6-timestampRef/1e6)<(histLeft+ToF))) // read new event and go to the begining of while loop, if time greater than the traping time+ejecttime+Tof or if it is smaller than traping time + Tof
		{
      		EOF1=((data = faster_file_reader_next (reader)) == NULL) ;

      		if (!EOF1){
        		alias = faster_data_type_alias (data);  
        		label = faster_data_label      (data);
        		clockt = faster_data_clock_sec  (data);
      		}
			continue; 
		}

		firstEvent=false;
		charge = qdc_conv_q_mVns (q_tdc.q1); // in mV.ns
		eventCounter++;
//    if (q_tdc.q1_saturated) printf ("        saturated : q1");
//      printf ("\n");
	}

/***************************

CONDITION TO ACCEPT ONLY CERTAIN EVENT REGARDING THE DISTANCE FROM THE timestampRef

*************************** */

// fill histogram of current bunch (or total bunch)
if ((alias == QDC_TDC_X1_TYPE_ALIAS)&&(label==2)) // avoid to write counters and the trigger
	{
	bunchHist->Fill((timestamp-timestampRef)/1e6);
	bunchSumHist->Fill((timestamp-timestampRef)/1e6);
cout<<(timestamp-timestampRef)/1e6<<endl;
}

      EOF1=((data = faster_file_reader_next (reader)) == NULL) ;

      if (!EOF1){
        alias = faster_data_type_alias (data);  
        label = faster_data_label      (data);
        clockt = faster_data_clock_sec  (data);
      }



}// end while

myfile<<"Average number of events per trigger: "<<float(IonSum)/float(trigCounter)<<endl;
bunchHist->Write();
bunchSumHist->Write();
bunchStatHist->Write();

f->Write();
f->Close();
myfile.close();
}//end mail




     
