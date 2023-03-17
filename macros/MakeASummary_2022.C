#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TTree.h"
#include "TROOT.h"
#include "TFile.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TBranch.h"
#include "TChain.h"
#include "TUnixSystem.h"
#include "TCanvas.h"
#include "TGraph.h"


// F u n c t i o n   P r o t o t y p e   
// 

Double_t expFunc(Double_t*, Double_t*);//Function prototype 
Double_t PolyFuncFit(Double_t*, Double_t*);//Function prototype 
Double_t PolyFunccalc(Double_t, Double_t, Double_t , Double_t, Double_t );//Function prototype 
Double_t expFunccalc(double, double, double);//Function prototype 
Double_t SigmoidFunc( Double_t*, Double_t*);//Function prototype 
Double_t Sigmoidcalc(Double_t,Double_t,Double_t,Double_t);//Function prototype 
Double_t difcalc(Double_t, Double_t, Double_t, Double_t);//Function prototype 
std::vector<std::string> blacklist(const char*);
std::vector< std::pair<std::string,std::string> > dictionary(const char*);// function prototype
std::vector< std::pair<std::string, double> > WPchannel(const char*, const char*);


// ----- F U N C T I O N S  T O   B E   I M P L E M E N T E D   I N   T H E   F I T 
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Double_t expFunc(Double_t* _x, Double_t* _par){
  return  TMath::Exp(_par[0]+_x[0]*_par[1]);
}
Double_t expFunccalc(double hv, double A, double B){
   return TMath::Exp(A + hv*B);
}
Double_t PolyFuncFit(Double_t *hv, Double_t *par){
   Double_t clsz = par[0] + par[1]*(hv[0]) + par[2]*((2*hv[0]*hv[0])-1) + par[3]*(4*(hv[0]*hv[0]*hv[0]) - 3*(hv[0]));
   return clsz;
}

Double_t PolyFunccalc(Double_t wp, Double_t a, Double_t b, Double_t c, Double_t d){
   Double_t clsz = a + b*(wp)+ c*((2*wp*wp)-1)+ d*(4*(wp*wp*wp)-3*wp) ;
   return clsz;
}

Double_t SigmoidFunc( Double_t* _x, Double_t* _par ){
  Double_t effmax = _par[0];
  Double_t S = _par[1];
  Double_t HV50 = _par[2];
  return effmax / (1.0 + TMath::Exp( S *( _x[0] - HV50 ) ) );//
}

Double_t Sigmoidcalc(Double_t hv,Double_t emax,Double_t S,Double_t hv50 ){
  return emax / (1.0 + TMath::Exp( S *( hv - hv50 ) ) );
}

Double_t difcalc(Double_t hv, Double_t emax, Double_t S, Double_t hv50){
  return -emax*S*TMath::Exp(S*(hv-hv50))/((1.0 + TMath::Exp( S *( hv - hv50 ) ) )*(1.0 + TMath::Exp( S *( hv - hv50 ) ) )) ;
}


void MakeASummary_2022(const char* subdetect, bool BL=false){

  std::ifstream fResEff, fResCls;
  std::vector< std::pair<std::string,std::string> > map;
  vector<std::string> blacklistrollv;
  std::vector< std::pair<std::string, double> > channel_map; 
  bool match = false;


  map = dictionary(subdetect);
  std::cout << map.size()<< std::endl;
  
  if( strcmp(subdetect,"endcap") == 0 )channel_map = WPchannel("endcap", "../data/EndcapChambersName_WP.txt") ; 
  else  channel_map = WPchannel("barrel","../data/BarrelChambersName_WP.txt") ; 

  if (BL) blacklistrollv = blacklist("../data/blacklist_2018.txt");
 
 

  Char_t   RollName[38];
  Char_t   Id[10];
  Double_t WorkingPoint;
  Double_t slope50, emax, hv50, chi2, EffWP, clsWP, chi2cls;
  Double_t emaxErr, slopeErr, hv50Err;
  Double_t slope;
  Double_t WPch, EffWPch;  
  
  TString outputFile;
  if ( (strcmp(subdetect,"barrel")==0) && BL==true )outputFile  = "../summary/barrel_summary_2018BlackList.root";
  else if ((strcmp(subdetect,"endcap")==0) && BL==true ) outputFile = "../summary/endcap_summary_2018BlackList.root";
  else if ((strcmp(subdetect,"barrel")==0) && BL==false )outputFile  = "../summary/barrel_summary_2022.root";
  else if ((strcmp(subdetect,"endcap")==0) && BL==false )outputFile  = "../summary/endcap_summary_2022.root";
  else{}; 

  Int_t nevents = 0;
  TFile *file = new TFile(outputFile,"RECREATE");
  TTree *filtered = new TTree("filtered","summary");
  TTree *removed; 
  if (BL) removed = new TTree("removed","summary");
    

  //To save the filtered 
  filtered->Branch("RollName",&RollName,"RollName/C");
  filtered->Branch("Id",&Id,"Id/C");
  filtered->Branch("WorkingPoint",&WorkingPoint,"WorkingPoint/D");
  filtered->Branch("emax",&emax,"emax/D");
  filtered->Branch("slope",&slope,"slope/D");
  filtered->Branch("hv50",&hv50,"hv50/D");
  filtered->Branch("chi2",&chi2,"chi2/D");
  filtered->Branch("slope50",&slope50,"slope50/D");
  filtered->Branch("EffWP",&EffWP,"EffWP/D");
  filtered->Branch("clsWP",&clsWP,"clsWP/D");
  filtered->Branch("WPch",&WPch,"WPch/D");
  filtered->Branch("EffWPch",&EffWPch,"EffWPch/D");
  filtered->Branch("chi2cls",&chi2cls,"chi2cls/D");
  filtered->Branch("emaxErr",&emaxErr,"emaxErr/D");
  filtered->Branch("slopeErr",&slopeErr,"slopeErr/D");
  filtered->Branch("hv50Err",&hv50Err,"hv50Err/D");

   // Only if the Black list is applied, save only the blacklisted rolls 
  if (BL){
  removed->Branch("RollName",&RollName,"RollName/C");
  removed->Branch("Id",&Id,"Id/C");
  removed->Branch("WorkingPoint",&WorkingPoint,"WorkingPoint/D");
  removed->Branch("emax",&emax,"emax/D");
  removed->Branch("slope",&slope,"slope/D");
  removed->Branch("hv50",&hv50,"hv50/D");
  removed->Branch("chi2",&chi2,"chi2/D");
  removed->Branch("slope50",&slope50,"slope50/D");
  removed->Branch("EffWP",&EffWP,"EffWP/D");
  removed->Branch("clsWP",&clsWP,"clsWP/D");
  removed->Branch("WPch",&WPch,"WPch/D");
  removed->Branch("EffWPch",&EffWPch,"EffWPch/D");
  removed->Branch("chi2cls",&chi2cls,"chi2cls/D");
  removed->Branch("emaxErr",&emaxErr,"emaxErr/D");
  removed->Branch("slopeErr",&slopeErr,"slopeErr/D");
  removed->Branch("hv50Err",&hv50Err,"hv50Err/D");
  }else ; 


  std::string id_;
  Double_t FitResEff[10];
  Double_t FitResCls[7];
  vector<std::pair<std::string, double> >::const_iterator itchannel;  
  vector<std::pair<std::string,std::string> >::const_iterator itmap;  
  for ( itmap = map.begin() ;itmap != map.end(); itmap++  ){
  	 id_ =itmap->first;
         fResEff.open(("../results_scan2022/"+id_+"/fitData.txt").c_str());
         while (1){
          fResEff >>FitResEff[0]//wp
                  >>FitResEff[1]//slope50
                  >>FitResEff[2]//emax
                  >>FitResEff[3]//hv50
                  >>FitResEff[4]//chi2
                  >>FitResEff[5]//effwp
                  >>FitResEff[6]//emaxerr
                  >>FitResEff[7]//slopeerr
                  >>FitResEff[8]//hv50err
                  >>FitResEff[9];//slope
                  if (fResEff.eof())break;
                  }
          fResEff.close();
          fResCls.open(("../results_scan2022/"+id_+"/fitDataCls.txt").c_str());
          while (1){
          fResCls >>FitResCls[0]//a
                  >>FitResCls[1]//b
                  >>FitResCls[2]//c
                  >>FitResCls[3]//d
                  >>FitResCls[4]//chi2cls
                  >>FitResCls[5]//wp
                  >>FitResCls[6];//clswp
                  if (fResEff.eof())break;
          }
          fResCls.close();
	
  	  /*TF1 *f1 = new TF1("f1",SigmoidFunc, 8.5, 9.9  ,3);//range of the function and how many parameters  
 	  f1->SetParNames("emax","slope","hv50");
 	  f1->SetParameter(0, FitResEff[2]);
 	  f1->SetParameter(1, FitResEff[9]);//
 	  f1->SetParameter(2, FitResEff[3]);
	   */

          if (BL){
		  match=false;
                  for (int i=0; i <int(blacklistrollv.size()); i++)if(id_ == blacklistrollv.at(i))match=true; 
		  }
          
          strcpy(RollName, (itmap->first).c_str());
          strcpy(Id, (itmap->second).c_str());
	  WorkingPoint 	= FitResEff[0]; 
          slope50 	= FitResEff[1]; 
	  emax 		= FitResEff[2]; 
	  hv50 		= FitResEff[3]; 
	  chi2 		= FitResEff[4]; 
	  EffWP 	= FitResEff[5]; 
	  emaxErr 	= FitResEff[6]; 
	  slopeErr 	= FitResEff[7];
	  hv50Err 	= FitResEff[8];
	  slope 	= FitResEff[9];
	  clsWP 	= FitResCls[6];
	  chi2cls	= FitResCls[4]; 
          for (itchannel = channel_map.begin() ; itchannel != channel_map.end(); itchannel++ ){
                 if(strcmp(subdetect,"barrel")==0 ){
			if ( (itmap->first).find(itchannel->first) != std::string::npos){
			WPch = ((itchannel->second)/1000.0);
	  		EffWPch = Sigmoidcalc( WPch , emax, slope, hv50); 
			std::cout << itchannel->first << "  " <<itmap->first << "  " << WPch << "  " << EffWPch << "  " <<  std::endl;
			continue;   
		   }else continue; //std::cout << itchannel->first << "  " <<itmap->first << endl;;   			
	       }
	        else {
			std::string roll1 = (itchannel->first).substr(0,12);
			std::string roll2 = (itchannel->first).substr(13,12);
			//std::cout << roll2 << "  " << roll1 << endl;
			if ( (itmap->first).find(roll1) != std::string::npos || (itmap->first).find(roll2) != std::string::npos ){
			WPch = (itchannel->second/1000.0);
	  		EffWPch = Sigmoidcalc( WPch ,FitResEff[2], FitResEff[9],FitResEff[3]); 
			std::cout << itchannel->first << "  " <<itmap->first << "  " << WPch << "  " << EffWPch << "  " << std::endl; 
			continue;   
		   }else continue; //std::cout << itchannel->first << "  " <<itmap->first << endl;;   			
	      }

	  }
		    
          if (match){
		  removed->Fill();	
		  }
	  else {
		  nevents++; 
		  filtered->Fill();
		  } 
	}
        std::cout << "rolls->  " <<  nevents << std::endl; 
	filtered->Write("",TObject::kOverwrite);
	if (BL)removed->Write("",TObject::kOverwrite);
	file->Write("",TObject::kOverwrite);
	file->Close();
}


std::vector<std::string> blacklist(const char* filename){
  ifstream infile; 
  
  infile.open(filename); 
  std::vector<std::string> rolls; 
  while (1){
   	std::string name; 
   	infile >> name; 
   	if (infile.eof())break;
   	rolls.push_back(name); 
   }   
  
  return rolls;   
}
// To do a map of the rolls. Every roll has an id match a name   
std::vector< std::pair<std::string,std::string> > dictionary(const char* subd){
  ifstream f;
  Double_t id;
  std::string name;
  if (strcmp(subd,"endcap")==0)f.open("../data/detIdEndCapScan2022.txt");
  else f.open("../data/detIdBarrelScan2022.txt");
  std::vector< std::pair<std::string,std::string> > dictionary_;
  while (1) {
        f >> name
          >> id;
        if (f.eof()) break;
        stringstream s; std::string stid;
        s << UInt_t(id);
        stid=s.str();
        std::pair<std::string,std::string > pair_map;
        pair_map.first = name ;
        pair_map.second = stid;
        dictionary_.push_back(pair_map);
   }
   f.close();
   return dictionary_ ;
}

std::vector< std::pair<std::string, double> > WPchannel(const char* subd, const char* filename){

 ifstream infile; 
 std::pair<std::string, double > pair_map ; 
 std::vector< std::pair<std::string, double> > WPchannel_ ; 
 infile.open(filename); 
 while(infile.good() ){
 	std::string  roll1, roll2; 
	double HVvalues=0;      
 	if (strcmp(subd,"endcap")==0){ 
	   infile >> roll1 >> roll2 >> HVvalues; 
	   if (infile.eof()) break;
	   pair_map.first = roll1+" "+roll2; 
	   pair_map.second = HVvalues; 
	   std::cout << pair_map.first << std::endl; 
	  }	

	else { 
	   infile >> roll1  >> HVvalues; 
	   if (infile.eof()) break;
	   pair_map.first = roll1; 
	   pair_map.second = HVvalues; 
	}		
	WPchannel_.push_back(pair_map);
        	
   }
 infile.close();
 
 return WPchannel_; 

}

