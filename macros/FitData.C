#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
//Root 
#include <TMath.h>
#include <TTree.h>
#include "TROOT.h"
#include "TFile.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TBranch.h"
#include "TChain.h"
#include "TUnixSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

/*
class Chamber : public TObject {
public:
  std::string   id;    //! not persistent
  Float_t  *eff;     //[run]  
  Float_t  *efferr;      //[run] 
  Float_t  *numext;     //[run]  
  Float_t   *cls;    //[run]
  Chamber() {
    *eff=0;     
    *efferr=0;      
    *numext=0;     
    *cls=0;  
  }

  ClassDef(Chamber,1)
};
*/

//P r o t o t y p e 
//
Double_t expFunc(Double_t*, Double_t*);//Function prototype 
Double_t PolyFuncFit(Double_t*, Double_t*);//Function prototype 
Double_t PolyFunccalc(Double_t, Double_t, Double_t, Double_t , Double_t );//Function prototype 
Double_t expFunccalc(double, double, double);//Function prototype 
Double_t SigmoidFunc( Double_t*, Double_t*);//Function prototype 
Double_t Sigmoidcalc(double,double,double,double);//Function prototype 
Double_t difcalc(double, double, double, double);//Function prototype 
std::vector< std::pair<Int_t,Float_t> > hvEff(const char* );//Function prototype 
std::vector< std::pair<std::string,std::string> > dictionary(const char*);//Function prototype 
void rData(const char *);//Function prototype
Double_t FitDataFuncEff(const char*,const char*, int, float*,float*, float*, float*);//Function prototype
void FitDataFuncCls(const char*,const char*, int, float *,float *, float *, float*, Double_t);



// ----- F U N C T I O N S  T O   B E   I M P L E M E N T E D   I N   T H E   F I T 
// - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Double_t expFunc(Double_t* _x, Double_t* _par){
  return  TMath::Exp(_par[0]+_x[0]*_par[1]);
}
Double_t expFunccalc(double hv, double A, double B){
   return TMath::Exp(A + hv*B);
}
Double_t PolyFuncFit(Double_t *hv, Double_t *par){
   Double_t clsz = par[0] + par[1]*(hv[0])+ par[2]*((2*hv[0]*hv[0])-1) + par[3]*(4*(hv[0]*hv[0]*hv[0]) - 3*(hv[0]));
   return clsz;
}

Double_t PolyFunccalc(Double_t wp, Double_t a, Double_t b, Double_t c,  Double_t d){
   Double_t clsz = a + b*(wp)+ c*((2*wp*wp)-1)+ d*(4*(wp*wp*wp)-3*wp) ;
   return clsz;
}

Double_t SigmoidFunc( Double_t* _x, Double_t* _par ){
  Double_t effmax = _par[0];
  Double_t S = _par[1];
  Double_t HV50 = _par[2];
  return effmax / (1.0 + TMath::Exp( S *( _x[0] - HV50 ) ) );//
}

Double_t Sigmoidcalc(double hv,double emax,double S,double hv50 ){
  return emax / (1.0 + TMath::Exp( S *( hv - hv50 ) ) );
}

Double_t difcalc(double hv, double emax, double S, double hv50){
  return -emax*S*TMath::Exp(S*(hv-hv50))/((1.0 + TMath::Exp( S *( hv - hv50 ) ) )*(1.0 + TMath::Exp( S *( hv - hv50 ) ) )) ;
}


//----F U N C T I O N S   T O   D A T A   H A N D L I N G ---------------
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

// Get the number of  runs avalaible and their HV eff   
std::vector < std::pair<Int_t,Float_t> > hvEff(const char* subd){
	Int_t run=0;
	Float_t hv_b=0;
	Float_t hv_ec=0;
	vector<Float_t> hv;  
	ifstream f;
	std::vector< std::pair<Int_t,Float_t> > hvEff_;
	f.open("../data/hvEffective.txt");
	while (1) {
		f >> hv_b >> hv_ec;
		if (f.eof()) break;
		if (strcmp(subd,"barrel")==0)hv.push_back(hv_b);
		else hv.push_back(hv_ec);
	}
	f.close();
        hv.erase( unique( hv.begin(), hv.end() ), hv.end() );
	cout << hv.size() << endl;
	for (unsigned int i=0; i<hv.size();i++){
		std::pair<Int_t, Float_t > pair_hv;
		if ( hv.at(i)== 0 )continue;
		run++;
                pair_hv.first = run; 
	        pair_hv.second = hv.at(i);
		hvEff_.push_back(pair_hv); 
	}
        hv.clear();
	cout << "runs available HV EFF" << endl;
	return hvEff_ ;
}

// To do a map of the rolls. Every roll has an id match a name   
std::vector< std::pair<std::string,std::string> > dictionary(const char* subd){
  ifstream f;
  Double_t id;
  cout << "dictionary 1" << endl;
  std::string name;
  if (strcmp(subd,"endcap")==0)f.open("../data/detIdEndCapScan2022.txt");
  else f.open("../data/detIdBarrelScan2022.txt");
  std::vector< std::pair<std::string,std::string> > dictionary_;
  cout << "dictionary 2" << endl;  
  while (1) {
        f >> name >> id;
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
   cout << "Map of the rolls" << endl;
   return dictionary_ ;    
}

// R E A D  and  S A V E   T H E   I N P U T   R O O T   F I L E S   D A T A   
void rData(const char *subdetect){
   using namespace std;
   //Define how many runs are available and their high voltage, reading the hvEffective.txt file.  
   std::vector<Int_t> samples; 
   std::vector< std::pair<Int_t,Float_t> > hvscan;
   hvscan = hvEff(subdetect);
   std::cout << "Number of points were taken: " <<hvscan.size() << std::endl;
   for (vector< std::pair<Int_t,Float_t> >::const_iterator it = hvscan.begin() ;it != hvscan.end(); it++  ){
   	if (strcmp(subdetect, "endcap")==0)samples.push_back(-(it->first));
   	else if (strcmp(subdetect,"barrel")==0)samples.push_back(it->first);
   	else {std::cout << "put \"barrel\" or \"endcap\" only" << std::endl;exit(1);}
	std::cout << subdetect <<" High Voltage effective:  "<< it->first << " " << it->second << " [kV] "<<   std::endl;
   } 
  
   // G e t   t h e   D a t a
   //- Define the variables of the tree 
   Float_t fiducialCutExp;
   Float_t fiducialCutEff;
   Float_t fiducialCutEffErr;
   Float_t clustersize;
   Int_t detId;
   std::ofstream RollEff; 
   //std::string dirFile ="/eos/cms/store/group/dpg_rpc/comm_rpc/Run-II/data2017/HVscan_Offline/HVScan_2017May_selPt7_chi8" ;  
   //std::string dirFile =" /eos/cms/store/group/dpg_rpc/comm_rpc/Run-II/data2018/HVscan_April";  
  // std::string dirFile ="/eos/cms/store/group/dpg_rpc/comm_rpc/Run-II/data2018/HVscan_April_Corrected/secondStep";
  //std::string dirFile ="/eos/cms/store/group/dpg_rpc/comm_rpc/Run-II/data2017/HVscan_Offline/HVScan_2017September_selPt7_chi8";
   std::string dirFile = "/eos/cms/store/group/dpg_rpc/comm_rpc/Run-III/HVscan2022/Input_files/";
   for (std::vector<int>::const_iterator Sample = samples.begin(); Sample != samples.end(); Sample++) {	
                gROOT->Reset();
                Int_t Run=0;
                if (strcmp(subdetect,"endcap")==0)Run = -1*(*Sample);
                else Run=(*Sample);	
                std::stringstream ss; std::string strun; ss << Run; 
                strun = ss.str();  
                TChain *ch = new TChain("ch","");
		if (*Sample == -1) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P86_RPCMon2022.root/effTreeEndcap").c_str());
 		if (*Sample == -2) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P88_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -3) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P90_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -4) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P91_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -5) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P92_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -6) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P93_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -7) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P94_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -8) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P95_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -9) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P96_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -10) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P97_RPCMon2022.root/effTreeEndcap").c_str());
  		if (*Sample == -11) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P98_RPCMon2022.root/effTreeEndcap").c_str());
		if (*Sample == 1) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P86_RPCMon2022.root/effTreeBarrel").c_str());
 		if (*Sample == 2) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P88_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 3) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P90_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 4) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P91_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 5) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P92_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 6) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P93_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 7) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P94_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 8) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P95_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 9) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P96_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 10) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P97_RPCMon2022.root/effTreeBarrel").c_str());
  		if (*Sample == 11) ch->Add((dirFile+"/AnalyzeEfficiency_HVscan_P98_RPCMon2022.root/effTreeBarrel").c_str());
   		TTree *t = (TTree*)ch;
        	t->SetBranchAddress("detId",&detId);
		t->SetBranchAddress("fiducialCutEff",&fiducialCutEff);
		t->SetBranchAddress("fiducialCutEffErr",&fiducialCutEffErr);
		t->SetBranchAddress("fiducialCutExp",&fiducialCutExp);
		t->SetBranchAddress("clustersize",&clustersize);
        cout << "passed" << endl;
      	   	Int_t ientries = (Int_t)t->GetEntries();
                std::cout <<" number of chambers in the detector " << ientries << std::endl;
                //Save  the data in a txt File 
                if (strcmp(subdetect,"endcap")==0)RollEff.open(("../data/rollEff_"+strun+"_ec.txt").c_str());
                else RollEff.open(("../data/rollEff_"+strun+"_b.txt").c_str());
                for (Long64_t i=0; i<ientries; i++) {
 			t->GetEntry(i);
			RollEff << UInt_t(detId) 	 << " "  
		                << 100*fiducialCutEff 	 << " "
		                << 100*fiducialCutEffErr << " "
		        	<< fiducialCutExp 	 << " "
		                << clustersize		 << std::endl;
        	}

   	       	RollEff.close();
			cout << "closed RollEff" << endl; 
   }
   cout << " end rData " << endl;
   return; 
}
 Double_t FitDataFuncEff(const char* subdetect,const char* chamber, int points, float *hv,float *hverr, float *eff, float *efferr){
 if (points ==0.) {std::cout << "N O  D A T A "; return 0;} 
 std::cout << subdetect <<":  "<< chamber << std::endl; 

 std::cout << "Look here -1" << endl;


 //Defining the fit function and setting parameters 
 TF1 *f1 = new TF1("f1",SigmoidFunc, 8.4, 9.9999  ,3);//Domain of the function and how many parameters  
 f1->SetParNames("emax","slope","hv50");
 f1->SetParLimits(0, 0.0, 99.9999);//bound emax parameter 
 f1->SetParameter(0, 99.99);//setting at value eff 95%
 f1->SetParLimits(2, 0.0, 97.999);//bound emax parameter 
 //f1->SetParameter(0, 0.00);//setting at value eff 95%
 if (strcmp(subdetect,"endcap")==0){
      f1->SetParameter(1, -12);//Setting to achieve convergence 
      f1->SetParameter(2, 9.3);}
 else{
      f1->SetParameter(1, -8.5);//Setting to achieve convergence
      f1->SetParameter(2, 8.9);}

 std::cout << "Look here 0" << endl;

 // Being made the fit on the efficiency vs HV distribution  
 TGraphErrors *hveff = new TGraphErrors(points, hv, eff, hverr, efferr);
 hveff->Fit(f1,"RN");
 // Get the results of the Fit
 Double_t wp=0., effwp=0., effknee=0.0001, knee=0., effkneemin=0., slope50=0.,slopeerr=0., hv50err=0., chi2=0.;
 Double_t resolution=0.01;
 chi2=(f1->GetChisquare())/points; 
 
  knee=f1->GetParameter(2);//start at HV50  
 
 effkneemin = (f1->GetParameter(0))*0.95; //f1->GetParameter(0) = emax
 
 std::cout << "Look here 1" << endl;

 if ((f1->GetParameter(0)>0)&&(f1->GetParameter(1)<0)){
    while (effknee < effkneemin){
           knee=knee+resolution;
           // reaching  95% efficiency  
           effknee=Sigmoidcalc(knee,f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
    }
 }
 else; 
 
  
 if ((f1->GetParameter(0))*100.0 < 0.0001  ) {knee = 0.0; wp = 0.00; }
 else {
 	if (strcmp(subdetect,"endcap")==0 )wp = knee +0.120;
 	else wp = knee +0.100;
	}
 effwp = Sigmoidcalc(wp,f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2));
 slope50 = difcalc(f1->GetParameter(2),f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2)); 
 
 std::cout << "Look here 2" << endl;


 //S A V E   T H E   R E S U L T S   I N   A   T X T   F I L E  
 ofstream fitResEff;
 std::string chamber_ = chamber;  
 gSystem->mkdir(("../results_scan2022/"+chamber_).c_str());
 
 fitResEff.open(("../results_scan2022/"+chamber_+"/fitData.txt").c_str());  
 cout << wp << " "
           <<slope50<<" "
           <<f1->GetParameter(0)   /*emax*/
           <<" "<<f1->GetParameter(2) /*hv50*/
           <<" "<<chi2<<" "   /*chi2 eff*/
           <<effwp<<" "       /*eff wp*/
           <<f1->GetParError(0)<<" " /* emax err */
           <<f1->GetParError(1)<<" "  /*slopeerr*/
           <<f1->GetParError(2)<<" "  /*hv50 err*/
           <<f1->GetParameter(1)<<endl;  /*slope*/

 fitResEff << wp <<" "
           <<slope50<<" "
           <<f1->GetParameter(0)   /*emax*/
           <<" "<<f1->GetParameter(2) /*hv50*/
           <<" "<<chi2<<" "   /*chi2 eff*/
           <<effwp<<" "       /*eff wp*/
           <<f1->GetParError(0)<<" " /* emax err */
           <<f1->GetParError(1)<<" "  /*slopeerr*/ 
           <<f1->GetParError(2)<<" "  /*hv50 err*/
           <<f1->GetParameter(1)<<endl;  /*slope*/
 fitResEff.close();
 
 std::cout << "Look here 3" << endl;
 
return wp; 
}
 
 void FitDataFuncCls(const char* subdetect,const char* chamber, int points, float *hv,float *hverr, float *cls, float *clserr, Double_t wp){
 gROOT->Reset();
 if (points ==0.) {std::cout << "N O  D A T A "; return;} 
 Double_t clswp=0.,clsknee=0, chi2cls=0.;  
 //Being f i t   o n   t h e   c l s   v s   H V   d i s t r i b u t i o n
 //
 TGraphErrors *hvcls = new TGraphErrors(points, hv, cls, hverr, clserr);
 // Defining the fit function and setting parameters  on the cls 
 TF1 *f2 = (TF1*) gROOT->GetFunction("chebyshev3");
 f2->SetParNames("a","b","c","d"); 
 f2->SetParLimits(0,-100,100);
 f2->SetParLimits(1,-100,100);
 f2->SetParLimits(2,-100,100);
 f2->SetParLimits(3,-100,100);
 f2->SetParameters(0.0,0.0,0.0,0.0); 
 
 //TF1 *f2 = new TF1("f2",PolyFuncFit, 8.5, 10.0 ,4);//Same model
 hvcls->Fit(f2,"RFN","",8.5,9.999);
 clswp = PolyFunccalc(wp,f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2),f2->GetParameter(3));
 if (strcmp(subdetect,"endcap")==0 )clsknee = PolyFunccalc(wp - 0.120,f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2),f2->GetParameter(3));
 else clsknee = PolyFunccalc(wp - 0.100,f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2),f2->GetParameter(3));
 chi2cls=(f2->GetChisquare())/(points-1);
 
 
 //S A V E   T H E   R E S U L T S   I N   A   T X T   F I L E  
 ofstream fitResCls;
 std::string chamber_ = chamber;  
 fitResCls.open(("../results_scan2022/"+chamber_+"/fitDataCls.txt").c_str());  
 fitResCls<<f2->GetParameter(0)<<" "
          <<f2->GetParameter(1)<<" "
          <<f2->GetParameter(2)<<" "
          <<f2->GetParameter(3)<<" "
          << chi2cls<<" "
          <<wp<<" "
          <<clswp<<" "<<endl;
 fitResCls.close();
 return;
}

void FitData(const char *subdetect){
   
   std::ifstream RollEff;
   std::ofstream runsData;
   std::vector< std::pair<Int_t,Float_t> > hvscan;
   hvscan = hvEff(subdetect);
   cout << "Look here 4" << endl;
   int run = int(hvscan.size()); 
   cout << "Look here 5" << endl;  
   rData(subdetect); 
   cout << "Look here 6" << endl;  
   std::string id_; 
   float hv;
   float eff;
   float err;
   float exp;
   float cls;
   const int RUN = run;
   float HV[RUN];
   float HVerr[RUN];
   float EFF[RUN];
   float ERR[RUN];
   float EXP[RUN];
   float CLS[RUN];
   float CLSerr[RUN];
   
   cout << "Look here 7" << endl;
   
   //The chamber has a name and id... Read detId.txt file and define a map  
   //map -> first == name and map-> second ==id 
   std::vector< std::pair<std::string,std::string> > map;
   cout << "Look here 8" << endl;
   map = dictionary(subdetect);
   std::cout << map.size() << std::endl;
   cout << "Look here 9" << endl;

   for (vector<std::pair<std::string,std::string> >::const_iterator itmap = map.begin() ;itmap != map.end(); itmap++  )
   {      
   
		for (vector<std::pair<Int_t,Float_t> >::const_iterator it = hvscan.begin() ;it != hvscan.end(); it++  ){  
               run = int(it->first);
               std::stringstream s; 
               std::string strR;
               s << run; strR = s.str();
               if (strcmp(subdetect,"barrel")==0) RollEff.open(("../data/rollEff_"+strR+"_b.txt").c_str());
               else  RollEff.open(("../data/rollEff_"+strR+"_ec.txt").c_str());
               while (1){
              		 RollEff >> id_ >> eff >> err >> exp >> cls;   
               		 if (RollEff.eof())break; 
               		 if ((itmap->second) != id_ )continue;
   			 HV[run-1] = it->second;
   			 HVerr[run-1] = 0.001;//Force an error in hv 
                         EFF[run-1]= eff; 
                	 ERR[run-1]=err; 
                	 EXP[run-1]=exp; 
                	 CLS[run-1]=cls; 
                	 if(ERR[run-1]==0.) ERR[run-1]=99.99;
                         if(EFF[run-1]==0.) EFF[run-1]=0.0001;
                         if(EXP[run-1]==0.) EXP[run-1]=1.;
                         CLSerr[run-1] = 1/sqrt((EXP[run-1])*(EFF[run-1]/100)); 
                         //CLSerr[run-1] = 0.3;//Force the error probably wrong 
                         }
                                  
            RollEff.close();
          } 
    
    id_ =itmap->first;
    const char *c = id_.c_str();
    Double_t wp;

	cout << "Chamber id " << id_ << endl; 

    wp = FitDataFuncEff(subdetect,c,RUN, HV, HVerr, EFF, ERR); 
    FitDataFuncCls(subdetect,c,RUN, HV, HVerr,CLS, CLSerr, wp); 
    
    //%%%%%% Only if you want to save the data to each roll
    
    gSystem->mkdir("../results_scan2022");
    gSystem->mkdir(("../results_scan2022/"+id_).c_str());
    runsData.open(("../results_scan2022/"+id_+"/runsData.txt").c_str());          
    cout << id_ << " " << " create directory" <<endl;
    for (int n=0;n<RUN;n++){
             runsData <<HV[n] << " "
             <<EFF[n]<< " "    
             <<ERR[n]<< " "
             <<EXP[n]<< " "   
             <<CLS[n]<<std::endl;
             } 
	runsData.close();
     //%%%%%%%%%%%%%%%%%%%%%%%%% 
	cout << "results" << endl;
     
    }
  exit(0); 
}

