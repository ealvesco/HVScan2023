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
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

//P r o t o t y p e 
//
Double_t expFunc(Double_t*, Double_t*);//Function prototype 
Double_t PolyFuncFit(Double_t*, Double_t*);//Function prototype 
Double_t PolyFunccalc(Double_t, Double_t, Double_t , Double_t, Double_t );//Function prototype 
Double_t expFunccalc(double, double, double);//Function prototype 
Double_t SigmoidFunc( Double_t*, Double_t*);//Function prototype 
Double_t Sigmoidcalc(Double_t,Double_t,Double_t,Double_t);//Function prototype 
Double_t difcalc(Double_t, Double_t, Double_t, Double_t);//Function prototype 
std::vector< std::pair<Int_t,Float_t> > hvEff(const char* );//Function prototype 
std::vector< std::pair<std::string,std::string> > dictionary(const char*);//Function prototype 
void DrawingEff(const char* ,const char* , int, Double_t *,Double_t *, Double_t *,Double_t *, Double_t *, Double_t *,Double_t *, Double_t *, const char*, const char*);
//void pngProducer(const char* );




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


//----F U N C T I O N S   T O   D A T A   H A N D L I N G ---------------
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
		if (strncmp(subd,"barrel",6)==0)hv.push_back(hv_b);
		else hv.push_back(hv_ec);
	}
	f.close();
        hv.erase( unique( hv.begin(), hv.end() ), hv.end() );
	for (unsigned int i=0; i<hv.size();i++){
		std::pair<Int_t, Float_t > pair_hv;
		if ( hv.at(i)== 0 )continue;
		run++;
                pair_hv.first = run; 
	        pair_hv.second = hv.at(i);
		hvEff_.push_back(pair_hv); 
	}
        hv.clear();
	return hvEff_ ;
}

// To do a map of the rolls. Every roll has an id match a name   
std::vector< std::pair<std::string,std::string> > dictionary(const char* subd){
  ifstream f;
  Double_t id;
  std::string name;
  if (strncmp(subd,"endcap",6)==0)f.open("../data/detIdEndCapScan2022.txt");
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

 
void DrawingEff(const char* subdetect,const char* chamber, int points, Double_t *hv,Double_t *hverr, Double_t *eff,Double_t *efferr,Double_t *fitparamEff, Double_t *cls, Double_t *clserr, Double_t* fitparamCls, const char* id, const char* name){

// E F F I C I E N C Y   P L O T S 

	if (points ==0.) {std::cout << "N O  D A T A "; return;} 
	std::cout << subdetect <<":  "<< chamber << std::endl; 
	std::string chamber_ = chamber; 
	std::string id_ = id; 
	std::string name_ = name;
	std::string canvasname1 = id_+"_eff"; 
	const int P=200;
	Float_t hvc[2];
	Float_t effc[2];
	Double_t xmax=9.9;
	Double_t xmin = 8.5;
	Double_t x[P];
	Double_t y[P];

 //Creating the graphics 
 
 	TGraphErrors *hveff = new TGraphErrors(points, hv, eff, hverr, efferr);// first graphics HV points  
    hveff->GetXaxis()->SetTitle("High voltage (kV)");

 	TF1 *f1 = new TF1("f1",SigmoidFunc, 8.5, 9.9  ,3);//range of the function and how many parameters  
 	f1->SetParNames("emax","slope","hv50");
	f1->SetParameter(0, fitparamEff[2]);
 	f1->SetParameter(1, fitparamEff[9]);//
 	f1->SetParameter(2, fitparamEff[3]);

 	Double_t knee;
 	Double_t effknee;
 	if (strncmp(subdetect,"endcap",6)==0 )knee = fitparamEff[0] - 0.120;  
 	else knee = fitparamEff[0] - 0.100;
 	effknee=Sigmoidcalc(knee, fitparamEff[2],fitparamEff[9], fitparamEff[3]);
 	hvc[0]=knee;
 	hvc[1]=fitparamEff[0];//wp
 	effc[0]=effknee;//eff at knee
 	effc[1]=fitparamEff[5];//effatwp 
	TGraph *hvpoints = new TGraph(2,hvc,effc);//second graphic, only two points 

 	for (int k=0; k<P; k++){
    	x[k]=xmin+k*(xmax-xmin)/(P-1);
      	y[k]=Sigmoidcalc(x[k],fitparamEff[2],fitparamEff[9],fitparamEff[3]);
  	}
 
 	TGraph *sigmoid = new TGraph(P,x,y);//third graphic, the fit curve. 
  

 	// Being made the fit on the efficiency vs HV distribution  
 	int W = 1800;
 	int H = 600; 
 	int H_ref = 1800; 
 	int W_ref = 600;

 	float T = 0.08*H_ref;
  	float B = 0.14*H_ref; 
  	float L = 0.14*W_ref;
  	float R = 0.03*W_ref;


  	gStyle->SetLineWidth(4);
  	gStyle->SetHistLineWidth(4); 
  	gStyle->SetFrameLineWidth(4); 

  
  	//TCanvas *c1 = new TCanvas(TString(canvasname1),"Sigmoid",200,10,600,400);
  	TCanvas* c1 = new TCanvas(TString(canvasname1),"Sigmoid",50,50,W,H);
  	c1->Divide(3,1);
	c1->cd(2);
    c1->SetFillColor(0);
  	c1->SetBorderMode(0);
  	c1->SetFrameFillStyle(0);
  	c1->SetFrameBorderMode(0);
  	c1->SetLeftMargin( L/W );
  	c1->SetRightMargin( R/W );
  	c1->SetTopMargin( T/H );
  	c1->SetBottomMargin( B/H );
  	c1->SetTickx(1);
  	c1->SetTicky(1); 
//  c1->SetGridx(1);
//  c1->SetGridy(1); 
 
 	hveff->SetLineColor(2);
 	hveff->SetMarkerStyle(20);
 	hveff->SetMarkerSize(2.0);
 	hveff->SetMinimum(-0.01);
 	hveff->SetMaximum(110);
 	TAxis *axis = hveff->GetXaxis();
 	axis->SetLimits(8.5,9.9);
 	hveff->SetTitle(("HV vs efficiency  " + chamber_).c_str());
 	hveff->GetXaxis()->SetTitle("High voltage (kV)");
 	hveff->GetYaxis()->SetTitle("Efficiency(%)"); 
 	hveff->GetXaxis()->SetTitleSize(0.04);
 	hveff->GetXaxis()->SetTitleOffset(1.1);
 	hveff->GetYaxis()->SetTitleSize(0.05);
 	hveff->GetYaxis()->SetTitleOffset(1.0);
 	hveff->Draw("AP");
 
 	hvpoints->SetMarkerStyle(28);
 	hvpoints->SetMarkerSize(3);
 	hvpoints->SetLineColor(4);
 	hvpoints->Draw("P");
   
 	sigmoid->SetLineColor(4);
 	sigmoid->SetLineWidth(3);
 	sigmoid->Draw("C");
 
 //S A V E   T H E   R E S U L T S   I N   A   P N G   F I L E  
 //gSystem->mkdir(("../results/"+chamber_).c_str());
 	//c1->SaveAs(("../results.04.08/"+name_+"_EFFvsHV.png").c_str());
 	//c1->SaveAs(("../results/"+name_+"_EFFvsHV.C").c_str());
 	//c1->Clear(); 
 

 //C L U S T E R   S I Z E   P L O T S
	c1->cd(3);
 	std::string canvasname2 = id_+"_cls";
 	const int P_cls=200;
 	Double_t hvc_cls[2];
 	Double_t clsc[2];
 	Double_t xmax_cls=9.9;
 	Double_t xmin_cls = 8.5;
 	Double_t x_cls[P];
 	Double_t y_cls[P];


 	TGraphErrors *hvcls = new TGraphErrors(points, hv, cls, hverr, clserr);
 // Defining the fit function and setting parameters  on the cls 
 	TF1 *f2 = (TF1*) gROOT->GetFunction("chebyshev3");//root6
 	f2->SetParNames("a","b","c","d"); 
 
 	f2->SetParameter(0, fitparamCls[0]);
 	f2->SetParameter(1, fitparamCls[1]);
 	f2->SetParameter(2, fitparamCls[2]);
 	f2->SetParameter(3, fitparamCls[3]);
 
 
 	Double_t knee_cls;
	Double_t clsknee;
 	
	if (strncmp(subdetect,"endcap",6)==0)knee_cls = fitparamCls[5] - 0.120;//set knee = wp - 0.120
 	else knee_cls = fitparamCls[5] - 0.100; //set knee = wp - 0.1
 	clsknee=PolyFunccalc(knee,fitparamCls[0],fitparamCls[1],fitparamCls[2], fitparamCls[3]);
 	hvc_cls[0]=knee_cls;
 	hvc_cls[1]=fitparamCls[5];//wp
 	clsc[0]=clsknee;//cls at knee
 	clsc[1]=fitparamCls[6];//clsatwp 
 	TGraph *hvpointsCls = new TGraph(2,hvc_cls,clsc);
 
 	for (int k=0; k<P; k++){
		x_cls[k]=xmin_cls+k*(xmax_cls-xmin_cls)/(P_cls-1);
		y_cls[k]=PolyFunccalc(x_cls[k],fitparamCls[0],fitparamCls[1],fitparamCls[2], fitparamCls[3]);
   	}
 	
	TGraph *poly = new TGraph(P_cls,x_cls,y_cls);//third graphic, the fit curve.  
 
 // Being made the fit on the efficiency vs HV distribution  
 // Being made the fit on the efficiency vs HV distribution  
 	/*int W = 800;
 	int H = 800; 
   	int H_ref = 800; 
 	int W_ref = 800;
 	float T = 0.08*H_ref;
 	float B = 0.14*H_ref; 
 	float L = 0.14*W_ref;
 	float R = 0.03*W_ref;*/



  	gStyle->SetLineWidth(4);
  	gStyle->SetHistLineWidth(4); 
  	gStyle->SetFrameLineWidth(4); 

   
  //TCanvas *c1 = new TCanvas(TString(canvasname1),"Sigmoid",200,10,600,400);
  /*	TCanvas* c2 = new TCanvas(TString(canvasname2),"Polynomial",50,50,W,H);
  	c2->SetFillColor(0);
  	c2->SetBorderMode(0);
 	c2->SetFrameFillStyle(0);
  	c2->SetFrameBorderMode(0);
  	c2->SetLeftMargin( L/W );
  	c2->SetRightMargin( R/W );
  	c2->SetTopMargin( T/H );
  	c2->SetBottomMargin( B/H );
  	c2->SetTickx(1);
  	c2->SetTicky(1);*/ 
  //TCanvas *c2 = new TCanvas(TString(canvasname2),"Polynom",200,10,600,400);
 	hvcls->SetLineColor(2);
 	hvcls->SetMarkerStyle(20);
 	hvcls->SetMarkerSize(2.0);
 	hvcls->SetMinimum(-0.01);
 	hvcls->SetMaximum(6);
 	TAxis *axis_cls = hvcls->GetXaxis();
 	axis_cls->SetLimits(8.5,9.9);
 	hvcls->SetTitle(("HV vs Cluster-Size" + chamber_).c_str());
 	hvcls->GetXaxis()->SetTitle("High Voltage (kV)");
 	hvcls->GetXaxis()->SetTitleSize(0.04);
 	hvcls->GetXaxis()->SetTitleOffset(1.0);
 	hvcls->GetYaxis()->SetTitleSize(0.04);
 	hvcls->GetYaxis()->SetTitleOffset(1.1);
 	hvcls->GetYaxis()->SetTitle("Cluster size"); 
 	hvcls->Draw("AP");
 	
 	hvpointsCls->SetMarkerStyle(28);
 	hvpointsCls->SetMarkerSize(3);
 	hvpointsCls->SetLineColor(4);
 	hvpointsCls->Draw("P");
 	  
 	poly->SetLineColor(4);
 	poly->SetLineWidth(3);
 	poly->Draw("C");
 
// P U T   A   S U M M A R Y   O F   E F F   F I T 

	c1->cd(1);

	TPaveText* Text = new TPaveText(0.5,0.92,0.9,0.96,"nbNDC");

	//char wp_text;
    //sprintf(wp_text,"%f",fitparamEff[0]);//wp
    
    double wp_text = fitparamEff[0];
	double slope50_text = fitparamEff[1];//slope50
    double emax_text = fitparamEff[2];//emax
    double hv50_text = fitparamEff[3];//hv50
    double chi2_text = fitparamEff[4];//chi2
	double clswp_text = fitparamCls[6];//clswp
    double effwp_text = fitparamEff[5];//effwp
    //double emaxer_text = FitResEff[6]//emaxerr
    //double slopeerr_text = FitResEff[7]//slopeerr
    //double hv50err_text = FitResEff[8]//hv50err
    //double slope_text = FitResEff[9];//slope*/

	const char *name_text = name_.c_str();

	cout << name_text << endl;

	TLatex chambername;
    chambername.SetTextSize(0.08);
    chambername.SetTextAlign(13);  //align at top
   	chambername.DrawLatex(.2,.92,name_text);

	TLatex parameters;

	chambername.SetTextSize(0.08);
    chambername.SetTextAlign(12);
   	chambername.DrawLatex(.2,.78,Form("wp = %g",wp_text));
    chambername.DrawLatex(.2,.68,Form("slope50 = %g",slope50_text));
    chambername.DrawLatex(.2,.58,Form("emax = %g",emax_text));
    chambername.DrawLatex(.2,.48,Form("hv50 = %g",hv50_text));
    chambername.DrawLatex(.2,.38,Form("chi2 = %g",chi2_text));
    chambername.DrawLatex(.2,.28,Form("clswp = %g",clswp_text));
    chambername.DrawLatex(.2,.18,Form("effwp = %g",effwp_text));
	 
 
	//TLatex *t = new TLatex(.5,.5,Form("wp = %g",wp_text));

 	//t->Draw();
	
 	//S A V E   T H E   R E S U L T S   I N   A   T X T   F I L E  
//gSystem->mkdir(("../results/"+chamber_).c_str());
 	c1->SaveAs(("../results_scan2022/"+name_+"_CLSvsHV.png").c_str());
 //c2->SaveAs(("../results/"+id_+"_CLSvsHV_Pt7chi8.C").c_str());
 	c1->Clear(); 
 	return; 
}

void pngProducer_split(const char* subdetect){
   
   std::ifstream RollEff;
   std::ifstream fResEff, fResCls;
   std::vector< std::pair<Int_t,Float_t> > hvscan;
   hvscan = hvEff(subdetect);
   cout << "Passed hveff" << endl;
   int run = int(hvscan.size());   
   //rData(subdetect); 
   
   std::string id_; 
   std::string name_; 
   Double_t hv;
   Double_t eff;
   Double_t err;
   Double_t exp;
   Double_t cls;
   const int RUN = run;
   Double_t HV[RUN];
   Double_t HVerr[RUN];
   Double_t EFF[RUN];
   Double_t ERR[RUN];
   Double_t EXP[RUN];
   Double_t CLS[RUN];
   Double_t CLSerr[RUN];
   Double_t FitResEff[10];
   Double_t FitResCls[7];
  
   cout << "Passed vars" << endl;

   //The chamber has a name and id... Read detId.txt file and define a map  
   //map -> first == name and map-> second ==id 
   std::vector< std::pair<std::string,std::string> > map;
   map = dictionary(subdetect);
   cout << "Passed dictionary" << endl;
   std::cout << map.size()<< std::endl;
   for (vector<std::pair<std::string,std::string> >::const_iterator itmap = map.begin() ;itmap != map.end(); itmap++  )
   {      
         //std::cout << itmap->first << " " <<itmap->second << std::endl;;
         for (vector<std::pair<Int_t,Float_t> >::const_iterator it = hvscan.begin() ;it != hvscan.end(); it++  ){  
               run = int(it->first);
               std::stringstream s; 
               std::string strR;
               s << run; strR = s.str();
               if (strncmp(subdetect,"barrel",6)==0) RollEff.open(("../data/rollEff_"+strR+"_b.txt").c_str());
               //if (subdetect=="barrel") RollEff.open(("../data/rollEff_"+strR+"_b.txt").c_str());
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
                	 if(ERR[run-1]==0.) ERR[run-1]=0.00001;
                         if(EFF[run-1]==0.) EFF[run-1]=0.001;
                         if(EXP[run-1]==0.) EXP[run-1]=1.;
                         CLSerr[run-1] = 1/sqrt((EXP[run-1])*(EFF[run-1]/100));; 
                         //CLSerr[run-1] = 0.3;//Force the error probably wrong 
                         }
          RollEff.close();
  
          } 
	  name_ =itmap->first;
      id_   = itmap->second;
	  const char *c = name_.c_str();
	  fResEff.open(("../results_scan2022/"+name_+"/fitData.txt").c_str());
      cout << "Passed loop" << endl;
      while (1){
		  cout << itmap->first << "  " << FitResEff[0] << endl;
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
		  cout << "inside loop" << endl;
		  if (fResEff.eof())break;
		  } 
	  fResEff.close();  
      cout << "Passed fResEff" << endl;

      fResCls.open(("../results_scan2022/"+name_+"/fitDataCls.txt").c_str());
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

      DrawingEff(subdetect, c, RUN, HV, HVerr, EFF, ERR, FitResEff, CLS, CLSerr, FitResCls, (id_).c_str(), (name_).c_str());
    }
  exit(0); 
}

