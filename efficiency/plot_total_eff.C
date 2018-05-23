#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
//#include "../../../CMS_lumi.C"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TCut.h"
#include "TChain.h"
#include "THStack.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TF1.h"



int xlo = 1;
int xhi = 2;


double getError(double A, double eA, double B, double eB);
const char *dir = "fig_eff";

const int nptbins=7;
double ptbin[nptbins]={0,2,4,6,8,10,160};
const int nptbin= sizeof(ptbin)/sizeof(double) - 1;


void make_hist(TH1D *&, Color_t , int );

void plot_total_eff(){
  
  //gROOT->LoadMacro("CMS_lumi.C");
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(1);
  
  TFile *outf= new TFile("all_plot.root","recreate");  
  
  TChain *ged = new TChain("diphoton");
  ged->Add("scged_diphoton_efficiency.root");
  
  
  const TCut purity = "pho_swissCrx_gen_reco_1 < 0.95 && fabs(pho_seedTime_gen_reco_1) < 3.6 && pho_swissCrx_gen_reco_2 < 0.95 && fabs(pho_seedTime_gen_reco_2) < 3.6 " ;

 const TCut crack = "fabs(phoEta_gen_reco_1)<2.4 && fabs(phoEta_gen_reco_2)<2.4 && (fabs(phoEta_gen_reco_1) < 1.4442 || fabs(phoEta_gen_reco_1) > 1.566) && (fabs(phoEta_gen_reco_2) < 1.4442 || fabs(phoEta_gen_reco_2) > 1.566) ";
  
  const TCut excl = "ngenTrk ==0 && nEle ==0 && ngsfEle==0";
  const TCut calocut = "leadingEmEnergy_EB< 0.55 && leadingEmEnergy_EE < 3.16 && leadingHadEnergy_HB < 2.0 && leadingHadEnergy_HE < 3.0 && leadingHadEnergy_HF_Plus < 4.85 && leadingHadEnergy_HF_Minus < 4.12 ";  

 
  //const TCut ged_cut     =  purity && crack && excl && calocut && "vSum_Pt_doubleEG2 < 1 && pho_aco_doubleEG2 < 0.01 && vSum_M_doubleEG2 > 5" ;

    const TCut ged_cut     =  purity && crack && "vSum_Pt_doubleEG2 < 1 && pho_aco_doubleEG2 < 0.01 && vSum_M_doubleEG2 > 5" ;

  //// ged histogram
  
  TH1D* hnum_ged          = new TH1D("hnum_ged" ,"",8,0,20);
  TH1D* hden_ged          = new TH1D("hden_ged" ,"",8,0,20);
  

  ged->Project(hnum_ged->GetName(),        "vSum_M_doubleEG2"  ,ged_cut);
  ged->Project(hden_ged->GetName(),        "vSum_gg_M_accp", "vSum_gg_M_accp > 5 && vSum_gg_Pt_accp < 1 && gen_aco < 0.01");
  
    
  double genError;
  double int_num = hnum_ged->IntegralAndError(1,9,genError);
  double int_num_err = genError;
  
  double denError;
  double int_den = hden_ged->IntegralAndError(1,9,denError);
  double int_den_err = denError;
  
  
  cout << "Integral:" << int_num << "  error num:" << int_num_err << "  den int:" << int_den << "  error" << int_den_err; 
  cout << "  efficiency:" << int_num/int_den << "  error: " << getError(int_num,int_num_err,int_den,int_den_err) <<  endl;  
  
  
  TGraphAsymmErrors* gr = new TGraphAsymmErrors();
  gr->Divide( hnum_ged, hden_ged );
  
  
  TLegend *leg2=new TLegend(0.69,0.27,0.96,0.39,NULL,"brNDC");
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(24);
  leg2->AddEntry(gr,"GED ","pl");
  
  
  
    
  int W = 700;
  int H = 600;
  
  float T = 0.08;
  float B = 0.14; 
  float L = 0.14;
  float R = 0.04;


  //// Diphoton invmass......................
  TCanvas* c1 = new TCanvas("c1","Reco efficiency",50,50,W,H);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin( L );
  c1->SetRightMargin( R );
  c1->SetTopMargin( T );
  c1->SetBottomMargin( B );
  c1->SetTickx(0);
  c1->SetTicky(0);
  
  c1->cd();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gr->SetMaximum(1.2);
  gr->SetMinimum(0.0);
  gr->SetMarkerColor(kRed);
  gr->SetLineColor(kRed);
  gr->GetXaxis()->SetTitle("m(#gamma#gamma) (GeV)");
  gr->GetYaxis()->SetTitle("Eff ^{#gamma#gamma}");
  gr->Draw("APE");
  
  leg2->Draw();
  //CMS_lumi( c1, 104, 10,lumi_PbPb2015 );
  c1->Print("total_diphoton_eff.gif");
  c1->Print("total_diphoton_eff.pdf");


  new TCanvas(); 
  hden_ged->Draw("p");
  hnum_ged->SetMarkerColor(kRed);
  hnum_ged->Draw("psames");





  // mass scale 

  TH2D* hmass_scale  = new TH2D("hmass_scale","",8,0,20,100,0.8,1.1);   
  ged->Project(hmass_scale->GetName(),     "vSum_M_doubleEG2/vSum_gg_M_accp:vSum_gg_M_accp", ged_cut && "vSum_gg_M_accp > 5 && vSum_gg_Pt_accp < 1 && gen_aco < 0.01");
  TProfile *hprofile_mass =  hmass_scale->ProfileX("hprofile_mass", 0, 100);  // raw barel profile 

  
  TCanvas* c2 = new TCanvas("c1","Reco efficiency",50,50,W,H);
  c2->Divide(1,2);
  c2->cd(1);
  hmass_scale->Draw("colz");
  c2->cd(2);
  hprofile_mass->Draw("p");

/*
  TH1D* hscale  = new TH1D("hscale","",80,0.8,1.2);   
  ged->Project(hscale->GetName(),     "vSum_M_doubleEG2/vSum_gg_M_accp",ged_cut);
  hscale->Draw("p");  */


  }

void make_hist(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.05, "XYZ");
 
  hist->SetMarkerColor(kcolor);
  hist->SetFillColor(kcolor);
  hist->SetMarkerStyle(kstyle);
  
 
}

//Ratio Error
double getError(double A, double eA, double B, double eB){
  double f=A/B;
  double fA=eA/A;
  double fB=eB/B;
  double eR=  f*sqrt( (fA*fA + fB*fB )) ;
  return eR;
}
