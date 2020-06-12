/*
 * macro to analyze and save di-Photon, di-Electron and di-Muon spectrum
 */
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCut.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>

using namespace std;

const int MAXJETS = 500;
const float cutPt = 10;
const float cutEta = 1.4791;
const float cutDeltaR = 0.15;
const long MAXTREESIZE = 200000000000; // set maximum tree size from 10 GB to 100 GB, so that the code does not switch to a new file after 10 GB7

const double eleMass = 0.000511;
const double muMass  = 0.105658;
//float wt  = 0.0083;
float wt  = 1;

const bool looseIsolation  = true;
#define PI 3.141592653589
const float awayRange = PI * 7./8.;

//void mycheck( );
Double_t getDR  ( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
Double_t getDETA(Double_t eta1, Double_t eta2);
Double_t getX(Double_t Pt, Double_t Phi);
Double_t getY(Double_t Pt, Double_t Phi);
Double_t getZ(Double_t Pt, Double_t Eta);
Double_t cosTheta12(Double_t Eta1, Double_t Phi1, Double_t Eta2, Double_t Phi2);
Double_t getInvMass(Double_t Energy1, Double_t Eta1, Double_t Phi1, Double_t Energy2, Double_t Eta2, Double_t Phi2);

//void check_data_driven_id_eff( std::string infile_Forest = "filelist_28.txt",
void check_data_driven_reco_eff( std::string infile_Forest = "data.txt",
			       std::string out = "test.root") 
{
  bool ismc = 1; 
  TH1::SetDefaultSumw2();
  
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  instr_Forest>>filename_Forest;
  cout << "opening file " << filename_Forest << endl;
  TFile *fin = TFile::Open(filename_Forest.c_str());
  
  TTree *treeHLT   = (TTree*)fin->Get("hltanalysis/HltTree");
  TTree *tree      = (TTree*)fin->Get("ggHiNtuplizer/EventTree");
  TTree *treePixel = (TTree*)fin->Get("pixel/PixelTree");
  
  TH1D* htrigger_Ereco_before_trigg = new TH1D("htrigger_Ereco_before_trigg","",20,0,10);
  TH1D* htrigger_Ereco_after_trigg  = new TH1D("htrigger_Ereco_after_trigg","",20,0,10);
  TH1D* htrigger_eff_vs_Ereco   = new TH1D("htrigger_eff_vs_Ereco","",20,0,10);
  
  TH1D* hsingle_pt   = new TH1D("hsingle_pt","",20,0,10);
  TH1D* hsingle_eta  = new TH1D("hsingle_eta","",6,-2.4,2.4);
  TH1D* hsingle_phi  = new TH1D("hsingle_phi","",8,-4,4);
  
  TH1D* hdpt   = new TH1D("hdpt","",10,0,5);
  TH1D* hdeta  = new TH1D("hdeta","",100,0,0.1);
  TH1D* hdphi  = new TH1D("hdphi","",100,0,0.5);
  TH1D* haco   = new TH1D("haco","",12,0,0.06);
  
  TH1D* hdouble_pt    = new TH1D("hdouble_pt","",12,0,2);
  TH1D* hdouble_rap   = new TH1D("hdouble_rap","",6,-2.4,2.4);
  TH1D* hdouble_mass  = new TH1D("hdouble_mass","",10,0,30);
  
  /// HLT tree ////////////////////////////////////////////////////////
  treeHLT->SetBranchStatus("*",0);     // disable all branches
  treeHLT->SetBranchStatus("HLT_HIUPCL1SingleEG5*",1);     // enable photon branches
  treeHLT->SetBranchStatus("HLT_HIUPCL1DoubleEG2*",1);     // enable photon branches
  
  Int_t HLT_HIUPCL1DoubleEG2NotHF2_v1 ;   Int_t HLT_HIUPCL1SingleEG5NotHF2_v1;
  
  treeHLT->SetBranchAddress("HLT_HIUPCL1SingleEG5NotHF2_v1",&HLT_HIUPCL1SingleEG5NotHF2_v1);     // enable photon branches
  treeHLT->SetBranchAddress("HLT_HIUPCL1DoubleEG2NotHF2_v1",&HLT_HIUPCL1DoubleEG2NotHF2_v1);     // enable photon branches
  
  // event information
  UInt_t run, lumis;
  ULong64_t event;
  tree->SetBranchStatus("*",0);     // disable all branches
  tree->SetBranchStatus("run",1);
  tree->SetBranchStatus("event",1);
  tree->SetBranchStatus("lumis",1);
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("lumis", &lumis);
  
  // vertex info
  Int_t nVtx;
  std::vector<float>* xVtx=0;
  std::vector<float>* yVtx=0;
  std::vector<float>* zVtx=0;
  std::vector<int>*   isFakeVtx=0;
  
  tree->SetBranchStatus("nVtx",1);     
  tree->SetBranchAddress("nVtx",&nVtx);
  tree->SetBranchStatus("xVtx",1);     
  tree->SetBranchStatus("yVtx",1);     
  tree->SetBranchStatus("zVtx",1);  
  tree->SetBranchStatus("isFakeVtx",1);  
  
  tree->SetBranchAddress("xVtx",&xVtx);
  tree->SetBranchAddress("yVtx",&yVtx);
  tree->SetBranchAddress("zVtx",&zVtx);
  tree->SetBranchAddress("isFakeVtx",&isFakeVtx);
  
  
  Int_t nMC;
  std::vector<float>* mcE=0;
  std::vector<float>* mcEt=0;
  std::vector<float>* mcPt=0;
  std::vector<float>* mcP=0;
  std::vector<float>* mcEta=0;
  std::vector<float>* mcPhi=0;
  
  tree->SetBranchStatus("nMC",1);     // enable photon branches
  tree->SetBranchStatus("mc*",1);     // enable photon branches
  tree->SetBranchAddress("nMC",&nMC);
  tree->SetBranchAddress("mcE",&mcE);
  tree->SetBranchAddress("mcEt",&mcEt);
  tree->SetBranchAddress("mcPt",&mcPt);
  tree->SetBranchAddress("mcP",&mcP);
  tree->SetBranchAddress("mcEta",&mcEta);
  tree->SetBranchAddress("mcPhi",&mcPhi);
  
  // GSF electron tracks
   // GED RECO electrons
  Int_t nEle;
  std::vector<int>*   eleCharge=0;
  std::vector<float>* eleP=0;
  std::vector<float>* elePt=0;
  std::vector<float>* eleEta=0;
  std::vector<float>* elePhi=0;
 
   
  tree->SetBranchStatus("nEle",1);     // enable electron branches
  tree->SetBranchStatus("ele*",1);     // enable electron branches
  tree->SetBranchAddress("nEle",&nEle);
  tree->SetBranchAddress("eleCharge",&eleCharge);
  tree->SetBranchAddress("eleP",&eleP);
  tree->SetBranchAddress("elePt",&elePt);
  tree->SetBranchAddress("eleEta",&eleEta);
  tree->SetBranchAddress("elePhi",&elePhi);
 

  // RECO photons
  Int_t nPho;
  std::vector<float>* phoE=0;
  std::vector<float>* phoEt=0;
  std::vector<float>* phoEta=0;
  std::vector<float>* phoPhi=0;
  

  tree->SetBranchStatus("nPho",1);     // enable photon branches
  tree->SetBranchAddress("nPho",&nPho);
  tree->SetBranchStatus("pho*",1);     // enable photon branches
  tree->SetBranchAddress("phoE",&phoE);
  tree->SetBranchAddress("phoEt",&phoEt); 
  tree->SetBranchAddress("phoEta",&phoEta);
  tree->SetBranchAddress("phoPhi",&phoPhi);


 
  // RECO muons
  Int_t nMu;
  
  tree->SetBranchStatus("nMu",1);     // enable muon branches
  tree->SetBranchAddress("nMu",&nMu);
  

  /// General track collection/////////
  Int_t ngenTrk;
  std::vector<float>* gentrkPt=0;
  std::vector<float>* gentrkP=0;
  std::vector<float>* gentrkEta=0;
  std::vector<float>* gentrkPhi=0;
  std::vector<int>*   gentrkcharge=0;
  std::vector<float>* gentrkchi2=0;
  std::vector<float>* gentrknormchi2=0;
  std::vector<int>*   gentrkValidHits=0;
  std::vector<int>*   gentrkMissHits=0;


  tree->SetBranchStatus("ngenTrk",1);     // enable electron branches
  tree->SetBranchAddress("ngenTrk",&ngenTrk);
  tree->SetBranchStatus("gen*",1);     // enable electron branches
  tree->SetBranchAddress("gentrkPt",&gentrkPt);
  tree->SetBranchAddress("gentrkP",&gentrkP);
  tree->SetBranchAddress("gentrkEta",&gentrkEta);
  tree->SetBranchAddress("gentrkPhi",&gentrkPhi);
  tree->SetBranchAddress("gentrkcharge",&gentrkcharge);
  tree->SetBranchAddress("gentrkchi2",&gentrkchi2);
  tree->SetBranchAddress("gentrknormchi2",&gentrknormchi2);
  tree->SetBranchAddress("gentrkValidHits",&gentrkValidHits);
  tree->SetBranchAddress("gentrkMissHits",&gentrkMissHits);


  int pix_hit1;
  int pix_hit2;
  int pix_hit3;
  int pix_hit4;
  int pix_hit5;


  float phodpt;
  float phodeta;
  float phodphi;
  float aco;
  
  
   TFile *output;
   output = new TFile(Form("%s",out.c_str()),"recreate");

  // output tree variables

  TTree *outputTreeHLT    = treeHLT->CloneTree(0);
  outputTreeHLT->SetName("hltTree");
 
 
 
  TTree *outputTreediPho  = new TTree("diphoton","");

  
  outputTreeHLT->SetMaxTreeSize(MAXTREESIZE);
  outputTreediPho->SetMaxTreeSize(MAXTREESIZE);
  
  
  ////**************************** write output branches      *************************/
  
  outputTreediPho->Branch("run",&run);
  outputTreediPho->Branch("event",&event);
  outputTreediPho->Branch("lumis",&lumis);
  
  outputTreediPho->Branch("nVtx",&nVtx);
  outputTreediPho->Branch("xVtx",&xVtx);
  outputTreediPho->Branch("yVtx",&yVtx);
  outputTreediPho->Branch("zVtx",&zVtx);
  outputTreediPho->Branch("isFakeVtx",&isFakeVtx);
  
  std::vector<float> ele_matchedPt;
  std::vector<float> ele_matchedP;
  std::vector<float> ele_matchedEta;
  std::vector<float> ele_matchedPhi;
  std::vector<int>   ele_matchedCharge;
 
  outputTreediPho->Branch("nEle",&nEle);
  outputTreediPho->Branch("ele_matchedPt",&ele_matchedPt);
  outputTreediPho->Branch("ele_matchedP",&ele_matchedP);
  outputTreediPho->Branch("ele_matchedEta",&ele_matchedEta);
  outputTreediPho->Branch("ele_matchedPhi",&ele_matchedPhi);
  outputTreediPho->Branch("ele_matchedCharge",&ele_matchedCharge);

  std::vector<float> trk_matchedPt;
  std::vector<float> trk_matchedP;
  std::vector<float> trk_matchedEta;
  std::vector<float> trk_matchedPhi;
  std::vector<int>   trk_matchedCharge;

  outputTreediPho->Branch("trk_matchedPt",&trk_matchedPt);
  outputTreediPho->Branch("trk_matchedP",&trk_matchedP);
  outputTreediPho->Branch("trk_matchedEta",&trk_matchedEta);
  outputTreediPho->Branch("trk_matchedPhi",&trk_matchedPhi);
  outputTreediPho->Branch("trk_matchedCharge",&trk_matchedCharge);


  std::vector<float> trk_unmatchedPt;
  std::vector<float> trk_unmatchedP;
  std::vector<float> trk_unmatchedEta;
  std::vector<float> trk_unmatchedPhi;
  std::vector<int>   trk_unmatchedCharge;

  outputTreediPho->Branch("trk_unmatchedPt",&trk_unmatchedPt);
  outputTreediPho->Branch("trk_unmatchedP",&trk_unmatchedP);
  outputTreediPho->Branch("trk_unmatchedEta",&trk_unmatchedEta);
  outputTreediPho->Branch("trk_unmatchedPhi",&trk_unmatchedPhi);
  outputTreediPho->Branch("trk_unmatchedCharge",&trk_unmatchedCharge);
 
  


  // write tower info ////
  
 
  outputTreediPho->Branch("ngenTrk",&ngenTrk);


  std::vector<float> phoEt_1;
  std::vector<float> phoEta_1;
  std::vector<float> phoPhi_1;
 

  outputTreediPho->Branch("phoEt_1",&phoEt_1);
  outputTreediPho->Branch("phoEta_1",&phoEta_1);
  outputTreediPho->Branch("phoPhi_1",&phoPhi_1);
  
  std::vector<float> dr_1;
  std::vector<float> dr_2;
 
  outputTreediPho->Branch("dr_1",&dr_1);
  outputTreediPho->Branch("dr_2",&dr_2);
  
  Long64_t entries = tree->GetEntries();
  Long64_t entriesAnalyzed = 0;   Long64_t id = 0;   Long64_t excl = 0; Long64_t ptcut = 0; Long64_t myloop = 0;
  std::cout << "entries         = " << entries << std::endl;
  std::cout<< "Loop : ggHiNtuplizer/EventTree" <<std::endl;
  
  for (Long64_t j_entry=0; j_entry<entries; ++j_entry){
    //for (Long64_t j_entry= 0; j_entry< 50000; ++j_entry){
    
    //if (j_entry % 10000 == 0)  {
    //  std::cout << "current entry = " <<j_entry<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j_entry/entries*100<<" %"<<std::endl;
    //}
    
    treeHLT->GetEntry(j_entry);
    tree->GetEntry(j_entry);
    treePixel->GetEntry(j_entry);
    
    if(HLT_HIUPCL1DoubleEG2NotHF2_v1==1 || HLT_HIUPCL1SingleEG5NotHF2_v1==1){   entriesAnalyzed++; 
      
      outputTreeHLT->Fill(); 

      // dielectron block
      
      ele_matchedPt.clear();
      ele_matchedP.clear();
      ele_matchedEta.clear();
      ele_matchedPhi.clear();
      ele_matchedCharge.clear();

      trk_matchedPt.clear();
      trk_matchedP.clear();
      trk_matchedEta.clear();
      trk_matchedPhi.clear();
      trk_matchedCharge.clear();


      trk_unmatchedPt.clear();
      trk_unmatchedP.clear();
      trk_unmatchedEta.clear();
      trk_unmatchedPhi.clear();
      trk_unmatchedCharge.clear();
     
      phoEt_1.clear();
      phoEta_1.clear();
      phoPhi_1.clear();
     
     
      dr_1.clear();       
      dr_2.clear();       
  

      //cout << "entry:" << j_entry << "  coming in loop"<< "  gentrk:" << ngenTrk << endl;
      if(ngenTrk ==2 && nEle ==1 && nPho==1 ){
	
	for(int i=0; i<ngenTrk; ++i){ //cout << " In first track loop " << endl;
	  for (int j=i+1; j<ngenTrk; ++j){
	    
	    if (gentrkcharge->at(i) == gentrkcharge->at(j)) continue;  
           	    
	    double dr1 = getDR(gentrkEta->at(i),gentrkEta->at(i),eleEta->at(0),elePhi->at(0));
	    double dr2 = getDR(gentrkEta->at(j),gentrkEta->at(j),eleEta->at(0),elePhi->at(0));
	    
	    dr_1.push_back(dr1);       
	    dr_2.push_back(dr2);      
	    
            if(dr1 < 1.0 && dr2 > 1.0 && gentrkPt->at(i) > 3.0 &&  gentrkPt->at(j)< 2.0){
	     
		cout << "Entry: "<< j_entry << "  Track 1 pt:" << gentrkPt->at(i) << "  Track 2 pt:"<< gentrkPt->at(j);
		cout << "  Elec pt:" << elePt->at(0) << "  photon et:" << phoEt->at(0)  ;
		cout <<"   dr1:"<< dr1 << "  dr2:" << dr2 << endl; 
	        myloop++;
		
		phoEt_1.push_back(phoEt->at(0));
		phoEta_1.push_back(phoEta->at(0));
		phoPhi_1.push_back(phoPhi->at(0));
		
		ele_matchedPt.push_back(elePt->at(0));
		ele_matchedP.push_back(eleP->at(0));
		ele_matchedEta.push_back(eleEta->at(0));
		ele_matchedPhi.push_back(elePhi->at(0));
		ele_matchedCharge.push_back(eleCharge->at(0));
		
		
		trk_matchedPt.push_back(gentrkPt->at(i));
		trk_matchedP.push_back(gentrkP->at(i));
		trk_matchedEta.push_back(gentrkEta->at(i));
		trk_matchedPhi.push_back(gentrkPhi->at(i));
		
		
		trk_unmatchedPt.push_back(gentrkPt->at(j));
		trk_unmatchedP.push_back(gentrkP->at(j));
		trk_unmatchedEta.push_back(gentrkEta->at(j));
		trk_unmatchedPhi.push_back(gentrkPhi->at(j));
		
	     
            }

	    if(dr2 < 1.0 && dr1 > 1.0 && gentrkPt->at(j) > 3.0 &&  gentrkPt->at(i)< 2.0){
	   
		cout << "Entry: "<< j_entry << "  Track 1 pt:" << gentrkPt->at(i) << "  Track 2 pt:"<< gentrkPt->at(j);
		cout << "  Elec pt:" << elePt->at(0) << "  photon et:" << phoEt->at(0)  ;
		cout <<"   dr1:"<< dr1 << "  dr2:" << dr2 << endl; 
	        myloop++;

		phoEt_1.push_back(phoEt->at(0));
		phoEta_1.push_back(phoEta->at(0));
		phoPhi_1.push_back(phoPhi->at(0));
		
		ele_matchedPt.push_back(elePt->at(0));
		ele_matchedP.push_back(eleP->at(0));
		ele_matchedEta.push_back(eleEta->at(0));
		ele_matchedPhi.push_back(elePhi->at(0));
		ele_matchedCharge.push_back(eleCharge->at(0));
		
		
		trk_matchedPt.push_back(gentrkPt->at(j));
		trk_matchedP.push_back(gentrkP->at(j));
		trk_matchedEta.push_back(gentrkEta->at(j));
		trk_matchedPhi.push_back(gentrkPhi->at(j));
		
		
		trk_unmatchedPt.push_back(gentrkPt->at(i));
		trk_unmatchedP.push_back(gentrkP->at(i));
		trk_unmatchedEta.push_back(gentrkEta->at(i));
		trk_unmatchedPhi.push_back(gentrkPhi->at(i));
	   
            }
	  } //trk j
	} // trk i

      } // cuts on photon , electron and tracks
      outputTreediPho->Fill();
      
    } //trigger
  } //event loop
  std::cout<<  "Loop ENDED : ggHiNtuplizer/EventTree" <<std::endl;
  std::cout << "entries            = " << entries << std::endl;
  std::cout << "In my loop         = " << myloop << std::endl;
  std::cout << "In three body loop = " << id << std::endl;
  
  output->cd();
  hsingle_pt->Write();
  hsingle_eta->Write();
  hsingle_phi->Write();
  hdouble_pt->Write();
  hdouble_rap->Write();
  hdouble_mass->Write();
  hdpt->Write();
  hdeta->Write();
  hdphi->Write();
  haco->Write();
  output->Write();
  output->Close();
}

int main(int argc, char *argv[])
{ 
  check_data_driven_reco_eff(argv[1]);
  //mycheck_ele_tree(argv[3],argv[4]);
  return 0;
}

Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
  Double_t theDphi = getDPHI( phi1, phi2);
  Double_t theDeta = eta1 - eta2;
  return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;
  
  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;
  
  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }
  
  return TMath::Abs(dphi);
  //return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
  return TMath::Abs(eta1 - eta2);
}
Double_t getX(Double_t Pt, Double_t Phi){
  return Pt*TMath::Cos(Phi);
}
Double_t getY(Double_t Pt, Double_t Phi){
  return Pt*TMath::Sin(Phi);
}
Double_t getZ(Double_t Pt, Double_t Eta){
  return Pt*sinh(Eta);
}
Double_t cosTheta12(Double_t Eta1, Double_t Phi1, Double_t Eta2, Double_t Phi2) {
    return ((cos(Phi1 - Phi2) + sinh(Eta1) * sinh(Eta2)) / (cosh(Eta1) * cosh(Eta2)));
}
Double_t getInvMass(Double_t Energy1, Double_t Eta1, Double_t Phi1, Double_t Energy2, Double_t Eta2, Double_t Phi2) {
  return (sqrt(2 * Energy1 * Energy2 * (1 - cosTheta12(Eta1, Phi1, Eta2, Phi2))));
}

