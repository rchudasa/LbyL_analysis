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

void check_data_driven_id_eff( std::string infile_Forest = "data.txt",
			 std::string out = "test.root") 
{

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
  
  TH1D* hsingle_pt   = new TH1D("hsingle_pt","",30,0,15);
  TH1D* hsingle_eta  = new TH1D("hsingle_eta","",6,-2.4,2.4);
  TH1D* hsingle_phi  = new TH1D("hsingle_phi","",8,-4,4);
  
  TH1D* hdpt   = new TH1D("hdpt","",10,0,5);
  TH1D* hdeta  = new TH1D("hdeta","",6,0,2.4);
  TH1D* hdphi  = new TH1D("hdphi","",20,2.65,3.15);
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
  std::vector<float>* mcP=0;
  std::vector<float>* mcPt=0;
  std::vector<float>* mcEta=0;
  std::vector<float>* mcPhi=0;
  
  tree->SetBranchStatus("nMC",1);     // enable photon branches
  tree->SetBranchStatus("mc*",1);     // enable photon branches
  tree->SetBranchAddress("nMC",&nMC);
  tree->SetBranchAddress("mcE",&mcE);
  tree->SetBranchAddress("mcEt",&mcEt);
  tree->SetBranchAddress("mcP",&mcP);
  tree->SetBranchAddress("mcPt",&mcPt);
  tree->SetBranchAddress("mcEta",&mcEta);
  tree->SetBranchAddress("mcPhi",&mcPhi);
 
  // GED RECO electrons
  Int_t nEle;
  std::vector<int>*   eleCharge=0;
  std::vector<float>* eleEn=0;
  std::vector<float>* eleP=0;
  std::vector<float>* elePt=0;
  std::vector<float>* eleEta=0;
  std::vector<float>* elePhi=0;
  std::vector<float>* eleHoverE=0;
  std::vector<float>* eleSigmaIEtaIEta=0;
  std::vector<float>* eleSigmaIEtaIEta_2012=0;
  std::vector<float>* eleSigmaIPhiIPhi=0;
  std::vector<float>* eleEoverP=0;
  std::vector<float>* eledEtaAtVtx=0;
  std::vector<float>* eledPhiAtVtx=0;
  std::vector<float>* eledEtaSeedAtVtx=0;
  std::vector<float>* eleD0=0;
  std::vector<float>* eleDz=0;
  std::vector<float>* eleTrkPt=0;
  std::vector<int>*   eleMissHits=0;
  std::vector<float>* eleTrackIso=0;
  std::vector<float>* eleHCalIso=0;
  std::vector<float>* eleECalIso=0;
  std::vector<float>* eleECalDriven=0;
  std::vector<float>* eleSCEn=0;
  std::vector<float>* eleSCRawEn=0;
  std::vector<float>* eleSCEta=0;
  std::vector<float>* eleSCPhi=0;
  std::vector<float>* eleBrem=0;
  std::vector<int>*   NClusters=0;

   
  tree->SetBranchStatus("nEle",1);     // enable electron branches
  tree->SetBranchStatus("ele*",1);     // enable electron branches
  tree->SetBranchAddress("nEle",&nEle);
  tree->SetBranchAddress("eleEn",&eleEn);
  tree->SetBranchAddress("eleCharge",&eleCharge);
  tree->SetBranchAddress("eleP",&eleP);
  tree->SetBranchAddress("elePt",&elePt);
  tree->SetBranchAddress("eleEta",&eleEta);
  tree->SetBranchAddress("elePhi",&elePhi);
  tree->SetBranchAddress("eleHoverE",&eleHoverE);
  tree->SetBranchAddress("eleSigmaIEtaIEta",&eleSigmaIEtaIEta);
  tree->SetBranchAddress("eleSigmaIEtaIEta_2012",&eleSigmaIEtaIEta_2012);
  tree->SetBranchAddress("eleSigmaIPhiIPhi",&eleSigmaIPhiIPhi);
  tree->SetBranchAddress("eleEoverP",&eleEoverP);
  tree->SetBranchAddress("eledEtaAtVtx",&eledEtaAtVtx);
  tree->SetBranchAddress("eledPhiAtVtx",&eledPhiAtVtx);
  tree->SetBranchAddress("eledEtaSeedAtVtx",&eledEtaSeedAtVtx);
  tree->SetBranchAddress("eleD0",&eleD0);
  tree->SetBranchAddress("eleDz",&eleDz);
  tree->SetBranchAddress("eleTrkPt",&eleTrkPt);
  tree->SetBranchAddress("eleMissHits",&eleMissHits);
  tree->SetBranchAddress("eleTrackIso",&eleTrackIso);
  tree->SetBranchAddress("eleHCalIso",&eleHCalIso);
  tree->SetBranchAddress("eleECalIso",&eleECalIso);
  tree->SetBranchAddress("eleECalDriven",&eleECalDriven);
  tree->SetBranchAddress("eleSCEn",&eleSCEn);
  tree->SetBranchAddress("eleSCRawEn",&eleSCRawEn);
  tree->SetBranchAddress("eleSCEta",&eleSCEta);
  tree->SetBranchAddress("eleSCPhi",&eleSCPhi);
  tree->SetBranchAddress("eleBrem",&eleBrem);
  tree->SetBranchAddress("NClusters",&NClusters);

  // RECO muons
  Int_t nMu;  
  tree->SetBranchStatus("nMu",1);     // enable muon branches
  tree->SetBranchAddress("nMu",&nMu);

  // calo tower
  Int_t nTower;
  std::vector<float>* CaloTower_hadE=0;
  std::vector<float>* CaloTower_emE=0;
  std::vector<float>* CaloTower_e=0;
  std::vector<float>* CaloTower_et=0;
  std::vector<float>* CaloTower_eta=0;
  std::vector<float>* CaloTower_phi=0;

  tree->SetBranchStatus("nTower",1);     // enable photon branches
  tree->SetBranchAddress("nTower",&nTower);
  tree->SetBranchStatus("Calo*",1);
  tree->SetBranchAddress("CaloTower_hadE",&CaloTower_hadE);
  tree->SetBranchAddress("CaloTower_emE",&CaloTower_emE);
  tree->SetBranchAddress("CaloTower_e",&CaloTower_e);
  tree->SetBranchAddress("CaloTower_et",&CaloTower_et);
  tree->SetBranchAddress("CaloTower_eta",&CaloTower_eta);
  tree->SetBranchAddress("CaloTower_phi",&CaloTower_phi);


  // pixel info  
  Int_t nhits1, nhits2, nhits3, nhits4, nhits5, nEv;
  treePixel->SetBranchStatus("nEv",1);     // enable photon branches
  treePixel->SetBranchAddress("nEv",&nEv);     // enable photon branches
  treePixel->SetBranchStatus("nhits*",1);     // enable pixel branches
  treePixel->SetBranchAddress("nhits1",&nhits1);     
  treePixel->SetBranchAddress("nhits2",&nhits2);     
  treePixel->SetBranchAddress("nhits3",&nhits3);   
  treePixel->SetBranchAddress("nhits4",&nhits4);     
  treePixel->SetBranchAddress("nhits5",&nhits5);  
  
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



  Int_t nPho;
  std::vector<float>* phoE=0;
  std::vector<float>* phoEt=0;
  std::vector<float>* phoEta=0;
  std::vector<float>* phoPhi=0;
  std::vector<float>* phoSCE=0;
  std::vector<float>* phoSCEt=0;
  std::vector<float>* phoSCEta=0;
  std::vector<float>* phoSCPhi=0;
  std::vector<float>* pho_ecalClusterIsoR4=0;
  std::vector<float>* pho_hcalRechitIsoR4=0;
  std::vector<float>* pho_trackIsoR4PtCut20=0;
  std::vector<float>* phoR9=0;
  std::vector<float>* phoHoverE=0;
  std::vector<float>* phoSigmaIEtaIEta=0;
  std::vector<float>* phoSigmaIEtaIEta_2012=0;
  std::vector<float>* phoSCRawE=0;
  std::vector<float>* phoE5x5=0;
  std::vector<float>* phoSCEtaWidth=0;
  std::vector<float>* phoSCPhiWidth=0;
  std::vector<float>* pho_swissCrx=0;
  std::vector<float>* pho_seedTime=0;
  std::vector<int>* phohasPixelSeed=0;
  std::vector<int>* phopassConversionVeto=0;

  tree->SetBranchStatus("nPho",1);     // enable photon branches
  tree->SetBranchStatus("pho*",1);     // enable photon branches
  tree->SetBranchAddress("nPho",&nPho);
  tree->SetBranchAddress("phoE",&phoE);
  tree->SetBranchAddress("phoEt",&phoEt);
  tree->SetBranchAddress("phoEta",&phoEta);
  tree->SetBranchAddress("phoPhi",&phoPhi);
  tree->SetBranchAddress("phoSCE",&phoSCE);
  tree->SetBranchAddress("phoSCEt",&phoSCEt);
  tree->SetBranchAddress("phoSCEta",&phoSCEta);
  tree->SetBranchAddress("phoSCPhi",&phoSCPhi);
  tree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
  tree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
  tree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
  tree->SetBranchAddress("phoR9",&phoR9);
  tree->SetBranchAddress("phoHoverE",&phoHoverE);
  tree->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
  tree->SetBranchAddress("phoSigmaIEtaIEta_2012",&phoSigmaIEtaIEta_2012);
  tree->SetBranchAddress("phoSCRawE",&phoSCRawE);
  tree->SetBranchAddress("phoE5x5",&phoE5x5);
  tree->SetBranchAddress("phoSCEtaWidth",&phoSCEtaWidth);
  tree->SetBranchAddress("phoSCPhiWidth",&phoSCPhiWidth);
  tree->SetBranchAddress("pho_swissCrx",&pho_swissCrx);
  tree->SetBranchAddress("pho_seedTime",&pho_seedTime);
  tree->SetBranchAddress("phohasPixelSeed",&phohasPixelSeed);
  tree->SetBranchAddress("phopassConversionVeto",&phopassConversionVeto);

  
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

  TTree *outputTreePho    = tree->CloneTree(0);
  outputTreePho->SetName("photons");
  outputTreePho->SetTitle("Event data + photons");
       
 
  TTree *outputTreediPho  = new TTree("diphoton","");
  TTree *outputTreediEle  = new TTree("dielectron","");
 
  outputTreeHLT->SetMaxTreeSize(MAXTREESIZE);
  outputTreediPho->SetMaxTreeSize(MAXTREESIZE);
  outputTreediEle->SetMaxTreeSize(MAXTREESIZE);

  ////**************************** write output branches      *************************/
  
  outputTreediPho->Branch("run",&run);
  outputTreediPho->Branch("event",&event);
  outputTreediPho->Branch("lumis",&lumis);

  outputTreediPho->Branch("nVtx",&nVtx);
  outputTreediPho->Branch("xVtx",&xVtx);
  outputTreediPho->Branch("yVtx",&yVtx);
  outputTreediPho->Branch("zVtx",&zVtx);
  outputTreediPho->Branch("isFakeVtx",&isFakeVtx);
  
  std::vector<int>   eleCharge_1;
  std::vector<float> eleP_1;
  std::vector<float> elePt_1;
  std::vector<float> eleEta_1;
  std::vector<float> elePhi_1;
  std::vector<float> eleHoverE_1;
  std::vector<float> eleSigmaIEtaIEta_1;
  std::vector<float> eleSigmaIEtaIEta_2012_1;
  std::vector<float> eleSigmaIPhiIPhi_1;
  std::vector<float> eleEoverP_1;
  std::vector<float> eledEtaAtVtx_1;
  std::vector<float> eledPhiAtVtx_1;
  std::vector<float> eledEtaSeedAtVtx_1;
  std::vector<float> eleD0_1;
  std::vector<float> eleDz_1;
  std::vector<int>   eleMissHits_1;
  std::vector<float> eleTrackIso_1;
  std::vector<float> eleHCalIso_1;
  std::vector<float> eleECalIso_1;
  std::vector<float> eleECalDriven_1;
  std::vector<float> eleBrem_1;
  
   
  std::vector<int>   eleCharge_2;
  std::vector<float> eleP_2;
  std::vector<float> elePt_2;
  std::vector<float> eleEta_2;
  std::vector<float> elePhi_2;
  std::vector<float> eleHoverE_2;
  std::vector<float> eleSigmaIEtaIEta_2;
  std::vector<float> eleSigmaIEtaIEta_2012_2;
  std::vector<float> eleSigmaIPhiIPhi_2;
  std::vector<float> eleEoverP_2;
  std::vector<float> eledEtaAtVtx_2;
  std::vector<float> eledPhiAtVtx_2;
  std::vector<float> eledEtaSeedAtVtx_2;
  std::vector<float> eleD0_2;
  std::vector<float> eleDz_2;
  std::vector<int>   eleMissHits_2;
  std::vector<float> eleTrackIso_2;
  std::vector<float> eleHCalIso_2;
  std::vector<float> eleECalIso_2;
  std::vector<float> eleECalDriven_2;
  std::vector<float> eleBrem_2;
 
  
  std::vector<float> vSum_ee_M;
  std::vector<float> vSum_ee_Energy;
  std::vector<float> vSum_ee_Pt;
  std::vector<float> vSum_ee_Eta;
  std::vector<float> vSum_ee_Phi;
  std::vector<float> vSum_ee_Rapidity;
  
  outputTreediPho->Branch("nEle",&nEle);
  outputTreediPho->Branch("nPho",&nPho);
  outputTreediPho->Branch("eleP_1",&eleP_1);
  outputTreediPho->Branch("elePt_1",&elePt_1);
  outputTreediPho->Branch("eleEta_1",&eleEta_1);
  outputTreediPho->Branch("elePhi_1",&elePhi_1);

  outputTreediPho->Branch("eleP_2",&eleP_2);
  outputTreediPho->Branch("elePt_2",&elePt_2);
  outputTreediPho->Branch("eleEta_2",&eleEta_2);
  outputTreediPho->Branch("elePhi_2",&elePhi_2);
  
  
  outputTreediPho->Branch("vSum_M",&vSum_ee_M);
  outputTreediPho->Branch("vSum_Energy",&vSum_ee_Energy);
  outputTreediPho->Branch("vSum_Pt",&vSum_ee_Pt);
  outputTreediPho->Branch("vSum_Eta",&vSum_ee_Eta);
  outputTreediPho->Branch("vSum_Phi",&vSum_ee_Phi);
  outputTreediPho->Branch("vSum_Rapidity",&vSum_ee_Rapidity);


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
  

  //EventMatcher* em = new EventMatcher();
  //Long64_t duplicateEntries = 0;
  
  Long64_t entries = tree->GetEntries();
  Long64_t entriesAnalyzed = 0;   Long64_t id = 0;   Long64_t excl = 0; Long64_t ptcut = 0; Long64_t myloop = 0;
  std::cout << "entries         = " << entries << std::endl;
  std::cout<< "Loop : ggHiNtuplizer/EventTree" <<std::endl;
  
  //for (Long64_t j_entry= 0; j_entry< 2000; ++j_entry){
  for (Long64_t j_entry=0; j_entry < entries; ++j_entry){
    
    //if (j_entry % 10000 == 0)  {
    //  std::cout << "current entry = " <<j_entry<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)j_entry/entries*100<<" %"<<std::endl;
    //}
    
    treeHLT->GetEntry(j_entry);
    tree->GetEntry(j_entry);
    treePixel->GetEntry(j_entry);
    
    if(HLT_HIUPCL1DoubleEG2NotHF2_v1==1 || HLT_HIUPCL1SingleEG5NotHF2_v1==1){   entriesAnalyzed++; 
      
      outputTreeHLT->Fill(); 
      
      // dielectron block
      eleCharge_1.clear();
      eleP_1.clear();
      elePt_1.clear();
      eleEta_1.clear();
      elePhi_1.clear();
     

      eleCharge_2.clear();
      eleP_2.clear();
      elePt_2.clear();
      eleEta_2.clear();
      elePhi_2.clear();
    
      
      vSum_ee_M.clear();
      vSum_ee_Energy.clear();
      vSum_ee_Pt.clear();
      vSum_ee_Eta.clear();
      vSum_ee_Phi.clear();
      vSum_ee_Rapidity.clear();

      phoEt_1.clear();
      phoEta_1.clear();
      phoPhi_1.clear();
            
      dr_1.clear();       
      dr_2.clear();       
  
      
    /*  if(nPho ==1 && ngenTrk == 2){
	
	for(int i=0; i<ngenTrk; ++i){
	  for (int j=i+1; j<ngenTrk; ++j){
	    
	    if (gentrkcharge->at(i) == gentrkcharge->at(j)) continue;  
	    cout << "Entry: "<< j_entry << "  nEle:" << nEle << "  muon:"<< nMu << "  Track 1 pt:" << gentrkPt->at(i) << "  Track 2 pt:"<< gentrkPt->at(j);//<< endl;
	    
	    for(int k=0; k<nPho; ++k){
	      cout <<"   Photon pt:"<< phoEt->at(k);// << endl; 
	      
	      double dr1 = getDR(gentrkEta->at(i),gentrkEta->at(i),phoEta->at(k),phoPhi->at(k));
	      double dr2 = getDR(gentrkEta->at(j),gentrkEta->at(j),phoEta->at(k),phoPhi->at(k));
	      cout << "  dr1:" << dr1 << "  dr2:"<< dr2 << endl;
	      
	      
	      TLorentzVector v1, v2, v3, vSum;
	      v1.SetPtEtaPhiM( gentrkPt->at(i), gentrkEta->at(i), gentrkPhi->at(i), eleMass);
	      v2.SetPtEtaPhiM( gentrkPt->at(j), gentrkEta->at(j), gentrkPhi->at(j), eleMass);
	      v3.SetPtEtaPhiE( phoEt->at(k), phoEta->at(k), phoPhi->at(k), phoE->at(k));
	      vSum = v1+v2+v3;
	      
	      if(vSum.Pt() < 1){cout <<"   three body-pt:" << vSum.Pt() << endl; id++; }
	      
	    }// npho
	  } //trk j
	} // trk i
	myloop++;
      }*/


      
      if(nPho ==1 && nEle == 2){
	
	for(int i=0; i<nEle; ++i){
	  for (int j=i+1; j<nEle; ++j){
	    
	    if (eleCharge->at(i) == eleCharge->at(j)) continue;  
	    cout << "Entry: "<< j_entry << "  nEle:" << nEle << "  muon:"<< nMu << "  Track 1 pt:" << elePt->at(i) << "  Track 2 pt:"<< elePt->at(j);//<< endl;
	    
	    for(int k=0; k<nPho; ++k){
	      cout <<"   Photon pt:"<< phoSCEt->at(k);// << endl; 
	      
	      double dr1 = getDR(eleEta->at(i),eleEta->at(i),phoSCEta->at(k),phoSCPhi->at(k));
	      double dr2 = getDR(eleEta->at(j),eleEta->at(j),phoSCEta->at(k),phoSCPhi->at(k));
	      cout << "  dr1:" << dr1 << "  dr2:"<< dr2 << endl;
	      
	      dr_1.push_back(dr1);       
	      dr_2.push_back(dr2);       
	      
	      TLorentzVector v1, v2, v3, vSum;
	      v1.SetPtEtaPhiM( elePt->at(i), eleEta->at(i), elePhi->at(i), eleMass);
	      v2.SetPtEtaPhiM( elePt->at(j), eleEta->at(j), elePhi->at(j), eleMass);
	      v3.SetPtEtaPhiE( phoSCEt->at(k), phoSCEta->at(k), phoSCPhi->at(k), phoSCE->at(k));
	      vSum = v1+v2+v3;
	      
	      eleP_1.push_back(eleP->at(i));
	      elePt_1.push_back(elePt->at(i));
	      eleEta_1.push_back(eleEta->at(i));
	      elePhi_1.push_back(elePhi->at(i));
	      
	      eleP_2.push_back(eleP->at(j));
	      elePt_2.push_back(elePt->at(j));
	      eleEta_2.push_back(eleEta->at(j));
	      elePhi_2.push_back(elePhi->at(j));
	      
	      vSum_ee_M.push_back(vSum.M());
	      vSum_ee_Energy.push_back(vSum.Energy());
	      vSum_ee_Pt.push_back(vSum.Pt());
	      vSum_ee_Eta.push_back(vSum.Eta());
	      vSum_ee_Phi.push_back(vSum.Phi());
	      vSum_ee_Rapidity.push_back(vSum.Rapidity());
	      
	      phoEt_1.push_back(phoSCEt->at(k));
	      phoEta_1.push_back(phoSCEta->at(k));
	      phoPhi_1.push_back(phoSCPhi->at(k));
	      
	      
	      if(vSum.Pt() < 1){cout <<"   three body-pt:" << vSum.Pt() << endl; id++; }
	      
	    }// npho
	  } //trk j
	} // trk i
	myloop++;
      }
      
      //outputTreePho->Fill();
      outputTreediPho->Fill();
      
    }
  }
  std::cout<<  "Loop ENDED : ggHiNtuplizer/EventTree" <<std::endl;
  std::cout << "entries            = " << entries << std::endl;
  std::cout << "In my loop         = " << myloop << std::endl;
  std::cout << "In three body loop = " << id << std::endl;
  //std::cout << "entriesAnalyzed  after trigger  = " << entriesAnalyzed << std::endl;
  //std::cout << "after 2 photon exclusive cut    = " << excl << std::endl;
  //std::cout << "for reconstruction number       = " << myloop << std::endl;
  //std::cout << "photons after applying ieta id  = " << id << std::endl;
  // std::cout << "outputTreeHLT->GetEntries()   = " << outputTreeHLT->GetEntries() << std::endl;
  //std::cout << "outputTreediPho->GetEntries()    = " << outputTreediPho->GetEntries() << std::endl;
  //std::cout << "outputTreediEle->GetEntries()    = " << outputTreediEle->GetEntries() << std::endl;
  
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
  check_data_driven_id_eff(argv[1]);
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


