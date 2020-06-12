// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"

#include "TNtuple.h"
#include "TRegexp.h"

using namespace std;
using namespace edm;

//
// class declaration
//

class TriggerObjectAnalyzer : public edm::EDAnalyzer {
public:
  explicit TriggerObjectAnalyzer(const edm::ParameterSet&);
  ~TriggerObjectAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  
  std::string   processName_;
  std::vector<std::string>   triggerNames_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;
  
  edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  std::string sl1name1 ;
  std::string sl1name2 ;
  
  HLTConfigProvider hltConfig_;
  
  unsigned int triggerIndex_;
  unsigned int moduleIndex_;
  string moduleLabel_;
  vector<string> moduleLabels_;
  
  edm::Service<TFileService> fs;
  vector<TTree*> nt_;
  TTree* mytree;
  int verbose_;
  
  std::map<std::string, bool> triggerInMenu;
  
  vector<double> id[500];
  vector<double> pt[500];
  vector<double> eta[500];
  vector<double> phi[500];
  vector<double> mass[500];

  vector<double> l1_doubleEG2_notHF2_pt_;
  vector<double> l1_doubleEG2_notHF2_eta_;
  vector<double> l1_doubleEG2_notHF2_phi_;
  vector<double> l1_singleEG5_notHF2_pt_;
  vector<double> l1_singleEG5_notHF2_eta_;
  vector<double> l1_singleEG5_notHF2_phi_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerObjectAnalyzer::TriggerObjectAnalyzer(const edm::ParameterSet& ps):
  processName_(ps.getParameter<std::string>("processName")),
  triggerNames_(ps.getParameter<std::vector<std::string> >("triggerNames")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent")),
  triggerResultsToken_(consumes<edm::TriggerResults>(triggerResultsTag_)),
  triggerEventToken_(consumes<trigger::TriggerEvent>(triggerEventTag_)),
  sl1name1(ps.getUntrackedParameter<string>("triggerL1Name1")),
  sl1name2(ps.getUntrackedParameter<string>("triggerL1Name2"))
{
  //now do what ever initialization is needed
  nt_.reserve(triggerNames_.size());
  for(unsigned int isize=0; isize<triggerNames_.size(); isize++){
    nt_[isize] = fs->make<TTree>(triggerNames_.at(isize).c_str(),Form("trigger %d",isize));
  }
   mytree = fs->make<TTree>("mytree", "Store L1 HF objects");
  
  mytree->Branch("l1_doubleEG2_notHF2_pt",      &l1_doubleEG2_notHF2_pt_);
  mytree->Branch("l1_doubleEG2_notHF2_eta",     &l1_doubleEG2_notHF2_eta_);
  mytree->Branch("l1_doubleEG2_notHF2_phi",     &l1_doubleEG2_notHF2_phi_);
  mytree->Branch("l1_singleEG5_notHF2_pt",      &l1_singleEG5_notHF2_pt_);
  mytree->Branch("l1_singleEG5_notHF2_eta",     &l1_singleEG5_notHF2_eta_);
  mytree->Branch("l1_singleEG5_notHF2_phi",     &l1_singleEG5_notHF2_phi_);  

  verbose_ = 0;
}


TriggerObjectAnalyzer::~TriggerObjectAnalyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerObjectAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  l1_doubleEG2_notHF2_pt_.clear();
  l1_doubleEG2_notHF2_eta_.clear();
  l1_doubleEG2_notHF2_phi_.clear();
  l1_singleEG5_notHF2_pt_.clear();
  l1_singleEG5_notHF2_eta_.clear();
  l1_singleEG5_notHF2_phi_.clear();

  if(hltConfig_.size() > 0){
    
    //float id = -99,pt=-99,eta=-99,phi=-99,mass=-99;
    
    using namespace edm;
    iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
    iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
    
    bool fl1 = false;
    bool fl2 = false;
    //hltConfig_.dump("Triggers");   
    ///cout <<" trigger size:" << triggerNames_.size() << " trigger result handle size:"<< triggerResultsHandle_->size() << endl; 
    
    for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){
      std::map<std::string,bool>::iterator inMenu = triggerInMenu.find(triggerNames_[itrig]);
      if (inMenu==triggerInMenu.end()){ continue; }
      
      triggerIndex_ = hltConfig_.triggerIndex(triggerNames_[itrig]);
      const unsigned int mIndex = triggerResultsHandle_->index(triggerIndex_);
      
      //cout << "even:" << iEvent.id() << "  Trigger index:"<<triggerIndex_ << "  mindex:"<< mIndex << "  trigger name:" << triggerNames_[itrig] << endl;
      
      // Results from TriggerEvent product - Attention: must look only for
      // modules actually run in this path for this event!
      
      for (unsigned int j=0; j<=mIndex; ++j) {
	// check whether the module is packed up in TriggerEvent product
	
	
	string trigFilterIndex = hltConfig_.moduleLabels(triggerIndex_).at(j); //this is simple to put into a loop to get all triggers...
        const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(trigFilterIndex,"",processName_)));
       // cout << "Filter index: "<<filterIndex << "  name:"<< trigFilterIndex << "  Event handle:" << triggerEventHandle_->sizeFilters() << endl;
	
	string moduleLabel = "";
        fl1 = false;
        fl2 = false;

      

	//for(unsigned int l = 0; l < trigFilterIndex.size(); l++){
	 //   moduleLabel = trigFilterIndex.at(l); 
         //   cout <<" modulelabel: " << moduleLabel << " sl1:" << sl1name1 << endl;
	 //   if(moduleLabel.find(sl1name1) != std::string::npos){ fl1 = true; cout << "event:" << iEvent.id() << "  fl1 true" << endl;}
	//}
	
	if(sl1name1 == trigFilterIndex){ fl1 = true;} //cout << "event:" << iEvent.id() << "  fl1 true" << endl;}
	if(sl1name2 == trigFilterIndex){ fl2 = true;} //cout << "event:" << iEvent.id() << "  fl2 true" << endl;}
	if(fl1 == true ){
	  const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl1name1,"",processName_)));
	  if (filterIndex < triggerEventHandle_->sizeFilters() ){
	    const trigger::Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
	    const trigger::Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
	    const unsigned int nI(VIDS.size());
	    const unsigned int nK(KEYS.size());
	    assert(nI==nK);
	    const unsigned int n2(max(nI,nK));
	    const trigger::TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
	    for (unsigned int i=0; i!=n2; ++i) { //cout << "coming in loop" << endl;
	      const trigger::TriggerObject& TO(TOC[KEYS[i]]);
	      l1_doubleEG2_notHF2_pt_  .push_back(TO.pt());
	      l1_doubleEG2_notHF2_eta_ .push_back(TO.eta());
	      l1_doubleEG2_notHF2_phi_ .push_back(TO.phi());
	     

	    }
	  } // filer size
	  //if(l1_doubleEG2_notHF2_pt_ == -1 ) l1_doubleEG2_notHF2_pt_ = 0; // this shouldn't happen...
	}
        //cout <<"fl1 bool:" << fl1 << "  fl2 bool:" << fl2 << endl;
	if(fl2 == true ){ //cout <<" In EG5 loop"<<endl;
	  const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(sl1name2,"",processName_)));
	  if (filterIndex < triggerEventHandle_->sizeFilters() ){ //cout <<" pass filter indexd loop"<<endl;
	    const trigger::Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
	    const trigger::Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
	    const unsigned int nI(VIDS.size());
	    const unsigned int nK(KEYS.size());
	    assert(nI==nK);  //cout <<" assert"<<endl;
	    const unsigned int n2(max(nI,nK)); 
	    const trigger::TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());  //cout <<" TOC  and n2:"<< n2 << endl;
	    for (unsigned int j=0; j!=n2; ++j) {  //cout <<" In n2 loop"<<endl;
	      const trigger::TriggerObject& TO(TOC[KEYS[j]]);
	      l1_singleEG5_notHF2_pt_  .push_back(TO.pt());
	      l1_singleEG5_notHF2_eta_ .push_back(TO.eta());
	      l1_singleEG5_notHF2_phi_ .push_back(TO.phi());
	      
             // cout << "single EG5 pT:"<< TO.pt() << endl;

	    }
	  } // filer size
	  //if(l1_singleEG5_notHF2_pt_ == -1 ) l1_singleEG5_notHF2_pt_ = 0; // this shouldn't happen...
	}
	
	
	if (filterIndex<triggerEventHandle_->sizeFilters()) {
      
	  //cout <<"  Filter index: "<<filterIndex << "  name:"<< trigFilterIndex << endl;
	  const trigger::Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
	  const trigger::Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
	  const unsigned int nI(VIDS.size());
	  const unsigned int nK(KEYS.size());
	  assert(nI==nK);
	  const unsigned int n(max(nI,nK));
	  
	  //cout <<" vidn:" << VIDS[n] << endl;
	  const trigger::TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
	  for (unsigned int i=0; i!=n; ++i) {
            //cout <<" vidi:" << VIDS[i] << endl;
	    const trigger::TriggerObject& TO(TOC[KEYS[i]]);
	    ///cout <<"trigger object id:"<<TO.id()<<endl;
	    //This check prevents grabbing the L1 trigger object (VIDS < 0), and finds the max trigger pt within all trigger collections
	    if(VIDS[i]>0){ // && pt<TO.pt()){
	      if(verbose_){
		cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
		     << TO.id() << " " << TO.pt() << " " << TO.et() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
		     << endl;
	      }
	      id[itrig].push_back(TO.id());
	      pt[itrig].push_back(TO.pt());
	      eta[itrig].push_back(TO.eta());
	      phi[itrig].push_back(TO.phi());
	      mass[itrig].push_back(TO.mass());
	    }
	  }
	}
      }
    }
    
    //nt_[0]->Fill(id,pt,eta,phi,mass);
  }
  for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){
    nt_[itrig]->Fill();
    id[itrig].clear();
    pt[itrig].clear();
    eta[itrig].clear();
    phi[itrig].clear();
    mass[itrig].clear();
  }
   mytree->Fill(); // fill even if no tri-lepton
}


// ------------ method called once each job just before starting event loop  ------------
void
TriggerObjectAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggerObjectAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
TriggerObjectAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      std::vector<std::string> activeHLTPathsInThisEvent = hltConfig_.triggerNames();
      
      triggerInMenu.clear(); 
      for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){
	for (std::vector<std::string>::const_iterator iHLT = activeHLTPathsInThisEvent.begin(); iHLT != activeHLTPathsInThisEvent.end(); ++iHLT){
	  //matching with regexp filter name. More than 1 matching filter is allowed so trig versioning is transparent to analyzer
	  if (TString(*iHLT).Contains(TRegexp(TString(triggerNames_[itrig])))){
	    triggerInMenu[*iHLT] = true;
	    triggerNames_[itrig] = TString(*iHLT);
	  }
	}
      }
      for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){
        std::map<std::string,bool>::iterator inMenu = triggerInMenu.find(triggerNames_[itrig]);
	if (inMenu==triggerInMenu.end()) {
	  cout << "<HLT Object Analyzer> Warning! Trigger " << triggerNames_[itrig] << " not found in HLTMenu. Skipping..." << endl;
        }
      }
      if(verbose_){
	hltConfig_.dump("ProcessName");
	hltConfig_.dump("GlobalTag");
	hltConfig_.dump("TableName");
	hltConfig_.dump("Streams");
	hltConfig_.dump("Datasets");
	hltConfig_.dump("PrescaleTable");
	hltConfig_.dump("ProcessPSet");
      }
    }
  } else {
    cout << "HLTObjectAnalyzer::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }
  
  for(unsigned int itrig=0; itrig<triggerNames_.size(); itrig++){
    nt_[itrig]->Branch("TriggerObjID",&(id[itrig]));
    nt_[itrig]->Branch("pt",&(pt[itrig]));
    nt_[itrig]->Branch("eta",&(eta[itrig]));
    nt_[itrig]->Branch("phi",&(phi[itrig]));
    nt_[itrig]->Branch("mass",&(mass[itrig]));
  }
}

// ------------ method called when ending the processing of a run  ------------
void
TriggerObjectAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TriggerObjectAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TriggerObjectAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerObjectAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerObjectAnalyzer);
