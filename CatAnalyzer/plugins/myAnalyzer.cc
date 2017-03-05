#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "CATTools/CatAnalyzer/interface/analysisUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TTree.h"
#include "TH1D.h"

using namespace std;
using namespace cat;

class myAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  explicit myAnalyzer(const edm::ParameterSet&);
  ~myAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override{};
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};
  
  edm::EDGetTokenT<int>  nGoodVertexToken_;
  edm::EDGetTokenT<float> puweightToken_, puweightToken_up_, puweightToken_dn_;
  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> elecToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;  

  TTree* ttree_;
  TH1D * h_nevents;
  float b_met, b_weight, b_puweight, b_puweight_up, b_puweight_dn;

  TLorentzVector b_genlep1; int b_genlep1_pid;
  TLorentzVector b_genlep2; int b_genlep2_pid;
  TLorentzVector b_gendilep;
  
  TLorentzVector b_lep1; int b_lep1_pid;
  TLorentzVector b_lep2; int b_lep2_pid;
  TLorentzVector b_dilep;
};
//
// constructors and destructor
//
myAnalyzer::myAnalyzer(const edm::ParameterSet& iConfig)
{
  nGoodVertexToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("nGoodVertex"));
  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  puweightToken_up_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_up"));
  puweightToken_dn_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_dn"));
  muonToken_ = consumes<cat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<cat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<cat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
  vtxToken_  = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  mcLabel_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  h_nevents = fs->make<TH1D>("nevents","nevents",1,0,1);       
  ttree_ = fs->make<TTree>("nom","nom");
    
  ttree_->Branch("met", &b_met, "met/F");
  ttree_->Branch("lep1", "TLorentzVector", &b_lep1);
  ttree_->Branch("lep1_pid", &b_lep1_pid, "lep1_pid/I");    
  ttree_->Branch("lep2", "TLorentzVector", &b_lep2);
  ttree_->Branch("lep2_pid", &b_lep2_pid, "lep2_pid/I");    
  ttree_->Branch("dilep", "TLorentzVector", &b_dilep);
    
  ttree_->Branch("weight", &b_weight, "weight/F");
  ttree_->Branch("puweight", &b_puweight, "puweight/F");
  ttree_->Branch("puweight_up", &b_puweight_up, "puweight_up/F");
  ttree_->Branch("puweight_dn", &b_puweight_dn, "puweight_dn/F");
  
  ttree_->Branch("genlep1", "TLorentzVector", &b_genlep1);
  ttree_->Branch("genlep1_pid", &b_genlep1_pid, "genlep1_pid/I");    
  ttree_->Branch("genlep2", "TLorentzVector", &b_genlep2);
  ttree_->Branch("genlep2_pid", &b_genlep2_pid, "genlep2_pid/I");    
  ttree_->Branch("gendilep", "TLorentzVector", &b_gendilep);
}

myAnalyzer::~myAnalyzer(){}

void myAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{    
  const bool runOnMC = !iEvent.isRealData();
  if (runOnMC){
    edm::Handle<float> puweightHandle;
    iEvent.getByToken(puweightToken_, puweightHandle);
    b_puweight = *puweightHandle;
    
    edm::Handle<float> puweightHandle_up;
    iEvent.getByToken(puweightToken_up_, puweightHandle_up);
    b_puweight_up = *puweightHandle_up;
    
    edm::Handle<float> puweightHandle_dn;
    iEvent.getByToken(puweightToken_dn_, puweightHandle_dn);
    b_puweight_dn = *puweightHandle_dn;
    
    b_weight = b_puweight;
    
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(mcLabel_,genParticles);
    bool bosonSample = false;
    TLorentzVector genLep1;
    TLorentzVector genLep2;
    for (const reco::GenParticle & g : *genParticles){
      if (abs(g.pdgId())!=13){
	continue;
      }
      bool isfromBoson = false;
      for (unsigned int i = 0; i < g.numberOfMothers(); ++i){
	//In case of pdgId() = 23, indicate Z-boson. if it's 25, that becomes higgs.
	if (g.mother(i)->pdgId() == 23 || g.mother(i)->pdgId() == 25){
	  isfromBoson = true;
	  bosonSample = true;
	}
      }
      if (isfromBoson){
	if (g.charge() > 0) genLep1.SetPtEtaPhiM(g.pt(), g.eta(), g.phi(), g.mass());
	else genLep2.SetPtEtaPhiM(g.pt(), g.eta(), g.phi(), g.mass());
      }
    }
    if (bosonSample){
      b_genlep1 = genLep1;
      b_genlep2 = genLep2;
      b_genlep1_pid = -13;
      b_genlep2_pid = 13;
      b_gendilep = b_genlep1 + b_genlep2;
    }    
  }
  
  h_nevents->Fill(0.5,b_weight);

  ////////////////////////////////////////////////////////////////////////////////
  // filters
  ////////////////////////////////////////////////////////////////////////////////
  // trigger
  ////////////////////////////////////////////////////////////////////////////////
  // vertex
  ////////////////////////////////////////////////////////////////////////////////
  // get physics objects
  edm::Handle<cat::MuonCollection> muons;          iEvent.getByToken(muonToken_, muons);
  edm::Handle<cat::ElectronCollection> electrons;  iEvent.getByToken(elecToken_, electrons);
  edm::Handle<cat::JetCollection> jets;            iEvent.getByToken(jetToken_, jets);
  edm::Handle<cat::METCollection> mets;            iEvent.getByToken(metToken_, mets);
  ////////////////////////////////////////////////////////////////////////////////
  // muon selection
  cat::MuonCollection selectedMuons;
  for (auto& m : *muons) {
    cat::Muon mu(m);
    if (mu.pt() < 20.) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    if (!mu.isTightMuon()) continue;
    if (mu.relIso(0.4) > 0.15) continue;
    printf("muon with pt %4.1f, POG loose id %d, tight id %d\n", mu.pt(), mu.isLooseMuon(), mu.isTightMuon());
    selectedMuons.push_back(mu);
  }
  if (selectedMuons.size() > 1){
    b_lep1 = selectedMuons[0].tlv();
    b_lep2 = selectedMuons[1].tlv();
  }
  ttree_->Fill();
}

//define this as a plug-in
DEFINE_FWK_MODULE(myAnalyzer);
