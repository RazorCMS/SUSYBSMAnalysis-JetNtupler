// -*- C++ -*-
// Class:      JetNtupler
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Caltech razor team
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#include "JetNtupler.h"
//------ Constructors and destructor ------//
JetNtupler::JetNtupler(const edm::ParameterSet& iConfig):
  //get inputs from config file
  isData_(iConfig.getParameter<bool> ("isData")),
  isFourJet_(iConfig.getParameter<bool> ("isFourJet")),
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  readGenVertexTime_(iConfig.getParameter<bool> ("readGenVertexTime")),
  triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),
  eleHLTFilterNamesFile_(iConfig.getParameter<string> ("eleHLTFilterNamesFile")),
  muonHLTFilterNamesFile_(iConfig.getParameter<string> ("muonHLTFilterNamesFile")),
  photonHLTFilterNamesFile_(iConfig.getParameter<string> ("photonHLTFilterNamesFile")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tracksTag_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  PFClustersToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfClusters"))),
  //genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  //genParticlesToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  hepMCToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMC"))),
  //triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  //triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
  metToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
//  metNoHFToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF"))),
  metPuppiToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsPuppi"))),
  metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  //hbheNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheNoiseFilter"))),
  //hbheTightNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheTightNoiseFilter"))),
  //hbheIsoNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheIsoNoiseFilter"))),
  //badChargedCandidateFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
  //badMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"))),
//  lheRunInfoTag_(iConfig.getParameter<edm::InputTag>("lheInfo")),
//  lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(lheRunInfoTag_)),
//  lheInfoToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator",""))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
  //hcalNoiseInfoToken_(consumes<HcalNoiseSummary>(iConfig.getParameter<edm::InputTag>("hcalNoiseInfo"))),
  secondaryVerticesToken_(consumes<vector<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
  rhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))),
  rhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll"))),
  rhoFastjetAllCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo"))),
  rhoFastjetCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo"))),
  rhoFastjetCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp"))),
  rhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
  esRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("esRecHits"))),
  ebeeClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("ebeeClusters"))),
  esClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("esClusters"))),
  conversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
  singleLegConversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("singleLegConversions"))),
  gedGsfElectronCoresToken_(consumes<vector<reco::GsfElectronCore> >(iConfig.getParameter<edm::InputTag>("gedGsfElectronCores"))),
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores")))
  //superClustersToken_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("superClusters"))),
  //  lostTracksToken_(consumes<vector<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks")))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;

  //set up output tree
  JetTree = fs->make<TTree>("Jets", "selected AOD information");
  //JetTree = new TTree("Jets", "selected AOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);
  //*****************************************************************************************
  //Read in HLT Trigger Path List from config file
  //*****************************************************************************************
  for (int i = 0; i<NTriggersMAX; ++i) triggerPathNames[i] = "";
  ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;
  if (myfile.is_open())
  {
    std::string line;
    int index;
    std::string hltpathname;

    while(myfile>>index>>hltpathname)
    {
      if (index < NTriggersMAX)
      {
        triggerPathNames[index] = hltpathname;
      }
    }
    myfile.close();
  }
  else
  {
    std::cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str() << "\n";
  }

  if(enableTriggerInfo_)
  {
    std::cout << "\n";
    std::cout << "****************** Trigger Paths Defined For Razor Ntuple ******************\n";
    for (int i = 0; i<NTriggersMAX; ++i)
    {
      if (triggerPathNames[i] != "") std::cout << "Trigger " << i << " " << triggerPathNames[i] << "\n";
    }
    std::cout << "****************************************************************************\n";
    std::cout << "\n";
  }
  if(readGenVertexTime_) genParticles_t0_Token_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genParticles_t0"));
  /*
  fJetPhotonRecHitEta = new std::vector<float>; fJetPhotonRecHitEta->clear();
  fJetPhotonRecHitPhi = new std::vector<float>; fJetPhotonRecHitPhi->clear();
  fJetPhotonRecHitE = new std::vector<float>; fJetPhotonRecHitE->clear();
  fJetPhotonRecHitTime = new std::vector<float>; fJetPhotonRecHitTime->clear();
*/
}

JetNtupler::~JetNtupler()
{
};

//------ Enable the desired set of branches ------//
void JetNtupler::setBranches(){


  
  JetTree->Branch("isData", &isData, "isData/O");
  JetTree->Branch("isFourJet", &isFourJet, "isFourJet/O");
  JetTree->Branch("runNum", &runNum, "runNum/i");
  JetTree->Branch("lumiNum", &lumiNum, "lumiNum/i");
  JetTree->Branch("eventNum", &eventNum, "eventNum/i");
  JetTree->Branch("pvX", &pvX, "pvX/F");
  JetTree->Branch("pvY", &pvY, "pvY/F");
  JetTree->Branch("pvZ", &pvZ, "pvZ/F");
  JetTree->Branch("nPV", &nPV, "nPV/I");
  JetTree->Branch("Rho", &Rho, "Rho/F");
  JetTree->Branch("nPU", &nPU, "nPU/I");
  JetTree->Branch("nPUmean", &nPUmean, "nPUmean/F");

  JetTree->Branch("nJets", &nJets,"nJets/I");
  JetTree->Branch("jetE", jetE,"jetE[nJets]/F");
  JetTree->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  JetTree->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  JetTree->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  JetTree->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  JetTree->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  JetTree->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  JetTree->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  JetTree->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  JetTree->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  JetTree->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  JetTree->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  JetTree->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  JetTree->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  JetTree->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  JetTree->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  JetTree->Branch("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction[nJets]/F");
  JetTree->Branch("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction[nJets]/F");
  JetTree->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F");
  JetTree->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F");
  JetTree->Branch("jetMatchedGenPt", jetMatchedGenPt,"jetMatchedGenPt[nJets]/F");
  JetTree->Branch("jetMatchedGenEta", jetMatchedGenEta,"jetMatchedGenEta[nJets]/F");
  JetTree->Branch("jetMatchedGenPhi", jetMatchedGenPhi,"jetMatchedGenPhi[nJets]/F");
  JetTree->Branch("jetMatchedGenMass", jetMatchedGenMass, "jetMatchedGenMass[nJets]/F");
  JetTree->Branch("jet_n_rechits", jet_n_rechits, "jet_n_rechits[nJets]/I");
  JetTree->Branch("jet_rechits_E", jet_rechits_E, "jet_rechits_E[nJets][1000]/F");
  JetTree->Branch("jet_rechits_T", jet_rechits_T, "jet_rechits_T[nJets][1000]/F");
  JetTree->Branch("jet_rechit_E_Ecut3", jet_rechit_E_Ecut3, "jet_rechit_E_Ecut3[nJets]/F");
  JetTree->Branch("jet_rechit_T_Ecut3", jet_rechit_T_Ecut3, "jet_rechit_T_Ecut3[nJets]/F");
  JetTree->Branch("jet_rechit_E_Ecut4", jet_rechit_E_Ecut4, "jet_rechit_E_Ecut4[nJets]/F");
  JetTree->Branch("jet_rechit_T_Ecut4", jet_rechit_T_Ecut4, "jet_rechit_T_Ecut4[nJets]/F");
  JetTree->Branch("jet_rechit_E_Ecut2", jet_rechit_E_Ecut2, "jet_rechit_E_Ecut2[nJets]/F");
  JetTree->Branch("jet_rechit_T_Ecut2", jet_rechit_T_Ecut2, "jet_rechit_T_Ecut2[nJets]/F");
  JetTree->Branch("jet_rechit_E_Ecut1p5", jet_rechit_E_Ecut1p5, "jet_rechit_E_Ecut1p5[nJets]/F");
  JetTree->Branch("jet_rechit_T_Ecut1p5", jet_rechit_T_Ecut1p5, "jet_rechit_T_Ecut1p5[nJets]/F");
  JetTree->Branch("jet_rechit_E_Ecut1", jet_rechit_E_Ecut1, "jet_rechit_E_Ecut1[nJets]/F");
  JetTree->Branch("jet_rechit_T_Ecut1", jet_rechit_T_Ecut1, "jet_rechit_T_Ecut1[nJets]/F");
  JetTree->Branch("jet_rechit_E_Ecut0p5", jet_rechit_E_Ecut0p5, "jet_rechit_E_Ecut0p5[nJets]/F");
  JetTree->Branch("jet_rechit_T_Ecut0p5", jet_rechit_T_Ecut0p5, "jet_rechit_T_Ecut0p5[nJets]/F");
  JetTree->Branch("jet_rechit_E", jet_rechit_E, "jet_rechit_E[nJets]/F");
  JetTree->Branch("jet_rechit_T", jet_rechit_T, "jet_rechit_T[nJets]/F");
  
  JetTree->Branch("jet_pv_rechits_T", jet_rechits_T, "jet_rechits_T[nJets][1000]/F");
  JetTree->Branch("jet_pv_rechit_T_Ecut3", jet_rechit_T_Ecut3, "jet_rechit_T_Ecut3[nJets]/F");
  JetTree->Branch("jet_pv_rechit_T_Ecut4", jet_rechit_T_Ecut4, "jet_rechit_T_Ecut4[nJets]/F");
  JetTree->Branch("jet_pv_rechit_T_Ecut2", jet_rechit_T_Ecut2, "jet_rechit_T_Ecut2[nJets]/F");
  JetTree->Branch("jet_pv_rechit_T_Ecut1p5", jet_rechit_T_Ecut1p5, "jet_rechit_T_Ecut1p5[nJets]/F");
  JetTree->Branch("jet_pv_rechit_T_Ecut1", jet_rechit_T_Ecut1, "jet_rechit_T_Ecut1[nJets]/F");
  JetTree->Branch("jet_pv_rechit_T_Ecut0p5", jet_rechit_T_Ecut0p5, "jet_rechit_T_Ecut0p5[nJets]/F");
  JetTree->Branch("jet_pv_rechit_T", jet_rechit_T, "jet_rechit_T[nJets]/F");

  JetTree->Branch("nPhotons", &fJetNPhotons,"nPhotons/I");
  JetTree->Branch("phoPt", fJetPhotonPt,"phoPt[nPhotons]/F");
  JetTree->Branch("phoEta", fJetPhotonEta,"phoEta[nPhotons]/F");
  JetTree->Branch("phoPhi", fJetPhotonPhi,"phoPhi[nPhotons]/F");
  JetTree->Branch("phoSeedRecHitEta", fJetPhotonSeedRecHitEta, "phoSeedRecHitEta[nPhotons]/F");
  JetTree->Branch("phoSeedRecHitPhi", fJetPhotonSeedRecHitPhi, "phoSeedRecHitPhi[nPhotons]/F");
  JetTree->Branch("phoSeedRecHitE", fJetPhotonSeedRecHitE, "phoSeedRecHitE[nPhotons]/F");
  JetTree->Branch("phoSeedRecHitT", fJetPhotonSeedRecHitTime, "phoSeedRecHitT[nPhotons]/F");

  // JetTree->Branch("fJetPhotonRecHitEta", "std::vector<float>",&fJetPhotonRecHitEta);
  // JetTree->Branch("fJetPhotonRecHitPhi", "std::vector<float>",&fJetPhotonRecHitPhi);
  // JetTree->Branch("fJetPhotonRecHitE", "std::vector<float>",&fJetPhotonRecHitE);
  // JetTree->Branch("fJetPhotonRecHitTime", "std::vector<float>",&fJetPhotonRecHitTime);

  cout << "BRANCHES\n";
  //enableMCBranches();
  enableGenParticleBranches();
  if (enableTriggerInfo_) enableTriggerBranches();
};


void JetNtupler::enableMCBranches(){
  JetTree->Branch("nGenJets", &nGenJets, "nGenJets/I");
  JetTree->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  JetTree->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  JetTree->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  JetTree->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  JetTree->Branch("genMetPt", &genMetPt, "genMetPt/F");
  JetTree->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
  JetTree->Branch("genVertexX", &genVertexX, "genVertexX/F");
  JetTree->Branch("genVertexY", &genVertexY, "genVertexY/F");
  JetTree->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  JetTree->Branch("genVertexT", &genVertexT, "genVertexT/F");
  JetTree->Branch("genWeight", &genWeight, "genWeight/F");
  JetTree->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  JetTree->Branch("genQScale", &genQScale, "genQScale/F");
  JetTree->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  JetTree->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
  /*scaleWeights = new std::vector<float>; scaleWeights->clear();
  pdfWeights = new std::vector<float>; pdfWeights->clear();
  alphasWeights = new std::vector<float>; alphasWeights->clear();
  if (isFastsim_) {
    JetTree->Branch("lheComments", "std::string",&lheComments);
  }
  JetTree->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  JetTree->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  JetTree->Branch("alphasWeights", "std::vector<float>",&alphasWeights);
  */
}

void JetNtupler::enableTriggerBranches()
{
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  JetTree->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
  //JetTree->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str());
  //JetTree->Branch("HLTMR", &HLTMR, "HLTMR/F");
  //JetTree->Branch("HLTRSQ", &HLTRSQ, "HLTRSQ/F");
};

void JetNtupler::enableGenParticleBranches(){
  JetTree->Branch("gLLP_prod_vertex_x", gLLP_prod_vertex_x, "gLLP_prod_vertex_x[2]/F");
  JetTree->Branch("gLLP_prod_vertex_y", gLLP_prod_vertex_y, "gLLP_prod_vertex_y[2]/F");
  JetTree->Branch("gLLP_prod_vertex_z", gLLP_prod_vertex_z, "gLLP_prod_vertex_z[2]/F");
  JetTree->Branch("gLLP_decay_vertex_x", gLLP_decay_vertex_x, "gLLP_decay_vertex_x[2]/F");
  JetTree->Branch("gLLP_decay_vertex_y", gLLP_decay_vertex_y, "gLLP_decay_vertex_y[2]/F");
  JetTree->Branch("gLLP_decay_vertex_z", gLLP_decay_vertex_z, "gLLP_decay_vertex_z[2]/F");
  JetTree->Branch("gLLP_beta", gLLP_beta, "gLLP_beta[2]/F");
  JetTree->Branch("gLLP_travel_time", gLLP_travel_time, "gLLP_travel_time[2]/F");

  JetTree->Branch("gLLP_daughter_travel_time", gLLP_daughter_travel_time, "gLLP_daughter_travel_time[4]/F");
  JetTree->Branch("gLLP_daughter_pt", gLLP_daughter_pt, "gLLP_daughter_pt[4]/F");
  JetTree->Branch("gLLP_daughter_eta", gLLP_daughter_eta, "gLLP_daughter_eta[4]/F");
  JetTree->Branch("gLLP_daughter_phi", gLLP_daughter_phi, "gLLP_daughter_phi[4]/F");
  JetTree->Branch("gLLP_daughter_e", gLLP_daughter_e, "gLLP_daughter_e[4]/F");
  JetTree->Branch("gLLP_daughter_match_jet_index", gLLP_daughter_match_jet_index, "gLLP_daughter_match_jet_index[4]/i");
  JetTree->Branch("gLLP_min_delta_r_match_jet", gLLP_min_delta_r_match_jet, "gLLP_min_delta_r_match_jet[4]/F");

  JetTree->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  JetTree->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  JetTree->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  JetTree->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  JetTree->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  JetTree->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  JetTree->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  JetTree->Branch("gParticlePx", gParticlePx, "gParticlePx[nGenParticle]/F");
  JetTree->Branch("gParticlePy", gParticlePy, "gParticlePy[nGenParticle]/F");
  JetTree->Branch("gParticlePz", gParticlePz, "gParticlePz[nGenParticle]/F");
  JetTree->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  JetTree->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
  JetTree->Branch("gParticleDecayVertexX", gParticleDecayVertexX, "gParticleDecayVertexX[nGenParticle]/F");
  JetTree->Branch("gParticleDecayVertexY", gParticleDecayVertexY, "gParticleDecayVertexY[nGenParticle]/F");
  JetTree->Branch("gParticleDecayVertexZ", gParticleDecayVertexZ, "gParticleDecayVertexZ[nGenParticle]/F");
}



//------ Load the miniAOD objects and reset tree variables for each event ------//
void JetNtupler::loadEvent(const edm::Event& iEvent){
  //load all miniAOD objects for the current event
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(hepMCToken_, hepMC);
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(metFilterBitsToken_, metFilterBits);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(tracksTag_,tracks);
  iEvent.getByToken(PFCandsToken_, pfCands);
  iEvent.getByToken(PFClustersToken_, pfClusters);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(tausToken_, taus);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsPuppiToken_, jetsPuppi);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  iEvent.getByToken(metToken_, mets);
  //iEvent.getByToken(metNoHFToken_, metsNoHF);
  iEvent.getByToken(metPuppiToken_, metsPuppi);
//  iEvent.getByToken(hcalNoiseInfoToken_,hcalNoiseInfo);
  iEvent.getByToken(secondaryVerticesToken_,secondaryVertices);
  iEvent.getByToken(rhoAllToken_,rhoAll);
  iEvent.getByToken(rhoFastjetAllToken_,rhoFastjetAll);
  iEvent.getByToken(rhoFastjetAllCaloToken_,rhoFastjetAllCalo);
  iEvent.getByToken(rhoFastjetCentralCaloToken_,rhoFastjetCentralCalo);
  iEvent.getByToken(rhoFastjetCentralChargedPileUpToken_,rhoFastjetCentralChargedPileUp);
  iEvent.getByToken(rhoFastjetCentralNeutralToken_,rhoFastjetCentralNeutral);
  iEvent.getByToken(beamSpotToken_,beamSpot);
  iEvent.getByToken(ebRecHitsToken_,ebRecHits);
  iEvent.getByToken(eeRecHitsToken_,eeRecHits);
  iEvent.getByToken(esRecHitsToken_,esRecHits);
  iEvent.getByToken(ebeeClustersToken_,ebeeClusters);
  iEvent.getByToken(esClustersToken_,esClusters);
  iEvent.getByToken(conversionsToken_,conversions);
  iEvent.getByToken(singleLegConversionsToken_,singleLegConversions);
  iEvent.getByToken(gedGsfElectronCoresToken_,gedGsfElectronCores);
  iEvent.getByToken(gedPhotonCoresToken_, gedPhotonCores);
//  iEvent.getByToken(superClustersToken_,superClusters);
//  iEvent.getByToken(lostTracksToken_,lostTracks);
//  iEvent.getByToken(hbheNoiseFilterToken_, hbheNoiseFilter);
//  iEvent.getByToken(hbheTightNoiseFilterToken_, hbheTightNoiseFilter);
//  iEvent.getByToken(hbheIsoNoiseFilterToken_, hbheIsoNoiseFilter);
  //iEvent.getByToken(badChargedCandidateFilterToken_, badChargedCandidateFilter);
  //iEvent.getByToken(badMuonFilterToken_, badMuonFilter);
  if(readGenVertexTime_) iEvent.getByToken(genParticles_t0_Token_,genParticles_t0);

  if (useGen_) {
//    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);

    //for Spring16 fastsim, this has been changed and removed
//    if (!isFastsim_) iEvent.getByToken(lheInfoToken_, lheInfo);

    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }


}

//called by the loadEvent() method
void JetNtupler::resetBranches(){
    //reset tree variables
    //Event
    eventNum = 0;
    lumiNum = 0;
    runNum = 0;
    pvX = -99.0;
    pvY = -99.0;
    pvZ = -99.0;
    nPV = -1;
    Rho = -99.0;
    nPUmean = -1;
    nPU = -1;

    //Photon
    fJetNPhotons = 0;
    for (int i=0; i< OBJECTARRAYSIZE; i++) {
      fJetPhotonPt[i] = 0.0;
      fJetPhotonEta[i] = 0.0;
      fJetPhotonPhi[i] = 0.0;
      fJetPhotonSeedRecHitE[i]      = -99.0;
      fJetPhotonSeedRecHitEta[i]      = -99.0;
      fJetPhotonSeedRecHitPhi[i]      = -99.0;
      fJetPhotonSeedRecHitTime[i]      = -99.0;
    }

/*
    fJetPhotonRecHitE->clear();
    fJetPhotonRecHitEta->clear();
    fJetPhotonRecHitPhi->clear();
    fJetPhotonRecHitTime->clear();
*/
    //Jet
    nJets = 0;
    for ( int i = 0; i < 1000; i++)
    {
      jetE[i] = 0.0;
      jetPt[i] = 0.0;
      jetEta[i] = 0.0;
      jetPhi[i] = 0.0;
      jetCISV[i] = 0.0;
      jetMass[i] =  -99.0;
      jetJetArea[i] = -99.0;
      jetPileupE[i] = -99.0;
      jetPileupId[i] = -99.0;
      jetPileupIdFlag[i] = -1;
      jetPassIDLoose[i] = false;
      jetPassIDTight[i] = false;
      jetPassMuFrac[i] = false;
      jetPassEleFrac[i] = false;
      jetPartonFlavor[i] = 0;
      jetHadronFlavor[i] = 0;
      jetChargedEMEnergyFraction[i] = -99.0;
      jetNeutralEMEnergyFraction[i] = -99.0;
      jetChargedHadronEnergyFraction[i] = -99.0;
      jetNeutralHadronEnergyFraction[i] = -99.0;
      jetMatchedGenPt[i] = 0.0;
      jetMatchedGenEta[i] = 0.0;
      jetMatchedGenPhi[i] = 0.0;
      jetMatchedGenMass[i] = 0.0;
      jetMatchedGenTime[i] = 0.0;
      jet_n_rechits[i] = 0;
      jet_rechit_E[i] = 0.0;
      jet_rechit_T[i] = 0.0;
      jet_rechit_E_Ecut4[i] = 0.0; //energy with a 2 GeV cut
      jet_rechit_T_Ecut4[i] = 0.0;
      jet_rechit_E_Ecut2[i] = 0.0; //energy with a 2 GeV cut
      jet_rechit_T_Ecut2[i] = 0.0;
      jet_rechit_E_Ecut1p5[i] = 0.0; //energy with a 2 GeV cut
      jet_rechit_T_Ecut1p5[i] = 0.0;
      jet_rechit_E_Ecut1[i] = 0.0; //energy with a 2 GeV cut
      jet_rechit_T_Ecut1[i] = 0.0;
      jet_rechit_E_Ecut0p5[i] = 0.0; //energy with a 2 GeV cut
      jet_rechit_T_Ecut0p5[i] = 0.0;
      

      jet_pv_rechit_T[i] = 0.0;
      jet_pv_rechit_T_Ecut4[i] = 0.0;
      jet_pv_rechit_T_Ecut2[i] = 0.0;
      jet_pv_rechit_T_Ecut1p5[i] = 0.0;
      jet_pv_rechit_T_Ecut1[i] = 0.0;
      jet_pv_rechit_T_Ecut0p5[i] = 0.0;

      for(int j =0;j<1000;j++)
      {
        jet_rechits_E[i][j] = -666.;
        jet_rechits_T[i][j] = -666.;
	jet_pv_rechits_T[i][j] = -666.;

      }
   }

    for ( int i = 0; i < 2; i++ )
    {
      gLLP_prod_vertex_x[i] = -666.;
      gLLP_prod_vertex_y[i] = -666.;
      gLLP_prod_vertex_z[i] = -666.;
      gLLP_decay_vertex_x[i] = -666.;
      gLLP_decay_vertex_y[i] = -666.;
      gLLP_decay_vertex_z[i] = -666.;
      gLLP_beta[i] = -666.;
      gLLP_travel_time[i] = -666.;
    }

    for ( int i = 0; i < 4; i++ )
    {
      gLLP_daughter_pt[i] = -666.;
      gLLP_daughter_eta[i] = -666.;
      gLLP_daughter_phi[i] = -666.;
      gLLP_daughter_e[i] = -666.;
      gLLP_daughter_travel_time[i] = -666.;
      gLLP_daughter_match_jet_index[i] = 666;
      gLLP_min_delta_r_match_jet[i] = -666.;
    }
}

//------ Methods to fill tree variables ------//




//------ Method called for each run ------//

void JetNtupler::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {


}


//------ Method called for each lumi block ------//
void JetNtupler::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {

}


//------ Method called for each event ------//

void JetNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  //initialize
  loadEvent(iEvent); //loads objects and resets tree branches
  NEvents->Fill(0); //increment event count
  //resetting output tree branches
  resetBranches();


  //*************************************
  //Fill Event-Level Info
  //*************************************

  //store basic event info
  isData = isData_;
  isFourJet = isFourJet_;
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();

 //select the primary vertex, if any
  nPV = 0;
  myPV = &(vertices->front());

  bool foundPV = false;
  for(unsigned int i = 0; i < vertices->size(); i++)
  {
    if(vertices->at(i).isValid() && !vertices->at(i).isFake())
    {
      if (!foundPV)
      {
        myPV = &(vertices->at(i));
        foundPV = true;
      }
      nPV++;
    }
  }

  pvX = myPV->x();
  pvY = myPV->y();
  pvZ = myPV->z();

  //get rho
  Rho = *rhoFastjetAll;

  //Fill Pileup info
  if (!isData)
  {
    for(const PileupSummaryInfo &pu : *puInfo)
    {
      if ( pu.getBunchCrossing() == 0)
      {
        nPU = pu.getPU_NumInteractions();
        nPUmean = pu.getTrueNumInteractions();
      }
    }
  }

  int i_jet = 0;
  for (const reco::PFJet &j : *jets)
  {
    //resetBranches();
    if (j.pt() < 10) continue;
    if (fabs(j.eta()) > 2.4) continue;
    //*************************************
    //Fill Jet-Level Info
    //*************************************
    jetE[i_jet] = j.energy();
    jetPt[i_jet] = j.pt();
    jetEta[i_jet] = j.eta();
    jetPhi[i_jet] = j.phi();
    jetMass[i_jet] = j.mass();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(jetPt[i_jet], jetEta[i_jet], jetPhi[i_jet], jetE[i_jet]);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    jetJetArea[i_jet] = j.jetArea();
    jetPileupE[i_jet] = j.pileup();

    jetPileupIdFlag[i_jet] = 0;
    jetPassIDLoose[i_jet] = passJetID(&j, 0);
    jetPassIDTight[i_jet] = passJetID(&j, 1);
    jetPassMuFrac[i_jet]  = ( j.muonEnergyFraction() < 0.80 );
    jetPassEleFrac[i_jet]  = ( j.electronEnergyFraction() < 0.90 );


    // if (useGen_) {
    //   jetPartonFlavor = j.partonFlavour();
    //   jetHadronFlavor = j.hadronFlavour();
    // }

    jetChargedEMEnergyFraction[i_jet] = j.chargedEmEnergyFraction();
    jetNeutralEMEnergyFraction[i_jet] = j.neutralEmEnergyFraction();
    jetChargedHadronEnergyFraction[i_jet] = j.chargedHadronEnergyFraction();
    jetNeutralHadronEnergyFraction[i_jet] = j.neutralHadronEnergyFraction();


    //*************************************
    //find photons inside the jet
    //*************************************
    for (const reco::Photon &pho : *photons) {
      //cout << "Nphoton: " << fJetNPhotons << "\n";

      if (!(deltaR(pho.eta(), pho.phi() , j.eta(), j.phi()) < 0.5)) continue;


      fJetPhotonPt[fJetNPhotons]  = pho.pt();
      fJetPhotonEta[fJetNPhotons] = pho.eta(); //correct this for the vertex
      fJetPhotonPhi[fJetNPhotons] = pho.phi(); //correct this for the vertex

      fJetPhotonSeedRecHitE[fJetNPhotons]      = pho.superCluster()->seed()->x();
      fJetPhotonSeedRecHitEta[fJetNPhotons]      = pho.superCluster()->seed()->y();
      fJetPhotonSeedRecHitPhi[fJetNPhotons]      = pho.superCluster()->seed()->z();
      fJetPhotonSeedRecHitTime[fJetNPhotons]      = pho.superCluster()->seed()->energy();

      // //get time coordinate for the seed
      // for (const reco::PFCluster &pfcluster : *pfClusters) {
      // 	if(pfcluster.seed() == pho.superCluster()->seed()->seed())
      // 	  {
      // 	    pho_superClusterSeedT[fJetNPhotons] = pfcluster.time();
      // 	    pho_pfClusterSeedE[fJetNPhotons]      = pfcluster.energy();
      // 	  }
      // }

      //*************************************
      //fill all rechits inside photons
      //*************************************

      fJetNPhotons++;

    }
    //***************************
    //Find RecHits Inside the Jet
    //***************************
    // geometry (from ECAL ELF)

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    const CaloSubdetectorGeometry *barrelGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    const CaloSubdetectorGeometry *endcapGeometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    double ecal_radius = 129.0;
    int n_matched_rechits = 0;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
      {
        if ( recHit->checkFlag(0) )
        {
          const DetId recHitId = recHit->detid();
          const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
          if ( deltaR(jetEta[i_jet], jetPhi[i_jet], recHitPos.eta(), recHitPos.phi())  < 0.4)
          {
            jet_rechit_E[i_jet] += recHit->energy();
            jet_rechit_T[i_jet] += recHit->time()*recHit->energy();
            jet_rechits_E[i_jet][n_matched_rechits] = recHit->energy();
	    jet_rechits_T[i_jet][n_matched_rechits] = recHit->time();
	    double rechit_x = ecal_radius * cos(recHitPos.phi());
	    double rechit_y = ecal_radius * sin(recHitPos.phi());
	    double rechit_z = ecal_radius * sinh(recHitPos.eta());
	    double photon_pv_travel_time = (1./30) * sqrt(pow(pvX-rechit_x,2)+pow(pvY-rechit_y,2)+pow(pvZ-rechit_z,2));
	    jet_pv_rechits_T[i_jet][n_matched_rechits] = recHit->time()+(1./30)*ecal_radius*cosh(recHitPos.eta()) - photon_pv_travel_time;
	    jet_pv_rechit_T[i_jet] += recHit->energy()*jet_pv_rechits_T[i_jet][n_matched_rechits]; 
	  
            if (recHit->energy() > 0.5)
	    {
		jet_rechit_E_Ecut0p5[i_jet] += recHit->energy();
		jet_rechit_T_Ecut0p5[i_jet] += recHit->time()*recHit->energy();
                jet_pv_rechit_T_Ecut0p5[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

	    }
            if (recHit->energy() > 1.0)
	    {
		jet_rechit_E_Ecut1[i_jet] += recHit->energy();
		jet_rechit_T_Ecut1[i_jet] += recHit->time()*recHit->energy();
                jet_pv_rechit_T_Ecut1[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

	    }
            if (recHit->energy() > 1.5)
	    {
		jet_rechit_E_Ecut1p5[i_jet] += recHit->energy();
		jet_rechit_T_Ecut1p5[i_jet] += recHit->time()*recHit->energy();
                jet_pv_rechit_T_Ecut1p5[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

	    }
            if (recHit->energy() > 2.0)
	    {
		jet_rechit_E_Ecut2[i_jet] += recHit->energy();
		jet_rechit_T_Ecut2[i_jet] += recHit->time()*recHit->energy();
                jet_pv_rechit_T_Ecut2[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

	    }
	    if (recHit->energy() > 3.0)
            {
                jet_rechit_E_Ecut3[i_jet] += recHit->energy();
                jet_rechit_T_Ecut3[i_jet] += recHit->time()*recHit->energy();
                jet_pv_rechit_T_Ecut3[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

            }

	    if (recHit->energy() > 4.0)
            {
                jet_rechit_E_Ecut4[i_jet] += recHit->energy();
                jet_rechit_T_Ecut4[i_jet] += recHit->time()*recHit->energy();
                jet_pv_rechit_T_Ecut4[i_jet] += jet_pv_rechits_T[i_jet][n_matched_rechits] *recHit->energy();

            }
	    n_matched_rechits++;
        
          //std::cout << recHitPos.eta() << std::endl;
            }
        //std::cout << recHitId << std::endl;
        }
    }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    jet_n_rechits[i_jet] = n_matched_rechits;
    jet_rechit_T[i_jet] = jet_rechit_T[i_jet]/jet_rechit_E[i_jet];
    jet_rechit_T_Ecut4[i_jet] = jet_rechit_T_Ecut4[i_jet]/jet_rechit_E_Ecut4[i_jet];
    jet_rechit_T_Ecut3[i_jet] = jet_rechit_T_Ecut3[i_jet]/jet_rechit_E_Ecut3[i_jet];
    jet_rechit_T_Ecut2[i_jet] = jet_rechit_T_Ecut2[i_jet]/jet_rechit_E_Ecut2[i_jet];
    jet_rechit_T_Ecut1p5[i_jet] = jet_rechit_T_Ecut1p5[i_jet]/jet_rechit_E_Ecut1p5[i_jet];
    jet_rechit_T_Ecut1[i_jet] = jet_rechit_T_Ecut1[i_jet]/jet_rechit_E_Ecut1[i_jet];
    jet_rechit_T_Ecut0p5[i_jet] = jet_rechit_T_Ecut0p5[i_jet]/jet_rechit_E_Ecut0p5[i_jet]; //incrementing jet counter
    jet_pv_rechit_T[i_jet] = jet_pv_rechit_T[i_jet]/jet_rechit_E[i_jet];
    jet_pv_rechit_T_Ecut4[i_jet] = jet_pv_rechit_T_Ecut4[i_jet]/jet_rechit_E_Ecut4[i_jet];
    jet_pv_rechit_T_Ecut3[i_jet] = jet_pv_rechit_T_Ecut3[i_jet]/jet_rechit_E_Ecut3[i_jet];
    jet_pv_rechit_T_Ecut2[i_jet] = jet_pv_rechit_T_Ecut2[i_jet]/jet_rechit_E_Ecut2[i_jet];
    jet_pv_rechit_T_Ecut1p5[i_jet] = jet_pv_rechit_T_Ecut1p5[i_jet]/jet_rechit_E_Ecut1p5[i_jet];
    jet_pv_rechit_T_Ecut1[i_jet] = jet_pv_rechit_T_Ecut1[i_jet]/jet_rechit_E_Ecut1[i_jet];
    jet_pv_rechit_T_Ecut0p5[i_jet] = jet_pv_rechit_T_Ecut0p5[i_jet]/jet_rechit_E_Ecut0p5[i_jet]; //incrementing jet counter
    nJets++;
    i_jet++;
  } //loop over jets

  //MC AND GEN LEVEL INFO
  fillGenParticles();
  //fillMC();

  if ( enableTriggerInfo_ ) fillTrigger( iEvent );
  JetTree->Fill();
}

//------ Method called once each job just before starting event loop ------//
void JetNtupler::beginJob(){
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//
void JetNtupler::endJob(){
}


bool JetNtupler::passJetID( const reco::PFJet *jet, int cutLevel) {
  bool result = false;

  double NHF = jet->neutralHadronEnergyFraction();
  double NEMF = jet->neutralEmEnergyFraction();
  int NumConst = jet->chargedMultiplicity() + jet->neutralMultiplicity() ;
  double CHF = jet->chargedHadronEnergyFraction();
  double MUF = jet->muonEnergyFraction();
  double CEMF = jet->chargedEmEnergyFraction();
  int NumNeutralParticles =jet->neutralMultiplicity();
  int CHM = jet->chargedMultiplicity();

  //Loose
  if (cutLevel == 0) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  //Tight
  else if (cutLevel == 1) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  //Tight Lep Veto
  else if (cutLevel == 2) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 && MUF < 0.8 ) result = true;
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;
    }
  }

  return result;
}

double JetNtupler::deltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
  {
    dphi -= TMath::TwoPi();
  }
  while (dphi <= -TMath::Pi())
  {
    dphi += TMath::TwoPi();
  }
  return dphi;
};

double JetNtupler::deltaR(double eta1, double phi1, double eta2, double phi2)
{
double dphi = deltaPhi(phi1,phi2);
double deta = eta1 - eta2;
return sqrt( dphi*dphi + deta*deta);
};


bool JetNtupler::fillMC()
{
  for(const reco::GenJet &j : *genJets)
  {
    genJetE[nGenJets] = j.energy();
    genJetPt[nGenJets] = j.pt();
    genJetEta[nGenJets] = j.eta();
    genJetPhi[nGenJets] = j.phi();
    nGenJets++;
  }

  const pat::MET &Met = mets->front();
  genMetPt = Met.genMET()->pt();
  genMetPhi = Met.genMET()->phi();

  bool foundGenVertex = false;
  for(size_t i=0; i<genParticles->size();i++)
  {
    if (!foundGenVertex)
    {
      for (unsigned int j=0; j<(*genParticles)[i].numberOfDaughters(); ++j)
      {
        const reco::Candidate *dau = (*genParticles)[i].daughter(j);
        if (dau)
        {
          genVertexX = dau->vx();
          genVertexY = dau->vy();
          genVertexZ = dau->vz();
          // if(readGenVertexTime_) genVertexT = *genParticles_t0;
          foundGenVertex = true;
          break;
        }
      }
    }
  }

  genWeight = genInfo->weight();
  genSignalProcessID = genInfo->signalProcessID();
  genQScale = genInfo->qScale();
  genAlphaQCD = genInfo->alphaQCD();
  genAlphaQED = genInfo->alphaQED();

    /*
    if (isFastsim_) {

      //get lhe weights for systematic uncertainties:
      double nomlheweight = genInfo->weights()[0];

      //fill scale variation weights
      if (genInfo->weights().size()>=10) {
	for (unsigned int iwgt=1; iwgt<10; ++iwgt) {
	  //normalize to
	  double wgtval = genInfo->weights()[iwgt]*genWeight/genInfo->weights()[1];
	  scaleWeights->push_back(wgtval);
	}
      }

      //fill pdf variation weights
      if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(genInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {

	//fill pdf variation weights after converting with mc2hessian transformation
	std::array<double, 100> inpdfweights;
	for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
	  inpdfweights[ipdf] = genInfo->weights()[iwgt]/genInfo->weights()[firstPdfWeight-1];
	}

	std::array<double, 60> outpdfweights;
	pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());

	for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
	  double wgtval = outpdfweights[iwgt]*genWeight;
	  pdfWeights->push_back(wgtval);
	}

	//fill alpha_s variation weights
	if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(genInfo->weights().size())) {
	  for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	    double wgtval = genInfo->weights()[iwgt]*genWeight/nomlheweight;
	    alphasWeights->push_back(wgtval);
	  }
	}

      }
    } else {

      if (lheInfo.isValid() && lheInfo->weights().size()>0) {

	double nomlheweight = lheInfo->weights()[0].wgt;

	//fill scale variation weights
	if (lheInfo->weights().size()>=9) {
	  for (unsigned int iwgt=0; iwgt<9; ++iwgt) {
	    double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	    scaleWeights->push_back(wgtval);
	  }
	}

	//fill pdf variation weights
	if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(lheInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {

	  //fill pdf variation weights after converting with mc2hessian transformation
	  std::array<double, 100> inpdfweights;
	  for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
	    inpdfweights[ipdf] = lheInfo->weights()[iwgt].wgt/nomlheweight;
	  }

	  std::array<double, 60> outpdfweights;
	  pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());

	  for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
	    double wgtval = outpdfweights[iwgt]*genWeight;
	    pdfWeights->push_back(wgtval);
	  }

	  //fill alpha_s variation weights
	  if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(lheInfo->weights().size())) {
	    for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	      double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	      alphasWeights->push_back(wgtval);
	    }
	  }
	}
      }
    }

    //fill sum of weights histograms
    sumWeights->Fill(0.,genWeight);

    for (unsigned int iwgt=0; iwgt<scaleWeights->size(); ++iwgt) {
      sumScaleWeights->Fill(double(iwgt),(*scaleWeights)[iwgt]);
    }
    for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
      sumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);
    }
    for (unsigned int iwgt=0; iwgt<alphasWeights->size(); ++iwgt) {
      sumAlphasWeights->Fill(double(iwgt),(*alphasWeights)[iwgt]);
    }
*/
    return true;
};

bool JetNtupler::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
  //Fills selected gen particles
  for(size_t i=0; i<genParticles->size();i++)
  {
    if( (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 && ( (*genParticles)[i].status() < 30 ))
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21 && (*genParticles)[i].status() < 30)
       || (abs((*genParticles)[i].pdgId()) >= 22 && abs((*genParticles)[i].pdgId()) <= 25 && ( (*genParticles)[i].status() < 30))
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039) )
       {
         prunedV.push_back(&(*genParticles)[i]);
       }
     }

  //Total number of gen particles
  nGenParticle = prunedV.size();

  //Look for mother particle and Fill gen variables
  for(unsigned int i = 0; i < prunedV.size(); i++)
  {
    gParticleId[i] = prunedV[i]->pdgId();
    gParticleStatus[i] = prunedV[i]->status();
    gParticleE[i] = prunedV[i]->energy();
    gParticlePt[i] = prunedV[i]->pt();
    gParticlePx[i] = prunedV[i]->px();
    gParticlePy[i] = prunedV[i]->py();
    gParticlePz[i] = prunedV[i]->pz();
    gParticleEta[i] = prunedV[i]->eta();
    gParticlePhi[i] = prunedV[i]->phi();
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;

    //For Neutralinos we try to find the decay vertex locaton.
    //Algorithm: Find the first daughter particle that is not a neutralino,
    //and call that the daughter. get the creation vertex of that daughter.
    if ( (gParticleId[i] == 1000022 && gParticleStatus[i] == 22) )
    {
      const reco::Candidate *dau = 0;
      bool foundDaughter = false;
      bool noDaughter = false;
      const reco::Candidate *tmpParticle = prunedV[i];

      while (!foundDaughter && !noDaughter)
      {
        if (tmpParticle->numberOfDaughters() > 0)
        {
          dau = tmpParticle->daughter(0);
          if (dau && dau->pdgId() != 1000022)
          {
            foundDaughter = true;
          } else
          {
            tmpParticle = dau;
          }
        }
        else
        {
          noDaughter = true;
        }
      }

      if (foundDaughter)
      {
        gParticleDecayVertexX[i] = dau->vx();
        gParticleDecayVertexY[i] = dau->vy();
        gParticleDecayVertexZ[i] = dau->vz();
      }
    }


    if(prunedV[i]->numberOfMothers() > 0)
    {
      //find the ID of the first mother that has a different ID than the particle itself
      const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      if (firstMotherWithDifferentID)
      {
        gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
      }

      //find the mother and keep going up the mother chain if the ID's are the same
      const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
      for(unsigned int j = 0; j < prunedV.size(); j++)
      {
        if(prunedV[j] == originalMotherWithSameID)
        {
          gParticleMotherIndex[i] = j;
          break;
        }
      }
    }
    else
    {
      gParticleMotherIndex[i] = -1;
    }
    //LLPs vextex
    if (isFourJet)
    {
      if ( (gParticleId[i] == 35 || gParticleId[i] == 36) && gParticleStatus[i] == 22 )
      {
        if (gParticleId[i] == 35)
        {
          gLLP_prod_vertex_x[0] = prunedV[i]->vx();
          gLLP_prod_vertex_y[0] = prunedV[i]->vy();
          gLLP_prod_vertex_z[0] = prunedV[i]->vz();
        //std::cout << prunedV[i]->tau() << std::endl;
        }   
        else if (gParticleId[i] == 36)
        {
          gLLP_prod_vertex_x[1] = prunedV[i]->vx();
          gLLP_prod_vertex_y[1] = prunedV[i]->vy();
          gLLP_prod_vertex_z[1] = prunedV[i]->vz();
        }

        const reco::Candidate *dau = 0;
        bool foundDaughter = false;
        bool noDaughter = false;
        const reco::Candidate *tmpParticle = prunedV[i];

      //std::cout << "tmpParticle->numberOfDaughters(): " << tmpParticle->numberOfDaughters() << std::endl;

        while (!foundDaughter && !noDaughter)
        {
          if (tmpParticle->numberOfDaughters() > 0)
          {
            dau = tmpParticle->daughter(0);
            if (dau && (dau->pdgId() != 35 && dau->pdgId() != 36))
            {
              foundDaughter = true;
            } else
            {
              tmpParticle = dau;
            }
          }
          else
          {
            noDaughter = true;
          }
        }

        if (foundDaughter)
        {
        //gParticleDecayVertexX[i] = dau->vx();
        //gParticleDecayVertexY[i] = dau->vy();
        //gParticleDecayVertexZ[i] = dau->vz();
          if (gParticleId[i] == 35)
          {
            gLLP_decay_vertex_x[0] = dau->vx();
            gLLP_decay_vertex_y[0] = dau->vy();
            gLLP_decay_vertex_z[0] = dau->vz();
            gLLP_beta[0] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
            gLLP_travel_time[0] = sqrt(pow(gLLP_decay_vertex_x[0]-gLLP_prod_vertex_x[0],2)
                                    +pow(gLLP_decay_vertex_y[0]-gLLP_prod_vertex_y[0],2)
                                    +pow(gLLP_decay_vertex_z[0]-gLLP_prod_vertex_z[0],2))/(30. * gLLP_beta[0]);//1/30 is to convert cm to ns
            double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
            double ecal_radius = 129.0;
          //gLLP_decays_px[0]

          //std::cout << tmpParticle->pdgId() << " number of daughters: " << tmpParticle->numberOfDaughters() << std::endl;
          /*
          First two LLP daughters belong to LLP->pdgID()=35
          */

            for (unsigned int id = 0; id < tmpParticle->numberOfDaughters(); id++ )
            {
            //std::cout << "====================" << std::endl;
            //std::cout << " -> "<< tmpParticle->daughter(id)->pdgId() << std::endl;
              if( id > 1 ) break;
              TLorentzVector tmp;
              tmp.SetPxPyPzE(tmpParticle->daughter(id)->px(), tmpParticle->daughter(id)->py(), tmpParticle->daughter(id)->pz(), tmpParticle->daughter(id)->energy());
              gLLP_daughter_pt[id] = tmp.Pt();
              gLLP_daughter_eta[id] = tmp.Eta();
              gLLP_daughter_phi[id] = tmp.Phi();
              gLLP_daughter_e[id]  = tmp.E();


              gLLP_daughter_travel_time[id] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            //Calculate dt from generation point to ECAL face
              double x_ecal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id];
              double y_ecal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id];
              double z_ecal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id];
              if( fabs(z_ecal) < 271.6561246934 )
              {
    	        double photon_travel_time = (1./30) * sqrt(pow(ecal_radius,2)+pow((gLLP_decay_vertex_z[0] + (ecal_radius-radius) * sinh(tmp.Eta())),2));
                gLLP_daughter_travel_time[id] = gLLP_daughter_travel_time[id] - photon_travel_time;          
	      //std::cout << "(x,y,z) @ ecal = (" << x_ecal << "," << y_ecal << "," << z_ecal << ")" << std::endl;
              //std::cout << "extrapolated r = " << sqrt(pow(x_ecal,2)+pow(y_ecal,2)) << std::endl;
              }
              else
              {
                gLLP_daughter_travel_time[id] = -666;
              }

	      double min_delta_r = 666.;
	      unsigned int match_jet_index = 666;
	      for ( int i_jet = 0; i_jet < nJets; i_jet++ )
	      {
		double current_delta_r = deltaR(tmp.Eta(), tmp.Phi() , jetEta[i_jet], jetPhi[i_jet]);
	      //std::cout << i_jet << " current dR = " << current_delta_r << std::endl;
	      if ( current_delta_r < min_delta_r )
	      {
		min_delta_r = current_delta_r;
		match_jet_index = i_jet;
		//std::cout << i_jet << " min dR = " << min_delta_r << std::endl;
	      }
	    }//end matching to jets
	    if ( min_delta_r < 0.3 )
	    {
	      gLLP_daughter_match_jet_index[id] = match_jet_index;
	      gLLP_min_delta_r_match_jet[id] = min_delta_r;
	      //std::cout << "min dR = " << min_delta_r << " matched to jet index " << match_jet_index << std::endl;
	    }
	  }
	}
	  else if (gParticleId[i] == 36)
	  {
	    gLLP_decay_vertex_x[1] = dau->vx();
	    gLLP_decay_vertex_y[1] = dau->vy();
	    gLLP_decay_vertex_z[1] = dau->vz();
	    gLLP_beta[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
	    gLLP_travel_time[1] = sqrt(pow(gLLP_decay_vertex_x[1]-gLLP_prod_vertex_x[1],2)
				      +pow(gLLP_decay_vertex_y[1]-gLLP_prod_vertex_y[1],2)
				      +pow(gLLP_decay_vertex_z[1]-gLLP_prod_vertex_z[1],2))/(30. * gLLP_beta[1]);//1/30 is to convert cm to ns
	    double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
	    double ecal_radius = 129.0;
	    /*
	    Second two LLP daughters belong to LLP->pdgID()=36
	    */
	    for (unsigned int id = 0; id < tmpParticle->numberOfDaughters(); id++ )
	    {
	      //std::cout << " -> "<< tmpParticle->daughter(id)->pdgId() << std::endl;
	      if( id > 1 ) break;
	      TLorentzVector tmp;
	      tmp.SetPxPyPzE(tmpParticle->daughter(id)->px(), tmpParticle->daughter(id)->py(), tmpParticle->daughter(id)->pz(), tmpParticle->daughter(id)->energy());
	      gLLP_daughter_pt[id+2] = tmp.Pt();
	      gLLP_daughter_eta[id+2] = tmp.Eta();
	      gLLP_daughter_phi[id+2] = tmp.Phi();
	      gLLP_daughter_e[id+2]  = tmp.E();
	      //gLLP_daughter_travel_time[id+2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E()) - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
	      gLLP_daughter_travel_time[id+2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns

	      //Calculate dt from generation point to ECAL face
	      double x_ecal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[id+2];
	      double y_ecal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[id+2];
	      double z_ecal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[id+2];
	      if( fabs(z_ecal) < 271.6561246934 )
	      {
		double photon_travel_time = (1./30) * sqrt(pow(ecal_radius,2)+pow((gLLP_decay_vertex_z[1] + (ecal_radius-radius) * sinh(tmp.Eta())),2));
		gLLP_daughter_travel_time[id+2] = gLLP_daughter_travel_time[id+2] - photon_travel_time;         //std::cout << "(x,y,z) @ ecal = (" << x_ecal << "," << y_ecal << "," << z_ecal << ")" << std::endl;
		//std::cout << "extrapolated r = " << sqrt(pow(x_ecal,2)+pow(y_ecal,2)) << std::endl;
	      }
	      else
	      {
		gLLP_daughter_travel_time[id+2] = -666;
	      }

	      double min_delta_r = 666;
	      unsigned int match_jet_index = 666;
	      for ( int i_jet = 0; i_jet < nJets; i_jet++ )
	      {
		double current_delta_r = deltaR(tmp.Eta(), tmp.Phi() , jetEta[i_jet], jetPhi[i_jet]);
		if ( current_delta_r < min_delta_r )
		{
		  min_delta_r = current_delta_r;
		  match_jet_index = i_jet;
		}
	      }//end matching to jets
	      if ( min_delta_r < 0.3 )
	      {
		gLLP_daughter_match_jet_index[id+2] = match_jet_index;
		gLLP_min_delta_r_match_jet[id+2] = min_delta_r;
	      }
	    }
	  }
	}
      }
    }
    else
    {
      if (abs(gParticleId[i]) == 5 || abs(gParticleId[i]) == 12)
      {
	if (gParticleMotherId[i] == 9000006 || gParticleMotherId[i] == 9000007)
	{
	  if (gParticleId[i] == 5)
	  {
	    const reco::Candidate *tmpParticle = prunedV[i];
	    TLorentzVector tmp;
	    tmp.SetPxPyPzE(tmpParticle->px(), tmpParticle->py(), tmpParticle->pz(), tmpParticle->energy());
	    gLLP_daughter_pt[0] = tmp.Pt();
	    
	    gLLP_decay_vertex_x[0] = tmpParticle->vx();
            gLLP_decay_vertex_y[0] = tmpParticle->vy();
            gLLP_decay_vertex_z[0] = tmpParticle->vz();
	    double ecal_radius = 129.0;
	    double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
	    gLLP_daughter_eta[0] = tmp.Eta();
	    gLLP_daughter_phi[0] = tmp.Phi();
	    gLLP_daughter_e[0]  = tmp.E();    
	    
	    gLLP_daughter_travel_time[0] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
	    double x_ecal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[0];
	    double y_ecal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[0];
	    double z_ecal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[0];
	    if( fabs(z_ecal) < 271.6561246934 )
	    {
	      double photon_travel_time = (1./30) * sqrt(pow(ecal_radius,2)+pow((gLLP_decay_vertex_z[0] + (ecal_radius-radius) * sinh(tmp.Eta())),2));
	      gLLP_daughter_travel_time[0] = gLLP_daughter_travel_time[0] - photon_travel_time;
	    }
	    else
	    {
	      gLLP_daughter_travel_time[0] = -666;
	    }
	    double min_delta_r = 666;
	    unsigned int match_jet_index = 666;
	    for ( int i_jet = 0; i_jet < nJets; i_jet++ )
	    {
	      double current_delta_r = deltaR(tmp.Eta(), tmp.Phi() , jetEta[i_jet], jetPhi[i_jet]);
	      if ( current_delta_r < min_delta_r )
	      {
		min_delta_r = current_delta_r;
		match_jet_index = i_jet;
	      }
	    }//end matching to jets
	    if ( min_delta_r < 0.3 )
	    {
	      gLLP_daughter_match_jet_index[0] = match_jet_index;
	      gLLP_min_delta_r_match_jet[0] = min_delta_r;
	    }
	  }
	  else if (gParticleId[i] == -5)
	  {
	    const reco::Candidate *tmpParticle = prunedV[i];
            TLorentzVector tmp;
            tmp.SetPxPyPzE(tmpParticle->px(), tmpParticle->py(), tmpParticle->pz(), tmpParticle->energy());
            gLLP_daughter_pt[1] = tmp.Pt();
	    gLLP_daughter_eta[1] = tmp.Eta();
            gLLP_daughter_phi[1] = tmp.Phi();
            gLLP_daughter_e[1]  = tmp.E();
	    double radius = sqrt( pow(gLLP_decay_vertex_x[0],2) + pow(gLLP_decay_vertex_y[0],2) );
	    double ecal_radius = 129.0;
	    gLLP_daughter_travel_time[1] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            double x_ecal = gLLP_decay_vertex_x[0] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[1];
            double y_ecal = gLLP_decay_vertex_y[0] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[1];
            double z_ecal = gLLP_decay_vertex_z[0] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[1];
            if( fabs(z_ecal) < 271.6561246934 )
            {
              double photon_travel_time = (1./30) * sqrt(pow(ecal_radius,2)+pow((gLLP_decay_vertex_z[0] + (ecal_radius-radius) * sinh(tmp.Eta())),2));
              gLLP_daughter_travel_time[1] = gLLP_daughter_travel_time[1] - photon_travel_time;
            }
            else
            {
              gLLP_daughter_travel_time[1] = -666;
            }
            double min_delta_r = 666;
            unsigned int match_jet_index = 666;
            for ( int i_jet = 0; i_jet < nJets; i_jet++ )
            {
              double current_delta_r = deltaR(tmp.Eta(), tmp.Phi() , jetEta[i_jet], jetPhi[i_jet]);
              if ( current_delta_r < min_delta_r )
              {
                min_delta_r = current_delta_r;
                match_jet_index = i_jet;
              }
            }//end matching to jets
            if ( min_delta_r < 0.3 )
            {
              gLLP_daughter_match_jet_index[1] = match_jet_index;
              gLLP_min_delta_r_match_jet[1] = min_delta_r;
            }
	  }
	  else if (gParticleId[i] == 12)
          {
            const reco::Candidate *tmpParticle = prunedV[i];
            TLorentzVector tmp;
            tmp.SetPxPyPzE(tmpParticle->px(), tmpParticle->py(), tmpParticle->pz(), tmpParticle->energy());
            gLLP_daughter_pt[2] = tmp.Pt();

            gLLP_decay_vertex_x[1] = tmpParticle->vx();
            gLLP_decay_vertex_y[1] = tmpParticle->vy();
            gLLP_decay_vertex_z[1] = tmpParticle->vz();
            double ecal_radius = 129.0;
	    double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
            gLLP_daughter_eta[2] = tmp.Eta();
            gLLP_daughter_phi[2] = tmp.Phi();
            gLLP_daughter_e[2]  = tmp.E();

            gLLP_daughter_travel_time[2] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            double x_ecal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[2];
            double y_ecal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[2];
            double z_ecal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[2];
	    if( fabs(z_ecal) < 271.6561246934 )
            {
              double photon_travel_time = (1./30) * sqrt(pow(ecal_radius,2)+pow((gLLP_decay_vertex_z[1] + (ecal_radius-radius) * sinh(tmp.Eta())),2));
              gLLP_daughter_travel_time[2] = gLLP_daughter_travel_time[2] - photon_travel_time;
            }
            else
            {
              gLLP_daughter_travel_time[2] = -666;
            }
	  }
	  else if (gParticleId[i] == -12)
          {
            const reco::Candidate *tmpParticle = prunedV[i];
            TLorentzVector tmp;
            tmp.SetPxPyPzE(tmpParticle->px(), tmpParticle->py(), tmpParticle->pz(), tmpParticle->energy());
            gLLP_daughter_pt[3] = tmp.Pt();
            double ecal_radius = 129.0;
	    double radius = sqrt( pow(gLLP_decay_vertex_x[1],2) + pow(gLLP_decay_vertex_y[1],2) );
            gLLP_daughter_eta[3] = tmp.Eta();
            gLLP_daughter_phi[3] = tmp.Phi();
            gLLP_daughter_e[3]  = tmp.E();

            gLLP_daughter_travel_time[3] = (1./30.)*(ecal_radius-radius)/(tmp.Pt()/tmp.E());// - (1./30.) * ecal_radius * cosh(tmp.Eta());//1/30 is to convert cm to ns
            double x_ecal = gLLP_decay_vertex_x[1] + 30. * (tmp.Px()/tmp.E())*gLLP_daughter_travel_time[3];
            double y_ecal = gLLP_decay_vertex_y[1] + 30. * (tmp.Py()/tmp.E())*gLLP_daughter_travel_time[3];
            double z_ecal = gLLP_decay_vertex_z[1] + 30. * (tmp.Pz()/tmp.E())*gLLP_daughter_travel_time[3];
            if( fabs(z_ecal) < 271.6561246934 )
            {
              double photon_travel_time = (1./30) * sqrt(pow(ecal_radius,2)+pow((gLLP_decay_vertex_z[1] + (ecal_radius-radius) * sinh(tmp.Eta())),2));
              gLLP_daughter_travel_time[3] = gLLP_daughter_travel_time[3] - photon_travel_time;
            }
            else
            {
              gLLP_daughter_travel_time[3] = -666;
            }
          }
	}

      }	





    }


  }

  return true;
};



bool JetNtupler::fillTrigger(const edm::Event& iEvent)
{

  //fill trigger information
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;
  //********************************************************************
  //Option to save all HLT path names in the ntuple per event
  //Expensive option in terms of ntuple size
  //********************************************************************
  nameHLT->clear();
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
    string hltPathNameReq = "HLT_";
    //if (triggerBits->accept(i))
    if ((names.triggerName(i)).find(hltPathNameReq) != string::npos) nameHLT->push_back(names.triggerName(i));
    /*
    std::cout << "Trigger " << names.triggerName(i) <<
    ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    if ((names.triggerName(i)).find(hltPathNameReq) != string::npos && triggerBits->accept(i)) std::cout << "Trigger " << names.triggerName(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    */
  }
  //std::cout << "n triggers: " <<  nameHLT->size() << std::endl;
  //std::cout << "====================" << std::endl;
  //for ( unsigned int i = 0; i < nameHLT->size(); i++ )
  //{
  //  std::cout << i << " -> " << nameHLT->at(i) << std::endl;
  //}
  //********************************************************************
  // Save trigger decisions in array of booleans
  //********************************************************************

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i)
  {
    string hltPathNameReq = "HLT_";
    if ((names.triggerName(i)).find(hltPathNameReq) == string::npos) continue;
    if ((names.triggerName(i)).find_last_of("_") == string::npos) continue;
    int lastUnderscorePos = (names.triggerName(i)).find_last_of("_");
    string hltPathNameWithoutVersionNumber = (names.triggerName(i)).substr(0,lastUnderscorePos);

    for (unsigned int j = 0; j < NTriggersMAX; ++j)
    {
      if (triggerPathNames[j] == "") continue;
      if (hltPathNameWithoutVersionNumber == triggerPathNames[j])
      {
        triggerDecision[j] = triggerBits->accept(i);
        //triggerHLTPrescale[j] = triggerPrescales->getPrescaleForIndex(i);
      }
    }
  }

  //********************************************************************
  // Print Trigger Objects
  //********************************************************************
/*
  for (pat::TriggerObjectStandAlone trigObject : *triggerObjects)
  {
    //cout << "triggerObj: " << trigObject.pt() << " " << trigObject.eta() << " " << trigObject.phi() << "\n";
    //bool foundRazor = false;
    //Need to unpack the filter labels before checking
    trigObject.unpackFilterLabels(iEvent, *triggerBits);
    for(int j=0; j<int(trigObject.filterLabels().size());j++)
    {
      //if ((trigObject.filterLabels())[j] == "hltRsqMR200Rsq0p0196MR100Calo") foundRazor = true;
      // trigObject.unpackPathNames(names);
      // cout << "filter: " << (trigObject.pathNames())[j] << " " << (trigObject.filterLabels())[j] << "\n";
      //cout << "filter: " << (trigObject.filterLabels())[j] << "\n";
    }
  }
*/
//define this as a plug-in
  return true;
};
DEFINE_FWK_MODULE(JetNtupler);
