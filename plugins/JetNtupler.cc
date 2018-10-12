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
  useGen_(iConfig.getParameter<bool> ("useGen")),
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
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

  /*
  fJetPhotonRecHitEta = new std::vector<float>; fJetPhotonRecHitEta->clear();
  fJetPhotonRecHitPhi = new std::vector<float>; fJetPhotonRecHitPhi->clear();
  fJetPhotonRecHitE = new std::vector<float>; fJetPhotonRecHitE->clear();
  fJetPhotonRecHitTime = new std::vector<float>; fJetPhotonRecHitTime->clear();
*/
}

JetNtupler::~JetNtupler()
{
}

//------ Enable the desired set of branches ------//
void JetNtupler::setBranches(){


  JetTree->Branch("isData", &isData, "isData/O");
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


  JetTree->Branch("jetE", &jetE,"jetE/F");
  JetTree->Branch("jetPt", &jetPt,"jetPt/F");
  JetTree->Branch("jetEta", &jetEta,"jetEta/F");
  JetTree->Branch("jetPhi", &jetPhi,"jetPhi/F");
  JetTree->Branch("jetCISV", &jetCISV,"jetCISV/F");
  JetTree->Branch("jetMass", &jetMass, "jetMass/F");
  JetTree->Branch("jetJetArea", &jetJetArea, "jetJetArea/F");
  JetTree->Branch("jetPileupE", &jetPileupE, "jetPileupE/F");
  JetTree->Branch("jetPileupId", &jetPileupId, "jetPileupId/F");
  JetTree->Branch("jetPileupIdFlag", &jetPileupIdFlag, "jetPileupIdFlag/I");
  JetTree->Branch("jetPassIDLoose", &jetPassIDLoose, "jetPassIDLoose/O");
  JetTree->Branch("jetPassIDTight", &jetPassIDTight, "jetPassIDTight/O");
  JetTree->Branch("jetPassMuFrac", &jetPassMuFrac, "jetPassMuFrac/O");
  JetTree->Branch("jetPassEleFrac", &jetPassEleFrac, "jetPassEleFrac/O");
  JetTree->Branch("jetPartonFlavor", &jetPartonFlavor, "jetPartonFlavor/I");
  JetTree->Branch("jetHadronFlavor", &jetHadronFlavor, "jetHadronFlavor/I");
  JetTree->Branch("jetChargedEMEnergyFraction", &jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction/F");
  JetTree->Branch("jetNeutralEMEnergyFraction", &jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction/F");
  JetTree->Branch("jetChargedHadronEnergyFraction", &jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction/F");
  JetTree->Branch("jetNeutralHadronEnergyFraction", &jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction/F");
  JetTree->Branch("jetMatchedGenPt", &jetMatchedGenPt,"jetMatchedGenPt/F");
  JetTree->Branch("jetMatchedGenEta", &jetMatchedGenEta,"jetMatchedGenEta/F");
  JetTree->Branch("jetMatchedGenPhi", &jetMatchedGenPhi,"jetMatchedGenPhi/F");
  JetTree->Branch("jetMatchedGenMass", &jetMatchedGenMass, "jetMatchedGenMass/F");
  JetTree->Branch("jet_n_rechits", &jet_n_rechits, "jet_n_rechits/I");
  JetTree->Branch("jet_rechit_E", &jet_rechit_E, "jet_rechit_E/F");
  JetTree->Branch("jet_rechit_T", &jet_rechit_T, "jet_rechit_T/F");

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

void JetNtupler::enableGenParticleBranches(){
  JetTree->Branch("gLLP_prod_vertex_x", gLLP_prod_vertex_x, "gLLP_prod_vertex_x[2]/F");
  JetTree->Branch("gLLP_prod_vertex_y", gLLP_prod_vertex_y, "gLLP_prod_vertex_y[2]/F");
  JetTree->Branch("gLLP_prod_vertex_z", gLLP_prod_vertex_z, "gLLP_prod_vertex_z[2]/F");
  JetTree->Branch("gLLP_decay_vertex_x", gLLP_decay_vertex_x, "gLLP_decay_vertex_x[2]/F");
  JetTree->Branch("gLLP_decay_vertex_y", gLLP_decay_vertex_y, "gLLP_decay_vertex_y[2]/F");
  JetTree->Branch("gLLP_decay_vertex_z", gLLP_decay_vertex_z, "gLLP_decay_vertex_z[2]/F");
  JetTree->Branch("gLLP_beta", gLLP_beta, "gLLP_beta[2]/F");
  JetTree->Branch("gLLP_decays_px", gLLP_decays_px, "gLLP_decays_px[4]/F");
  JetTree->Branch("gLLP_decays_py", gLLP_decays_py, "gLLP_decays_py[4]/F");
  JetTree->Branch("gLLP_decays_pz", gLLP_decays_pz, "gLLP_decays_pz[4]/F");
  JetTree->Branch("gLLP_decays_e", gLLP_decays_e, "gLLP_decays_e[4]/F");
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
//  iEvent.getByToken(triggerObjectsToken_, triggerObjects);
//  iEvent.getByToken(triggerPrescalesToken_, triggerPrescales);
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
    jetE = 0.0;
    jetPt = 0.0;
    jetEta = 0.0;
    jetPhi = 0.0;
    jetCISV = 0.0;
    jetMass =  -99.0;
    jetJetArea = -99.0;
    jetPileupE = -99.0;
    jetPileupId = -99.0;
    jetPileupIdFlag = -1;
    jetPassIDLoose = false;
    jetPassIDTight = false;
    jetPassMuFrac = false;
    jetPassEleFrac = false;
    jetPartonFlavor = 0;
    jetHadronFlavor = 0;
    jetChargedEMEnergyFraction = -99.0;
    jetNeutralEMEnergyFraction = -99.0;
    jetChargedHadronEnergyFraction = -99.0;
    jetNeutralHadronEnergyFraction = -99.0;
    jetMatchedGenPt = 0.0;
    jetMatchedGenEta = 0.0;
    jetMatchedGenPhi = 0.0;
    jetMatchedGenMass = 0.0;
    jetMatchedGenTime = 0.0;
    jet_n_rechits = 0;
    jet_rechit_E = 0.0;
    jet_rechit_T = 0.0;
    for ( int i = 0; i < 2; i++ )
    {
      gLLP_prod_vertex_x[i] = -666.;
      gLLP_prod_vertex_y[i] = -666.;
      gLLP_prod_vertex_z[i] = -666.;
      gLLP_decay_vertex_x[i] = -666.;
      gLLP_decay_vertex_y[i] = -666.;
      gLLP_decay_vertex_z[i] = -666.;
      gLLP_beta[i] = -666.;
    }

    for ( int i = 0; i < 4; i++ )
    {
      gLLP_decays_px[i] = -666.;
      gLLP_decays_py[i] = -666.;
      gLLP_decays_pz[i] = -666.;
      gLLP_decays_e[i] = -666.;
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

  for (const reco::PFJet &j : *jets) {
    resetBranches();
    if (j.pt() < 10) continue;
    if (fabs(j.eta()) > 1.4) continue;

    //*************************************
    //Fill Event-Level Info
    //*************************************

    //store basic event info
    isData = isData_;
    runNum = iEvent.id().run();
    lumiNum = iEvent.luminosityBlock();
    eventNum = iEvent.id().event();

   //select the primary vertex, if any
    nPV = 0;
    myPV = &(vertices->front());

    bool foundPV = false;
    for(unsigned int i = 0; i < vertices->size(); i++){
      if(vertices->at(i).isValid() && !vertices->at(i).isFake()){
	if (!foundPV) {
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
    if (!isData) {
      for(const PileupSummaryInfo &pu : *puInfo){
	if ( pu.getBunchCrossing() == 0) {
	  nPU = pu.getPU_NumInteractions();
	  nPUmean = pu.getTrueNumInteractions();
	}
      }
    }


    //*************************************
    //Fill Jet-Level Info
    //*************************************
    jetE = j.energy();
    jetPt = j.pt();
    jetEta = j.eta();
    jetPhi = j.phi();
    jetMass = j.mass();

    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiE(jetPt, jetEta, jetPhi, jetE);
    //jetCISV = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

    jetJetArea = j.jetArea();
    jetPileupE = j.pileup();

    jetPileupIdFlag = 0;
    jetPassIDLoose = passJetID(&j, 0);
    jetPassIDTight = passJetID(&j, 1);
    jetPassMuFrac  = ( j.muonEnergyFraction() < 0.80 );
    jetPassEleFrac  = ( j.electronEnergyFraction() < 0.90 );


    // if (useGen_) {
    //   jetPartonFlavor = j.partonFlavour();
    //   jetHadronFlavor = j.hadronFlavour();
    // }

    jetChargedEMEnergyFraction = j.chargedEmEnergyFraction();
    jetNeutralEMEnergyFraction = j.neutralEmEnergyFraction();
    jetChargedHadronEnergyFraction = j.chargedHadronEnergyFraction();
    jetNeutralHadronEnergyFraction = j.neutralHadronEnergyFraction();


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

    int n_matched_rechits = 0;
    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
      {
        if ( recHit->checkFlag(0) )
        {
          const DetId recHitId = recHit->detid();
          const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
          if ( deltaR(jetEta, jetPhi, recHitPos.eta(), recHitPos.phi())  < 0.4)
          {
            n_matched_rechits++;
            jet_rechit_E += recHit->energy();
            jet_rechit_T += recHit->time()*recHit->energy();
          }
          //std::cout << recHitPos.eta() << std::endl;
        }
        //std::cout << recHitId << std::endl;
      }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";
    //std::cout << "n: " << n_matched_rechits << std::endl;
    jet_n_rechits = n_matched_rechits;
    jet_rechit_T = jet_rechit_T/jet_rechit_E;

    //MC AND GEN LEVEL INFO
    fillGenParticles();
    //fillMC();
    JetTree->Fill();
  } //loop over jets

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

    //LLPs vextex
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
          //gLLP_decays_px[0]
        }
        else if (gParticleId[i] == 36)
        {
          gLLP_decay_vertex_x[1] = dau->vx();
          gLLP_decay_vertex_y[1] = dau->vy();
          gLLP_decay_vertex_z[1] = dau->vz();
          gLLP_beta[1] = sqrt(gParticlePx[i]*gParticlePx[i]+gParticlePy[i]*gParticlePy[i]+gParticlePz[i]*gParticlePz[i])/gParticleE[i];
        }
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
  }

  return true;
};

//define this as a plug-in
DEFINE_FWK_MODULE(JetNtupler);
