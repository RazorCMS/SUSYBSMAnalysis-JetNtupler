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
  //JetTree = fs->make<TTree>("Jets", "selected miniAOD information");
  JetTree = new TTree("Jets", "selected AOD information");
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

    for (EcalRecHitCollection::const_iterator recHit = ebRecHits->begin(); recHit != ebRecHits->end(); ++recHit)
      {
        if ( recHit->checkFlag(0) )
        {
          const DetId recHitId = recHit->detid();
          const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
          std::cout << recHitPos.eta() << std::endl;
        }
        //std::cout << recHitId << std::endl;
      }
    //cout << "Last Nphoton: " << fJetNPhotons << "\n";

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


//define this as a plug-in
DEFINE_FWK_MODULE(JetNtupler);
