import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("jetNtupler")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'/store/mc/RunIISummer17DRPremix/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10-v2/00000/FED6864A-CEAB-E711-81CF-008CFAE45308.root',
	#'/store/mc/RunIISummer17DRPremix/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/92X_upgrade2017_realistic_v10-v2/00000/F661A06D-F0AB-E711-AB08-0242AC110004.root'
	'/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/143C8F1C-D3B0-E711-87D6-FA163EA92854.root',
       'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/26B84519-D3B0-E711-B2CD-FA163EC714FC.root',
        'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/FCCC0FEA-D4B0-E711-BCEB-0CC47AD98BF0.root'
	#'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/422E4951-82AC-E711-9642-90B11C2CC76A.root',
	#'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/DCAA63C0-C3AC-E711-89BC-0CC47A13D3A8.root'

#	'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/5C922B1E-ECAB-E711-BB3C-FA163E1A8CC6.root',
#	'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/7CDC5245-ECAC-E711-8846-02163E013DE2.root',
 #       'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/8250B221-13AD-E711-A78D-FA163EAD94F3.root',
#	 'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/D40EFDBC-12AD-E711-8345-0CC47AA53D5A.root',
#	 'file:/mnt/hadoop/store/mc/RunIISummer17DRStdmix/XXTo4J_M-500_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/A6489F83-BCAC-E711-9F03-FA163EE923D3.root'


 #'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000_DR-AODSIM_CaltechT2/190122_173529/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_1.root',
  #      'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000_DR-AODSIM_CaltechT2/190122_173529/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_2.root',
   #     'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000_DR-AODSIM_CaltechT2/190121_173529/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_3.root',
    #    'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m50_pl100_ev10000_DR-AODSIM_CaltechT2/190121_173529/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_4.root'
#       'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000_DR-AODSIM_CaltechT2/190121_105636/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_1.root',
 #     'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000_DR-AODSIM_CaltechT2/190121_105636/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_2.root',
  #    'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000_DR-AODSIM_CaltechT2/190121_105636/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_3.root',
    #  'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m10_pl1000_ev10000_DR-AODSIM_CaltechT2/190121_105636/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_4.root'
        #'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m30_pl1000_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m30_pl1000_ev10000_DR-AODSIM_CaltechT2/190120_072245/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_1.root',
        
	#'file:/mnt/hadoop/store/user/christiw/ppTohToSS1SS2_SS1Tobb_SS2Toveve_MC_prod/ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m30_pl1000_ev10000/crab_CMSSW_8_0_31_ppTohToSS1SS2_SS1Tobb_SS2Toveve_run_m30_pl1000_ev10000_DR-AODSIM_CaltechT2/190120_072245/0000/ppTohToSS1SS2_SS1Tobb_SS2Toveve_m_10_pl_10_ev_10000_step2_2.root' 
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(201) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string("jetNtuple_M-500_CTau-1000mm.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#


#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'
process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'

#------ If we add any inputs beyond standard miniAOD event content, import them here ------#

#process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
#process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
#process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#process.BadPFMuonFilter.taggingMode = cms.bool(True)

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('JetNtupler',
    isData = cms.bool(False),
    useGen = cms.bool(True),
    isFastsim = cms.bool(False),
    enableTriggerInfo = cms.bool(True),
    readGenVertexTime = cms.bool(True),
    isFourJet = cms.bool(True), #false means glueball model, true means four-jet model 
    genParticles_t0 = cms.InputTag("genParticles", "t0", ""),
    triggerPathNamesFile = cms.string("SUSYBSMAnalysis/JetNtupler/data/TriggerNames_LLP_V1.dat"),
    eleHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorElectronHLTFilterNames.dat"),
    muonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorMuonHLTFilterNames.dat"),
    photonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorPhotonHLTFilterNames.dat"),

    vertices = cms.InputTag("offlinePrimaryVerticesWithBS"),  # for non-timing case
    muons = cms.InputTag("muons"),
    electrons = cms.InputTag("gedGsfElectrons"),
    taus = cms.InputTag("hpsPFTauProducer"),
    photons = cms.InputTag("gedPhotons"),
    jets = cms.InputTag("ak4PFJetsCHS"),
    jetsPuppi = cms.InputTag("ak4PFJets"),
    jetsAK8 = cms.InputTag("ak8PFJetsCHS"),
    mets = cms.InputTag("pfMet"),
    #metsNoHF = cms.InputTag("pfMet30"),
    metsPuppi = cms.InputTag("pfMet"),
    pfCands = cms.InputTag("particleFlow","","RECO"),

    #packedPfCands = cms.InputTag("packedPFCandidates"),

    genParticles = cms.InputTag("genParticles"),

    #packedGenParticles = cms.InputTag("packedGenParticles"),
    #prunedGenParticles = cms.InputTag("prunedGenParticles"),
    genJets = cms.InputTag("ak4GenJets"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    #triggerBits = cms.InputTag("TriggerResults","","RECO"),
    hepMC = cms.InputTag("generatorSmeared", "", "SIM"),

    #triggerPrescales = cms.InputTag("patTrigger"),
    #triggerObjects = cms.InputTag("selectedPatTrigger"),

    metFilterBits = cms.InputTag("TriggerResults", "", "RECO"),

    #hbheNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
    #hbheTightNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    #hbheIsoNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),

    #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
    #BadMuonFilter = cms.InputTag("BadPFMuonFilter",""),

    #lheInfo = cms.InputTag("externalLHEProducer", "", ""),
    genInfo = cms.InputTag("generator", "", "SIM"),

    tracks = cms.InputTag("generalTracks", "", "RECO"),
    #trackTime = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    #trackTimeReso = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),

    puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    #hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    secondaryVertices = cms.InputTag("inclusiveSecondaryVertices", "", "RECO"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),

    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    pfClusters = cms.InputTag("particleFlowClusterECAL","","RECO"),
    ebRecHits = cms.InputTag("reducedEcalRecHitsEB", "","RECO"),
    #ebRecHits = cms.InputTag("EcalRecHit", "reducedEcalRecHitsEB", "RECO"),
    eeRecHits  = cms.InputTag("reducedEcalRecHitsEE", "","RECO"),
    esRecHits = cms.InputTag("reducedEcalRecHitsES", "","RECO"),
    #ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "RECO"),
    ebeeClusters = cms.InputTag("particleFlowEGamma", "EBEEClusters", "RECO"),
    esClusters = cms.InputTag("particleFlowEGamma", "ESClusters", "RECO"),
    #conversions = cms.InputTag("reducedEgamma", "reducedConversions", "RECO"),
    conversions = cms.InputTag("allConversions", "", "RECO"),

    #singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "RECO"),
    singleLegConversions = cms.InputTag("particleFlowEGamma", "", "RECO"),

    gedGsfElectronCores = cms.InputTag("gedGsfElectronCores", "", "RECO"),
    gedPhotonCores = cms.InputTag("gedPhotonCore", "", "RECO"),
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    #lostTracks = cms.InputTag("lostTracks", "", "RECO")
)

#run
process.p = cms.Path( #process.HBHENoiseFilterResultProducer*
                      #process.BadChargedCandidateFilter*
                      #process.BadPFMuonFilter*
                      process.ntuples)
