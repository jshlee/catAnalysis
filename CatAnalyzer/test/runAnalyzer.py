import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("myAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),)
process.source.fileNames = cms.untracked.vstring("file:/pnfs/user/CAT/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/v8-0-3_RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/161219_090504/0000/catTuple_1.root")
    
process.cattree = cms.EDAnalyzer("myAnalyzer",
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag('lumiMask'),
    puweight = cms.InputTag('pileupWeight'),
    puweight_up = cms.InputTag('pileupWeight',"up"),
    puweight_dn = cms.InputTag('pileupWeight',"dn"),
    vertices = cms.InputTag("catVertex"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag('catMETs'),
    mcLabel = cms.InputTag("prunedGenParticles"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree.root")
)

process.p = cms.Path(process.cattree)
#process.MessageLogger.cerr.FwkReport.reportEvery = 50000
