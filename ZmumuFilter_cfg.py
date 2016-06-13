import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("opendata")

goodJSON = 'Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')

import FWCore.Utilities.FileUtils as FileUtils
from FWCore.MessageLogger.MessageLogger_cfi import *

doubleMuFiles = FileUtils.loadListFromFile('CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10001_file_index.txt')
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*doubleMuFiles)
                            )

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
process.source.lumisToProcess.extend(myLumis)

process.load('ZmumuFilter.ZmumuFilter.ZmumuFilter_cfi')

process.ZmumuFilter.csvFileName = cms.string('Zmumu_Run2011A.csv')
process.ZmumuFilter.minMuonPt = cms.double(20.0)
process.ZmumuFilter.maxMuonEta = cms.double(2.1)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500000))

process.mypath = cms.Path(process.ZmumuFilter)
process.schedule = cms.Schedule(process.mypath)
