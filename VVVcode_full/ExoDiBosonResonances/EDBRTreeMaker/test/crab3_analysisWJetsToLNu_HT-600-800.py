from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'WJetsToLNu_HT-600To800'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt','Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt','Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt','Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt','Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt','Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt','Summer16_23Sep2016V3_MC_L1FastJet_AK8PFPuppi.txt','Summer16_23Sep2016V3_MC_L2Relative_AK8PFPuppi.txt','Summer16_23Sep2016V3_MC_L3Absolute_AK8PFPuppi.txt','Summer16_23Sep2016V3_MC_L1FastJet_AK4PFPuppi.txt','Summer16_23Sep2016V3_MC_L2Relative_AK4PFPuppi.txt','Summer16_23Sep2016V3_MC_L3Absolute_AK4PFPuppi.txt']
#config.JobType.inputFiles = ['PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt','PHYS14_25_V2_All_L2Relative_AK4PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK4PFchs.txt','PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt','PHYS14_25_V2_All_L2Relative_AK8PFchs.txt','PHYS14_25_V2_All_L3Absolute_AK8PFchs.txt']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
#config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDataset = '/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =2
config.Data.totalUnits = -1
config.Data.publication = False
name = 'ZZZ'
steam_dir = 'xulyu'
config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/' + steam_dir + '/' + name + '/'
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'WJetsToLNu_HT-600To800'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
