#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom3.h>
#include <iostream>
#include <fstream>

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BarePhotons.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"

#include "MitAnalysisRunII/macros/80x/factors.h"
#include "MitAnalysisRunII/macros/80x/BTagCalibrationStandalone.cc"
#include "MitZHAnalysis/macros/80x/zhMVA.h"

bool       isMINIAOD               = true;
int        whichSkim               = 4;
bool       usePureMC               = true; 
bool       useEMFromData           = false;
const bool useDYPT                 = true;
double     mcPrescale              = 1.;
bool       verbose                 = true;
enum selType                     { ZHGSEL,   BTAGSEL,	WWSEL,   PRESEL,  ZLLSEL,   ZLGSEL,   WZSEL,   ZZSEL,   ZZLOOSESEL, nSelTypes};
TString selTypeName[nSelTypes]=  {"ZHGSEL", "BTAGSEL", "WWSEL", "PRESEL","ZLLSEL", "ZLGSEL", "WZSEL", "ZZSEL", "ZZLOOSESEL"};
enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};
enum categoryType                 { DATA=0, EM, DY, WZ, ZZ, VVV, ZH, ggZH};
const TString typeLepSel = "medium";
const double bTagCuts[1] = {0.8484}; // 0.5426/0.8484/0.9535 (check BTagCalibration2Reader!)
const double sf_el_gamma[2] = {1.13, 1.25};

void zhgAnalysis(
 unsigned int nJetsType = 1,
 Int_t typeSel = 3,
 Int_t plotModel = 0,
 bool isMIT = true,
 double ptGMIN = 25,
 double metMIN = 100
 ){

  TString subFolder = "";
  if(ptGMIN != 25 || metMIN != 100){
    subFolder = Form("pt%d_met%d",(int)ptGMIN,(int)metMIN);
  }

  system(Form("mkdir -p MitZHAnalysis/datacards_zhg%s",subFolder.Data()));
  system(Form("mkdir -p MitZHAnalysis/plots_zhg%s",subFolder.Data()));

  Int_t period = 1;
  // File instances on EOS
  TString filesPathDA   = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/met_";
  TString filesPathMC   = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/met_";
  TString filesPathMC2  = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/met_";
  TString filesPathDMMC = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/";
  if(whichSkim == 0){
  TString filesPathDA   = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/";
  TString filesPathMC   = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/";
  TString filesPathMC2  = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/";
  TString filesPathDMMC = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/";
  }
  // File instances on T3 hadoop
  if(isMIT){
    filesPathDA   = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/data/met_";
    filesPathMC	  = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/met_";
    filesPathMC2  = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/met_";
    filesPathDMMC = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/dm/";
  }
  Double_t lumi = 35.9;
  TString processTag = "";

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev, signalName_;
  vector<Int_t> infilecatv, signalIndex_;  

  TString puPath = "";
  TString triggerSuffix = "*";
  if(isMINIAOD) triggerSuffix = "";

  puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";

  // Data files
  if(isMINIAOD) {
    infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016D.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016E.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016F.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016G.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016H.root",filesPathDA.Data())); infilecatv.push_back(0);
  } else {
  }

  // Monte carlo backgrounds
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV.root",filesPathMC.Data()));					      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg.root",filesPathMC2.Data()));					      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));infilecatv.push_back(1);

  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root",filesPathMC.Data()));	              infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8.root",filesPathMC.Data()));              infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                            infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data())); 	      infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.root",filesPathMC.Data()));   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));           infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infilecatv.push_back(1);

  if(useDYPT==false){
  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));      infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	    infilecatv.push_back(2);
  }
  else {
  infilenamev.push_back(Form("%sDYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));  infilecatv.push_back(2);
  }
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(2);

  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8.root",filesPathMC.Data()));  				      infilecatv.push_back(4);

  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8.root",filesPathMC.Data()));					      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	 	      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(5);
  infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                                  infilecatv.push_back(5);

  for(int ifile=0; ifile<(int)infilenamev.size(); ifile++) {
    signalIndex_.push_back(-1); // Populate vector of signal indices with -1 for the non-MC-signal files
  }

  // Monte Carlo signals
  { // Model 0: standard model Higgs (125) with glu-glu
    int mH=125;
    signalName_.push_back("sm_powheg");
    infilenamev.push_back(Form("%sZH_ZToLL_HToInvG_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH)); infilecatv.push_back(6); signalIndex_.push_back(0);
    signalName_.push_back("sm_amcnlo");
    infilenamev.push_back(Form("%sZH_ZToLL_HToInvG_M%d_13TeV_amcnlo_pythia8.root",filesPathDMMC.Data(),mH)); infilecatv.push_back(6); signalIndex_.push_back(1);
  }  // Models 1 thru 8: standard-model-like Higgs mass points without glu-glu (8 models)

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}
  
  //signalName_.clear();infilenamev.clear();infilecatv.clear();signalIndex_.clear();int i = 0;
  //signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
  //infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathMC.Data(),600)); infilecatv.push_back(5);

  int nSigModels=signalName_.size();
  
  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff2016_m");sprintf(momMName,"CMS_scale2016_m");
  sprintf(effEName,"CMS_eff2016_e");sprintf(momEName,"CMS_scale2016_e");
  if     (typeSel == 3 && nJetsType == 0) {sprintf(finalStateName,"ll0j");}
  else if(typeSel == 1 && nJetsType == 0) {sprintf(finalStateName,"mm0j");}
  else if(typeSel == 2 && nJetsType == 0) {sprintf(finalStateName,"ee0j");}
  else if(typeSel == 3 && nJetsType == 1) {sprintf(finalStateName,"ll1j");}
  else if(typeSel == 1 && nJetsType == 1) {sprintf(finalStateName,"mm1j");}
  else if(typeSel == 2 && nJetsType == 1) {sprintf(finalStateName,"ee1j");}
  else if(typeSel == 3 && nJetsType == 2) {sprintf(finalStateName,"ll2j");}
  else if(typeSel == 1 && nJetsType == 2) {sprintf(finalStateName,"mm2j");}
  else if(typeSel == 2 && nJetsType == 2) {sprintf(finalStateName,"ee2j");}
  else if(typeSel == 0 && nJetsType == 0) {sprintf(finalStateName,"em0j");}
  else if(typeSel == 0 && nJetsType == 1) {sprintf(finalStateName,"em1j");}
  else if(typeSel == 0 && nJetsType == 2) {sprintf(finalStateName,"em2j");}
  else {printf("Wrong lSel/nJetsType: %d/%d\n",typeSel,nJetsType); assert(0); return;}

  double denBTagging[5][5][3],jetEpsBtagMEDIUM[5][5][3];
  double numBTaggingMEDIUM[5][5][3];
  for(int i0=0; i0<5; i0++) {
    for(int i1=0; i1<5; i1++) {
      for(int i2=0; i2<3; i2++) {
        denBTagging[i0][i1][i2] = 0.0;
        numBTaggingMEDIUM[i0][i1][i2] = 0.0;
	if     (i2==BTagEntry::FLAV_B)    jetEpsBtagMEDIUM[i0][i1][i2] = jetEpsBtagBMEDIUM[i0][i1];
	else if(i2==BTagEntry::FLAV_C)    jetEpsBtagMEDIUM[i0][i1][i2] = jetEpsBtagCMEDIUM[i0][i1];
	else if(i2==BTagEntry::FLAV_UDSG) jetEpsBtagMEDIUM[i0][i1][i2] = jetEpsBtagLMEDIUM[i0][i1];
      }
    }
  }

  //Float_t fMVACut[4][4];
  //InitializeJetIdCuts(fMVACut);
  
  BTagCalibration2 *btagCalib = new BTagCalibration2("csvv2","MitAnalysisRunII/data/80x/CSVv2_Moriond17_B_H.csv");
  BTagCalibration2Reader btagReaderBCMEDIUM(btagCalib,BTagEntry::OP_MEDIUM,"comb","central");
  BTagCalibration2Reader btagReaderLMEDIUM(btagCalib,BTagEntry::OP_MEDIUM,"incl","central");
  BTagCalibration2Reader btagReaderBCMEDIUMUP(btagCalib,BTagEntry::OP_MEDIUM,"comb","up");
  BTagCalibration2Reader btagReaderLMEDIUMUP(btagCalib,BTagEntry::OP_MEDIUM,"incl","up");
  BTagCalibration2Reader btagReaderBCMEDIUMDOWN(btagCalib,BTagEntry::OP_MEDIUM,"comb","down");
  BTagCalibration2Reader btagReaderLMEDIUMDOWN(btagCalib,BTagEntry::OP_MEDIUM,"incl","down");

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;

  TFile *fTrackElectronReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDeltrksf= (TH2D*)(fTrackElectronReco_SF->Get("scalefactors_Reco_Electron")); assert(fhDeltrksf); fhDeltrksf->SetDirectory(0);
  delete fTrackElectronReco_SF;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_Medium_Electron"));
  TH2D *fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF->SetDirectory(0);
  delete fElSF;

  TFile *fElVeryTightSF = TFile::Open(Form("MitAnalysisRunII/data/80x/veryTightSF_37ifb.root"));
  TH1D *fhDVeryTightSF = (TH1D*)(fElVeryTightSF->Get("veryTightSF"));
  assert(fhDVeryTightSF);
  fhDVeryTightSF->SetDirectory(0);
  delete fElVeryTightSF;

  TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/Tracking_EfficienciesAndSF_BCDEFGH.root"));
  TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("ratio_eff_eta3_dr030e030_corr")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
  delete fTrackMuonReco_SF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("scalefactors_TightId_Muon")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  //TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonID_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  delete fMuSF;

  TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("scalefactors_Iso_MuonTightId")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  //TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonIso_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  delete fMuIsoSF;

  TFile *fZHEwkCorr = TFile::Open(Form("MitAnalysisRunII/data/80x/Zll_nloEWK_weight_unnormalized.root"));
  TH1D *fhDZHEwkCorr     = (TH1D*)(fZHEwkCorr->Get("SignalWeight_nloEWK_rebin"));      assert(fhDZHEwkCorr);     fhDZHEwkCorr    ->SetDirectory(0);
  TH1D *fhDZHEwkCorrUp   = (TH1D*)(fZHEwkCorr->Get("SignalWeight_nloEWK_up_rebin"));   assert(fhDZHEwkCorrUp);   fhDZHEwkCorrUp  ->SetDirectory(0);
  TH1D *fhDZHEwkCorrDown = (TH1D*)(fZHEwkCorr->Get("SignalWeight_nloEWK_down_rebin")); assert(fhDZHEwkCorrDown); fhDZHEwkCorrDown->SetDirectory(0);
  delete fZHEwkCorr;

  TFile *fPhotonSF = TFile::Open(Form("MitAnalysisRunII/data/80x/photon_scalefactors_37ifb.root"));
  TH2D *fhDPhotonSF       = (TH2D*)(fPhotonSF->Get("EGamma_SF2D")); assert(fhDPhotonSF); fhDPhotonSF->SetDirectory(0);
  TH2D *fhDElectronVetoSF = (TH2D*)(fPhotonSF->Get("Scaling_Factors_HasPix_R9 Inclusive")); assert(fhDElectronVetoSF); fhDElectronVetoSF->SetDirectory(0);
  delete fPhotonSF;

  TString ECMsb  = "13TeV2016";
  
  // MVA variable types:
  // 0: MT(g-MET)
  // 1: MT(g-MET) vs. eta

  const int MVAVarType = 0; const int nBinMVA = 5; Double_t xbins[nBinMVA+1] = {0, 75, 100, 125, 150, 175}; TString addChan = "";
  //const int MVAVarType = 1; const int nBinMVA = 6; Double_t xbins[nBinMVA+1] = {0, 75, 120, 175, 275, 320, 375}; TString addChan = "";

  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  TFile *fVJetsKfactorFile = TFile::Open(Form("MitAnalysisRunII/data/80x/kfactors_vjets.root"));
  TH1D *fhDVjetsNum = (TH1D*)(fVJetsKfactorFile->Get("EWKcorr/Z"));
  TH1D *fhDVjetsDen = (TH1D*)(fVJetsKfactorFile->Get("ZJets_LO/inv_pt"));
  assert(fhDVjetsNum);
  assert(fhDVjetsDen);
  fhDVjetsNum->SetDirectory(0);
  fhDVjetsDen->SetDirectory(0);
  delete fVJetsKfactorFile;

  const int numberCuts = 12;
  TH1D* histoZHSEL[4];
  histoZHSEL[0] = new TH1D("histoZHSEL_0", "histoZHSEL_0", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[1] = new TH1D("histoZHSEL_1", "histoZHSEL_1", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[2] = new TH1D("histoZHSEL_2", "histoZHSEL_2", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[3] = new TH1D("histoZHSEL_3", "histoZHSEL_3", numberCuts+1, -0.5, numberCuts+0.5);
  TString cutName[numberCuts+1] = {"ptl>25/20","3rd lepton veto","ptll>60","Z mass","npho>=1","btag-veto","tauVeto","Njets","MET>100","dPhi(ZG-MET)>2.5","|ptllg-MET|/ptll<0.4","dPhiJetMet>0.5","MTGMET<200"};

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 40;
  const int histBins = 8;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "...DY", "...WZ", "....ZZ", "...VVV", "....ZH", "..ggZH"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot ==  0) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =1000.0;}
    else if(thePlot ==  1) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot ==  2) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot ==  3) {nBinPlot = 200; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot ==  4) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot ==  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot ==  6) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot ==  7) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot ==  8) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot ==  9) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot == 10) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot =  2.5;}
    else if(thePlot == 11) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot == 12) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 13) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 14) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 15) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot == 16) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot == 17) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot == 18) {nBinPlot = 500; xminPlot =-0.5; xmaxPlot = 499.5;}
    else if(thePlot == 19) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot == 20) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   2.5;}
    else if(thePlot == 21) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot == 22) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot == 23) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot == 24) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot == 25) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot == 26) {nBinPlot =  60; xminPlot = 0.0; xmaxPlot =   3.0;}
    else if(thePlot == 27) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot == 28) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 29) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 30) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 250.0;}
    else if(thePlot == 31) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot == 32) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot == 33) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot == 34) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot = 3.5;}
    else if(thePlot == allPlots-2)          {nBinPlot =  numberCuts+1; xminPlot =-0.5; xmaxPlot =  numberCuts+0.5;}
    TH1D* histos;
    if(thePlot != allPlots-1) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    else                      histos = new TH1D("histos", "histos", nBinMVA, xbins);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  TH1D *histo_Data     = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_Zjets    = (TH1D*) histoMVA->Clone("histo_Zjets");	 
  TH1D *histo_VVV      = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZ       = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZ       = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_EM       = (TH1D*) histoMVA->Clone("histo_EM");	 
  TH1D *histo_ggZH_hinv= (TH1D*) histoMVA->Clone("histo_ggZH_hinv"); 
  TH1D *histo_ZjetsNoW    = (TH1D*) histoMVA->Clone("histo_Zjets");	 
  TH1D *histo_VVVNoW      = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZNoW       = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZNoW       = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_EMNoW       = (TH1D*) histoMVA->Clone("histo_EM");	 
  TH1D *histo_ggZH_hinvNoW= (TH1D*) histoMVA->Clone("histo_ggZH_hinv"); 
  TH1D *histo_ZH_hinv[nSigModels];
  TH1D *histo_ZH_hinvNoW[nSigModels];
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histo_ZH_hinv[nModel]     = (TH1D*) histoMVA->Clone(Form("histo_ZH_hinv_%s",   signalName_[nModel].Data())); 
    histo_ZH_hinvNoW[nModel]  = (TH1D*) histoMVA->Clone(Form("histo_ZH_hinvNoW_%s",signalName_[nModel].Data())); 
  }

  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nSigModels];
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]      = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sUp"  ,finalStateName, signalName_[nModel].Data(), ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sUp"  ,finalStateName, signalName_[nModel].Data(), ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]    = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sDown",finalStateName, signalName_[nModel].Data(), ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sDown",finalStateName, signalName_[nModel].Data(), ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]->Sumw2();
  }
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingUp     = new TH1D( Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingDown   = new TH1D( Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp           = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown         = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_EM_CMS_MVAEMStatBoundingUp           = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingUp  ->Sumw2();
  TH1D* histo_EM_CMS_MVAEMStatBoundingDown         = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp  = new TH1D( Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown= new TH1D( Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown->Sumw2();

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  TH1D* histo_ZH_hinv_CMS_QCDScaleBounding[nSigModels][6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ggZH_hinv_CMS_QCDScaleBounding[6];
  for(int nb=0; nb<6; nb++){
    for(int nModel=0; nModel<nSigModels; nModel++) {
      histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nb]   = new TH1D(Form("histo_ZH_hinv_%s_QCDScale_f%d", signalName_[nModel].Data(), nb), Form("histo_ZH_hinv_%s_QCDScale_f%d", signalName_[nModel].Data(), nb),nBinMVA, xbins); histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nb]->Sumw2();
    }
    histo_VVV_CMS_QCDScaleBounding[nb]       = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding[nb]	     = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),      Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_WZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding[nb]	     = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ggZH_hinv_CMS_QCDScaleBounding[nb] = new TH1D(Form("histo_ggZH_hinv_QCDScale_f%d",nb), Form("histo_ggZH_hinv_QCDScale_f%d",nb),nBinMVA, xbins); histo_ggZH_hinv_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_ZH_hinv_CMS_PDFBounding[nSigModels][100];
  TH1D* histo_VVV_CMS_PDFBounding[100];
  TH1D* histo_WZ_CMS_PDFBounding[100];
  TH1D* histo_ZZ_CMS_PDFBounding[100];
  TH1D* histo_ggZH_hinv_CMS_PDFBounding[100];
  for(int nb=0; nb<100; nb++) {
    for(int nModel=0; nModel<nSigModels; nModel++) {
      histo_ZH_hinv_CMS_PDFBounding[nModel][nb]   = new TH1D(Form("histo_ZH_hinv_%s_PDF_f%d", signalName_[nModel].Data(), nb), Form("histo_ZH_hinv_%s_PDF_f%d", signalName_[nModel].Data(), nb),nBinMVA, xbins); histo_ZH_hinv_CMS_PDFBounding[nModel][nb]->Sumw2();
    }
    histo_VVV_CMS_PDFBounding[nb]       = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WZ_CMS_PDFBounding[nb]        = new TH1D(Form("histo_WZ_PDF_f%d",nb),      Form("histo_WZ_PDF_f%d",nb),     nBinMVA, xbins); histo_WZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ZZ_CMS_PDFBounding[nb]        = new TH1D(Form("histo_ZZ_PDF_f%d",nb),      Form("histo_ZZ_PDF_f%d",nb),     nBinMVA, xbins); histo_ZZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ggZH_hinv_CMS_PDFBounding[nb] = new TH1D(Form("histo_ggZH_hinv_PDF_f%d",nb), Form("histo_ggZH_hinv_PDF_f%d",nb),nBinMVA, xbins); histo_ggZH_hinv_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nSigModels][nBinMVA];
  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nSigModels][nBinMVA];
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];
  TH1D* histo_EM_CMS_MVAEMStatBoundingBinUp[nBinMVA];
  TH1D* histo_EM_CMS_MVAEMStatBoundingBinDown[nBinMVA];
  TH1D* histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[nBinMVA];
  TH1D* histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++) {
    for(int nModel=0; nModel<nSigModels; nModel++) { 
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][nb]        = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dUp"       ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dUp"       ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb),nBinMVA, xbins); 
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][nb]      = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dDown"     ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dDown"     ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb),nBinMVA, xbins); 
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][nb]   ->Sumw2();
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][nb] ->Sumw2();
    }
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]       = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"      ,finalStateName,  ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"      ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]     = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown"    ,finalStateName,  ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown"    ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]            = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]            = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"        ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"        ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]                = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]            = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]                = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]            = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]                = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]            = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[nb]    = new TH1D(Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dUp"   ,finalStateName,  ECMsb.Data(),nb), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dUp"   ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[nb]  = new TH1D(Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dDown" ,finalStateName,  ECMsb.Data(),nb), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dDown" ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]  ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	      ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]	  ->Sumw2();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[nb] ->Sumw2();
  }

  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	   = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	   = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingUp    	   = new TH1D( Form("histo_WZ_%sUp",effMName)  , Form("histo_WZ_%sUp",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingDown  	   = new TH1D( Form("histo_WZ_%sDown",effMName), Form("histo_WZ_%sDown",effMName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp    	   = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown  	   = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingUp   = new TH1D( Form("histo_ggZH_hinv_%sUp",effMName)  , Form("histo_ggZH_hinv_%sUp",effMName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingDown = new TH1D( Form("histo_ggZH_hinv_%sDown",effMName), Form("histo_ggZH_hinv_%sDown",effMName), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nSigModels]; 
  TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	   = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg    	   = new TH1D( Form("histo_WZ_%sAvg",effMName)  , Form("histo_WZ_%sAvg",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    	   = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg  = new TH1D( Form("histo_ggZH_hinv_%sAvg",effMName)  , Form("histo_ggZH_hinv_%sAvg",effMName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nSigModels];

  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	   = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	   = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingUp    	   = new TH1D( Form("histo_WZ_%sUp",effEName)  , Form("histo_WZ_%sUp",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingDown  	   = new TH1D( Form("histo_WZ_%sDown",effEName), Form("histo_WZ_%sDown",effEName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp    	   = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown  	   = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingUp   = new TH1D( Form("histo_ggZH_hinv_%sUp",effEName)  , Form("histo_ggZH_hinv_%sUp",effEName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingDown = new TH1D( Form("histo_ggZH_hinv_%sDown",effEName), Form("histo_ggZH_hinv_%sDown",effEName), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	   = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg    	   = new TH1D( Form("histo_WZ_%sAvg",effEName)  , Form("histo_WZ_%sAvg",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    	   = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg  = new TH1D( Form("histo_ggZH_hinv_%sAvg",effEName)  , Form("histo_ggZH_hinv_%sAvg",effEName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nSigModels];

  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown"), nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAMETBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_scale_metUp")  , Form("histo_ggZH_hinv_CMS_scale_metUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAMETBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_scale_metDown"), Form("histo_ggZH_hinv_CMS_scale_metDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAMETBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_MVAMETBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_MVAJESBoundingUp      	= new TH1D( Form("histo_VVV_CMS_eff_b_2016Up")  , Form("histo_VVV_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown    	= new TH1D( Form("histo_VVV_CMS_eff_b_2016Down"), Form("histo_VVV_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_WZ_CMS_eff_b_2016Up")  , Form("histo_WZ_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_WZ_CMS_eff_b_2016Down"), Form("histo_WZ_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_ZZ_CMS_eff_b_2016Up")  , Form("histo_ZZ_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_ZZ_CMS_eff_b_2016Down"), Form("histo_ZZ_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_eff_b_2016Up")  , Form("histo_ggZH_hinv_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_eff_b_2016Down"), Form("histo_ggZH_hinv_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_MVABTAGBoundingUp      	= new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVABTAGBoundingDown    	= new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVABTAGBoundingUp       	= new TH1D( Form("histo_WZ_CMS_scale_jUp")  , Form("histo_WZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVABTAGBoundingDown     	= new TH1D( Form("histo_WZ_CMS_scale_jDown"), Form("histo_WZ_CMS_scale_jDown"), nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVABTAGBoundingUp       	= new TH1D( Form("histo_ZZ_CMS_scale_jUp")  , Form("histo_ZZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVABTAGBoundingDown     	= new TH1D( Form("histo_ZZ_CMS_scale_jDown"), Form("histo_ZZ_CMS_scale_jDown"), nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVABTAGBoundingUp   = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jUp")  , Form("histo_ggZH_hinv_CMS_scale_jUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVABTAGBoundingDown = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jDown"), Form("histo_ggZH_hinv_CMS_scale_jDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVABTAGBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_MVABTAGBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_PUBoundingUp           	= new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown         	= new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingUp            	= new TH1D( Form("histo_WZ_CMS_puUp")  , Form("histo_WZ_CMS_puUp")  , nBinMVA, xbins); histo_WZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingDown          	= new TH1D( Form("histo_WZ_CMS_puDown"), Form("histo_WZ_CMS_puDown"), nBinMVA, xbins); histo_WZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingUp            	= new TH1D( Form("histo_ZZ_CMS_puUp")  , Form("histo_ZZ_CMS_puUp")  , nBinMVA, xbins); histo_ZZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingDown  	        = new TH1D( Form("histo_ZZ_CMS_puDown"), Form("histo_ZZ_CMS_puDown"), nBinMVA, xbins); histo_ZZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_PUBoundingUp        = new TH1D( Form("histo_ggZH_hinv_CMS_puUp")  , Form("histo_ggZH_hinv_CMS_puUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_PUBoundingDown      = new TH1D( Form("histo_ggZH_hinv_CMS_puDown"), Form("histo_ggZH_hinv_CMS_puDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_PUBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_PUBoundingDown[nSigModels];

  TH1D* histo_ZH_hinv_CMS_EWKCorrUp[nSigModels];
  TH1D* histo_ZH_hinv_CMS_EWKCorrDown[nSigModels];
  TH1D* histo_WZ_CMS_EWKCorrUp    	        = new TH1D( Form("histo_WZ_EWKCorrUp")  , Form("histo_WZ_EWKCorrUp")  , nBinMVA, xbins); histo_WZ_CMS_EWKCorrUp  ->Sumw2();
  TH1D* histo_WZ_CMS_EWKCorrDown                = new TH1D( Form("histo_WZ_EWKCorrDown"), Form("histo_WZ_EWKCorrDown"), nBinMVA, xbins); histo_WZ_CMS_EWKCorrDown->Sumw2();
  TH1D* histo_ZZ_CMS_EWKCorrUp                  = new TH1D( Form("histo_ZZ_EWKCorrUp")  , Form("histo_ZZ_EWKCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_EWKCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_EWKCorrDown                = new TH1D( Form("histo_ZZ_EWKCorrDown"), Form("histo_ZZ_EWKCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_EWKCorrDown->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrUp                   = new TH1D( Form("histo_ZZ_ggCorrUp")  , Form("histo_ZZ_ggCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_ggCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrDown                 = new TH1D( Form("histo_ZZ_ggCorrDown"), Form("histo_ZZ_ggCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_ggCorrDown->Sumw2();
  TH1D* histo_Zjets_CMS_ZjetsSystUp    	        = new TH1D( Form("histo_Zjets_ZjetsSystUp")  , Form("histo_Zjets_ZjetsSystUp")  , nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_ZjetsSystDown           = new TH1D( Form("histo_Zjets_ZjetsSystDown"), Form("histo_Zjets_ZjetsSystDown"), nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystDown->Sumw2();

  for(int nModel=0; nModel<nSigModels; nModel++) { 
    histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]          = new TH1D( Form("histo_ZH_hinv_%s_%sUp",   signalName_[nModel].Data(), effMName), Form("histo_ZH_hinv_%s_%sUp",  signalName_[nModel].Data(), effMName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]        = new TH1D( Form("histo_ZH_hinv_%s_%sDown", signalName_[nModel].Data(), effMName), Form("histo_ZH_hinv_%s_%sDown",signalName_[nModel].Data(), effMName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffMBoundingAvg [nModel]        = new TH1D( Form("histo_ZH_hinv_%s_%sAvg",  	   signalName_[nModel].Data(), effMName), Form("histo_ZH_hinv_%s_%sAvg" ,	   signalName_[nModel].Data(), effMName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingUp [nModel]         = new TH1D( Form("histo_ZH_hinv_%s_%sUp",		   signalName_[nModel].Data(), effEName), Form("histo_ZH_hinv_%s_%sUp"  ,	   signalName_[nModel].Data(), effEName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingDown [nModel]       = new TH1D( Form("histo_ZH_hinv_%s_%sDown", 	   signalName_[nModel].Data(), effEName), Form("histo_ZH_hinv_%s_%sDown",	   signalName_[nModel].Data(), effEName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingAvg [nModel]        = new TH1D( Form("histo_ZH_hinv_%s_%sAvg",  	   signalName_[nModel].Data(), effEName), Form("histo_ZH_hinv_%s_%sAvg" ,	   signalName_[nModel].Data(), effEName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAMETBoundingUp [nModel]             = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_metUp"  , signalName_[nModel].Data()), 	  Form("histo_ZH_hinv_%s_CMS_scale_metUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAMETBoundingDown [nModel]           = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_metDown", signalName_[nModel].Data()), 	  Form("histo_ZH_hinv_%s_CMS_scale_metDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVAJESBoundingUp [nModel]             = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_jUp"	 , signalName_[nModel].Data()), 	  Form("histo_ZH_hinv_%s_CMS_scale_jUp"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAJESBoundingDown [nModel]           = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_jDown"  , signalName_[nModel].Data()), 	  Form("histo_ZH_hinv_%s_CMS_scale_jDown"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVABTAGBoundingUp [nModel]            = new TH1D( Form("histo_ZH_hinv_%s_CMS_eff_b_2016Up"	, signalName_[nModel].Data()),  Form("histo_ZH_hinv_%s_CMS_eff_b_2016Up"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVABTAGBoundingDown [nModel]          = new TH1D( Form("histo_ZH_hinv_%s_CMS_eff_b_2016Down"  , signalName_[nModel].Data()), 	 Form("histo_ZH_hinv_%s_CMS_eff_b_2016Down"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_PUBoundingUp [nModel]                 = new TH1D( Form("histo_ZH_hinv_%s_CMS_puUp"	 , signalName_[nModel].Data()), 	  Form("histo_ZH_hinv_%s_CMS_puUp"	   , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_PUBoundingDown [nModel]               = new TH1D( Form("histo_ZH_hinv_%s_CMS_puDown"	 , signalName_[nModel].Data()), 	  Form("histo_ZH_hinv_%s_CMS_puDown"	   , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_EWKCorrUp[nModel]                     = new TH1D( Form("histo_ZH_hinv_%s_%sUp",   signalName_[nModel].Data(), "CMS_EWKCorr"), Form("histo_ZH_hinv_%s_%sUp",  signalName_[nModel].Data(), "CMS_EWKCorr"), nBinMVA, xbins); histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_EWKCorrDown[nModel]                   = new TH1D( Form("histo_ZH_hinv_%s_%sDown", signalName_[nModel].Data(), "CMS_EWKCorr"), Form("histo_ZH_hinv_%s_%sDown",signalName_[nModel].Data(), "CMS_EWKCorr"), nBinMVA, xbins); histo_ZH_hinv_CMS_EWKCorrDown[nModel]->Sumw2();

  }

  double bgdDecay[nSigModels][nSelTypes*4][histBins],weiDecay[nSigModels][nSelTypes*4][histBins];
  for(int nModel=0; nModel<nSigModels; nModel++) { for(unsigned int i=0; i<nSelTypes*4; i++) { for(int j=0; j<histBins; j++) {       
    bgdDecay[nModel][i][j] = 0.0; weiDecay[nModel][i][j] = 0.0; 
  }}}

  unsigned int numberOfLeptons = 2;
  TString signalName="";
  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());

    TFile *the_input_file = TFile::Open(infilenamev[ifile].Data());
    int nModel = (infilecatv[ifile]==6 || infilecatv[ifile]==7) ? signalIndex_[ifile] : -1;
    if(nModel>=0) signalName=signalName_[nModel];
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    //TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file->FindObjectAny("SelBit_tree");
    TTree *the_PDF_tree   = (TTree*)the_input_file->FindObjectAny("pdfReweight");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    eventJets.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BareTaus eventTaus;
    eventTaus.SetExtend();
    eventTaus.setBranchAddresses(the_input_tree);

    BarePhotons eventPhotons;
    eventPhotons.SetExtend();
    eventPhotons.setBranchAddresses(the_input_tree);

    BareMet eventMet;
    eventMet.SetExtend();
    eventMet.setBranchAddresses(the_input_tree);

    BareTrigger eventTrigger;
    eventTrigger.setBranchAddresses(the_input_tree);

    BareVertex eventVertex;
    eventVertex.setBranchAddresses(the_input_tree);

    BareMonteCarlo eventMonteCarlo;
    eventMonteCarlo.SetExtend();
    eventMonteCarlo.setBranchAddresses(the_input_tree);

    TNamed *triggerNames = (TNamed*)the_input_file->FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(ifile == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }
    else {
    }

    int initPDFTag = -1;

    bool errorMsgQCDscale = false;
    unsigned int selBit_= 0;
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    histoZHSEL[0]->Scale(0.0);
    histoZHSEL[1]->Scale(0.0);
    histoZHSEL[2]->Scale(0.0);
    histoZHSEL[3]->Scale(0.0);
    double countGenPhotons[3] = {0, 0, 0};
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    if(the_input_tree->GetEntries() != the_SelBit_tree->GetEntries()) {printf("BIG SKIMMING FAILURE\n"); return;}
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;
      the_input_tree->GetEntry(i);

      Bool_t passFilter[4] = {kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      if(infilecatv[ifile] != 999) {
	for (int nt = 0; nt <(int)numtokens; nt++) {
          if((*eventTrigger.triggerFired)[nt] == 0) continue;
          if((strcmp(tokens[nt],Form("HLT_Ele25_eta2p1_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Ele27_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_IsoMu24_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_IsoTkMu24_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Ele27_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Ele30_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Ele35_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_IsoMu22_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_IsoTkMu22_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu45_eta2p1_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu50_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	     (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0)
             ) passFilter[1] = kTRUE;
	}
      } else { passFilter[1] = kTRUE;}

      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;

      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {

        //if(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() <= 10) continue;

        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepBaseline )== BareLeptons::LepBaseline ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepVeto	  )== BareLeptons::LepVeto     ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake	  )== BareLeptons::LepFake     ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoft	  )== BareLeptons::LepSoft     ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepLoose    )== BareLeptons::LepLoose    ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepMedium   )== BareLeptons::LepMedium   ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepTight    )== BareLeptons::LepTight    ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepMediumIP )== BareLeptons::LepMediumIP ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepTightIP  )== BareLeptons::LepTightIP  ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP   )== BareLeptons::LepSoftIP   ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepVetoIso  )== BareLeptons::LepVetoIso  ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepLooseIso )== BareLeptons::LepLooseIso ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepMediumIso)== BareLeptons::LepMediumIso){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepTightIso )== BareLeptons::LepTightIso ){idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}

        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        //else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        //else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepLoose)  == BareLeptons::LepLoose ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(selectIdIsoCut("veto",TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}

      }
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()>=numberOfLeptons) passFilter[2] = kTRUE;
      
      if(passFilter[2] == kFALSE) continue; 

      if(idLep.size()==goodIsTight) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 25 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 20) continue;

      double dPhiLepMETMin = 999.;
      int signQ = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);

      vector<int> idB,idC;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if     (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idB.push_back(ngen0);
        else if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idC.push_back(ngen0);
      }

      vector<int> idPho,idPhoIsLep;
      for(int npho=0; npho<eventPhotons.p4->GetEntriesFast(); npho++) {
        if(TMath::Abs(((TLorentzVector*)(*eventPhotons.p4)[npho])->Eta()) >= 2.5) continue;
	if(((TLorentzVector*)(*eventPhotons.p4)[npho])->Pt() <= ptGMIN) continue;
	//if((double)(*eventPhotons.r9)[npho] <= 0.9) continue;
	//if((double)(*eventPhotons.sieie)[npho] >= 0.011) continue;
        bool isGoodPhoton = ((int)(*eventPhotons.selBits)[npho] & BarePhotons::PhoMedium)        == BarePhotons::PhoMedium &&
	                    ((int)(*eventPhotons.selBits)[npho] & BarePhotons::PhoElectronVeto)  == BarePhotons::PhoElectronVeto &&
			    ((int)(*eventPhotons.selBits)[npho] & BarePhotons::PhoPixelSeedVeto) == BarePhotons::PhoPixelSeedVeto;
        if(isGoodPhoton == false) continue;

        bool isRecoLepton = false;
	for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventPhotons.p4)[npho])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.1)
	    {isRecoLepton = true; break;}
        }
	if(isRecoLepton == false) idPho.push_back(npho);
	else                      idPhoIsLep.push_back(npho);
      }

      double minMassll = 999.0;
      double minMassZ = 999.0;
      int type3l = 0;
      int tagZ[4] = {-1,-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
	  if((int)(*eventLeptons.pdgId)[idLep[nl0]] * (int)(*eventLeptons.pdgId)[idLep[nl1]] > 0) continue;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl1])) ) ));

	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl1]]) &&
	     TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ-91.1876)) {
	     minMassZ = dilepAux.M();tagZ[0]=nl0;tagZ[1]=nl1;
	     if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 13) type3l = 0;
	     else                                                         type3l = 1;
	  }

	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }
      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if((int)nl0==tagZ[0]||(int)nl0==tagZ[1]) continue;
        if     (tagZ[2] == -1) tagZ[2] = nl0;
        else if(tagZ[3] == -1) tagZ[3] = nl0;
      }
      if(tagZ[2] != -1){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]) == 13) type3l += 0;
        else                                                             type3l += 2;
      }

      if(tagZ[0] == -1 && tagZ[1] != -1) printf("NOOOOOO1\n");
      if(tagZ[0] != -1 && tagZ[1] == -1) printf("NOOOOOO2\n");
      if(tagZ[0] == -1 || tagZ[1] == -1) {tagZ[0] = 0; tagZ[1] = 1; tagZ[2] = -1; tagZ[3] = -1;}

      bool tight3rdLepId = true;
      if(tagZ[2] != -1 && idTight[tagZ[2]] == 1 && idLep.size() == 3){
        tight3rdLepId = selectIdIsoCut("default",TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Eta()),(double)(*eventLeptons.iso)[idLep[tagZ[2]]],(int)(*eventLeptons.selBits)[idLep[tagZ[2]]],(double)(*eventLeptons.mva)[idLep[tagZ[2]]]);
      }
      if(tight3rdLepId == false) continue;

      bool isZBosonIn = false;
      TLorentzVector dilep (( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) ) )); 
      TLorentzVector dilepg(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) ) ));
      double massLG[2] = {0.,0.}; double dPhiLG[2] = {0.,0.}; double dPhiGMET = 0; double mTGMET = 0.0;
      if(idPho.size() >= 1 && idLep.size() == 2){
        TLorentzVector theG = ( *(TLorentzVector*)(eventPhotons.p4->At(idPho[0])));
        dilepg = dilepg + theG;
	massLG[0] = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + theG ).M();      
	massLG[1] = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) ) + theG ).M();
	dPhiLG[0] = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaPhi(theG));
	dPhiLG[1] = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->DeltaPhi(theG));
	dPhiGMET = TMath::Abs(theG.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
        mTGMET = TMath::Sqrt(2.0*theG.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(dPhiGMET)));
      }
      else if(idPho.size() == 0 && idLep.size() == 3 && tagZ[2] != -1){       
        TLorentzVector theG = ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[2]])) ); 
        dilepg = dilepg + theG;
	massLG[0] = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + theG ).M();      
	massLG[1] = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) ) + theG ).M();
	dPhiLG[0] = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaPhi(theG));
	dPhiLG[1] = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->DeltaPhi(theG));
	dPhiGMET = TMath::Abs(theG.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
        mTGMET = TMath::Sqrt(2.0*theG.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(dPhiGMET)));
      }
      else if(idPho.size() == 0 && idLep.size() == 4 && tagZ[2] != -1 && tagZ[3] != -1){     
        TLorentzVector theG = (( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[2]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[3]])) ) ));
	if((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]+(int)(*eventLeptons.pdgId)[idLep[tagZ[3]]] == 0 &&
	   theG.M() > 76.1876 && theG.M() < 106.1876) isZBosonIn = true;
        dilepg = dilepg + theG;
	massLG[0] = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + theG ).M();      
	massLG[1] = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) ) + theG ).M();
	dPhiLG[0] = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaPhi(theG));
	dPhiLG[1] = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->DeltaPhi(theG));
	dPhiGMET = TMath::Abs(theG.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
        mTGMET = TMath::Sqrt(2.0*theG.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(dPhiGMET)));
      }
      else if(idPho.size() == 0 && idLep.size() == 2){ }
      else if(idPho.size() >= 1 && idLep.size() == 3){ }
      else if(idPho.size() >= 1 && idLep.size() == 4){ }
      else if(idPho.size() == 0 && idLep.size() == 3 && tagZ[2] == -1){ }
      else if(idPho.size() == 0 && idLep.size() == 4 && tagZ[2] == -1 && tagZ[3] == -1){ }
      else if(idLep.size() >= 5){ }
      else {
        printf("PROB: %d %d %d %d\n",(int)idPho.size(),(int)idLep.size(),tagZ[2],tagZ[3]);
      }

      vector<int> idJet,idJetUp,idJetDown,idBJet;
      double total_bjet_probMEDIUM[2] = {1,1};double total_bjet_probMEDIUMUP[2] = {1,1};double total_bjet_probMEDIUMDOWN[2] = {1,1};
      bool isBtag = kFALSE;
      double sumPtJets = 0.0;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      TLorentzVector dilepJet = dilep;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 20) continue;
        //bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;

       if(((int)(*eventJets.selBits)[nj] & BareJets::JetLoose)!= BareJets::JetLoose) continue;

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.4) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

	if(infilecatv[ifile] != 0){
          BTagEntry::JetFlavor jetFlavor = BTagEntry::FLAV_UDSG;
	  for(unsigned int ng=0; ng<idB.size(); ng++) {
	    if (((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[idB[ng]])) < 0.4) {jetFlavor = BTagEntry::FLAV_B; break;}
	  }
	  if(jetFlavor == BTagEntry::FLAV_UDSG){
	    for(unsigned int ng=0; ng<idC.size(); ng++) {
	      if (((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[idC[ng]])) < 0.4) {jetFlavor = BTagEntry::FLAV_C; break;}
	    }
          }

          int nJPt = 0;
	  if     (((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) nJPt = 0;
	  else if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 40) nJPt = 1;
	  else if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 60) nJPt = 2;
	  else if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 80) nJPt = 3;
          else                                                       nJPt = 4;
          int nJEta = 0;
	  if     (TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 0.5) nJEta = 0;
	  else if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 1.0) nJEta = 1;
	  else if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 1.5) nJEta = 2;
	  else if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 2.0) nJEta = 3;
          else                                                                     nJEta = 4;
          denBTagging[nJEta][nJPt][jetFlavor]++;
          if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]) numBTaggingMEDIUM[nJEta][nJPt][jetFlavor]++;

          double bjet_SFMEDIUM = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFMEDIUM = btagReaderLMEDIUM.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFMEDIUM = btagReaderBCMEDIUM.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFMEDIUM == 0) bjet_SFMEDIUM = 1;
          double bjet_SFMEDIUMUP = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFMEDIUMUP = btagReaderLMEDIUMUP.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFMEDIUMUP = btagReaderBCMEDIUMUP.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFMEDIUMUP == 0) bjet_SFMEDIUMUP = 1;
          double bjet_SFMEDIUMDOWN = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFMEDIUMDOWN = btagReaderLMEDIUMDOWN.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFMEDIUMDOWN = btagReaderBCMEDIUMDOWN.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFMEDIUMDOWN == 0) bjet_SFMEDIUMDOWN = 1;

	  if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]){
	    total_bjet_probMEDIUM[0] = total_bjet_probMEDIUM[0] * jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor];
	    total_bjet_probMEDIUM[1] = total_bjet_probMEDIUM[1] * jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor] * bjet_SFMEDIUM;
	    total_bjet_probMEDIUMUP[0] = total_bjet_probMEDIUMUP[0] * jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor];
	    total_bjet_probMEDIUMUP[1] = total_bjet_probMEDIUMUP[1] * jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor] * bjet_SFMEDIUMUP;
	    total_bjet_probMEDIUMDOWN[0] = total_bjet_probMEDIUMDOWN[0] * jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor];
	    total_bjet_probMEDIUMDOWN[1] = total_bjet_probMEDIUMDOWN[1] * jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor] * bjet_SFMEDIUMDOWN;
	  } else {
	    total_bjet_probMEDIUM[0] = total_bjet_probMEDIUM[0] * (1.0 - jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor]);
	    total_bjet_probMEDIUM[1] = total_bjet_probMEDIUM[1] * (1.0 - jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor] * bjet_SFMEDIUM);
	    total_bjet_probMEDIUMUP[0] = total_bjet_probMEDIUMUP[0] * (1.0 - jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor]);
	    total_bjet_probMEDIUMUP[1] = total_bjet_probMEDIUMUP[1] * (1.0 - jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor] * bjet_SFMEDIUMUP);
	    total_bjet_probMEDIUMDOWN[0] = total_bjet_probMEDIUMDOWN[0] * (1.0 - jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor]);
	    total_bjet_probMEDIUMDOWN[1] = total_bjet_probMEDIUMDOWN[1] * (1.0 - jetEpsBtagMEDIUM[nJEta][nJPt][jetFlavor] * bjet_SFMEDIUMDOWN);
	  }

        }

        Bool_t isPhoton = kFALSE;
        if(idPho.size() >= 1){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])) < 0.3) isPhoton = kTRUE;
	}
	if(isPhoton == kTRUE) continue;

        if(dPhiJetMET   == -1 && ((TLorentzVector*)(*eventJets.p4)[nj])->Pt()> 30) {
          dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
        }

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20) {
           sumPtJets = sumPtJets + ((TLorentzVector*)(*eventJets.p4)[nj])->Pt();
	   if ((float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];
           if ((float)(*eventJets.bDiscr)[nj] > bTagCuts[0]) idBJet.push_back(nj);
        }
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()      > 30) {idJet.push_back(nj); dilepJet = dilepJet + ( *(TLorentzVector*)(eventJets.p4->At(nj)) );}
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*1.03 > 30) idJetUp.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*0.97 > 30) idJetDown.push_back(nj);
      }

      int numberGoodTaus = 0;
      for(int ntau=0; ntau<eventTaus.p4->GetEntriesFast(); ntau++) {
        if(((TLorentzVector*)(*eventTaus.p4)[ntau])->Pt() <= 18.0 ||
           TMath::Abs(((TLorentzVector*)(*eventTaus.p4)[ntau])->Eta()) >= 2.3) continue;
        bool isElMu = false;
        for(unsigned nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventTaus.p4)[ntau])) < 0.3) {
            isElMu = true;
            break;
          }
        }
        if(isElMu == false &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFinding	 ) == BareTaus::TauDecayModeFinding &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFindingNewDMs) == BareTaus::TauDecayModeFindingNewDMs &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::byVLooseIsolationMVArun2v1DBnewDMwLT) == BareTaus::byVLooseIsolationMVArun2v1DBnewDMwLT){
           //(double)(*eventTaus.iso)[ntau] < 5.0){
          numberGoodTaus++;
        }
      }

      TLorentzVector theCaloMET;
      theCaloMET.SetPx((double)eventMet.CaloMet->Px());
      theCaloMET.SetPy((double)eventMet.CaloMet->Py());
      for(unsigned int nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]])==13) {
	  theCaloMET.SetPx(theCaloMET.Px()-((TLorentzVector*)(*eventLeptons.p4)[nl])->Px());
	  theCaloMET.SetPy(theCaloMET.Py()-((TLorentzVector*)(*eventLeptons.p4)[nl])->Py());
	}
      }

      // Determine flavor of the pair (0 means e-mu pair, 1 means mu-mu, 2 means e-e)
      int typePair = 0;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==13) typePair = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==11) typePair = 2;

      // Calculate a lot of physics quantities used for the rectangular selection
      double deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));

      double dPhiDiLepGMET = TMath::Abs(dilepg.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));

      TVector2 metv(((TLorentzVector*)(*eventMet.p4)[0])->Px(), ((TLorentzVector*)(*eventMet.p4)[0])->Py());
      TVector2 dilv(dilep.Px(), dilep.Py());
      TVector2 utv = -1.*(metv+dilv);
      double phiv = utv.DeltaPhi(dilv);
      double the_upara = TMath::Abs(utv.Mod()*TMath::Cos(phiv))/dilep.Pt();
      
      bool passZMass = minMassll > 4.0 && dilep.M() > 76.1876 && dilep.M() < 106.1876;
      bool passNjets = idJet.size() <= nJetsType;

      double ptSkimCut = 60.0; if(whichSkim == 0) ptSkimCut = 0.0;
      bool passMETMin  = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > ptSkimCut;
      bool passPTLL    = dilep.Pt() > ptSkimCut;
      bool passMET     = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > metMIN;

      double ptFracG = TMath::Abs(dilepg.Pt()-((TLorentzVector*)(*eventMet.p4)[0])->Pt())/dilepg.Pt();
      bool passPTFracG   = ptFracG < 0.4;
      bool passDPhiZGMET = dPhiDiLepGMET > 2.5;

      bool passBtagVeto  = bDiscrMax < bTagCuts[0];
      bool pass3rdLVeto  = idLep.size() == numberOfLeptons && TMath::Abs(signQ) == 0;
      bool pass3rdLSel   = idLep.size() == 3 && TMath::Abs(signQ) == 1 && idPho.size() == 0 && tagZ[2] != -1;
      bool pass4thLSel   = idLep.size() == 4 && TMath::Abs(signQ) == 0 && idPho.size() == 0 && tagZ[2] != -1 && tagZ[3] != -1 && isZBosonIn == true;
      double dphill = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]));
      double detall = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Eta()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Eta());
      double drll = sqrt(dphill*dphill+detall*detall);

      bool passZMassSB = (dilep.M() > 110.0 && dilep.M() < 200.0);

      bool passDPhiJetMET   = dPhiJetMET == -1 || dPhiJetMET >= 0.5;
      bool passTauVeto      = numberGoodTaus == 0;
      bool passPhotonSel    = idPho.size() >= 1;
      bool passLepPhotonSel = idPhoIsLep.size() == 1;
      bool passMT           = mTGMET < 175.0;

      //0            1                2            3          4              5                6               7           8               9                 10           11
      bool passNMinusOne[12] = {
                     passPhotonSel && passNjets && passMET && passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
        passZMass &&                  passNjets && passMET && passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
        passZMass && passPhotonSel &&              passMET && passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
        passZMass && passPhotonSel && passNjets &&            passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
        passZMass && passPhotonSel && passNjets && passMET &&                passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
        passZMass && passPhotonSel && passNjets && passMET && passPTFracG &&                  passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
        passZMass && passPhotonSel && passNjets && passMET && passPTFracG && passDPhiZGMET &&                 passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
        passZMass && passPhotonSel && passNjets && passMET && passPTFracG && passDPhiZGMET && passBtagVeto &&             pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,
	passZMass && passPhotonSel && passNjets && passMET && passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL &&                 passDPhiJetMET && passTauVeto && passMT,
	passZMass && passPhotonSel && passNjets && passMET && passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto &&                   passTauVeto && passMT,
	passZMass && passPhotonSel && passNjets && passMET && passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET                && passMT,
	passZMass && passPhotonSel && passNjets && passMET && passPTFracG && passDPhiZGMET && passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto
                               };

      bool passAllCuts[nSelTypes] = {                   
                         passZMass && passPhotonSel && passNjets && passMET    && passPTFracG && passDPhiZGMET &&  passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,  // ZHGSEL
                         passZMass && passPhotonSel && passNjets && passMET    && passPTFracG && passDPhiZGMET && !passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,  // BTAGSEL
      passZMassSB    && !passZMass && passPhotonSel && passNjets && passMET    && passPTFracG && passDPhiZGMET &&  passBtagVeto && passPTLL && pass3rdLVeto && passDPhiJetMET && passTauVeto && passMT,  // WWSEL
                         passZMass && passPhotonSel &&              passMET    &&                                                  passPTLL && pass3rdLVeto &&                   passTauVeto          ,  // PRESEL
                         passZMass &&                               passMETMin &&                                  passBtagVeto && passPTLL && pass3rdLVeto &&		         passTauVeto          ,  // ZLLSEL
                         passZMass && passLepPhotonSel &&           passMETMin &&                                  passBtagVeto && passPTLL && pass3rdLVeto &&		         passTauVeto          ,  // ZLGSEL
                         passZMass &&                  passNjets && passMET    && passPTFracG && passDPhiZGMET &&  passBtagVeto && passPTLL && pass3rdLSel  && passDPhiJetMET &&                passMT,  // WZSEL
                         passZMass &&                  passNjets && passMETMin && passPTFracG && passDPhiZGMET &&  passBtagVeto && passPTLL && pass4thLSel  && passDPhiJetMET &&                passMT,  // ZZSEL
                         passZMass &&                               passMETMin &&                passDPhiZGMET &&                  passPTLL && pass4thLSel                                               // ZZLOOSESEL
                                    };

     bool passEvolFilter[numberCuts] = {pass3rdLVeto,passPTLL,passZMass,passPhotonSel,passBtagVeto,passTauVeto,passNjets,passMET,passDPhiZGMET,passPTFracG,passDPhiJetMET,passMT};

     int sumEvol = 0;
     bool totalSel = kTRUE;
     for(int isel=0; isel<numberCuts; isel++) {
       totalSel = totalSel && passEvolFilter[isel];
       if(totalSel == kTRUE) sumEvol++;
     }
     // Syst cuts for MVA
     bool passSystCuts[nSystTypes] = {
          passNMinusOne[2] && idJetUp.size()   <= nJetsType,
          passNMinusOne[2] && idJetDown.size() <= nJetsType,
	  passNMinusOne[3] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() > metMIN,
	  passNMinusOne[3] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > metMIN
     };
      
      // begin event weighting
      vector<int>wBoson;
      vector<int>zBoson;
      vector<bool> isGenDupl;double bosonPtMin = 1000000000; bool isBosonFound = false;vector<bool> isNeuDupl;vector<int> idGenPho;vector<int> idGenLep;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) wBoson.push_back(ngen0);
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23) zBoson.push_back(ngen0);
        if((TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23||TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) &&
	   ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() < bosonPtMin) {bosonPtMin = ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt(); isBosonFound = true;}
        // begin neutrinos
        isNeuDupl.push_back(0);
	bool isGoodNFlags = ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::PromptFinalState) == BareMonteCarlo::PromptFinalState ||
            		   ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::DirectPromptTauDecayProductFinalState) == BareMonteCarlo::DirectPromptTauDecayProductFinalState;
        isGoodNFlags = isGoodNFlags && (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 12 || TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 14 || TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 16);
        if(isGoodNFlags == false) isNeuDupl[ngen0] = 1;

	// begin leptons	
        isGenDupl.push_back(0);
	bool isGoodFlags = ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::PromptFinalState) == BareMonteCarlo::PromptFinalState ||
            		   ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::DirectPromptTauDecayProductFinalState) == BareMonteCarlo::DirectPromptTauDecayProductFinalState;
        isGoodFlags = isGoodFlags && (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 11 || TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 13);
        if(isGoodFlags == false) {isGenDupl[ngen0] = 1; idGenLep.push_back(ngen0);}

	// begin photons	
	bool isGoodPhFlags = ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::PromptFinalState) == BareMonteCarlo::PromptFinalState ||
            		     ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::DirectPromptTauDecayProductFinalState) == BareMonteCarlo::DirectPromptTauDecayProductFinalState;
        isGoodPhFlags = isGoodPhFlags && (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 22);
        if(isGoodPhFlags == true) idGenPho.push_back(ngen0);
      }
      if(isBosonFound==false) bosonPtMin = 0;
      int numberGoodGenLep[3] = {0,0,0};
      TLorentzVector the_rhoP4(0,0,0,0);
      for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
        if(isNeuDupl[ngen] == 0 || isGenDupl[ngen] == 0) {
	  the_rhoP4 = the_rhoP4 + *(TLorentzVector*)(*eventMonteCarlo.p4)[ngen];
	}
        if(isNeuDupl[ngen] == 0) numberGoodGenLep[0]++;
	if(isGenDupl[ngen] == 1) continue;
	numberGoodGenLep[1]++;
	if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt() <= 10 ||
	   TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()) >= 2.5) continue;
	numberGoodGenLep[2]++;
      }
      vector<int> isGenLep; unsigned int goodIsGenLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenLepton = false;
        for(unsigned int ngenlep=0; ngenlep<idGenLep.size(); ngenlep++) {
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[idGenLep[ngenlep]]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[idGenLep[ngenlep]])) < 0.1) {
	    isGenLepton = true;
	    break;
	  }
	}
	if(isGenLepton == true) {isGenLep.push_back(1); goodIsGenLep++;}
	else                    {isGenLep.push_back(0);}
      }

      int isGenPho = 0;
      bool isGenPhoton = false;
      if(idPho.size() >= 1){
        for(unsigned int ngenpho=0; ngenpho<idGenPho.size(); ngenpho++) {
          if(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[idGenPho[ngenpho]])) < 0.1) {
	    isGenPhoton = true;
	    break;
	  }
	}
	if(isGenPhoton == true) isGenPho = 1;
	if(isGenPho == 0){
          bool isGenLepton = false;
          for(unsigned int ngenlep=0; ngenlep<idGenLep.size(); ngenlep++) {
            if(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[idGenLep[ngenlep]])) < 0.1) {
	      isGenLepton = true;
	      break;
	    }
	  }
	  if(isGenLepton == true) isGenPho = 2;
        }
      }

      double photonSF = 1.0;
      if(isGenPho == 2) {
        if(TMath::Abs(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Eta()) < 1.5) photonSF = sf_el_gamma[0]; 
	else                                                                         photonSF = sf_el_gamma[1]; 
	if(passAllCuts[ZHGSEL]) countGenPhotons[2]++;
      }
      else if(isGenPho == 1) {
        photonSF = effhDPhotonScaleFactor(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt(), ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Eta(), "medium", fhDPhotonSF, fhDElectronVetoSF);
        if(passAllCuts[ZHGSEL]) countGenPhotons[1]++;
      }
      else if(idPho.size() >= 1) {if(passAllCuts[ZHGSEL]) countGenPhotons[0]++;}

      // trigger efficiency
      double trigEff = 1.0;
      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight     = 1.0; if(infilecatv[ifile] != 0) puWeight     = nPUScaleFactor(fhDPU    , (double)eventMonteCarlo.puTrueInt);
      double puWeightUp   = 1.0; if(infilecatv[ifile] != 0) puWeightUp   = nPUScaleFactor(fhDPUUp  , (double)eventMonteCarlo.puTrueInt);
      double puWeightDown = 1.0; if(infilecatv[ifile] != 0) puWeightDown = nPUScaleFactor(fhDPUDown, (double)eventMonteCarlo.puTrueInt);
      // lepton efficiency
      double effSF = 1.0;
      if(infilecatv[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,fhDVeryTightSF,true);
        }
      }

      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff*photonSF;
      //printf("totalWeight: %f * %f * %f * %f * %f * %f * %f * %f = %f\n",mcWeight,theLumi,puWeight,effSF,fakeSF,theMCPrescale,trigEff,photonSF,totalWeight);

      // Btag scale factor
      totalWeight = totalWeight * total_bjet_probMEDIUM[1]/total_bjet_probMEDIUM[0];

      double btagCorr[2] = {(total_bjet_probMEDIUMUP[1]  /total_bjet_probMEDIUMUP[0]  )/(total_bjet_probMEDIUM[1]/total_bjet_probMEDIUM[0]),
                            (total_bjet_probMEDIUMDOWN[1]/total_bjet_probMEDIUMDOWN[0])/(total_bjet_probMEDIUM[1]/total_bjet_probMEDIUM[0])};

      if(theCategory == 2 && idPho.size() > 0 && idGenPho.size() > 0 && infilenamev[ifile].Contains("ZGTo2LG") == kFALSE){
        if(isGenPhoton == true) totalWeight = 0;
      }

      if(totalWeight == 0) continue;

      // ZH EWK correction (only for SM case)
      if(theCategory == 6 && zBoson.size() >= 1 && nModel == 0) {
        totalWeight = totalWeight * weightZHEWKCorr(fhDZHEwkCorr, ((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->Pt());
      }
      else if(theCategory == 6 && zBoson.size() == 0 && nModel == 0) {
        printf("zBoson = 0\n");
      }

      // ZZ
      double the_rho = 0.0; if(the_rhoP4.P() > 0) the_rho = the_rhoP4.Pt()/the_rhoP4.P();
      double theZZCorr[2] {1,1};
      if(theCategory == 4 && infilenamev[ifile].Contains("GluGlu") == kFALSE) {
	theZZCorr[0] = weightEWKCorr(bosonPtMin,1);

        //float GENdPhiZZ = 5;
	//if(zBoson.size() >= 2) GENdPhiZZ = TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->DeltaPhi(*((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[1]])));
	//theZZCorr[1] = kfactor_qqZZ_qcd_dPhi(GENdPhiZZ);
        float GENmZZ = 0.0;
	if(zBoson.size() >= 2) GENmZZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zBoson[1])) ) ).M();
	theZZCorr[1] = kfactor_qqZZ_qcd_M(GENmZZ);
        //float GENptZZ = 0.0;
	//if(zBoson.size() >= 2) GENptZZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zBoson[1])) ) ).Pt();
	//theZZCorr[1] = kfactor_qqZZ_qcd_M(GENptZZ);

        totalWeight = totalWeight * (theZZCorr[0]*theZZCorr[1]);
      }

      // WZ
      if(theCategory == 3 && wBoson.size() >= 1 && zBoson.size() >= 1) {
        float GENmWZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(wBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zBoson[0])) ) ).M();

        totalWeight = totalWeight * weightEWKWZCorr(GENmWZ);
        //printf("Possible to perform WZ EWK correction: %d %d\n",(int)wBoson.size(),(int)zBoson.size());
      }
      //else if(theCategory == 3) {
      //  printf("No possible to perform WZ EWK correction: %d %d\n",(int)wBoson.size(),(int)zBoson.size());
      //}
      // end event weighting
      //totalWeight = 1;

      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[ZHGSEL]) sumEventsProcess[ifile] += totalWeight;

      for(int nl=0; nl <=sumEvol; nl++) histo[allPlots-2][theCategory]->Fill((double)nl,totalWeight);
      for(int nl=0; nl <=sumEvol; nl++) histoZHSEL[typePair ]         ->Fill((double)nl,totalWeight);
      if(typePair == 1 || typePair == 2)
      for(int nl=0; nl <=sumEvol; nl++) histoZHSEL[3]                 ->Fill((double)nl,totalWeight);

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i] && (theCategory == 6 || theCategory == 7)) {
          bgdDecay[nModel][i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[nModel][i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        } else if(passAllCuts[i]) { for(int mModel=0; mModel<nSigModels; mModel++) { 
          bgdDecay[mModel][i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[mModel][i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        }}
      }

      if((typeSel == typePair) || (typeSel == 3 && (typePair == 1 || typePair == 2))) {
	if(nModel<0 || nModel==plotModel){
	  for(int thePlot=0; thePlot<allPlots-2; thePlot++){
	    double theVar = 0.0;
	    bool makePlot = false;
	    if     (thePlot ==  0 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(mtW,999.999);}
	    else if(thePlot ==  1 && passNMinusOne[0])       {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.M()-91.1876),99.999);}
	    else if(thePlot ==  2 && passNMinusOne[2])       {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	    else if(thePlot ==  3 && passNMinusOne[3])       {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),399.999);}
	    else if(thePlot ==  4 && passNMinusOne[4])       {makePlot = true;theVar = TMath::Min(ptFracG,0.999);}
	    else if(thePlot ==  5 && passNMinusOne[5])       {makePlot = true;theVar = dPhiDiLepGMET;}
	    else if(thePlot ==  6 && passNMinusOne[6])       {makePlot = true;theVar = TMath::Max(TMath::Min(bDiscrMax,0.999),0.001);}
	    else if(thePlot ==  7 && passNMinusOne[7])       {makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
	    else if(thePlot ==  8 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Pt(),199.999);}
	    else if(thePlot ==  9 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Pt(),199.999);}
	    else if(thePlot == 10 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(TMath::Abs(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Eta()),2.499);}
	    else if(thePlot == 11 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	    else if(thePlot == 12 && passNMinusOne[9])       {makePlot = true;theVar = dPhiJetMET;}
	    else if(thePlot == 13 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = dPhiLepMETMin;}
	    else if(thePlot == 14 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = dphill;}
	    else if(thePlot == 15 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),399.999);}
	    else if(thePlot == 16 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min(ptFracG,0.999);}
	    else if(thePlot == 17 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
	    else if(thePlot == 18 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = (double)(numberGoodGenLep[0]+10*numberGoodGenLep[1]+100*numberGoodGenLep[2]);}
	    else if(thePlot == 19 && passNMinusOne[10])      {makePlot = true;theVar = TMath::Min((double)numberGoodTaus,3.499);}
	    else if(thePlot == 20 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.Eta()),2.499);}
	    else if(thePlot == 21 && passNMinusOne[1])       {makePlot = true;theVar = idPho.size();}
	    else if(thePlot == 22 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(the_rho,0.999);}
	    else if(thePlot == 23 && passAllCuts[PRESEL])    {makePlot = true;theVar = idJet.size();}
	    else if(thePlot == 24 && passAllCuts[PRESEL])    {makePlot = true;theVar = idBJet.size();}
	    else if(thePlot == 25 && passNMinusOne[8])       {makePlot = true;theVar = idLep.size();;}
	    else if(thePlot == 26 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(drll,2.999);}
	    else if(thePlot == 27 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(TMath::Min(massLG[0],massLG[1]),199.999);}
	    else if(thePlot == 28 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min(dPhiLG[0],dPhiLG[1]);}
	    else if(thePlot == 29 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = dPhiGMET;}
	    else if(thePlot == 30 && passNMinusOne[11])      {makePlot = true;theVar = TMath::Min(mTGMET,249.999);}
	    else if(thePlot == 31 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min(mTGMET,249.999);}
	    else if(thePlot == 32 && passAllCuts[WZSEL])     {makePlot = true;theVar = mTGMET;}
	    else if(thePlot == 33 && passAllCuts[ZZLOOSESEL]){makePlot = true;theVar = mTGMET;}
	    else if(thePlot == 34 && passAllCuts[WZSEL])     {makePlot = true;theVar = type3l;}

	    if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
	  }
        }
      }

      if((typeSel == typePair) || (typeSel == 3 && (typePair == 1 || typePair == 2))) {
	double MVAVar = TMath::Min(mTGMET,xbins[nBinMVA]-0.0001);
	if(MVAVarType == 1 && idPho.size() >= 1) {
	  MVAVar = mTGMET;
	  if(TMath::Abs(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Eta()) > 1.5) MVAVar = MVAVar + 200.0;
	}
	double MVAVarMETSyst[2] = {MVAVar, MVAVar};	

	// Avoid QCD scale weights that are anomalous high
	double maxQCDscale = (TMath::Abs((double)eventMonteCarlo.r1f2)+TMath::Abs((double)eventMonteCarlo.r1f5)+TMath::Abs((double)eventMonteCarlo.r2f1)+
                              TMath::Abs((double)eventMonteCarlo.r2f2)+TMath::Abs((double)eventMonteCarlo.r5f1)+TMath::Abs((double)eventMonteCarlo.r5f5))/6.0;
        double PDFAvg = 0.0;
        if(infilecatv[ifile] != 0 && passAllCuts[ZHGSEL]){
          if(initPDFTag != -1)
          for(int npdf=0; npdf<100; npdf++) PDFAvg = PDFAvg + TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]);
          PDFAvg = PDFAvg/100.0;
        }

        if     (theCategory == 0){
	  if(passAllCuts[ZHGSEL]) histo_Data->Fill(MVAVar,totalWeight);
	  if((passAllCuts[ZHGSEL] || passAllCuts[ZZSEL] || passAllCuts[ZZLOOSESEL]) && verbose) {
            if      (typePair==1) printf("mm data event: ");
            else if (typePair==2) printf("ee data event: " );
            else if (typePair==3) printf("em data event: " );
            printf("runnumber %d lumisection %d eventnumber %lld (%d/%d/%d) ptll %f MET %f njets %d l1_pt %f l1_eta %f l2_pt %f l2_eta %f balance %f mtg: %f\n",
              eventEvent.runNum,
              eventEvent.lumiNum,
              eventEvent.eventNum,
	      passAllCuts[ZHGSEL],passAllCuts[ZZSEL],passAllCuts[ZZLOOSESEL],
              dilep.Pt(),
              ((TLorentzVector*)(*eventMet.p4)[0])->Pt(),
              (int)idJet.size(),
              ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Pt(), ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Eta(),
              ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Pt(), ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Eta(),
              ptFracG,mTGMET
            );
          }
        }
        else if(theCategory == 1){
	  if(passAllCuts[ZHGSEL]) histo_EM   ->Fill(MVAVar,totalWeight);
	  if(passAllCuts[ZHGSEL]) histo_EMNoW->Fill(MVAVar,1.);
        }
        else if(theCategory == 2){
	  if(passAllCuts[ZHGSEL]) histo_Zjets   ->Fill(MVAVar,totalWeight);
	  if(passAllCuts[ZHGSEL]) histo_ZjetsNoW->Fill(MVAVar,1.);
        }
        else if(theCategory == 3){
	  if(passAllCuts[ZHGSEL]) {
	     histo_WZ              ->Fill(MVAVar,totalWeight);
	     histo_WZNoW           ->Fill(MVAVar,1.);
	     histo_WZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*1.02);
	     histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_WZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_WZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_WZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_WZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
	}
        else if(theCategory == 4){
	  if(passAllCuts[ZHGSEL]) {
	     histo_ZZ              ->Fill(MVAVar,totalWeight);
	     histo_ZZNoW           ->Fill(MVAVar,1.);
	     if(infilenamev[ifile].Contains("GluGlu") == kFALSE) {
	       if(the_rho <= 0.3) histo_ZZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*(1.0+TMath::Abs((theZZCorr[0]-1)*(15.99/9.89-1))));
	       else               histo_ZZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*(1.0+TMath::Abs((theZZCorr[0]-1)               )));
	     } else {
               histo_ZZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight);
	     }
	     if(infilenamev[ifile].Contains("GluGlu") == kFALSE) histo_ZZ_CMS_ggCorrUp->Fill(MVAVar,totalWeight);
	     else                                                histo_ZZ_CMS_ggCorrUp->Fill(MVAVar,totalWeight*1.30);
	     histo_ZZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_ZZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ZZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_ZZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
          }
          if(passSystCuts[JESUP])  histo_ZZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 5){
	  if(passAllCuts[ZHGSEL]) {
	     histo_VVV   ->Fill(MVAVar,totalWeight);
	     histo_VVVNoW->Fill(MVAVar,1.);
	     histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_VVV_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_VVV_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
          }
          if(passSystCuts[JESUP])  histo_VVV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VVV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 6){
	  if(passAllCuts[ZHGSEL]) {
	     histo_ZH_hinv[nModel]   ->Fill(MVAVar,totalWeight);
	     histo_ZH_hinvNoW[nModel]->Fill(MVAVar,1.);
	     if(zBoson.size() >= 1 && nModel == 0) {
	       histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->Fill(MVAVar,totalWeight*weightZHEWKCorr(fhDZHEwkCorrUp,   ((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->Pt())/weightZHEWKCorr(fhDZHEwkCorr, ((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->Pt()));
	       histo_ZH_hinv_CMS_EWKCorrDown[nModel]->Fill(MVAVar,totalWeight*weightZHEWKCorr(fhDZHEwkCorrDown, ((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->Pt())/weightZHEWKCorr(fhDZHEwkCorr, ((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->Pt()));
             } else {
	       histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->Fill(MVAVar,totalWeight);
	       histo_ZH_hinv_CMS_EWKCorrDown[nModel]->Fill(MVAVar,totalWeight);
	     }
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_ZH_hinv_CMS_PDFBounding[nModel][npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_ZH_hinv_CMS_PDFBounding[nModel][npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingAvg [nModel]->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingAvg [nModel]->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingUp  [nModel]->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingUp  [nModel]->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Fill(MVAVar,totalWeight*0.98);
             histo_ZH_hinv_CMS_PUBoundingUp  [nModel]->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZH_hinv_CMS_PUBoundingDown[nModel]->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_ZH_hinv_CMS_MVAJESBoundingUp  [nModel]->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZH_hinv_CMS_MVAMETBoundingUp  [nModel]->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 7){
	  if(passAllCuts[ZHGSEL]) {
	     histo_ggZH_hinv   ->Fill(MVAVar,totalWeight);
	     histo_ggZH_hinvNoW->Fill(MVAVar,1.0);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_ggZH_hinv_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_ggZH_hinv_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_ggZH_hinv_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ggZH_hinv_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ggZH_hinv_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_ggZH_hinv_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ggZH_hinv_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ggZH_hinv_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      }
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);
    if(countGenPhotons[0]+countGenPhotons[1]+countGenPhotons[2] > 0)
    printf("GenPhotons fake/real/electron(%f) = %f %f %f\n",
    countGenPhotons[0]+countGenPhotons[1]+countGenPhotons[2],
    countGenPhotons[0]/(countGenPhotons[0]+countGenPhotons[1]+countGenPhotons[2]),
    countGenPhotons[1]/(countGenPhotons[0]+countGenPhotons[1]+countGenPhotons[2]),
    countGenPhotons[2]/(countGenPhotons[0]+countGenPhotons[1]+countGenPhotons[2]));
    for(int nc=0; nc<numberCuts+1; nc++){
      printf("(%20s): %10.2f %10.2f %10.2f %10.2f\n",cutName[nc].Data(),histoZHSEL[0]->GetBinContent(nc+1),histoZHSEL[1]->GetBinContent(nc+1),histoZHSEL[2]->GetBinContent(nc+1),histoZHSEL[3]->GetBinContent(nc+1));
    }
    the_input_file->Close();
  } // end of chain

  // "-1" to remove the Higgs contribution
  double sumEvents = 0;
  for(int np=1; np<histBins-1; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  //printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  //printf(" = %.3f\n",sumEvents);

  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Printing yields and statistical uncertainties for all signal models\n\n");
  for(int nModel=0; nModel<nSigModels; nModel++) {
    printf("Model: %s\n", signalName_[nModel].Data()); 
    printf("                    em                      mm                      ee                      ll\n");
    printf("-----------------------------------------------------------------------------------------------------------\n");
    for(int ns=0; ns<nSelTypes; ns++) {
      printf("Selection: %s\n",selTypeName[ns].Data());
      double sumEventsType[4] = {0,0,0,0}; double sumEventsTypeE[4] = {0,0,0,0};
      for(int np=0; np<histBins; np++) {       
         bgdDecay[nModel][ns+nSelTypes*3][np] = bgdDecay[nModel][ns+nSelTypes*1][np] + bgdDecay[nModel][ns+nSelTypes*2][np];
         weiDecay[nModel][ns+nSelTypes*3][np] = weiDecay[nModel][ns+nSelTypes*1][np] + weiDecay[nModel][ns+nSelTypes*2][np];
         printf("(%6s): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
         processName[np].Data(),bgdDecay[nModel][ns+nSelTypes*0][np],sqrt(weiDecay[nModel][ns+nSelTypes*0][np]),bgdDecay[nModel][ns+nSelTypes*1][np],sqrt(weiDecay[nModel][ns+nSelTypes*1][np]),
                                bgdDecay[nModel][ns+nSelTypes*2][np],sqrt(weiDecay[nModel][ns+nSelTypes*2][np]),bgdDecay[nModel][ns+nSelTypes*3][np],sqrt(weiDecay[nModel][ns+nSelTypes*3][np]));
         if(np!=0 && np!=6 && np!=7){
           sumEventsType[0] = sumEventsType[0] + bgdDecay[nModel][ns+nSelTypes*0][np]; sumEventsTypeE[0] = sumEventsTypeE[0] + weiDecay[nModel][ns+nSelTypes*0][np];
           sumEventsType[1] = sumEventsType[1] + bgdDecay[nModel][ns+nSelTypes*1][np]; sumEventsTypeE[1] = sumEventsTypeE[1] + weiDecay[nModel][ns+nSelTypes*1][np];
           sumEventsType[2] = sumEventsType[2] + bgdDecay[nModel][ns+nSelTypes*2][np]; sumEventsTypeE[2] = sumEventsTypeE[2] + weiDecay[nModel][ns+nSelTypes*2][np];
           sumEventsType[3] = sumEventsType[3] + bgdDecay[nModel][ns+nSelTypes*3][np]; sumEventsTypeE[3] = sumEventsTypeE[3] + weiDecay[nModel][ns+nSelTypes*3][np];
         }
      }
      printf("(...bkg): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
             sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
  	     sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]));
    }
  }

  // This is the alternative method
  double kEoverM  = sqrt(bgdDecay[0][ZLLSEL+nSelTypes*2][0]/bgdDecay[0][ZLLSEL+nSelTypes*1][0]);
  double NemFact[3] = {(1.0/kEoverM)*0.5, (kEoverM)*0.5, (kEoverM+1.0/kEoverM)*0.5};

  // This is the default method
  double NemFact_FromMLLSB[3] = {NemFact[0], NemFact[1], NemFact[2]}; double NemFact_FromMLLSBE[3] = {1.0, 1.0, 1.0};
  if(bgdDecay[0][WWSEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWSEL+nSelTypes*1][0]) > 0) {
    NemFact_FromMLLSB[0]  = (bgdDecay[0][WWSEL+nSelTypes*1][0])/bgdDecay[0][WWSEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[0] = sqrt(1./bgdDecay[0][WWSEL+nSelTypes*1][0]+1./bgdDecay[0][WWSEL+nSelTypes*0][0])*NemFact_FromMLLSB[0];
  }
  if(bgdDecay[0][WWSEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWSEL+nSelTypes*2][0]) > 0) {
    NemFact_FromMLLSB[1]  = (bgdDecay[0][WWSEL+nSelTypes*2][0])/bgdDecay[0][WWSEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[1] = sqrt(1./bgdDecay[0][WWSEL+nSelTypes*2][0]+1./bgdDecay[0][WWSEL+nSelTypes*0][0])*NemFact_FromMLLSB[1];
  }
  if(bgdDecay[0][WWSEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWSEL+nSelTypes*3][0]) > 0) {
    NemFact_FromMLLSB[2]  = (bgdDecay[0][WWSEL+nSelTypes*3][0])/bgdDecay[0][WWSEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[2] = sqrt(1./bgdDecay[0][WWSEL+nSelTypes*3][0]+1./bgdDecay[0][WWSEL+nSelTypes*0][0])*NemFact_FromMLLSB[2];
  }
  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Computing the flavor k-factors used for the EM background\n\n");

  printf("(mm) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[0],NemFact_FromMLLSB[0],NemFact_FromMLLSBE[0]);
  printf("(ee) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[1],NemFact_FromMLLSB[1],NemFact_FromMLLSBE[1]);
  printf("(ll) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[2],NemFact_FromMLLSB[2],NemFact_FromMLLSBE[2]);

  // There uncertainties: closure test, different between default and alternative method, and data statistics
  // Only the second item is used as systematic uncertainty
  double EMSystTotal[3] = {1.0,1.0,1.0}; double EMSyst[2][3] = {0.0,0.0,0.0,0.0,0.0,0.0};

  EMSyst[0][0] = bgdDecay[0][ZHGSEL+nSelTypes*1][1]/bgdDecay[0][ZHGSEL+nSelTypes*0][1]/NemFact[0];
  if(EMSyst[0][0] < 1.0) EMSyst[0][0] = 1/EMSyst[0][0]; EMSyst[0][0] = EMSyst[0][0] - 1.0;

  EMSyst[0][1] = bgdDecay[0][ZHGSEL+nSelTypes*2][1]/bgdDecay[0][ZHGSEL+nSelTypes*0][1]/NemFact[1];
  if(EMSyst[0][1] < 1.0) EMSyst[0][1] = 1/EMSyst[0][1]; EMSyst[0][1] = EMSyst[0][1] - 1.0;

  EMSyst[0][2] = bgdDecay[0][ZHGSEL+nSelTypes*3][1]/bgdDecay[0][ZHGSEL+nSelTypes*0][1]/NemFact[2];
  if(EMSyst[0][2] < 1.0) EMSyst[0][2] = 1/EMSyst[0][2]; EMSyst[0][2] = EMSyst[0][2] - 1.0;

  EMSyst[1][0] = NemFact_FromMLLSB[0]/NemFact[0];
  if(EMSyst[1][0] < 1.0) EMSyst[1][0] = 1/EMSyst[1][0]; EMSyst[1][0] = EMSyst[1][0] - 1.0;

  EMSyst[1][1] = NemFact_FromMLLSB[1]/NemFact[1];
  if(EMSyst[1][1] < 1.0) EMSyst[1][1] = 1/EMSyst[1][1]; EMSyst[1][1] = EMSyst[1][1] - 1.0;

  EMSyst[1][2] = NemFact_FromMLLSB[2]/NemFact[2];
  if(EMSyst[1][2] < 1.0) EMSyst[1][2] = 1/EMSyst[1][2]; EMSyst[1][2] = EMSyst[1][2] - 1.0;

  if(bgdDecay[0][ZHGSEL][0] > 0) EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1/bgdDecay[0][ZHGSEL][0]);
  else                        EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1.0);

  if(bgdDecay[0][ZHGSEL][0] > 0) EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1/bgdDecay[0][ZHGSEL][0]);
  else                        EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1.0);

  if(bgdDecay[0][ZHGSEL][0] > 0) EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1/bgdDecay[0][ZHGSEL][0]);
  else                        EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1.0);

  double EMbkg = bgdDecay[0][ZHGSEL][DY]+bgdDecay[0][ZHGSEL][WZ]+bgdDecay[0][ZHGSEL][ZZ]+bgdDecay[0][ZHGSEL][VVV];

  printf("(mm) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][ZHGSEL+nSelTypes*1][1],sqrt(weiDecay[0][ZHGSEL+nSelTypes*1][4]),
	 bgdDecay[0][ZHGSEL][1]*NemFact[0]  ,sqrt(weiDecay[0][ZHGSEL][4])*NemFact[0],bgdDecay[0][ZHGSEL][0],EMbkg,EMSystTotal[0],EMSyst[0][0],EMSyst[1][0],sqrt(EMSystTotal[0]*EMSystTotal[0] - EMSyst[0][0]*EMSyst[0][0] - EMSyst[1][0]*EMSyst[1][0]));

  printf("(ee) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][ZHGSEL+nSelTypes*2][1],sqrt(weiDecay[0][ZHGSEL+nSelTypes*2][4]),
	 bgdDecay[0][ZHGSEL][1]*NemFact[1]  ,sqrt(weiDecay[0][ZHGSEL][4])*NemFact[1],bgdDecay[0][ZHGSEL][0],EMbkg,EMSystTotal[1],EMSyst[0][1],EMSyst[1][1],sqrt(EMSystTotal[1]*EMSystTotal[1] - EMSyst[0][1]*EMSyst[0][1] - EMSyst[1][1]*EMSyst[1][1]));

  printf("(ll) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][ZHGSEL+nSelTypes*3][1],sqrt(weiDecay[0][ZHGSEL+nSelTypes*3][4]),
	 bgdDecay[0][ZHGSEL][1]*NemFact[2]  ,sqrt(weiDecay[0][ZHGSEL][4])*NemFact[2],bgdDecay[0][ZHGSEL][0],EMbkg,EMSystTotal[2],EMSyst[0][2],EMSyst[1][2],sqrt(EMSystTotal[2]*EMSystTotal[2] - EMSyst[0][2]*EMSyst[0][2] - EMSyst[1][2]*EMSyst[1][2]));

  double systEM[2] = {1.0, 1.0};
  if(histo_EM->GetSumOfWeights() > 1) systEM[0] = 1. + 1./sqrt(histo_EM->GetSumOfWeights());
  else  			      systEM[0] = 2.;
  if(useEMFromData == true && histo_Data->GetSumOfWeights() > 0.0){
    EMbkg = bgdDecay[0][ZHGSEL][DY]+bgdDecay[0][ZHGSEL][WZ]+bgdDecay[0][ZHGSEL][ZZ]+bgdDecay[0][ZHGSEL][VVV];
    double EMNormFact[4] = {((bgdDecay[0][ZHGSEL][DATA]-EMbkg)*NemFact[0])/bgdDecay[0][ZHGSEL+nSelTypes*(1)][EM],
                            ((bgdDecay[0][ZHGSEL][DATA]-EMbkg)*NemFact[1])/bgdDecay[0][ZHGSEL+nSelTypes*(2)][EM],
                            ((bgdDecay[0][ZHGSEL][DATA]-EMbkg)*NemFact[2])/bgdDecay[0][ZHGSEL+nSelTypes*(3)][EM],
                            ((bgdDecay[0][ZHGSEL][DATA]-EMbkg)*1.00000000)/bgdDecay[0][ZHGSEL+nSelTypes*(0)][EM]};

    printf("EM(1): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][ZHGSEL+nSelTypes*1][EM],bgdDecay[0][ZHGSEL][DATA],EMbkg,NemFact[0],bgdDecay[0][ZHGSEL+nSelTypes*(1)][EM],bgdDecay[0][ZHGSEL+nSelTypes*1][EM]*EMNormFact[0],bgdDecay[0][ZHGSEL+nSelTypes*1][EM]*EMNormFact[0]*EMSystTotal[0]);
    printf("EM(2): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][ZHGSEL+nSelTypes*2][EM],bgdDecay[0][ZHGSEL][DATA],EMbkg,NemFact[1],bgdDecay[0][ZHGSEL+nSelTypes*(2)][EM],bgdDecay[0][ZHGSEL+nSelTypes*2][EM]*EMNormFact[1],bgdDecay[0][ZHGSEL+nSelTypes*2][EM]*EMNormFact[1]*EMSystTotal[1]);
    printf("EM(3): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][ZHGSEL+nSelTypes*3][EM],bgdDecay[0][ZHGSEL][DATA],EMbkg,NemFact[2],bgdDecay[0][ZHGSEL+nSelTypes*(3)][EM],bgdDecay[0][ZHGSEL+nSelTypes*3][EM]*EMNormFact[2],bgdDecay[0][ZHGSEL+nSelTypes*3][EM]*EMNormFact[2]*EMSystTotal[2]);
    printf("EM(4): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][ZHGSEL+nSelTypes*0][EM],bgdDecay[0][ZHGSEL][DATA],EMbkg,1.00000000,bgdDecay[0][ZHGSEL+nSelTypes*(0)][EM],bgdDecay[0][ZHGSEL+nSelTypes*0][EM]*EMNormFact[3],bgdDecay[0][ZHGSEL+nSelTypes*0][EM]*EMNormFact[3]*0.000000000000);

    //systEM[0] = 1. + EMSystTotal[typeSel-1];
    //systEM[0] = 1. + EMSyst[1][typeSel-1];
    systEM[0] = 1. + 0.20;
    systEM[1] = TMath::Max(bgdDecay[0][ZHGSEL][DATA],1.0);
    if(bgdDecay[0][ZHGSEL][DATA] > 0) {
      histo_EM->Scale(EMNormFact[typeSel-1]);
      histo_EM->SetBinContent(1,histo_EM->GetBinContent(1)*EMNormFact[3]/EMNormFact[typeSel-1]);
    }
  }
  
  // electron to photon scale factor
  {
  double sig_den = bgdDecay[0][ZLLSEL+nSelTypes*(2)][DY]+bgdDecay[0][ZLLSEL+nSelTypes*(2)][WZ]+bgdDecay[0][ZLLSEL+nSelTypes*(2)][ZZ]+bgdDecay[0][ZLLSEL+nSelTypes*(2)][VVV];
  double bck_den = bgdDecay[0][ZLLSEL+nSelTypes*(2)][EM];
  double dat_den = bgdDecay[0][ZLLSEL+nSelTypes*(2)][DATA];
  double sig_num = bgdDecay[0][ZLGSEL+nSelTypes*(2)][DY]+bgdDecay[0][ZLGSEL+nSelTypes*(2)][WZ]+bgdDecay[0][ZLGSEL+nSelTypes*(2)][ZZ]+bgdDecay[0][ZLGSEL+nSelTypes*(2)][VVV];
  double bck_num = bgdDecay[0][ZLGSEL+nSelTypes*(2)][EM];
  double dat_num = bgdDecay[0][ZLGSEL+nSelTypes*(2)][DATA];
  printf("Electron to photon SF: ((%f-%f)/%f)/((%f-%f)/%f) = %f/%f = %f\n",dat_num,bck_num,sig_num,dat_den,bck_den,sig_den,
                                 (dat_num-bck_num)/sig_num,  (dat_den-bck_den)/sig_den,
				((dat_num-bck_num)/sig_num)/((dat_den-bck_den)/sig_den));
  }

  double mean,up,diff;
  for(int i=1; i<=histo_ZH_hinv[0]->GetNbinsX(); i++){
    mean = histo_Zjets                ->GetBinContent(i);
    up   = histo_Zjets_CMS_ZjetsSystUp->GetBinContent(i);
    diff = mean-up;
    histo_Zjets_CMS_ZjetsSystDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));

    mean = histo_WZ              ->GetBinContent(i);
    up   = histo_WZ_CMS_EWKCorrUp->GetBinContent(i);
    diff = mean-up;
    histo_WZ_CMS_EWKCorrDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));

    mean = histo_ZZ              ->GetBinContent(i);
    up   = histo_ZZ_CMS_EWKCorrUp->GetBinContent(i);
    diff = mean-up;
    histo_ZZ_CMS_EWKCorrDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));

    mean = histo_ZZ             ->GetBinContent(i);
    up   = histo_ZZ_CMS_ggCorrUp->GetBinContent(i);
    diff = mean-up;
    histo_ZZ_CMS_ggCorrDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
  }
  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Computing diboson corrections\n\n");
  
  printf("EWK Corr: WZ(%f/%f/%f) ZZ(%f/%f/%f) ggZZ(%f/%f/%f) ZH(%f/%f/%f)\n",
                     histo_WZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_WZ->GetSumOfWeights(),histo_WZ_CMS_EWKCorrDown->GetSumOfWeights(),
                     histo_ZZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_EWKCorrDown->GetSumOfWeights(),
		     histo_ZZ_CMS_ggCorrUp ->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_ggCorrDown ->GetSumOfWeights(),
		     histo_ZH_hinv_CMS_EWKCorrUp[0]->GetSumOfWeights(),histo_ZH_hinv[0]->GetSumOfWeights(),histo_ZH_hinv_CMS_EWKCorrDown[0]->GetSumOfWeights());

  printf("QCD Corr: WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) ZH(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZH_hinv[0]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][0]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][1]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][2]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][3]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][4]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][5]->GetSumOfWeights());

  histo[allPlots-1][0]->Add(histo_Data);
  histo[allPlots-1][1]->Add(histo_EM);
  histo[allPlots-1][2]->Add(histo_Zjets);
  histo[allPlots-1][3]->Add(histo_WZ);
  histo[allPlots-1][4]->Add(histo_ZZ);
  histo[allPlots-1][5]->Add(histo_VVV);
  histo[allPlots-1][6]->Add(histo_ZH_hinv[plotModel]);
  histo[allPlots-1][7]->Add(histo_ggZH_hinv);
  
  double qcdScaleTotal[2] = {0.035, 0.231};
  double pdfTotal[2] = {0.016, 0.051};
  
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"MitZHAnalysis/plots_zhg%s/histo%szhg%s_nice_%s_%d.root",subFolder.Data(),addChan.Data(),finalStateName,signalName_[plotModel].Data(),thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }

  for(int i=1; i<=histo_ZH_hinv[0]->GetNbinsX(); i++) {
    double factorUp = +1.0; double factorDown = -1.0;
    histo_VVV_CMS_MVAVVVStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorUp  *histo_VVV	->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorDown*histo_VVV	->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	       ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorUp  *histo_WZ	->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorDown*histo_WZ	->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	       ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorUp  *histo_ZZ	->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorDown*histo_ZZ	->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingUp	       ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorUp  *histo_EM	->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorDown*histo_EM       ->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorUp  *histo_ggZH_hinv->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorDown*histo_ggZH_hinv->GetBinError(i),0.000001));

    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	         ->Add(histo_VVV      ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]         ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorUp  *histo_VVV      ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]         ->Add(histo_VVV      ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]       ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorDown*histo_VVV      ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	         ->Add(histo_WZ       ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]           ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorUp  *histo_WZ       ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	         ->Add(histo_WZ       ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]         ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorDown*histo_WZ       ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	         ->Add(histo_ZZ       ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]           ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorUp  *histo_ZZ       ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	         ->Add(histo_ZZ       ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]         ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorDown*histo_ZZ       ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]	         ->Add(histo_EM       ); histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]           ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorUp  *histo_EM       ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]           ->Add(histo_EM       ); histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]         ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorDown*histo_EM       ->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[i-1]    ->Add(histo_ggZH_hinv); histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorUp  *histo_ggZH_hinv->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[i-1]  ->Add(histo_ggZH_hinv); histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorDown*histo_ggZH_hinv->GetBinError(i),0.000001));

    for(int nModel=0; nModel<nSigModels; nModel++) { 
      histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]      ->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorUp  *histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
      histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]    ->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorDown*histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][i-1]  ->Add(histo_ZH_hinv[nModel]);
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][i-1]->Add(histo_ZH_hinv[nModel]);
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][i-1]  ->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorUp  *histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][i-1]->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorDown*histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
    }

  }

  double process_syst[nSigModels][7];
  double yield_processTypes[nSigModels][histBins+1];
  double stat_processTypes[nSigModels][histBins+1];
  double syst_processTypes[nSigModels][histBins+1];
  double syst_types_allBackground[28];
  char outputLimits[200];
  // Output the limits for all the models
  for(int nModel=0; nModel<nSigModels; nModel++) { 
    sprintf(outputLimits,"MitZHAnalysis/plots_zhg%s/zll%szhg%s_%s_input_%s.root",subFolder.Data(),addChan.Data(),finalStateName,signalName_[nModel].Data(),ECMsb.Data());
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
    
    if(verbose) {
      cout << histo_Data             ->GetSumOfWeights() << " ";
      cout << histo_ZH_hinv[nModel]  ->GetSumOfWeights() << " ";
      cout << histo_Zjets            ->GetSumOfWeights() << " ";
      cout << histo_VVV              ->GetSumOfWeights() << " ";
      cout << histo_WZ               ->GetSumOfWeights() << " ";
      cout << histo_ZZ               ->GetSumOfWeights() << " ";
      cout << histo_EM               ->GetSumOfWeights() << " ";
      cout << histo_ggZH_hinv        ->GetSumOfWeights() << " ";
      cout << endl;
      printf("uncertainties Stat\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]	->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]	 ->GetBinContent(i)/histo_ZH_hinv[nModel]  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]	->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]	 ->GetBinContent(i)/histo_ZH_hinv[nModel]  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_VVV	->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp 	 ->GetBinContent(i)/histo_VVV	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_VVV	->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown	 ->GetBinContent(i)/histo_VVV	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp		 ->GetBinContent(i)/histo_WZ       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown 	 ->GetBinContent(i)/histo_WZ	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp		 ->GetBinContent(i)/histo_ZZ       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown 	 ->GetBinContent(i)/histo_ZZ	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_EM	->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingUp		 ->GetBinContent(i)/histo_EM       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_EM	->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingDown 	 ->GetBinContent(i)/histo_EM	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp   ->GetBinContent(i)/histo_ggZH_hinv->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown ->GetBinContent(i)/histo_ggZH_hinv->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffM\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]  ->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffE\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]  ->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties MET\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]      ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]   ->GetBinContent(i)/histo_ZH_hinv[nModel]   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAMETBoundingDown   ->GetBinContent(i)/histo_ggZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties JES\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]    ->GetBinContent(i)/histo_ZH_hinv[nModel]     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]	 ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_ggZH_hinv     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAJESBoundingDown	 ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties BTAG\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]    ->GetBinContent(i)/histo_ZH_hinv[nModel]     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]	 ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVABTAGBoundingUp    ->GetBinContent(i)/histo_ggZH_hinv     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVABTAGBoundingDown	 ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties PU\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_PUBoundingUp[nModel]      ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_PUBoundingDown[nModel]   ->GetBinContent(i)/histo_ZH_hinv[nModel]   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_PUBoundingUp	  ->GetBinContent(i)/histo_ggZH_hinv	  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_PUBoundingDown   ->GetBinContent(i)/histo_ggZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties EWK Corr\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_EWKCorrUp  ->GetBinContent(i)/histo_WZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_EWKCorrDown->GetBinContent(i)/histo_WZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_EWKCorrDown->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ggCorrUp  ->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ggCorrDown->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->GetBinContent(i)/histo_ZH_hinv[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_EWKCorrDown[nModel]->GetBinContent(i)/histo_ZH_hinv[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties Zjets\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystUp  ->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystDown->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
    } 

    histo_Data             ->Write();
    histo_ZH_hinv[nModel]  ->Write();
    histo_Zjets            ->Write();
    histo_VVV              ->Write();
    histo_WZ               ->Write();
    histo_ZZ               ->Write();
    histo_EM               ->Write();
    histo_ggZH_hinv        ->Write();
    
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]   ->Write();
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel] ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingUp	            ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingDown            ->Write();
    histo_WZ_CMS_MVAWZStatBoundingUp	            ->Write();
    histo_WZ_CMS_MVAWZStatBoundingDown	            ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingUp	            ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingDown	            ->Write();
    histo_EM_CMS_MVAEMStatBoundingUp	            ->Write();
    histo_EM_CMS_MVAEMStatBoundingDown	            ->Write();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp       ->Write();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown     ->Write();
    
    histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]  ->Write();
    histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Write();
    histo_VVV_CMS_MVALepEffMBoundingUp	            ->Write();
    histo_VVV_CMS_MVALepEffMBoundingDown	    ->Write();
    histo_WZ_CMS_MVALepEffMBoundingUp	            ->Write();
    histo_WZ_CMS_MVALepEffMBoundingDown	            ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingUp	            ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingDown	            ->Write();
    histo_ggZH_hinv_CMS_MVALepEffMBoundingUp        ->Write();
    histo_ggZH_hinv_CMS_MVALepEffMBoundingDown      ->Write();
    
    histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]  ->Write();
    histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Write();
    histo_VVV_CMS_MVALepEffEBoundingUp	            ->Write();
    histo_VVV_CMS_MVALepEffEBoundingDown	    ->Write();
    histo_WZ_CMS_MVALepEffEBoundingUp	            ->Write();
    histo_WZ_CMS_MVALepEffEBoundingDown	            ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingUp	            ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingDown	            ->Write();
    histo_ggZH_hinv_CMS_MVALepEffEBoundingUp        ->Write();
    histo_ggZH_hinv_CMS_MVALepEffEBoundingDown      ->Write();
    
    histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]      ->Write();
    histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]    ->Write();
    histo_VVV_CMS_MVAMETBoundingUp	            ->Write();
    histo_VVV_CMS_MVAMETBoundingDown	            ->Write();
    histo_WZ_CMS_MVAMETBoundingUp	            ->Write();
    histo_WZ_CMS_MVAMETBoundingDown	            ->Write();
    histo_ZZ_CMS_MVAMETBoundingUp	            ->Write();
    histo_ZZ_CMS_MVAMETBoundingDown	            ->Write();
    histo_ggZH_hinv_CMS_MVAMETBoundingUp            ->Write();
    histo_ggZH_hinv_CMS_MVAMETBoundingDown          ->Write();
    
    histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]	    ->Write();
    histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]    ->Write(); 
    histo_VVV_CMS_MVAJESBoundingUp	            ->Write();
    histo_VVV_CMS_MVAJESBoundingDown	            ->Write();
    histo_WZ_CMS_MVAJESBoundingUp 	            ->Write();
    histo_WZ_CMS_MVAJESBoundingDown	            ->Write();
    histo_ZZ_CMS_MVAJESBoundingUp 	            ->Write();
    histo_ZZ_CMS_MVAJESBoundingDown	            ->Write();
    histo_ggZH_hinv_CMS_MVAJESBoundingUp	    ->Write();
    histo_ggZH_hinv_CMS_MVAJESBoundingDown          ->Write();
    
    histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]	    ->Write();
    histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]   ->Write(); 
    histo_VVV_CMS_MVABTAGBoundingUp	            ->Write();
    histo_VVV_CMS_MVABTAGBoundingDown	            ->Write();
    histo_WZ_CMS_MVABTAGBoundingUp 	            ->Write();
    histo_WZ_CMS_MVABTAGBoundingDown	            ->Write();
    histo_ZZ_CMS_MVABTAGBoundingUp 	            ->Write();
    histo_ZZ_CMS_MVABTAGBoundingDown	            ->Write();
    histo_ggZH_hinv_CMS_MVABTAGBoundingUp	    ->Write();
    histo_ggZH_hinv_CMS_MVABTAGBoundingDown         ->Write();
    
    histo_ZH_hinv_CMS_PUBoundingUp[nModel]          ->Write();
    histo_ZH_hinv_CMS_PUBoundingDown[nModel]        ->Write();
    histo_VVV_CMS_PUBoundingUp	                    ->Write();
    histo_VVV_CMS_PUBoundingDown	            ->Write();
    histo_WZ_CMS_PUBoundingUp	                    ->Write();
    histo_WZ_CMS_PUBoundingDown	                    ->Write();
    histo_ZZ_CMS_PUBoundingUp	                    ->Write();
    histo_ZZ_CMS_PUBoundingDown	                    ->Write();
    histo_ggZH_hinv_CMS_PUBoundingUp                ->Write();
    histo_ggZH_hinv_CMS_PUBoundingDown              ->Write();
    
    histo_WZ_CMS_EWKCorrUp	                    ->Write();
    histo_WZ_CMS_EWKCorrDown	                    ->Write();
    histo_ZZ_CMS_EWKCorrUp	                    ->Write();
    histo_ZZ_CMS_EWKCorrDown	                    ->Write();
    histo_ZZ_CMS_ggCorrUp	                    ->Write();
    histo_ZZ_CMS_ggCorrDown	                    ->Write();
    histo_Zjets_CMS_ZjetsSystUp	                    ->Write();
    histo_Zjets_CMS_ZjetsSystDown	            ->Write();
    outFileLimits->Close();
    double lumiE = 1.025;
    double systLepResE[5] = {1.01,1.01,1.01,1.01,1.01};
    double systLepResM[5] = {1.01,1.01,1.01,1.01,1.01};
    double syst_btag = 1.02;
    double syst_WZl[2] = {1.010, 1.003};
    if(nJetsType > 0) syst_WZl[1] = 1.012;
    
    
    for(int nb=1; nb<=nBinMVA; nb++){
      double nggZHEvt = histo_ggZH_hinv->GetBinContent(nb);
      if(nModel != 0) nggZHEvt = 0.0;
      // QCD study
      double systQCDScale[5] = {TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nModel][0]->GetBinContent(nb)-histo_ZH_hinv[nModel]->GetBinContent(nb)),
                                TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]            ->GetBinContent(nb)-histo_VVV            ->GetBinContent(nb)),
                                TMath::Abs(histo_WZ_CMS_QCDScaleBounding[0]             ->GetBinContent(nb)-histo_WZ             ->GetBinContent(nb)),
                                TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0]             ->GetBinContent(nb)-histo_ZZ             ->GetBinContent(nb)),
  			      TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[0]        ->GetBinContent(nb)-histo_ggZH_hinv      ->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nqcd]->GetBinContent(nb)  -histo_ZH_hinv[nModel]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nqcd]->GetBinContent(nb)  -histo_ZH_hinv[nModel]  ->GetBinContent(nb));
        if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)              -histo_VVV    ->GetBinContent(nb))         > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)      -histo_VVV      ->GetBinContent(nb));
        if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)               -histo_WZ     ->GetBinContent(nb))         > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_WZ       ->GetBinContent(nb));
        if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)               -histo_ZZ     ->GetBinContent(nb))         > systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_ZZ       ->GetBinContent(nb));
        if(TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)        -histo_ggZH_hinv->GetBinContent(nb))       > systQCDScale[4]) systQCDScale[4] = TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggZH_hinv->GetBinContent(nb));
      }                 
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_ZH_hinv[nModel]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histo_VVV->GetBinContent(nb)             > 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_VVV            ->GetBinContent(nb); else systQCDScale[1] = 1;
      if(histo_WZ->GetBinContent(nb)              > 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_WZ             ->GetBinContent(nb); else systQCDScale[2] = 1;
      if(histo_ZZ->GetBinContent(nb)              > 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_ZZ             ->GetBinContent(nb); else systQCDScale[3] = 1;
      if(histo_ggZH_hinv->GetBinContent(nb)       > 0) systQCDScale[4] = 1 + systQCDScale[4]/histo_ggZH_hinv      ->GetBinContent(nb); else systQCDScale[4] = 1;
      for(int ntype=0; ntype<5; ntype++) if(systQCDScale[ntype] < 0) systQCDScale[ntype] = 1.0;
      if(verbose) printf("QCDScale(%d): %f %f %f %f %f - ",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[4]);

      if(systQCDScale[0]-1 > qcdScaleTotal[0]) systQCDScale[0] = 1 + sqrt((systQCDScale[0]-1.0)*(systQCDScale[0]-1.0)-qcdScaleTotal[0]*qcdScaleTotal[0]); else systQCDScale[0] = 1.0;
      if(systQCDScale[4]-1 > qcdScaleTotal[1]) systQCDScale[4] = 1 + sqrt((systQCDScale[4]-1.0)*(systQCDScale[4]-1.0)-qcdScaleTotal[1]*qcdScaleTotal[1]); else systQCDScale[4] = 1.0;
      if(verbose) printf("%f/%f %f/%f\n",systQCDScale[0],qcdScaleTotal[0],systQCDScale[4],qcdScaleTotal[1]);

      // PDF study
      double systPDF[5];
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_ZH_hinv_CMS_PDFBounding[nModel][npdf]->GetBinContent(nb)-histo_ZH_hinv[nModel]->GetBinContent(nb))/histo_ZH_hinv[nModel]->GetBinContent(nb));
      systPDF[0] = 1.0+sqrt(TMath::Max(histo_Diff->GetRMS()*histo_Diff->GetRMS()-pdfTotal[0]*pdfTotal[0],0.0));
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
      systPDF[1] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_WZ_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WZ->GetBinContent(nb))/histo_WZ->GetBinContent(nb));
      systPDF[2] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ZZ->GetBinContent(nb))/histo_ZZ->GetBinContent(nb));
      systPDF[3] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) {
        double aux=0;
        if(histo_ggZH_hinv->GetBinContent(nb) > 0) histo_Diff->Fill((histo_ggZH_hinv_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ggZH_hinv->GetBinContent(nb))/histo_ggZH_hinv->GetBinContent(nb));
      }
      systPDF[4] = 1.0+sqrt(TMath::Max(histo_Diff->GetRMS()*histo_Diff->GetRMS()-pdfTotal[1]*pdfTotal[1],0.0));
      if(verbose) printf("PDF(%d): %f %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4]); // ZH, VVV, WZ, ZZ, ggZH
  
      double systLepEffM[5] = {1.0,1.0,1.0,1.0,1.0};
      if     (histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]     ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb);
      else if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->GetBinContent(nb);
      if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingUp         ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffMBoundingUp          ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffMBoundingDown        ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp          ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown        ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[4] = histo_ggZH_hinv_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[4] = histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
  
      double systLepEffE[5] = {1.0,1.0,1.0,1.0,1.0};
      if     (histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]     ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb);
      else if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->GetBinContent(nb);
      if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffEBoundingUp         ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
      else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffEBoundingUp          ->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
      else if(histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffEBoundingDown        ->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffEBoundingUp          ->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
      else if(histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffEBoundingDown        ->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb) > 0) systLepEffE[4] = histo_ggZH_hinv_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
      else if(histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[4] = histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
  
      double systMetUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systMetDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]     ->GetBinContent(nb) > 0) systMetUp  [0] = histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]   ->GetBinContent(nb) > 0) systMetDown[0] = histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_MVAMETBoundingUp         ->GetBinContent(nb) > 0) systMetUp  [1] = histo_VVV_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMetDown[1] = histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)        > 0 && histo_WZ_CMS_MVAMETBoundingUp          ->GetBinContent(nb) > 0) systMetUp  [2] = histo_WZ_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)        > 0 && histo_WZ_CMS_MVAMETBoundingDown        ->GetBinContent(nb) > 0) systMetDown[2] = histo_WZ_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_MVAMETBoundingUp          ->GetBinContent(nb) > 0) systMetUp  [3] = histo_ZZ_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_MVAMETBoundingDown        ->GetBinContent(nb) > 0) systMetDown[3] = histo_ZZ_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMetUp  [4] = histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMetDown[4] = histo_ggZH_hinv_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      for(int i=0; i<5; i++) if(systMetUp  [i] == 1) systMetUp  [i] = 0.998;
      for(int i=0; i<5; i++) if(systMetDown[i] == 1) systMetDown[i] = 1.002;
      for(int nmet=0; nmet<5; nmet++) if(systMetUp[nmet]   > 1.04) systMetUp[nmet]   = 1.04;
      for(int nmet=0; nmet<5; nmet++) if(systMetUp[nmet]   < 0.96) systMetUp[nmet]   = 0.96;
      for(int nmet=0; nmet<5; nmet++) if(systMetDown[nmet] > 1.04) systMetDown[nmet] = 1.04;
      for(int nmet=0; nmet<5; nmet++) if(systMetDown[nmet] < 0.96) systMetDown[nmet] = 0.96;
  
      double systJesUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systJesDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]	->GetBinContent(nb) > 0) systJesUp  [0] = histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]	->GetBinContent(nb) > 0) systJesDown[0] = histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVAJESBoundingUp 	->GetBinContent(nb) > 0) systJesUp  [1] = histo_VVV_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[1] = histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_MVAJESBoundingUp  	->GetBinContent(nb) > 0) systJesUp  [2] = histo_WZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[2] = histo_WZ_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVAJESBoundingUp  	->GetBinContent(nb) > 0) systJesUp  [3] = histo_ZZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[3] = histo_ZZ_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAJESBoundingUp	->GetBinContent(nb) > 0) systJesUp  [4] = histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJesDown[4] = histo_ggZH_hinv_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      for(int njes=0; njes<5; njes++) if(systJesUp[njes]   > 1.10) systJesUp[njes]   = 1.10;
      for(int njes=0; njes<5; njes++) if(systJesUp[njes]   < 0.90) systJesUp[njes]   = 0.90;
      for(int njes=0; njes<5; njes++) if(systJesDown[njes] > 1.10) systJesDown[njes] = 1.10;
      for(int njes=0; njes<5; njes++) if(systJesDown[njes] < 0.90) systJesDown[njes] = 0.90;

      double systBtagUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBtagDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]	->GetBinContent(nb) > 0) systBtagUp  [0] = histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]	->GetBinContent(nb) > 0) systBtagDown[0] = histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVABTAGBoundingUp 	   ->GetBinContent(nb) > 0) systBtagUp  [1] = histo_VVV_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVABTAGBoundingDown	   ->GetBinContent(nb) > 0) systBtagDown[1] = histo_VVV_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_MVABTAGBoundingUp  	   ->GetBinContent(nb) > 0) systBtagUp  [2] = histo_WZ_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_MVABTAGBoundingDown	   ->GetBinContent(nb) > 0) systBtagDown[2] = histo_WZ_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVABTAGBoundingUp  	   ->GetBinContent(nb) > 0) systBtagUp  [3] = histo_ZZ_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVABTAGBoundingDown	   ->GetBinContent(nb) > 0) systBtagDown[3] = histo_ZZ_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVABTAGBoundingUp   ->GetBinContent(nb) > 0) systBtagUp  [4] = histo_ggZH_hinv_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVABTAGBoundingDown ->GetBinContent(nb) > 0) systBtagDown[4] = histo_ggZH_hinv_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      for(int njes=0; njes<5; njes++) if(systBtagUp[njes]   > 1.10) systBtagUp[njes]   = 1.10;
      for(int njes=0; njes<5; njes++) if(systBtagUp[njes]   < 0.90) systBtagUp[njes]   = 0.90;
      for(int njes=0; njes<5; njes++) if(systBtagDown[njes] > 1.10) systBtagDown[njes] = 1.10;
      for(int njes=0; njes<5; njes++) if(systBtagDown[njes] < 0.90) systBtagDown[njes] = 0.90;
   
      double systPUUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systPUDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_PUBoundingUp[nModel]     ->GetBinContent(nb) > 0) systPUUp  [0] = histo_ZH_hinv_CMS_PUBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_PUBoundingDown[nModel]   ->GetBinContent(nb) > 0) systPUDown[0] = histo_ZH_hinv_CMS_PUBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_PUBoundingUp         ->GetBinContent(nb) > 0) systPUUp  [1] = histo_VVV_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_PUBoundingDown       ->GetBinContent(nb) > 0) systPUDown[1] = histo_VVV_CMS_PUBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)        > 0 && histo_WZ_CMS_PUBoundingUp          ->GetBinContent(nb) > 0) systPUUp  [2] = histo_WZ_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)        > 0 && histo_WZ_CMS_PUBoundingDown        ->GetBinContent(nb) > 0) systPUDown[2] = histo_WZ_CMS_PUBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_PUBoundingUp          ->GetBinContent(nb) > 0) systPUUp  [3] = histo_ZZ_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_PUBoundingDown        ->GetBinContent(nb) > 0) systPUDown[3] = histo_ZZ_CMS_PUBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPUUp  [4] = histo_ggZH_hinv_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPUDown[4] = histo_ggZH_hinv_CMS_PUBoundingDown->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      for(int npu=0; npu<5; npu++) if(systPUUp[npu]   > 1.02) systPUUp[npu]   = 1.02;
      for(int npu=0; npu<5; npu++) if(systPUUp[npu]   < 0.98) systPUUp[npu]   = 0.98;
      for(int npu=0; npu<5; npu++) if(systPUDown[npu] > 1.02) systPUDown[npu] = 1.02;
      for(int npu=0; npu<5; npu++) if(systPUDown[npu] < 0.98) systPUDown[npu] = 0.98;
   
      double systZjetsUp  [1] = {1.0};
      double systZjetsDown[1] = {1.0};
      if(histo_Zjets->GetBinContent(nb) > 0 && histo_Zjets_CMS_ZjetsSystUp   ->GetBinContent(nb) > 0) systZjetsUp  [0] = histo_Zjets_CMS_ZjetsSystUp  ->GetBinContent(nb)/histo_Zjets->GetBinContent(nb);
      if(histo_Zjets->GetBinContent(nb) > 0 && histo_Zjets_CMS_ZjetsSystDown ->GetBinContent(nb) > 0) systZjetsDown[0] = histo_Zjets_CMS_ZjetsSystDown->GetBinContent(nb)/histo_Zjets->GetBinContent(nb);
  
      double syst_EWKCorrUp[4]   = {1.0,1.0,1.0,1.0}; // WZ, ZZ, ggZZ, ZH
      double syst_EWKCorrDown[4] = {1.0,1.0,1.0,1.0};
      if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [0] = histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[0] = histo_WZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [1] = histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[1] = histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb) > 0) syst_EWKCorrUp  [2] = histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb) > 0) syst_EWKCorrDown[2] = histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0 &&  histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [3] = histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0 &&  histo_ZH_hinv_CMS_EWKCorrDown[nModel]->GetBinContent(nb) > 0) syst_EWKCorrDown[3] = histo_ZH_hinv_CMS_EWKCorrDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);

      char outputLimitsShape[200];                                            
      sprintf(outputLimitsShape,"MitZHAnalysis/datacards_zhg%s/histo_limits_zll%szhg%s_%s_shape_%s_bin%d.txt",subFolder.Data(),addChan.Data(),finalStateName,signalName_[nModel].Data(),ECMsb.Data(),nb-1);
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
      newcardShape << Form("bin zhg%2s%4s%d zhg%2s%4s%d zhg%2s%4s%d zhg%2s%4s%d zhg%2s%4s%d zhg%2s%4s%d zhg%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process ZH_hinv Zjets VVV WZ ZZ EM ggZH_hinv\n");
      newcardShape << Form("process 0 1 2 3 4 5 -1\n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",histo_ZH_hinv[nModel]->GetBinContent(nb),histo_Zjets->GetBinContent(nb),TMath::Max(histo_VVV->GetBinContent(nb),0.0),histo_WZ->GetBinContent(nb),histo_ZZ->GetBinContent(nb),histo_EM->GetBinContent(nb),nggZHEvt);
      newcardShape << Form("lumi_%4s                               lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE);		     
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3],systLepEffM[4]);
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3],systLepEffE[4]);
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3],systLepResM[4]);
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3],systLepResE[4]);
      newcardShape << Form("CMS_pu2016                             lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systPUUp[0],systPUDown[0],systPUUp[1],systPUDown[1],systPUUp[2],systPUDown[2],systPUUp[3],systPUDown[3],systPUUp[0],systPUDown[0]); // 0 --> 4
      newcardShape << Form("CMS_scale_met                          lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systMetUp[0],systMetDown[0],systMetUp[1],systMetDown[1],systMetUp[2],systMetDown[2],systMetUp[3],systMetDown[3],systMetUp[0],systMetDown[0]); // 0 --> 4
      newcardShape << Form("CMS_scale_j                            lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systJesUp[0],systJesDown[0],systJesUp[1],systJesDown[1],systJesUp[2],systJesDown[2],systJesUp[3],systJesDown[3],systJesUp[0],systJesDown[0]); // 0 --> 4		 
      newcardShape << Form("CMS_trigger2016                        lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",1.01,1.01,1.01,1.01,1.01);
      newcardShape << Form("CMS_eff_photon                         lnN  %7.5f %7.5f   -     -     -     -   %7.5f\n",1.02,1.02,1.02);
      newcardShape << Form("CMS_fake_el                            lnN    -     -   %7.5f %7.5f %7.5f   -     -\n",1.10,1.10,1.10);
      newcardShape << Form("UEPS			           lnN  1.030   -     -     -     -     -   1.030\n");
      newcardShape << Form("CMS_eff_b_2016                         lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systBtagUp[0],systBtagDown[0],systBtagUp[1],systBtagDown[1],systBtagUp[2],systBtagDown[2],systBtagUp[3],systBtagDown[3],systBtagUp[0],systBtagDown[0]); // 0 --> 4
      newcardShape << Form("pdf_qqbar_ACCEPT                       lnN  %7.5f   -   %7.5f %7.5f %7.5f   -     -  \n",TMath::Max(systPDF[0],1.01),TMath::Max(systPDF[1],1.01),TMath::Max(systPDF[2],1.01),TMath::Max(systPDF[3],1.01));
      if(systPDF[4] != 1.0)
      newcardShape << Form("pdf_gg_ACCEPT                          lnN    -     -     -     -     -     -   %7.5f\n",systPDF[4]);
      newcardShape << Form("pdf_qqbar                              lnN  %7.5f   -     -     -     -     -     -  \n",1.0+pdfTotal[0]);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0)
      newcardShape << Form("pdf_gg                                 lnN    -     -     -     -     -     -   %7.5f\n",1.0+pdfTotal[1]);
      if(systQCDScale[0] != 1.0)
      newcardShape << Form("QCDscale_VH_ACCEPT		         lnN  %7.5f   -     -     -     -     -     -  \n",systQCDScale[0]);  
      if(systQCDScale[4] != 1.0 && histo_ggZH_hinv->GetBinContent(nb) > 0)
      newcardShape << Form("QCDscale_ggVH_ACCEPT		         lnN    -     -     -     -     -     -   %7.5f\n",systQCDScale[4]);  
      newcardShape << Form("QCDscale_VH		                 lnN  %7.5f   -     -     -     -     -     -  \n",1.0+qcdScaleTotal[0]);  
      if(histo_ggZH_hinv->GetBinContent(nb) > 0)
      newcardShape << Form("QCDscale_ggVH		                 lnN    -     -     -     -     -     -   %7.5f\n",1.0+qcdScaleTotal[1]);  
      newcardShape << Form("QCDscale_VVV		                 lnN    -     -   %7.5f   -     -     -     -  \n",systQCDScale[1]);		
      newcardShape << Form("QCDscale_WZ		                   lnN    -     -     -   %7.5f      -    -      -  \n",systQCDScale[2]);
      newcardShape << Form("QCDscale_ZZ		                   lnN    -     -     -   -      %7.5f    -      -  \n",systQCDScale[3]);
      newcardShape << Form("CMS_zllhinv_ZH_EWKCorr                 lnN    %7.5f/%7.5f      -     -   -  -      -      -  \n",syst_EWKCorrUp[3],syst_EWKCorrDown[3]);		

      newcardShape << Form("CMS_zllhinv_WZ_EWKCorr                 lnN    -     -     -   %7.5f/%7.5f   -      -      -  \n",syst_EWKCorrUp[0],syst_EWKCorrDown[0]);		
      newcardShape << Form("CMS_zllhinv_ZZ_EWKCorr                 lnN    -     -     -     -    %7.5f/%7.5f   -      -  \n",syst_EWKCorrUp[1],syst_EWKCorrDown[1]);		

      newcardShape << Form("CMS_zllhinv_ggZZCorr                   lnN    -     -     -     -   %7.5f/%7.5f   -     -  \n",syst_EWKCorrUp[2],syst_EWKCorrDown[2]);		

      newcardShape << Form("CMS_zllhinv_WZ_lep2016                 lnN     -	 -     -   %7.5f   -	  -    -  \n",syst_WZl[0]);	    
      newcardShape << Form("CMS_zllhinv_WZ_tau2016                 lnN     -	 -     -   %7.5f   -	  -    -  \n",syst_WZl[1]);	    
      newcardShape << Form("CMS_zllhinv_ZLLNorm2016_%s_%s              lnN	-   %7.5f   -	  -     -     -     -  \n",finalStateName,ECMsb.Data(),2.0);	    
      newcardShape << Form("CMS_zllhinv_EMSyst2016_%s_%s               lnN	-     -     -	  -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),systEM[0]);	       

      //newcardShape << Form("CMS_zllhinv_EMNorm2016_%s_%s               lnU	-     -     -	  -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),2.0);      
  
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding2016_%s_Bin%d    lnN    %7.5f -      -	 -    -    -	 -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ZH_hinv[nModel]->GetBinError(nb)  /histo_ZH_hinv[nModel]->GetBinContent(nb),0.999));
  
      if(histo_Zjets->GetBinContent(nb)           > 0) newcardShape << Form("CMS_zllhinv%s_MVAZjetsStatBounding2016_%s_Bin%d  lnN	-  %7.5f   -	-    -    -	-  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_Zjets->GetBinError(nb)    /histo_Zjets->GetBinContent(nb),0.999));
  
      if(histo_VVV->GetBinContent(nb)             > 0) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding2016_%s_Bin%d    lnN	-   -  %7.5f   -    -	 -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_VVV->GetBinError(nb)	   /histo_VVV->GetBinContent(nb),0.999));
  
      if(histo_WZ->GetBinContent(nb)              > 0) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding2016_%s_Bin%d     lnN      -  -     -  %7.5f  -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WZ->GetBinError(nb)	/histo_WZ->GetBinContent(nb),0.999));
  
      if(histo_ZZ->GetBinContent(nb)	          > 0) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding2016_%s_Bin%d     lnN	  -  -     -	-  %7.5f  -	-  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ZZ->GetBinError(nb)	    /histo_ZZ->GetBinContent(nb),0.999));
  
      if(histo_EM->GetBinContent(nb)              > 0) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding2016_%s_Bin%d     lnN	  -  -     -	-    -  %7.5f	-  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_EM->GetBinError(nb)	    /histo_EM->GetBinContent(nb),0.999));
  
      if(histo_ggZH_hinv->GetBinContent(nb)       > 0) newcardShape << Form("CMS_zllhinv%s_MVAggZHStatBounding2016_%s_Bin%d   lnN      -    -     -    -    -    -   %7.5f\n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ggZH_hinv->GetBinError(nb)/histo_ggZH_hinv->GetBinContent(nb),0.999));
    }

  } // end loop over models
  
}
