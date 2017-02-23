#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"

#include "MitAnalysisRunII/macros/80x/factors.h"
#include "MitAnalysisRunII/macros/80x/BTagCalibrationStandalone.cc"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"
#include "MitZHAnalysis/macros/80x/zhMVA.h"

enum selType                    { WZSEL,  WHSEL,   WZLOOSESEL, nSelTypes};
TString selTypeName[nSelTypes]= {"WZSEL","WHSEL", "WZLOOSESEL"};

enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};

bool isMINIAOD = true;
int whichSkim = 2;
double mcPrescale = 1.0;
bool usePureMC = false;
bool isWZhinv = true;
bool useZZWZEWKUnc = true;
const double bTagCuts[1] = {0.9535}; // 0.5426/0.8484/0.9535 (check BTagCalibration2Reader!)

double overallZgSF[1] = {0.8};

void wzAnalysis(
 bool applyWHSel = false,
 double minMass =  76,
 double maxMass = 106,
 bool applyBtagging = true,
 TString typeLepSel = "medium",
 TString type3rdLepSel = "default",
 string the_BDT_weights="",
 string subdirectory="",
 bool isMIT =false
 ){
  if(subdirectory!="" && subdirectory.c_str()[0]!='/') subdirectory = "/"+subdirectory;
  system(("mkdir -p MitZHAnalysis/datacards"+subdirectory).c_str());
  system(("mkdir -p MitZHAnalysis/plots"+subdirectory).c_str());

  bool printMCEventList=false;
  bool useBDT=false;
  Int_t period = 1;

  // File instances on EOS
  TString filesPathDA = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/met_";
  TString filesPathMC  = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/met_";
  TString filesPathMC2 = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/met_";
  // File instances on T3 hadoop
  if(isMIT){
    filesPathDA   = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/data/met_";
    filesPathMC	  = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/met_";
    filesPathMC2  = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/met_";
  }
  Double_t lumi = 35.9;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  TString triggerSuffix = "*";
  if(isMINIAOD) triggerSuffix = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";

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

  if(usePureMC == true){
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV.root",filesPathMC.Data()));					      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));        infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg.root",filesPathMC2.Data()));					      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root",filesPathMC.Data()));   	      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8.root",filesPathMC.Data()));      	      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));		              infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8.root",filesPathMC.Data()));  				      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(1);
  }
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(2);

  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));			   infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8.root",filesPathMC.Data()));					   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		   infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                        infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                        infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                        infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));        infilecatv.push_back(5);
  infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                               infilecatv.push_back(5);

  if(applyWHSel){
  infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data())); 		    infilecatv.push_back(6);
  infilenamev.push_back(Form("%sttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.root",filesPathMC.Data()));   infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sHWplusJ_HToWW_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                         infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sHWminusJ_HToWW_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                        infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sHZJ_HToWW_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                             infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sWplusHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                        infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sWminusHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                       infilecatv.push_back(6);
  }
  }
  else {assert(0);}

  if(infilenamev.size() != infilecatv.size()) {printf("Failure %d != %d\n",(int)infilenamev.size(),(int)infilecatv.size()); assert(0); return;}

  //infilenamev.clear();infilecatv.clear();
  //infilenamev.push_back(Form("%sWZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8.root",filesPathMC.Data()));  infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));    infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZJJ_EWK_QCD_13TeV-madgraph-pythia8.root",filesPathMC.Data()));	       infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sWZJJ_EWK_13TeV-madgraph-pythia8.root",filesPathMC.Data()));	       infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sWZJJ_QCD_13TeV-madgraph-pythia8.root",filesPathMC.Data()));	       infilecatv.push_back(5);

  double denBTagging[5][5][3],jetEpsBtagTIGHT[5][5][3];
  double numBTaggingTIGHT[5][5][3];
  for(int i0=0; i0<5; i0++) {
    for(int i1=0; i1<5; i1++) {
      for(int i2=0; i2<3; i2++) {
        denBTagging[i0][i1][i2] = 0.0;
        numBTaggingTIGHT[i0][i1][i2] = 0.0;
	if     (i2==BTagEntry::FLAV_B)    jetEpsBtagTIGHT[i0][i1][i2] = jetEpsBtagBTIGHT[i0][i1];
	else if(i2==BTagEntry::FLAV_C)    jetEpsBtagTIGHT[i0][i1][i2] = jetEpsBtagCTIGHT[i0][i1];
	else if(i2==BTagEntry::FLAV_UDSG) jetEpsBtagTIGHT[i0][i1][i2] = jetEpsBtagLTIGHT[i0][i1];
      }
    }
  }

  //Float_t fMVACut[4][4];
  //InitializeJetIdCuts(fMVACut);
  
  BTagCalibration2 *btagCalib = new BTagCalibration2("csvv2","MitAnalysisRunII/data/80x/CSVv2_Moriond17_B_H.csv");
  BTagCalibration2Reader btagReaderBCTIGHT(btagCalib,BTagEntry::OP_TIGHT,"comb","central");
  BTagCalibration2Reader btagReaderLTIGHT(btagCalib,BTagEntry::OP_TIGHT,"incl","central");
  BTagCalibration2Reader btagReaderBCTIGHTUP(btagCalib,BTagEntry::OP_TIGHT,"comb","up");
  BTagCalibration2Reader btagReaderLTIGHTUP(btagCalib,BTagEntry::OP_TIGHT,"incl","up");
  BTagCalibration2Reader btagReaderBCTIGHTDOWN(btagCalib,BTagEntry::OP_TIGHT,"comb","down");
  BTagCalibration2Reader btagReaderLTIGHTDOWN(btagCalib,BTagEntry::OP_TIGHT,"incl","down");

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

  const int allStates = 6;
  double totalFakeDataCount[allStates][9];
  for(int i=0; i<allStates; i++) for(int j=0; j<9; j++) totalFakeDataCount[i][j] = 0;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 45;
  const int histBins = 7;
  TH1D* histo[allStates][allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "Zgamma", "....WZ", "...ZZ", "...VVV", ".Higgs"};

  const int numberCuts = 7;
  TH1D* histoWZSEL[6];
  histoWZSEL[0] = new TH1D("histoWZSEL_0", "histoWZSEL_0", numberCuts+1, -0.5, numberCuts+0.5);
  histoWZSEL[1] = new TH1D("histoWZSEL_1", "histoWZSEL_1", numberCuts+1, -0.5, numberCuts+0.5);
  histoWZSEL[2] = new TH1D("histoWZSEL_2", "histoWZSEL_2", numberCuts+1, -0.5, numberCuts+0.5);
  histoWZSEL[3] = new TH1D("histoWZSEL_3", "histoWZSEL_3", numberCuts+1, -0.5, numberCuts+0.5);
  histoWZSEL[4] = new TH1D("histoWZSEL_4", "histoWZSEL_4", numberCuts+1, -0.5, numberCuts+0.5);
  histoWZSEL[5] = new TH1D("histoWZSEL_5", "histoWZSEL_5", numberCuts+1, -0.5, numberCuts+0.5);
  TString cutName[numberCuts+1] = {"presel","minmass","met","massz","ptl3","mass3l","btag","add"};

  double bgdDecay[nSelTypes*allStates][histBins],weiDecay[nSelTypes*allStates][histBins];
  for(unsigned int i=0; i<nSelTypes*allStates; i++) {
    for(int j=0; j<histBins; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }

  TString ECMsb  = "13TeV2016";
  //const int MVAVarType = 1; const int nBinMVA = 8; Float_t xbins[nBinMVA+1] = {0, 1, 2, 3, 4, 5, 100, 200, 300};
  //if(isWZhinv) {
  //  xbins[0] = 0;   xbins[1] = 50;  xbins[2] = 100; xbins[3] = 125; xbins[4] = 150; 
  //  xbins[5] = 175; xbins[6] = 200; xbins[7] = 250; xbins[8] = 350;
  //}
  const int MVAVarType = 1; const int nBinMVA = 12; Float_t xbins[nBinMVA+1] = {0, 1, 2, 3, 4, 5, 100, 200, 300, 400, 500, 600, 700};
  if(isWZhinv) {
    xbins[0] = 0;   xbins[1] = 50;  xbins[2] = 100; xbins[3] = 125; xbins[4] = 150; 
    xbins[5] = 175; xbins[6] = 200; xbins[7] = 250; xbins[8] = 300; xbins[9] = 350;
    xbins[10]= 400; xbins[11]= 500; xbins[12]= 600;
  }
  //const int MVAVarType = 1; const int nBinMVA = 11; Float_t xbins[nBinMVA+1] = {0, 1, 2, 3, 4, 5, 100, 200, 300, 400, 500, 600};
  //if(isWZhinv) {
  //  xbins[0] = 0;   xbins[1] = 50;  xbins[2] = 100; xbins[3] = 160; xbins[4] = 240; 
  //  xbins[5] = 320; xbins[6] = 400; xbins[7] = 480; xbins[8] = 560; xbins[9] = 640;
  //  xbins[10]= 800; xbins[11]=1200;
  //}
  //const int MVAVarType = 3; const int nBinMVA = 18; Float_t xbins[nBinMVA+1] =  {-2, -1, 0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4}; TString addChan = "3";
  //const int MVAVarType = 3; const int nBinMVA = 11; Float_t xbins[nBinMVA+1] =  {-2, -1, .5, .6, .7, .74, .78, .80, .82, .84, .86, 1.}; TString addChan = "3";
  //const int MVAVarType = 1; const int nBinMVA = 21; Float_t xbins[nBinMVA+1] = {0, 1, 2, 3, 4, 5, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
  //if(isWZhinv) {
  //  xbins[0] = 0;   xbins[1] = 50;  xbins[2] = 100; xbins[3] = 150; xbins[4] = 200; 
  //  xbins[5] = 250; xbins[6] = 300; xbins[7] = 350; xbins[8] = 400; xbins[9] = 450;
  //  xbins[10]= 500; xbins[11]= 550; xbins[12]= 600; xbins[13]= 650; xbins[14]= 700;
  //  xbins[15]= 750; xbins[16]= 800; xbins[17]= 850; xbins[18]= 900; xbins[19]= 950;
  //  xbins[20]=1000; xbins[21]=1200;
  //}
  if(MVAVarType==3 || MVAVarType==4) useBDT=true;
  TMVA::Reader *reader; // =new TMVA::Reader();
  Float_t  mva_balance,
           mva_cos_theta_star_l1,
           mva_cos_theta_CS_l1,
           mva_delphi_ptll_MET,
           mva_delphi_ll,
           mva_delphi_jet_MET,
           mva_deltaR_ll,
           mva_etall,
           mva_etal1,
           mva_etal2,
           mva_MET,
           mva_mll_minus_mZ,
           mva_mTjetMET,
           mva_mTll,
           mva_mTl1MET,
           mva_mTl2MET,
           mva_ptll,
           mva_ptl1,
           mva_ptl2,
           mva_ptl1mptl2_over_ptll,
           mva_response,
           mva_weight;
  UChar_t  mva_njets,
           mva_ntaus;
  Bool_t   mva_btag_veto,
           mva_3lveto;
  if(useBDT) {
    reader=new TMVA::Reader();
    if(MVAVarType==3) {
      reader->AddVariable( "mva_balance"                       , &mva_balance            );
      reader->AddVariable( "mva_cos_theta_star_l1"             , &mva_cos_theta_star_l1  );
      reader->AddVariable( "TMath::Abs(mva_cos_theta_CS_l1)"   , &mva_cos_theta_CS_l1    );
      reader->AddVariable( "mva_delphi_ptll_MET"               , &mva_delphi_ptll_MET    );
      reader->AddVariable( "mva_delphi_ll"                     , &mva_delphi_ll          );
      reader->AddVariable( "mva_deltaR_ll"                     , &mva_deltaR_ll          );
      reader->AddVariable( "TMath::Abs(mva_etall)"             , &mva_etall              );
      reader->AddVariable( "TMath::Abs(mva_etal1)"             , &mva_etal1              );
      reader->AddVariable( "TMath::Abs(mva_etal2)"             , &mva_etal2              );
      reader->AddVariable( "mva_MET"                           , &mva_MET                );
      reader->AddVariable( "mva_mll_minus_mZ"                  , &mva_mll_minus_mZ       );
      reader->AddVariable( "mva_mTll"                          , &mva_mTll               );
      reader->AddVariable( "mva_mTl1MET"                       , &mva_mTl1MET            );
      reader->AddVariable( "mva_mTl2MET"                       , &mva_mTl2MET            );
      reader->AddVariable( "mva_ptll"                          , &mva_ptll               );
      reader->AddVariable( "mva_ptl1"                          , &mva_ptl1               );
      reader->AddVariable( "mva_ptl2"                          , &mva_ptl2               );
      reader->AddVariable( "mva_ptl1mptl2_over_ptll"           , &mva_ptl1mptl2_over_ptll);
    } else if(MVAVarType==4) {
      reader->AddVariable( "TMath::Abs(mva_cos_theta_CS_l1)"   , &mva_cos_theta_CS_l1    );
      reader->AddVariable( "mva_deltaR_ll"                     , &mva_deltaR_ll          );
      reader->AddVariable( "TMath::Abs(mva_etall)"             , &mva_etall              );
      reader->AddVariable( "TMath::Abs(mva_etal1)"             , &mva_etal1              );
      reader->AddVariable( "TMath::Abs(mva_etal2)"             , &mva_etal2              );
      reader->AddVariable( "mva_mll_minus_mZ"                  , &mva_mll_minus_mZ       );
      reader->AddVariable( "mva_ptll"                          , &mva_ptll               );
      reader->AddVariable( "mva_ptl1"                          , &mva_ptl1               );
      reader->AddVariable( "mva_ptl2"                          , &mva_ptl2               );
      reader->AddVariable( "mva_ptl1mptl2_over_ptll"           , &mva_ptl1mptl2_over_ptll);
    }
    reader->BookMVA("BDT", the_BDT_weights);
  }

  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  2) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 150.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   4.0;}
    else if(thePlot >=  8 && thePlot <= 11) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 12 && thePlot <= 12) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 13 && thePlot <= 13) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >= 14 && thePlot <= 14) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >= 15 && thePlot <= 15) {nBinPlot =   5; xminPlot =-0.5; xmaxPlot =   4.5;}
    else if(thePlot >= 16 && thePlot <= 17) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 18 && thePlot <= 18) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 20 && thePlot <= 20) {nBinPlot = 200; xminPlot =20.0; xmaxPlot = 220.0;}
    else if(thePlot >= 21 && thePlot <= 21) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   4.0;}
    else if(thePlot >= 22 && thePlot <= 22) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 23 && thePlot <= 24) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =1000.0;}
    else if(thePlot >= 25 && thePlot <= 25) {nBinPlot =  50; xminPlot =-0.5; xmaxPlot =  49.5;}
    else if(thePlot >= 26 && thePlot <= 31) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 32 && thePlot <= 32) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 33 && thePlot <= 36) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =2000.0;}
    else if(thePlot >= 39 && thePlot <= 39) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   8.0;}
    else if(thePlot >= 42 && thePlot <= 42) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 43 && thePlot <= 43) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   TMath::Pi();}
    TH1D* histos;
    if(thePlot != allPlots-1 && thePlot != 40 && thePlot != 41) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    else                                                        histos = new TH1D("histos", "histos", nBinMVA, xbins);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) {
      for(int j=0; j<allStates; j++) histo[j][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    }
    histos->Reset();histos->Clear();
  }

  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_Zg     = (TH1D*) histoMVA->Clone("histo_Zg"); 
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZ     = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZ     = (TH1D*) histoMVA->Clone("histo_ZZ");	 
  TH1D *histo_FakeM  = (TH1D*) histoMVA->Clone("histo_FakeM");	 
  TH1D *histo_FakeE  = (TH1D*) histoMVA->Clone("histo_FakeE");	 
  TH1D *histo_Higgs  = (TH1D*) histoMVA->Clone("histo_Higgs");	 

  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff2016_m");sprintf(momMName,"CMS_scale2016_m");
  sprintf(effEName,"CMS_eff2016_e");sprintf(momEName,"CMS_scale2016_e");
  sprintf(finalStateName,"3l");

  TH1D* histo_Zg_CMS_MVAZHStatBoundingUp           = new TH1D( Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAZHStatBoundingDown         = new TH1D( Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp           = new TH1D( Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown         = new TH1D( Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingUp     = new TH1D( Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->Sumw2();
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingDown   = new TH1D( Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingDown->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingUp     = new TH1D( Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingDown   = new TH1D( Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingUp     = new TH1D( Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingDown   = new TH1D( Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingDown->Sumw2();

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  TH1D* histo_Zg_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBounding[6];
  for(int nb=0; nb<6; nb++){
    histo_Zg_CMS_QCDScaleBounding[nb]      = new TH1D(Form("histo_Zg_QCDScale_f%d",nb),      Form("histo_Zg_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Zg_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]     = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),      Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_WZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_Higgs_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_Higgs_QCDScale_f%d",nb),   Form("histo_Higgs_QCDScale_f%d",nb),nBinMVA, xbins);   histo_Higgs_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_Zg_CMS_PDFBounding[100];
  TH1D* histo_VVV_CMS_PDFBounding[100];
  TH1D* histo_WZ_CMS_PDFBounding[100];
  TH1D* histo_ZZ_CMS_PDFBounding[100];
  TH1D* histo_Higgs_CMS_PDFBounding[100];
  for(int nb=0; nb<100; nb++){
    histo_Zg_CMS_PDFBounding[nb]      = new TH1D(Form("histo_Zg_PDF_f%d",nb),      Form("histo_Zg_PDF_f%d",nb),     nBinMVA, xbins); histo_Zg_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]     = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_WZ_PDF_f%d",nb),      Form("histo_WZ_PDF_f%d",nb),     nBinMVA, xbins); histo_WZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ZZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_ZZ_PDF_f%d",nb),      Form("histo_ZZ_PDF_f%d",nb),     nBinMVA, xbins); histo_ZZ_CMS_PDFBounding[nb]->Sumw2();
    histo_Higgs_CMS_PDFBounding[nb]   = new TH1D(Form("histo_Higgs_PDF_f%d",nb),   Form("histo_Higgs_PDF_f%d",nb),  nBinMVA, xbins); histo_Higgs_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_Zg_CMS_MVAZHStatBoundingBinUp[nBinMVA];
  TH1D* histo_Zg_CMS_MVAZHStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nBinMVA];
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nBinMVA];
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nBinMVA];
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_Zg_CMS_MVAZHStatBoundingBinUp[nb]         = new TH1D(Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_Zg_CMS_MVAZHStatBoundingBinDown[nb]       = new TH1D(Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingBinDown[nb] ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      ->Sumw2();
    histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nb]   = new TH1D(Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nb]       ->Sumw2();
    histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nb] = new TH1D(Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nb]	 ->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]   = new TH1D(Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]       ->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb] = new TH1D(Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb]	 ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]   = new TH1D(Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]      ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb] = new TH1D(Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]    ->Sumw2();
  }

  TH1D* histo_Zg_CMS_MVALepEffMBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effMName)  , Form("histo_Zg_%sUp",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffMBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effMName), Form("histo_Zg_%sDown",effMName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effMName)  , Form("histo_WZ_%sUp",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effMName), Form("histo_WZ_%sDown",effMName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_Higgs_%sUp",effMName)  , Form("histo_Higgs_%sUp",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_Higgs_%sDown",effMName), Form("histo_Higgs_%sDown",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffMBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effMName)  , Form("histo_Zg_%sAvg",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effMName)  , Form("histo_WZ_%sAvg",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_Higgs_%sAvg",effMName)  , Form("histo_Higgs_%sAvg",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffEBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effEName)  , Form("histo_Zg_%sUp",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffEBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effEName), Form("histo_Zg_%sDown",effEName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effEName)  , Form("histo_WZ_%sUp",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effEName), Form("histo_WZ_%sDown",effEName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_Higgs_%sUp",effEName)  , Form("histo_Higgs_%sUp",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_Higgs_%sDown",effEName), Form("histo_Higgs_%sDown",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffEBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effEName)  , Form("histo_Zg_%sAvg",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effEName)  , Form("histo_WZ_%sAvg",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_Higgs_%sAvg",effEName)  , Form("histo_Higgs_%sAvg",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_Zg_CMS_MVAMETBoundingUp           = new TH1D( Form("histo_Zg_CMS_scale_metUp")  , Form("histo_Zg_CMS_scale_metUp")  , nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAMETBoundingDown         = new TH1D( Form("histo_Zg_CMS_scale_metDown"), Form("histo_Zg_CMS_scale_metDown"), nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown"), nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_Higgs_CMS_scale_metUp")  , Form("histo_Higgs_CMS_scale_metUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_Higgs_CMS_scale_metDown"), Form("histo_Higgs_CMS_scale_metDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVAJESBoundingUp           = new TH1D( Form("histo_Zg_CMS_scale_jesUp")  , Form("histo_Zg_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Zg_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAJESBoundingDown         = new TH1D( Form("histo_Zg_CMS_scale_jesDown"), Form("histo_Zg_CMS_scale_jesDown"), nBinMVA, xbins); histo_Zg_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_jesUp")  , Form("histo_VVV_CMS_scale_jesUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_jesDown"), Form("histo_VVV_CMS_scale_jesDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_jesUp")  , Form("histo_WZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_jesDown"), Form("histo_WZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_jesUp")  , Form("histo_ZZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_jesDown"), Form("histo_ZZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingUp   	= new TH1D( Form("histo_Higgs_CMS_scale_jesUp")  , Form("histo_Higgs_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingDown 	= new TH1D( Form("histo_Higgs_CMS_scale_jesDown"), Form("histo_Higgs_CMS_scale_jesDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVABTAGBoundingUp          = new TH1D( Form("histo_Zg_CMS_eff_b_2016Up")  , Form("histo_Zg_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_Zg_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVABTAGBoundingDown        = new TH1D( Form("histo_Zg_CMS_eff_b_2016Down"), Form("histo_Zg_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_Zg_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVABTAGBoundingUp   	= new TH1D( Form("histo_VVV_CMS_eff_b_2016Up")  , Form("histo_VVV_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVABTAGBoundingDown 	= new TH1D( Form("histo_VVV_CMS_eff_b_2016Down"), Form("histo_VVV_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVABTAGBoundingUp    	= new TH1D( Form("histo_WZ_CMS_eff_b_2016Up")  , Form("histo_WZ_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVABTAGBoundingDown  	= new TH1D( Form("histo_WZ_CMS_eff_b_2016Down"), Form("histo_WZ_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVABTAGBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_eff_b_2016Up")  , Form("histo_ZZ_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVABTAGBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_eff_b_2016Down"), Form("histo_ZZ_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVABTAGBoundingUp   	= new TH1D( Form("histo_Higgs_CMS_eff_b_2016Up")  , Form("histo_Higgs_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_Higgs_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVABTAGBoundingDown 	= new TH1D( Form("histo_Higgs_CMS_eff_b_2016Down"), Form("histo_Higgs_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_Higgs_CMS_MVABTAGBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_BDTMuonScaleBoundingUp           = new TH1D( Form("histo_Zg_CMS_bdt_muonUp")  , Form("histo_Zg_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Zg_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_BDTMuonScaleBoundingDown         = new TH1D( Form("histo_Zg_CMS_bdt_muonDown"), Form("histo_Zg_CMS_scale_jesDown"), nBinMVA, xbins); histo_Zg_CMS_BDTMuonScaleBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_BDTMuonScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_muonUp")  , Form("histo_VVV_CMS_scale_jesUp")  , nBinMVA, xbins); histo_VVV_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMuonScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_muonDown"), Form("histo_VVV_CMS_scale_jesDown"), nBinMVA, xbins); histo_VVV_CMS_BDTMuonScaleBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_BDTMuonScaleBoundingUp    	= new TH1D( Form("histo_WZ_CMS_bdt_muonUp")  , Form("histo_WZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_WZ_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_BDTMuonScaleBoundingDown  	= new TH1D( Form("histo_WZ_CMS_bdt_muonDown"), Form("histo_WZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_WZ_CMS_BDTMuonScaleBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_muonUp")  , Form("histo_ZZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_ZZ_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_muonDown"), Form("histo_ZZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_ZZ_CMS_BDTMuonScaleBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMuonScaleBoundingUp   	= new TH1D( Form("histo_Higgs_CMS_bdt_muonUp")  , Form("histo_Higgs_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Higgs_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMuonScaleBoundingDown 	= new TH1D( Form("histo_Higgs_CMS_bdt_muonDown"), Form("histo_Higgs_CMS_scale_jesDown"), nBinMVA, xbins); histo_Higgs_CMS_BDTMuonScaleBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_BDTElectronScaleBoundingUp           = new TH1D( Form("histo_Zg_CMS_bdt_electronUp")  , Form("histo_Zg_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Zg_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_BDTElectronScaleBoundingDown         = new TH1D( Form("histo_Zg_CMS_bdt_electronDown"), Form("histo_Zg_CMS_scale_jesDown"), nBinMVA, xbins); histo_Zg_CMS_BDTElectronScaleBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_BDTElectronScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_electronUp")  , Form("histo_VVV_CMS_scale_jesUp")  , nBinMVA, xbins); histo_VVV_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTElectronScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_electronDown"), Form("histo_VVV_CMS_scale_jesDown"), nBinMVA, xbins); histo_VVV_CMS_BDTElectronScaleBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_BDTElectronScaleBoundingUp    	= new TH1D( Form("histo_WZ_CMS_bdt_electronUp")  , Form("histo_WZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_WZ_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_BDTElectronScaleBoundingDown  	= new TH1D( Form("histo_WZ_CMS_bdt_electronDown"), Form("histo_WZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_WZ_CMS_BDTElectronScaleBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_electronUp")  , Form("histo_ZZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_ZZ_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_electronDown"), Form("histo_ZZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_ZZ_CMS_BDTElectronScaleBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_BDTElectronScaleBoundingUp   	= new TH1D( Form("histo_Higgs_CMS_bdt_electronUp")  , Form("histo_Higgs_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Higgs_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_BDTElectronScaleBoundingDown 	= new TH1D( Form("histo_Higgs_CMS_bdt_electronDown"), Form("histo_Higgs_CMS_scale_jesDown"), nBinMVA, xbins); histo_Higgs_CMS_BDTElectronScaleBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_BDTMETScaleBoundingUp           = new TH1D( Form("histo_Zg_CMS_bdt_METUp")  , Form("histo_Zg_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Zg_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_BDTMETScaleBoundingDown         = new TH1D( Form("histo_Zg_CMS_bdt_METDown"), Form("histo_Zg_CMS_scale_jesDown"), nBinMVA, xbins); histo_Zg_CMS_BDTMETScaleBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_BDTMETScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_METUp")  , Form("histo_VVV_CMS_scale_jesUp")  , nBinMVA, xbins); histo_VVV_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMETScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_METDown"), Form("histo_VVV_CMS_scale_jesDown"), nBinMVA, xbins); histo_VVV_CMS_BDTMETScaleBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_BDTMETScaleBoundingUp    	= new TH1D( Form("histo_WZ_CMS_bdt_METUp")  , Form("histo_WZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_WZ_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_BDTMETScaleBoundingDown  	= new TH1D( Form("histo_WZ_CMS_bdt_METDown"), Form("histo_WZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_WZ_CMS_BDTMETScaleBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMETScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_METUp")  , Form("histo_ZZ_CMS_scale_jesUp")  , nBinMVA, xbins); histo_ZZ_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMETScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_METDown"), Form("histo_ZZ_CMS_scale_jesDown"), nBinMVA, xbins); histo_ZZ_CMS_BDTMETScaleBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMETScaleBoundingUp   	= new TH1D( Form("histo_Higgs_CMS_bdt_METUp")  , Form("histo_Higgs_CMS_scale_jesUp")  , nBinMVA, xbins); histo_Higgs_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMETScaleBoundingDown 	= new TH1D( Form("histo_Higgs_CMS_bdt_METDown"), Form("histo_Higgs_CMS_scale_jesDown"), nBinMVA, xbins); histo_Higgs_CMS_BDTMETScaleBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_PUBoundingUp               = new TH1D( Form("histo_Zg_CMS_puUp")  , Form("histo_Zg_CMS_puUp")  , nBinMVA, xbins); histo_Zg_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_PUBoundingDown             = new TH1D( Form("histo_Zg_CMS_puDown"), Form("histo_Zg_CMS_puDown"), nBinMVA, xbins); histo_Zg_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingUp           	= new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown         	= new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingUp            	= new TH1D( Form("histo_WZ_CMS_puUp")  , Form("histo_WZ_CMS_puUp")  , nBinMVA, xbins); histo_WZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingDown  	        = new TH1D( Form("histo_WZ_CMS_puDown"), Form("histo_WZ_CMS_puDown"), nBinMVA, xbins); histo_WZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingUp    	        = new TH1D( Form("histo_ZZ_CMS_puUp")  , Form("histo_ZZ_CMS_puUp")  , nBinMVA, xbins); histo_ZZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingDown  	        = new TH1D( Form("histo_ZZ_CMS_puDown"), Form("histo_ZZ_CMS_puDown"), nBinMVA, xbins); histo_ZZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingUp           	= new TH1D( Form("histo_Higgs_CMS_puUp")  , Form("histo_Higgs_CMS_puUp")  , nBinMVA, xbins); histo_Higgs_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingDown         	= new TH1D( Form("histo_Higgs_CMS_puDown"), Form("histo_Higgs_CMS_puDown"), nBinMVA, xbins); histo_Higgs_CMS_PUBoundingDown->Sumw2();

  char outputEventList[200];
  sprintf(outputEventList,"MitZHAnalysis/datacards%s/eventlist_wz_MC.csv",subdirectory.c_str());
  ofstream eventList;
  if(printMCEventList) {
    eventList.open(outputEventList);
    eventList << "eventNumber,njets,realMET,fakeMET,balance,dPhiZMET,ZpT,dPhiJetMET,pdgl1,pdgl2,pdgl3\n";
  }

  unsigned int numberOfLeptons = 3;

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());

    TFile *the_input_file = TFile::Open(infilenamev[ifile].Data());
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file->FindObjectAny("SelBit_tree");
    TTree *the_PDF_tree   = (TTree*)the_input_file->FindObjectAny("pdfReweight");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    eventJets.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BareTaus eventTaus;
    eventTaus.setBranchAddresses(the_input_tree);

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
    histoWZSEL[0]->Scale(0.0);
    histoWZSEL[1]->Scale(0.0);
    histoWZSEL[2]->Scale(0.0);
    histoWZSEL[3]->Scale(0.0);
    histoWZSEL[4]->Scale(0.0);
    histoWZSEL[5]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    if(the_input_tree->GetEntries() != the_SelBit_tree->GetEntries()) {printf("BIG SKIMMING FAILURE\n"); return;}
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;
      the_input_tree->GetEntry(i);

      Bool_t passFilter[12] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
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
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep); }
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
        //else if(selectIdIsoCut("veto",TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	//   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	//                                                                                               {idTight.push_back(0); idLep.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()==numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 25 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() <= 20) continue;

      double dPhiLepMETMin = 999.;
      int signQ = 0;
      double systTotLep[2] = {1.0, 1.0}; // m/e
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) systTotLep[0] = systTotLep[0] * 1.005;
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 11) systTotLep[1] = systTotLep[1] * 1.015;
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);

      passFilter[4] = TMath::Abs(signQ) == 1;
      if(passFilter[4] == kFALSE) continue;

      int nFakeCount = 0;
      double minMassll = 999.0;
      double minMassZ = 999.0;
      double mass3l = 0.0;
      double deltaRllMin = 999.0;
      int type3l = 0;
      int tagZ[3] = {-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
          double deltaRllAux = ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl1]]));
          if(deltaRllAux < deltaRllMin) deltaRllMin = deltaRllAux;

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
        tagZ[2] = nl0;
        break;
      }

      bool tight3rdLepId = true;
      if(idTight[tagZ[2]] == 1 && tagZ[0] != -1){
        tight3rdLepId = selectIdIsoCut(type3rdLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Eta()),(double)(*eventLeptons.iso)[idLep[tagZ[2]]],(int)(*eventLeptons.selBits)[idLep[tagZ[2]]],(double)(*eventLeptons.mva)[idLep[tagZ[2]]]);
      }
      else if(tagZ[0] == -1){
        if(idTight[0] == 1) tight3rdLepId = tight3rdLepId &&
	   selectIdIsoCut(type3rdLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()),(double)(*eventLeptons.iso)[idLep[0]],(int)(*eventLeptons.selBits)[idLep[0]],(double)(*eventLeptons.mva)[idLep[0]]);
        if(idTight[1] == 1) tight3rdLepId = tight3rdLepId &&
	   selectIdIsoCut(type3rdLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()),(double)(*eventLeptons.iso)[idLep[1]],(int)(*eventLeptons.selBits)[idLep[1]],(double)(*eventLeptons.mva)[idLep[1]]);
        if(idTight[2] == 1) tight3rdLepId = tight3rdLepId &&
	   selectIdIsoCut(type3rdLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()),(double)(*eventLeptons.iso)[idLep[2]],(int)(*eventLeptons.selBits)[idLep[2]],(double)(*eventLeptons.mva)[idLep[2]]);
      }

      if(tight3rdLepId == false) continue;

      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]) == 13) type3l += 0;
      else							       type3l += 2;
      TLorentzVector trilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) ));
      mass3l = trilep.M();
      if(tagZ[0] == -1) type3l = 4;
      // Determine flavor of the pair (0 means e-mu pair, 1 means mu-mu, 2 means e-e)
      int typePair = 0;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==13) typePair = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==11) typePair = 2;

      TLorentzVector theFakeMET, theFakeMETUp, theFakeMETDown;
      theFakeMET    .SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px()                    +((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px());
      theFakeMET    .SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py()                    +((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py());
      theFakeMETUp  .SetPx(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px());
      theFakeMETUp  .SetPy(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py());
      theFakeMETDown.SetPx(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px());
      theFakeMETDown.SetPy(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py());

      vector<int> idB,idC;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if     (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idB.push_back(ngen0);
        else if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idC.push_back(ngen0);
      }

      vector<int> idJet,idJetUp,idJetDown; double btagjet[2] = {0., 0.};
      double total_bjet_probTIGHT[2] = {1,1};double total_bjet_probTIGHTUP[2] = {1,1};double total_bjet_probTIGHTDOWN[2] = {1,1};
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      double theHT = 0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 20) continue;
        //bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1 && ((TLorentzVector*)(*eventJets.p4)[nj])->Pt()> 30) {
          if(isWZhinv==false) dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
	  else                dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(theFakeMET));
        }

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
          denBTagging[nJPt][nJEta][jetFlavor]++;
          if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]) numBTaggingTIGHT[nJPt][nJEta][jetFlavor]++;

          double bjet_SFTIGHT = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFTIGHT = btagReaderLTIGHT.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFTIGHT = btagReaderBCTIGHT.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFTIGHT == 0) bjet_SFTIGHT = 1;
          double bjet_SFTIGHTUP = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFTIGHTUP = btagReaderLTIGHTUP.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFTIGHTUP = btagReaderBCTIGHTUP.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFTIGHTUP == 0) bjet_SFTIGHTUP = 1;
          double bjet_SFTIGHTDOWN = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFTIGHTDOWN = btagReaderLTIGHTDOWN.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFTIGHTDOWN = btagReaderBCTIGHTDOWN.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFTIGHTDOWN == 0) bjet_SFTIGHTDOWN = 1;

	  if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]){
	    total_bjet_probTIGHT[0] = total_bjet_probTIGHT[0] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor];
	    total_bjet_probTIGHT[1] = total_bjet_probTIGHT[1] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHT;
	    total_bjet_probTIGHTUP[0] = total_bjet_probTIGHTUP[0] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor];
	    total_bjet_probTIGHTUP[1] = total_bjet_probTIGHTUP[1] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTUP;
	    total_bjet_probTIGHTDOWN[0] = total_bjet_probTIGHTDOWN[0] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor];
	    total_bjet_probTIGHTDOWN[1] = total_bjet_probTIGHTDOWN[1] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTDOWN;
	  } else {
	    total_bjet_probTIGHT[0] = total_bjet_probTIGHT[0] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor]);
	    total_bjet_probTIGHT[1] = total_bjet_probTIGHT[1] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHT);
	    total_bjet_probTIGHTUP[0] = total_bjet_probTIGHTUP[0] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor]);
	    total_bjet_probTIGHTUP[1] = total_bjet_probTIGHTUP[1] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTUP);
	    total_bjet_probTIGHTDOWN[0] = total_bjet_probTIGHTDOWN[0] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor]);
	    total_bjet_probTIGHTDOWN[1] = total_bjet_probTIGHTDOWN[1] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTDOWN);
	  }

        }

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()      > 30) idJet.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*1.03 > 30) idJetUp.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*0.97 > 30) idJetDown.push_back(nj);

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;
	
        if     ((float)(*eventJets.bDiscr)[nj] > btagjet[0]) {btagjet[1] = btagjet[0]; btagjet[0] = (float)(*eventJets.bDiscr)[nj];}
	else if((float)(*eventJets.bDiscr)[nj] > btagjet[1]) {btagjet[1] = (float)(*eventJets.bDiscr)[nj];}

        theHT = theHT + ((TLorentzVector*)(*eventJets.p4)[nj])->Pt();
      }

      passFilter[ 5] = minMassll > 4;
      passFilter[ 6] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30;
      passFilter[ 7] = minMassZ > minMass && minMassZ < maxMass && type3l != 4;
      passFilter[ 8] = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt() > 20;
      passFilter[ 9] = mass3l > 100;
      passFilter[10] = true; if(applyBtagging) passFilter[10] = bDiscrMax < bTagCuts[0];
      passFilter[11] = true;
      if(passFilter[7]==kTRUE && (tagZ[0] == tagZ[1] || tagZ[0] == tagZ[2] || tagZ[1] == tagZ[2])) {printf("ZPROBLEM!\n");assert(0);return;}

      bool passWHFilter[7] = {false, false, false, false, false, false, false};
      passWHFilter[0] = minMassll > 4 && minMassll < 80.0;
      passWHFilter[1] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30;
      passWHFilter[2] = (minMassZ < 71 || minMassZ > 111) || type3l == 4;
      passWHFilter[3] = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt() > 20;
      passWHFilter[4] = (mass3l < 76 || mass3l > 106) || type3l == 4;
      passWHFilter[5] =  bDiscrMax < 0.800 && idSoft.size() == 0;
      passWHFilter[6] = idJet.size() == 0 || (double)((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Pt() < 40;
      if(applyWHSel){
        passFilter[ 5] = passWHFilter[0];
        passFilter[ 6] = passWHFilter[1];
        passFilter[ 7] = passWHFilter[2];
        passFilter[ 8] = passWHFilter[3];
        passFilter[ 9] = passWHFilter[4];
        passFilter[10] = passWHFilter[5];
        passFilter[11] = passWHFilter[6];
      }

     bool passEvolFilter[numberCuts] = {passFilter[ 5],passFilter[ 6],passFilter[ 7],passFilter[ 8],passFilter[ 9],passFilter[10],passFilter[11]};

     int sumEvol = 0;
     bool totalSel = kTRUE;
     for(int isel=0; isel<numberCuts; isel++) {
       totalSel = totalSel && passEvolFilter[isel];
       if(totalSel == kTRUE) sumEvol++;
     }

      bool passNMinusOne[7] = {                 passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11],
                               passFilter[5] &&                  passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11],
                               passFilter[5] && passFilter[6] && 		  passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11],
                               passFilter[5] && passFilter[6] && passFilter[7] &&		   passFilter[9] && passFilter[10] && passFilter[11],
                               passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] 		 && passFilter[10] && passFilter[11],
			       passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9]                   && passFilter[11],
			       passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};

      bool passAllCuts[nSelTypes] = {passFilter[5]   && passFilter[6]   && passFilter[7]   && passFilter[8]   && passFilter[9]   && passFilter[10]  && passFilter[11],
                                     passWHFilter[0] && passWHFilter[1] && passWHFilter[2] && passWHFilter[3] && passWHFilter[4] && passWHFilter[5] && passWHFilter[6],
                                     passFilter[5]   && passFilter[6]   && passFilter[7]   && passFilter[8]   && passFilter[9]   && passFilter[10]  && passFilter[11]};

      bool controlSel[3] = {passFilter[5] && passFilter[6] && !passFilter[7] && passFilter[9] && bDiscrMax > 0.605,
                            passFilter[5] &&                  !passFilter[7] && passFilter[9] && bDiscrMax > 0.605,
			    passFilter[5] && !passFilter[6] && passFilter[10] && ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() > 20};

      bool passTTZSel[2] = {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() >= 4 && btagjet[0] > 0.560 && btagjet[1] > 0.560,
                            passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() == 3 && btagjet[0] > 0.800 && btagjet[1] > 0.560};

      bool passVBFLoose = passAllCuts[WZSEL] && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 250;
      bool passVBF = passAllCuts[0] && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 500 &&
		     TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) > 2.5;

      TLorentzVector dilepZ(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) ) ));
      double dPhiDiLepMET = TMath::Abs(dilepZ.DeltaPhi(theFakeMET));
      double ptFrac = TMath::Abs(dilepZ.Pt()-theFakeMET.Pt())/dilepZ.Pt();
      double dphill = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]));
      double detall = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Eta()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Eta());
      double drll = sqrt(dphill*dphill+detall*detall);
      bool useZHcuts = true;
      bool passDelphiLL   = useZHcuts ? (drll < 1.8) : (true);
      bool passNjets      = idJet.size() <= 1;
      bool passFakeMET    = theFakeMET.Pt() > 50;
      bool passPTFrac     = ptFrac < (useZHcuts ? 0.4 : 1.0);
      bool passDPhiZMET   = dPhiDiLepMET > (useZHcuts ? 2.6 : 2.0);
      bool passPTLL       = dilepZ.Pt() > 60;
      bool passDPhiJetMET = useZHcuts ? (dPhiJetMET == -1 || dPhiJetMET >= 0.5) : (true);

      if(isWZhinv) passAllCuts[WZSEL] = passAllCuts[WZSEL] && passNjets && passFakeMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL;
      if(printMCEventList && infilecatv[ifile] == 3 && theFakeMET.Pt() > 100 && passAllCuts[WZSEL])
        eventList << Form("%lld,%d,%f,%f,%f,%f,%f,%f,%d,%d,%d\n", eventEvent.eventNum, (int)idJet.size(), ((TLorentzVector*)(*eventMet.p4)[0])->Pt(), theFakeMET.Pt(), ptFrac, dPhiDiLepMET, dilepZ.Pt(), dPhiJetMET, (int)(*eventLeptons.pdgId)[idLep[tagZ[0]]], (int)(*eventLeptons.pdgId)[idLep[tagZ[1]]], (int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]);
      double deltaPhiLeptonMet = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtLN = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiLeptonMet)));

      double deltaPhiTriLeptonMet = TMath::Abs(trilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtEvent = TMath::Sqrt(2.0*trilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiTriLeptonMet)));
      
      // Evaluate nominal BDT value
      double bdt_value=-2;
      if(useBDT && passAllCuts[WZSEL]) {
        TLorentzVector lepton1 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]),
                       lepton2 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]),
                       MET     = theFakeMET,
                       jet1    = idJet.size() > 0 ? *((TLorentzVector*)(*eventJets.p4)[idJet[0]]) : TLorentzVector(0,0,0,0);
         bdt_value = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0,0,0,0);
      }

      bool passSystCuts[nSystTypes] = {
        passFilter[5]   && passFilter[6]   						     && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11],
        passFilter[5]   && passFilter[6]   						     && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11],
        passFilter[5] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11],
        passFilter[5] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11]
      };
      if(isWZhinv) {
        passSystCuts[JESUP]   = passFilter[5]   && passFilter[6]   && passFilter[7]   && passFilter[8]   && passFilter[9]   && passFilter[10]  && passFilter[11] && 
	                        idJetUp.size()   <= 1 && passFakeMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL;
        passSystCuts[JESDOWN] = passFilter[5]   && passFilter[6]   && passFilter[7]   && passFilter[8]   && passFilter[9]   && passFilter[10]  && passFilter[11] && 
	                        idJetDown.size() <= 1 && passFakeMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL;
        passSystCuts[METUP]   = passAllCuts[WZSEL];
        passSystCuts[METDOWN] = passAllCuts[WZSEL];
      }
      // Do the BDT nuisance evaluations
      double bdt_muonScaleDown=-1, bdt_muonScaleUp=-1, bdt_electronScaleDown=-1, bdt_electronScaleUp=-1, bdt_METScaleDown=-1, bdt_METScaleUp=-1, bdt_jetScaleDown=-1, bdt_jetScaleUp=-1;
      double MVAVar_muonScaleDown=-9999, MVAVar_muonScaleUp=-9999, MVAVar_electronScaleDown=-9999, MVAVar_electronScaleUp=-9999, MVAVar_METScaleDown=-9999, MVAVar_METScaleUp=-9999, MVAVar_jetScaleDown=-9999, MVAVar_jetScaleUp=-9999;
      if(useBDT && passAllCuts[WZSEL]) {
        TLorentzVector lepton1 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]),
                       lepton2 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]),
                       MET     = theFakeMET,
                       jet1    = idJet.size() > 0 ? *((TLorentzVector*)(*eventJets.p4)[idJet[0]]) : TLorentzVector(0,0,0,0);
        double lepton1_scale_variation, lepton2_scale_variation, neutrino1_scale_variation;
        // BDT variation with the muon scale variation (flat 1%)
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==13 ? 0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==13 ? 0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==13 ? 0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py());
         bdt_muonScaleUp = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_muonScaleUp = getMVAVar(MVAVarType, passAllCuts[WZSEL], typePair, MET.Pt(), 0, dilepZ.M(), bdt_muonScaleUp, xbins[nBinMVA]);
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==13 ? -0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==13 ? -0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==13 ? -0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py());
        bdt_muonScaleDown = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_muonScaleDown = getMVAVar(MVAVarType, passAllCuts[WZSEL], typePair, MET.Pt(), 0, dilepZ.M(), bdt_muonScaleUp, xbins[nBinMVA]);
        // BDT variation with the electron scale variation (flat 1%)
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==11 ? 0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==11 ? 0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==11 ? 0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py());
        bdt_electronScaleUp = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_electronScaleUp = getMVAVar(MVAVarType, passAllCuts[WZSEL], typePair, MET.Pt(), 0, dilepZ.M(), bdt_electronScaleUp, xbins[nBinMVA]);
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==11 ? -0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==11 ? -0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==11 ? -0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py());
        bdt_electronScaleDown = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_electronScaleDown = getMVAVar(MVAVarType, passAllCuts[WZSEL], typePair, MET.Pt(), 0, dilepZ.M(), bdt_electronScaleUp, xbins[nBinMVA]);
        // BDT variation with the MET scale variation (from Nero)
        MET = theFakeMETUp;
        bdt_METScaleUp = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, 0, 0);
	MVAVar_METScaleUp = getMVAVar(MVAVarType, passAllCuts[WZSEL], typePair, MET.Pt(), 0, dilepZ.M(), bdt_METScaleUp, xbins[nBinMVA]);
        MET = theFakeMETDown;
        bdt_METScaleDown = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, 0, 0);
	MVAVar_METScaleDown = getMVAVar(MVAVarType, passAllCuts[WZSEL], typePair, MET.Pt(), 0, dilepZ.M(), bdt_METScaleUp, xbins[nBinMVA]);
      }

      // begin event weighting
      vector<int>wBoson;
      vector<int>zBoson;
      vector<bool> isGenDupl;
      int numberQuarks[2] = {0,0};
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) wBoson.push_back(ngen0);
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23) zBoson.push_back(ngen0);
        if     (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[0]++;
        else if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[1]++;
        isGenDupl.push_back(0);
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) isGenDupl[ngen0] = 1;
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) continue;
        for(int ngen1=ngen0+1; ngen1<eventMonteCarlo.p4->GetEntriesFast(); ngen1++) {
	  if((int)(*eventMonteCarlo.pdgId)[ngen0] != (int)(*eventMonteCarlo.pdgId)[ngen1]) continue;
          if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen1])) < 0.02) {
	    isGenDupl[ngen0] = 1;
	    break;
	  }
        }
      }
      vector<int> isGenLep; unsigned int goodIsGenLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenLepton = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenLepton = true;
	    break;
	  }
	}
	if(isGenLepton == true) {isGenLep.push_back(1); goodIsGenLep++;}
	else                    {isGenLep.push_back(0);}
      }

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
	  if(tagZ[2] != (int)nl && tagZ[0] != -1)
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,true);
          else
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		type3rdLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,true);
        }
      }

      // fake rate
      nFakeCount = 0;
      unsigned int typeFakeLepton[2] = {0,0};
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     ((infilecatv[ifile] == 0 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add Z+jets from data
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    if(tagZ[0] == (int)nl) nFakeCount = nFakeCount + 1;
	    if(tagZ[1] == (int)nl) nFakeCount = nFakeCount + 2;
	    if(tagZ[2] == (int)nl) nFakeCount = nFakeCount + 4;
	    if(tagZ[2] != (int)nl && tagZ[0] != -1)
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    else
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,type3rdLepSel.Data());
	    theCategory = 1;
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) typeFakeLepton[0]++;
	    else                                                        typeFakeLepton[1]++;
	  }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-3) fakeSF = -1.0 * fakeSF; // triple fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF = +1.0 * fakeSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-3) fakeSF = +1.0 * fakeSF; // triple fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF = +1.0 * fakeSF; // single fake, data
	  /*
	  if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) {
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13){
	      if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 15) fakeSF = 2.6 * fakeSF;
	      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 20) fakeSF = 1.6 * fakeSF;
	      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 15) fakeSF = 3.7 * fakeSF;
	      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 20) fakeSF = 2.0 * fakeSF;
            }
          }
	  */
	}
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 2 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] == 2 && goodIsTight != idTight.size()){ // remove Z+gamma, fakeable objects
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 2){ // data or Z+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infilecatv[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }

      double addSpecialWeight = 1.0;
      if(infilecatv[ifile] == 2) addSpecialWeight = overallZgSF[0];

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*addSpecialWeight;

      // Btag scale factor
      totalWeight = totalWeight * total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0];

      double btagCorr[2] = {(total_bjet_probTIGHTUP[1]  /total_bjet_probTIGHTUP[0]  )/(total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0]),
                            (total_bjet_probTIGHTDOWN[1]/total_bjet_probTIGHTDOWN[0])/(total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0])};

      if(totalWeight == 0) continue;

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
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[0]) sumEventsProcess[ifile] += totalWeight;

      if(passAllCuts[0] && infilecatv[ifile] == 0) totalFakeDataCount[type3l][nFakeCount]++;

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[i+type3l*nSelTypes][theCategory] += totalWeight;
          weiDecay[i+type3l*nSelTypes][theCategory] += totalWeight*totalWeight;
        }
      }

      for(int nl=0; nl <=sumEvol; nl++) if(fakeSF == 1) histoWZSEL[type3l]->Fill((double)nl,totalWeight);
      for(int nl=0; nl <=sumEvol; nl++) if(fakeSF == 1) histoWZSEL[5]     ->Fill((double)nl,totalWeight);

      //if(theCategory == 0 && passAllCuts[WZSEL]){
      //  	printf("DATA %d %d %lld %f %f %f %d %f %f %f\n",eventEvent.runNum,eventEvent.lumiNum,eventEvent.eventNum,
      //    minMassll,((TLorentzVector*)(*eventMet.p4)[0])->Pt(),minMassZ,type3l,((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),mass3l,bDiscrMax);
      //}

      for(int thePlot=0; thePlot<allPlots; thePlot++){
	double theVar = 0.0;
	bool makePlot = false;
 	if     (thePlot ==  0 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot ==  1 && passNMinusOne[0])                       {makePlot = true;theVar = TMath::Min(minMassll,199.999);}
	else if(thePlot ==  2 && passNMinusOne[1])                       {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot ==  3 && passNMinusOne[2])                       {makePlot = true;theVar = TMath::Max(TMath::Min(minMassZ,149.999),50.001);}
	else if(thePlot ==  4 && passNMinusOne[3])                       {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot ==  5 && passNMinusOne[4])                       {makePlot = true;theVar = TMath::Min(mass3l,399.999);}
	else if(thePlot ==  6 && passNMinusOne[5])                       {makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot ==  7 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min(deltaRllMin,3.999);}
	else if(thePlot ==  8 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min(minPMET,199.999);}
	else if(thePlot ==  9 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Pt(),199.999);}
	else if(thePlot == 10 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Pt(),199.999);}
	else if(thePlot == 11 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 12 && passNMinusOne[6])                       {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 13 && passNMinusOne[6])                       {makePlot = true;if(idJet.size() >= 1) theVar = TMath::Min((double)((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Pt(),99.999);}
	else if(thePlot == 14 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	else if(thePlot == 15 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 16 && passAllCuts[WZSEL])                     {makePlot = true;theVar = dPhiJetMET*180/TMath::Pi();}
	else if(thePlot == 17 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 18 && controlSel[0])                          {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 19 && controlSel[1])                          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 20 && controlSel[2])                          {makePlot = true;theVar = TMath::Max(TMath::Min(minMassZ,219.999),20.001);}
	else if(thePlot == 21 && passAllCuts[WHSEL])                     {makePlot = true;theVar = TMath::Min(deltaRllMin,3.999);}
	else if(thePlot == 22 && passVBF)                                {makePlot = true;theVar = TMath::Min(mtEvent,499.999);}
	else if(thePlot == 23 && passTTZSel[0])                          {makePlot = true;theVar = TMath::Min(theHT,999.999);}
	else if(thePlot == 24 && passTTZSel[1])                          {makePlot = true;theVar = TMath::Min(theHT,999.999);}
	else if(thePlot == 25 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min((double)(numberQuarks[0]+10*numberQuarks[1]),49.499);}
	else if(thePlot == 26 && passNMinusOne[3])	                 {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot == 27 && passAllCuts[WZSEL])                     {makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot == 28 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 29 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 30 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 31 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 32 && controlSel[2])                          {makePlot = true;theVar = TMath::Max(TMath::Min(mass3l,249.999),50.001);}
	else if(thePlot == 33 && passNMinusOne[5] && numberQuarks[1] == 0                  ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 34 && passNMinusOne[5] && numberQuarks[1] == 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 35 && passNMinusOne[5] && numberQuarks[1]  > 0                  ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 36 && passNMinusOne[5] && numberQuarks[1]  > 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 37 && passVBFLoose)                           {makePlot = true;theVar = TMath::Min(mtEvent,499.999);}
	else if(thePlot == 38 && passVBFLoose)                           {makePlot = true;theVar = TMath::Min((( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M(),1999.999);}
	else if(thePlot == 39 && passVBFLoose)                           {makePlot = true;theVar = TMath::Min(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()),7.999);}
	else if(thePlot == 40 && passAllCuts[WZSEL])                     {makePlot = true;theVar = TMath::Min(theFakeMET.Pt(),xbins[nBinMVA]-0.001);}
	else if(thePlot == 41 && passAllCuts[WZSEL])                     {makePlot = true;theVar = TMath::Min(dilepZ.Pt(),xbins[nBinMVA]-0.001);}
	else if(thePlot == 42 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
	else if(thePlot == 43 && passAllCuts[WZLOOSESEL])                {makePlot = true;theVar = dPhiDiLepMET;}

	if(makePlot) histo[     5][thePlot][theCategory]->Fill(theVar,totalWeight);
	if(makePlot) histo[type3l][thePlot][theCategory]->Fill(theVar,totalWeight);
      }
      
      if(1) {
	double MVAVar = (double)type3l;
	double MVAVarMETSyst[2] = {(double)type3l, (double)type3l};
	if(isWZhinv){
          if(MVAVarType==1) {
  	    MVAVar = TMath::Min(theFakeMET.Pt(),xbins[nBinMVA]-0.001);
  	    MVAVarMETSyst[0] = TMath::Min(theFakeMETUp.Pt(),xbins[nBinMVA]-0.001);
  	    MVAVarMETSyst[1] = TMath::Min(theFakeMETDown.Pt(),xbins[nBinMVA]-0.001);
          } else if(MVAVarType==3) {
            MVAVar = getMVAVar(MVAVarType, passAllCuts[WZSEL], typePair, theFakeMET.Pt(), 0, dilepZ.M(), bdt_value, xbins[nBinMVA]);
            MVAVarMETSyst[0] = MVAVar; 
            MVAVarMETSyst[1] = MVAVar; 
          }
	}

	// Avoid QCD scale weights that are anomalous high
	double maxQCDscale = (TMath::Abs((double)eventMonteCarlo.r1f2)+TMath::Abs((double)eventMonteCarlo.r1f5)+TMath::Abs((double)eventMonteCarlo.r2f1)+
                              TMath::Abs((double)eventMonteCarlo.r2f2)+TMath::Abs((double)eventMonteCarlo.r5f1)+TMath::Abs((double)eventMonteCarlo.r5f5))/6.0;
        double PDFAvg = 0.0;
        if(infilecatv[ifile] != 0 && passAllCuts[WZSEL]){
          if(initPDFTag != -1)
          for(int npdf=0; npdf<100; npdf++) PDFAvg = PDFAvg + TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]);
          PDFAvg = PDFAvg/100.0;
        }

        if     (theCategory == 0){
	  if(passAllCuts[WZSEL]) histo_Data->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 1){
	  if(passAllCuts[WZSEL]) {
	    if(typeFakeLepton[0]+typeFakeLepton[1] != goodIsTight == idTight.size()) {printf("PROBLEMFake %d %d %d %d\n",typeFakeLepton[0],typeFakeLepton[1],goodIsTight,(int)idTight.size()); return;}
	    if     (typeFakeLepton[0] > typeFakeLepton[1]) histo_FakeM->Fill(MVAVar,totalWeight);
	    else if(typeFakeLepton[0] < typeFakeLepton[1]) histo_FakeE->Fill(MVAVar,totalWeight);
	    else {
	      histo_FakeM->Fill(MVAVar,totalWeight/2.);
	      histo_FakeE->Fill(MVAVar,totalWeight/2.);
	    }
	  }
        }
        else if(theCategory == 2){
	  if(passAllCuts[WZSEL]) {
             histo_Zg->Fill(MVAVar,totalWeight);

	     histo_Zg_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Zg_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Zg_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Zg_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Zg_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_Zg_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_Zg_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	     if(useBDT) {
               histo_Zg_CMS_BDTMuonScaleBoundingUp      ->Fill(MVAVar_muonScaleUp      , totalWeight);
               histo_Zg_CMS_BDTMuonScaleBoundingDown    ->Fill(MVAVar_muonScaleDown    , totalWeight);
               histo_Zg_CMS_BDTElectronScaleBoundingUp  ->Fill(MVAVar_electronScaleUp  , totalWeight);
               histo_Zg_CMS_BDTElectronScaleBoundingDown->Fill(MVAVar_electronScaleDown, totalWeight);
               histo_Zg_CMS_BDTMETScaleBoundingUp       ->Fill(MVAVar_METScaleUp       , totalWeight);
               histo_Zg_CMS_BDTMETScaleBoundingDown     ->Fill(MVAVar_METScaleDown     , totalWeight);
             }
	  }
          if(passSystCuts[JESUP])  histo_Zg_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_Zg_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_Zg_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_Zg_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 3){
	  if(passAllCuts[WZSEL]) {
	     histo_WZ->Fill(MVAVar,totalWeight);
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
             histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_WZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_WZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	     if(useBDT) {
               histo_WZ_CMS_BDTMuonScaleBoundingUp      ->Fill(MVAVar_muonScaleUp      , totalWeight);
               histo_WZ_CMS_BDTMuonScaleBoundingDown    ->Fill(MVAVar_muonScaleDown    , totalWeight);
               histo_WZ_CMS_BDTElectronScaleBoundingUp  ->Fill(MVAVar_electronScaleUp  , totalWeight);
               histo_WZ_CMS_BDTElectronScaleBoundingDown->Fill(MVAVar_electronScaleDown, totalWeight);
               histo_WZ_CMS_BDTMETScaleBoundingUp       ->Fill(MVAVar_METScaleUp       , totalWeight);
               histo_WZ_CMS_BDTMETScaleBoundingDown     ->Fill(MVAVar_METScaleDown     , totalWeight);
             }
	  }
          if(passSystCuts[JESUP])  histo_WZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
	}
        else if(theCategory == 4){
	  if(passAllCuts[WZSEL]) {
	     histo_ZZ->Fill(MVAVar,totalWeight);
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
             histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_ZZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ZZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_ZZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	     if(useBDT) {
               histo_ZZ_CMS_BDTMuonScaleBoundingUp      ->Fill(MVAVar_muonScaleUp      , totalWeight);
               histo_ZZ_CMS_BDTMuonScaleBoundingDown    ->Fill(MVAVar_muonScaleDown    , totalWeight);
               histo_ZZ_CMS_BDTElectronScaleBoundingUp  ->Fill(MVAVar_electronScaleUp  , totalWeight);
               histo_ZZ_CMS_BDTElectronScaleBoundingDown->Fill(MVAVar_electronScaleDown, totalWeight);
               histo_ZZ_CMS_BDTMETScaleBoundingUp       ->Fill(MVAVar_METScaleUp       , totalWeight);
               histo_ZZ_CMS_BDTMETScaleBoundingDown     ->Fill(MVAVar_METScaleDown     , totalWeight);
             }
          }
          if(passSystCuts[JESUP])  histo_ZZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 5){
	  if(passAllCuts[WZSEL]) {
	     histo_VVV->Fill(MVAVar,totalWeight);
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
             histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
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
	  if(passAllCuts[WZSEL]) {
	     histo_Higgs->Fill(MVAVar,totalWeight);
	     histo_Higgs_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Higgs_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Higgs_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Higgs_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Higgs_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_Higgs_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_Higgs_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	     if(useBDT) {
               histo_Higgs_CMS_BDTMuonScaleBoundingUp      ->Fill(MVAVar_muonScaleUp      , totalWeight);
               histo_Higgs_CMS_BDTMuonScaleBoundingDown    ->Fill(MVAVar_muonScaleDown    , totalWeight);
               histo_Higgs_CMS_BDTElectronScaleBoundingUp  ->Fill(MVAVar_electronScaleUp  , totalWeight);
               histo_Higgs_CMS_BDTElectronScaleBoundingDown->Fill(MVAVar_electronScaleDown, totalWeight);
               histo_Higgs_CMS_BDTMETScaleBoundingUp       ->Fill(MVAVar_METScaleUp       , totalWeight);
               histo_Higgs_CMS_BDTMETScaleBoundingDown     ->Fill(MVAVar_METScaleDown     , totalWeight);
             }
          }
          if(passSystCuts[JESUP])  histo_Higgs_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_Higgs_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_Higgs_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_Higgs_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      } // making data cards
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);
    for(int nc=0; nc<numberCuts+1; nc++){
      printf("(%20s): %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",cutName[nc].Data(),histoWZSEL[0]->GetBinContent(nc+1),histoWZSEL[1]->GetBinContent(nc+1),histoWZSEL[2]->GetBinContent(nc+1),histoWZSEL[3]->GetBinContent(nc+1),histoWZSEL[4]->GetBinContent(nc+1),histoWZSEL[5]->GetBinContent(nc+1));
    }
    the_input_file->Close();
  } // end of chain

  if(printMCEventList) eventList.close();
  histo[5][allPlots-1][0]->Add(histo_Data);
  histo[5][allPlots-1][1]->Add(histo_FakeM);
  histo[5][allPlots-1][1]->Add(histo_FakeE);
  histo[5][allPlots-1][2]->Add(histo_Zg);
  histo[5][allPlots-1][3]->Add(histo_WZ);
  histo[5][allPlots-1][4]->Add(histo_ZZ);
  histo[5][allPlots-1][5]->Add(histo_VVV);
  histo[5][allPlots-1][6]->Add(histo_Higgs);  

  double theWZSF = 1.0;
  if(1) {
    theWZSF = (histo_Data->GetSumOfWeights()-(histo_FakeM->GetSumOfWeights()+histo_FakeE->GetSumOfWeights()+
               histo_Zg->GetSumOfWeights()+histo_ZZ->GetSumOfWeights()+histo_VVV->GetSumOfWeights()))/histo_WZ->GetSumOfWeights(); 
    printf("SF: %f\n",theWZSF); histo[5][allPlots-1][3]->Scale(theWZSF);
  }

  printf("----------------------totalFakeDataCount--------------------------------\n");
  printf("      TTT    FTT    TFT    FFT    TTF    FTF    TFF    FFF\n");
  for(int ni=0; ni<4; ni++) {
    printf("(%d): ",ni);
    for(int nj=0; nj<9; nj++) printf("%6.1f ",totalFakeDataCount[ni][nj]);
    printf("\n");
  }

  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsType[5] = {0,0,0,0,0}; double sumEventsTypeE[5] = {0,0,0,0,0};
    for(int np=0; np<histBins; np++) {   
      double theSF = 1.0;
      if(ns==WZSEL && np == 3) theSF = theWZSF;
       printf("(%6s): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f -> %8.2f +/- %6.2f\n",
       processName[np].Data(),bgdDecay[ns+nSelTypes*0][np]*theSF,sqrt(weiDecay[ns+nSelTypes*0][np])*theSF,bgdDecay[ns+nSelTypes*1][np]*theSF,sqrt(weiDecay[ns+nSelTypes*1][np])*theSF,
                              bgdDecay[ns+nSelTypes*2][np]*theSF,sqrt(weiDecay[ns+nSelTypes*2][np])*theSF,bgdDecay[ns+nSelTypes*3][np]*theSF,sqrt(weiDecay[ns+nSelTypes*3][np])*theSF,
			      bgdDecay[ns+nSelTypes*4][np]*theSF,sqrt(weiDecay[ns+nSelTypes*4][np])*theSF,
                              (bgdDecay[ns+nSelTypes*0][np]+bgdDecay[ns+nSelTypes*1][np]+bgdDecay[ns+nSelTypes*2][np]+bgdDecay[ns+nSelTypes*3][np]+bgdDecay[ns+nSelTypes*4][np])*theSF,
                              sqrt(weiDecay[ns+nSelTypes*0][np]+weiDecay[ns+nSelTypes*1][np]+weiDecay[ns+nSelTypes*2][np]+weiDecay[ns+nSelTypes*3][np]+weiDecay[ns+nSelTypes*4][np])*theSF);
       if(np != 0 && np != 6){
         sumEventsType[0] = sumEventsType[0] + bgdDecay[ns+nSelTypes*0][np]*theSF; sumEventsTypeE[0] = sumEventsTypeE[0] + weiDecay[ns+nSelTypes*0][np]*theSF;
         sumEventsType[1] = sumEventsType[1] + bgdDecay[ns+nSelTypes*1][np]*theSF; sumEventsTypeE[1] = sumEventsTypeE[1] + weiDecay[ns+nSelTypes*1][np]*theSF;
         sumEventsType[2] = sumEventsType[2] + bgdDecay[ns+nSelTypes*2][np]*theSF; sumEventsTypeE[2] = sumEventsTypeE[2] + weiDecay[ns+nSelTypes*2][np]*theSF;
         sumEventsType[3] = sumEventsType[3] + bgdDecay[ns+nSelTypes*3][np]*theSF; sumEventsTypeE[3] = sumEventsTypeE[3] + weiDecay[ns+nSelTypes*3][np]*theSF;
         sumEventsType[4] = sumEventsType[4] + bgdDecay[ns+nSelTypes*4][np]*theSF; sumEventsTypeE[4] = sumEventsTypeE[4] + weiDecay[ns+nSelTypes*4][np]*theSF;
      }
    }
    printf("(...bkg): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f -> %8.2f +/- %6.2f\n",
           sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
	   sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]),
	   sumEventsType[4],sqrt(sumEventsTypeE[4]),
           sumEventsType[0]+sumEventsType[1]+sumEventsType[2]+sumEventsType[3]+sumEventsType[4],
           sqrt(sumEventsTypeE[0]+sumEventsTypeE[1]+sumEventsTypeE[2]+sumEventsTypeE[3]+sumEventsTypeE[4]));
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }

  double sumEvents = 0;
  for(int np=1; np<histBins; np++) sumEvents += histo[5][0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[5][0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[5][0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);
  double sumEventsType[5] = {0,0,0,0,0};
  double sumEventsTypeE[5] = {0,0,0,0,0};
  printf("                  all                 mmm                 eem                 mme                 eee\n");
  printf("-----------------------------------------------------------------------------------------------------------\n");
  for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=4; i++) {sumEvents = sumEvents + histo[5][15][np]->GetBinContent(i); sumEventsE = sumEventsE + histo[5][15][np]->GetBinError(i)*histo[5][15][np]->GetBinError(i);}
    if(np != 0 && np != 6){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                          sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE;
    sumEventsType[1] = sumEventsType[1] + histo[5][15][np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[5][15][np]->GetBinError(1) * histo[5][15][np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + histo[5][15][np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[5][15][np]->GetBinError(2) * histo[5][15][np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + histo[5][15][np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + histo[5][15][np]->GetBinError(3) * histo[5][15][np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + histo[5][15][np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + histo[5][15][np]->GetBinError(4) * histo[5][15][np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    processName[np].Data(),sumEvents,sumEventsE,
    histo[5][15][np]->GetBinContent(1),histo[5][15][np]->GetBinError(1),histo[5][15][np]->GetBinContent(2),histo[5][15][np]->GetBinError(2),
    histo[5][15][np]->GetBinContent(3),histo[5][15][np]->GetBinError(3),histo[5][15][np]->GetBinContent(4),histo[5][15][np]->GetBinError(4));
    if(np==0)
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    "   all",sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]),sumEventsType[4],sqrt(sumEventsTypeE[4]));
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    for(int j=0; j<allStates; j++){
      char output[200];
      sprintf(output,"MitZHAnalysis/plots%s/histowz_nice_%d_%d.root",subdirectory.c_str(), j,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[j][thePlot][np]->Write();
      outFilePlotsNote->Close();
    }
  }

  printf("QCD Corr: WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) ZH(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Zg->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  for(int i=1; i<=histo_Zg->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_Zg_CMS_MVAZHStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorUp  *histo_Zg->GetBinError(i),0.000001));
    histo_Zg_CMS_MVAZHStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorDown*histo_Zg->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorUp  *histo_WZ     ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorDown*histo_WZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorUp  *histo_ZZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorDown*histo_ZZ     ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorUp  *histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingDown->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorDown*histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorUp  *histo_FakeE  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingDown->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorDown*histo_FakeE  ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_Higgs	->GetBinContent(i)+factorUp  *histo_Higgs->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingDown->SetBinContent(i,TMath::Max(histo_Higgs	->GetBinContent(i)+factorDown*histo_Higgs->GetBinError(i),0.000001));

    histo_Zg_CMS_MVAZHStatBoundingBinUp[i-1]        ->Add(histo_Zg     ); histo_Zg_CMS_MVAZHStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorUp  *histo_Zg->GetBinError(i),0.000001));
    histo_Zg_CMS_MVAZHStatBoundingBinDown[i-1]      ->Add(histo_Zg     ); histo_Zg_CMS_MVAZHStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorDown*histo_Zg->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	    ->Add(histo_VVV    ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]     ->SetBinContent(i,TMath::Max(histo_VVV	->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]    ->Add(histo_VVV    ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]   ->SetBinContent(i,TMath::Max(histo_VVV	->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	    ->Add(histo_WZ     ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]       ->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorUp  *histo_WZ     ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	    ->Add(histo_WZ     ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]     ->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorDown*histo_WZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	    ->Add(histo_ZZ     ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]       ->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorUp  *histo_ZZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	    ->Add(histo_ZZ     ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]     ->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorDown*histo_ZZ     ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[i-1]  ->Add(histo_FakeM  ); histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorUp  *histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[i-1]->Add(histo_FakeM  ); histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorDown*histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]  ->Add(histo_FakeE  ); histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorUp  *histo_FakeE  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]->Add(histo_FakeE  ); histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorDown*histo_FakeE  ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1]  ->Add(histo_Higgs  ); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1]     ->SetBinContent(i,TMath::Max(histo_Higgs ->GetBinContent(i)+factorUp  *histo_Higgs    ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]->Add(histo_Higgs  ); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]   ->SetBinContent(i,TMath::Max(histo_Higgs ->GetBinContent(i)+factorDown*histo_Higgs    ->GetBinError(i),0.000001));
  }
  char outputLimits[200];
  sprintf(outputLimits,"MitZHAnalysis/plots%s/wz3l%2s_input_%4s.root",subdirectory.c_str(),finalStateName,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data   ->Write();
  histo_Zg     ->Write();
  histo_VVV    ->Write();
  histo_WZ     ->Write();
  histo_ZZ     ->Write();
  histo_FakeM  ->Write();
  histo_FakeE  ->Write();
  histo_Higgs  ->Write();
  cout << histo_Data   ->GetSumOfWeights() << " ";
  cout << histo_Zg     ->GetSumOfWeights() << " ";
  cout << histo_VVV    ->GetSumOfWeights() << " ";
  cout << histo_WZ     ->GetSumOfWeights() << " ";
  cout << histo_ZZ     ->GetSumOfWeights() << " ";
  cout << histo_FakeM  ->GetSumOfWeights() << " ";
  cout << histo_FakeE  ->GetSumOfWeights() << " ";
  cout << histo_Higgs  ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_Zg_CMS_MVAZHStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAZHStatBoundingUp   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVAZHStatBoundingDown      ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAZHStatBoundingDown ->GetBinContent(i)/histo_Zg	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp	     ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeM->GetBinContent(i)>0)printf("%5.1f ",histo_FakeM_CMS_MVAFakeMStatBoundingUp	     ->GetBinContent(i)/histo_FakeM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeM_CMS_MVAFakeMStatBoundingDown->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeM->GetBinContent(i)>0)printf("%5.1f ",histo_FakeM_CMS_MVAFakeMStatBoundingDown      ->GetBinContent(i)/histo_FakeM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingUp	     ->GetBinContent(i)/histo_FakeE   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeE_CMS_MVAFakeEStatBoundingDown->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingDown      ->GetBinContent(i)/histo_FakeE   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingUp	 ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingDown->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingDown	 ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffM\n");
  histo_Zg_CMS_MVALepEffMBoundingUp       ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVALepEffMBoundingDown     ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingUp	 ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffMBoundingDown  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingDown  ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffE\n");
  histo_Zg_CMS_MVALepEffEBoundingUp       ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVALepEffEBoundingDown     ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingUp	 ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffEBoundingDown  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingDown  ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties MET\n");
  histo_Zg_CMS_MVAMETBoundingUp           ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_Zg      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVAMETBoundingDown         ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAMETBoundingDown   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_Zg_CMS_MVAJESBoundingUp           ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_Zg      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVAJESBoundingDown         ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAJESBoundingDown   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingDown  ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties BTAG\n");
  histo_Zg_CMS_MVABTAGBoundingUp          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	     ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_Zg      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVABTAGBoundingDown        ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	     ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVABTAGBoundingDown   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVABTAGBoundingUp          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVABTAGBoundingUp          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVABTAGBoundingUp    ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVABTAGBoundingDown  ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties PU\n");
  histo_Zg_CMS_PUBoundingUp               ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_PUBoundingUp      ->GetBinContent(i)/histo_Zg      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_PUBoundingDown             ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_PUBoundingDown   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingUp   ->GetBinContent(i)/histo_Higgs->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_PUBoundingDown          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Higgs ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingDown ->GetBinContent(i)/histo_Higgs->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  outFileLimits->Close();

  double lumiE = 1.026;
  double systLepResE[5] = {1.01,1.01,1.01,1.01,1.01};
  double systLepResM[5] = {1.01,1.01,1.01,1.01,1.01};

  for(int nb=1; nb<=nBinMVA; nb++){
     // QCD study
    double systQCDScale[5] = {TMath::Abs(histo_Zg_CMS_QCDScaleBounding[0]   ->GetBinContent(nb)-histo_Zg   ->GetBinContent(nb)),
                              TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]  ->GetBinContent(nb)-histo_VVV  ->GetBinContent(nb)),
                              TMath::Abs(histo_WZ_CMS_QCDScaleBounding[0]   ->GetBinContent(nb)-histo_WZ   ->GetBinContent(nb)),
                              TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0]   ->GetBinContent(nb)-histo_ZZ   ->GetBinContent(nb)),
                              TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))};
    for(int nqcd=1; nqcd<6; nqcd++) {
      if(TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_Zg ->GetBinContent(nb))     > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd]   ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb));
      if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))     > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]  ->GetBinContent(nb)-histo_VVV->GetBinContent(nb));
      if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_WZ ->GetBinContent(nb))     > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]   ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb));
      if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_ZZ ->GetBinContent(nb))     > systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]   ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb));
      if(TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb)) > systQCDScale[4]) systQCDScale[4] = TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb));
    }                 
    if(histo_Zg ->GetBinContent(nb)   > 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_Zg   ->GetBinContent(nb); else systQCDScale[0] = 1;
    if(histo_VVV->GetBinContent(nb)   > 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_VVV  ->GetBinContent(nb); else systQCDScale[1] = 1;
    if(histo_WZ ->GetBinContent(nb)   > 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_WZ   ->GetBinContent(nb); else systQCDScale[2] = 1;
    if(histo_ZZ ->GetBinContent(nb)   > 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_ZZ   ->GetBinContent(nb); else systQCDScale[3] = 1;
    if(histo_Higgs->GetBinContent(nb) > 0) systQCDScale[4] = 1 + systQCDScale[4]/histo_Higgs->GetBinContent(nb); else systQCDScale[4] = 1;
    printf("QCDScale(%d): %f %f %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[4]);
    
    // PDF study
    double systPDF[5];
    histo_Diff->Reset();
    if(histo_Zg ->GetBinContent(nb) > 0) for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_Zg_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb))/histo_Zg ->GetBinContent(nb));
    systPDF[0] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_VVV->GetBinContent(nb) > 0) for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
    systPDF[1] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_WZ ->GetBinContent(nb) > 0) for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_WZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb))/histo_WZ ->GetBinContent(nb));
    systPDF[2] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_ZZ ->GetBinContent(nb) > 0) for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb))/histo_ZZ ->GetBinContent(nb));
    systPDF[3] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_Higgs->GetBinContent(nb) > 0) for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_Higgs_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))/histo_Higgs->GetBinContent(nb));
    systPDF[4] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    printf("PDF(%d): %f %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4]);

    double systLepEffM[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_Zg_CMS_MVALepEffMBoundingUp      ->GetBinContent(nb) > 0) systLepEffM[0] = histo_Zg_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_Zg_CMS_MVALepEffMBoundingDown    ->GetBinContent(nb) > 0) systLepEffM[0] = histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Zg_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)   > 0 && histo_VVV_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)   > 0 && histo_VVV_CMS_MVALepEffMBoundingDown   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_WZ_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_WZ_CMS_MVALepEffMBoundingDown	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[4] = histo_Higgs_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[4] = histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingDown->GetBinContent(nb);

    double systLepEffE[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_Zg_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb)   > 0) systLepEffE[0] = histo_Zg_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_Zg_CMS_MVALepEffEBoundingDown ->GetBinContent(nb)   > 0) systLepEffE[0] = histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Zg_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)   > 0 && histo_VVV_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb)   > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)   > 0 && histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb)   > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_WZ_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb)   > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_WZ_CMS_MVALepEffEBoundingDown ->GetBinContent(nb)   > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_ZZ_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb)   > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_ZZ_CMS_MVALepEffEBoundingDown ->GetBinContent(nb)   > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb) > 0) systLepEffE[4] = histo_Higgs_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingDown->GetBinContent(nb) > 0) systLepEffE[4] = histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingDown->GetBinContent(nb);

    double systMet[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_MVAMETBoundingUp      ->GetBinContent(nb) > 0) systMet[0] = histo_Zg_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    else if(histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_MVAMETBoundingDown    ->GetBinContent(nb) > 0) systMet[0] = histo_Zg->GetBinContent(nb)/histo_Zg_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)  > 0 && histo_VVV_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[1] = histo_VVV_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)  > 0 && histo_VVV_CMS_MVAMETBoundingDown   ->GetBinContent(nb) > 0) systMet[1] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[2] = histo_WZ_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_MVAMETBoundingDown    ->GetBinContent(nb) > 0) systMet[2] = histo_WZ->GetBinContent(nb)/histo_WZ_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[3] = histo_ZZ_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_MVAMETBoundingDown    ->GetBinContent(nb) > 0) systMet[3] = histo_ZZ->GetBinContent(nb)/histo_ZZ_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMet[4] = histo_Higgs_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMet[4] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVAMETBoundingDown->GetBinContent(nb);
    for(int i=0; i<5; i++) if(systMet[i] == 1) systMet[i] = 0.998;
    for(int nmet=0; nmet<5; nmet++) if(systMet[nmet]   > 1.10) systMet[nmet]   = 1.10;
    for(int nmet=0; nmet<5; nmet++) if(systMet[nmet]   < 0.90) systMet[nmet]   = 0.90;

    double systJes[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_MVAJESBoundingUp      ->GetBinContent(nb) > 0) systJes[0] = histo_Zg_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    else if(histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_MVAJESBoundingDown    ->GetBinContent(nb) > 0) systJes[0] = histo_Zg->GetBinContent(nb)/histo_Zg_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)  > 0 && histo_VVV_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[1] = histo_VVV_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)  > 0 && histo_VVV_CMS_MVAJESBoundingDown   ->GetBinContent(nb) > 0) systJes[1] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[2] = histo_WZ_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_MVAJESBoundingDown    ->GetBinContent(nb) > 0) systJes[2] = histo_WZ->GetBinContent(nb)/histo_WZ_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[3] = histo_ZZ_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_MVAJESBoundingDown    ->GetBinContent(nb) > 0) systJes[3] = histo_ZZ->GetBinContent(nb)/histo_ZZ_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJESBoundingUp   ->GetBinContent(nb) > 0) systJes[4] = histo_Higgs_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJes[4] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVAJESBoundingDown->GetBinContent(nb);
    for(int njes=0; njes<5; njes++) if(systJes[njes]   > 1.10) systJes[njes]   = 1.10;
    for(int njes=0; njes<5; njes++) if(systJes[njes]   < 0.90) systJes[njes]   = 0.90;

    double systBtag[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_MVABTAGBoundingUp      ->GetBinContent(nb) > 0) systBtag[0] = histo_Zg_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    else if(histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_MVABTAGBoundingDown    ->GetBinContent(nb) > 0) systBtag[0] = histo_Zg->GetBinContent(nb)/histo_Zg_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)  > 0 && histo_VVV_CMS_MVABTAGBoundingUp     ->GetBinContent(nb) > 0) systBtag[1] = histo_VVV_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)  > 0 && histo_VVV_CMS_MVABTAGBoundingDown   ->GetBinContent(nb) > 0) systBtag[1] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_MVABTAGBoundingUp	     ->GetBinContent(nb) > 0) systBtag[2] = histo_WZ_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_MVABTAGBoundingDown    ->GetBinContent(nb) > 0) systBtag[2] = histo_WZ->GetBinContent(nb)/histo_WZ_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_MVABTAGBoundingUp	     ->GetBinContent(nb) > 0) systBtag[3] = histo_ZZ_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_MVABTAGBoundingDown    ->GetBinContent(nb) > 0) systBtag[3] = histo_ZZ->GetBinContent(nb)/histo_ZZ_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVABTAGBoundingUp   ->GetBinContent(nb) > 0) systBtag[4] = histo_Higgs_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVABTAGBoundingDown ->GetBinContent(nb) > 0) systBtag[4] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    for(int nbtag=0; nbtag<5; nbtag++) if(systBtag[nbtag]   > 1.10) systBtag[nbtag]   = 1.10;
    for(int nbtag=0; nbtag<5; nbtag++) if(systBtag[nbtag]   < 0.90) systBtag[nbtag]   = 0.90;

    double systBDTMuonUp  [5] = {1.0,1.0,1.0,1.0,1.0};
    double systBDTMuonDown[5] = {1.0,1.0,1.0,1.0,1.0};
    if(histo_Zg->GetBinContent(nb)	  > 0 && histo_Zg_CMS_BDTMuonScaleBoundingUp 	->GetBinContent(nb) > 0) systBDTMuonUp  [0] = histo_Zg_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    if(histo_Zg->GetBinContent(nb)	  > 0 && histo_Zg_CMS_BDTMuonScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMuonDown[0] = histo_Zg_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMuonScaleBoundingUp 	->GetBinContent(nb) > 0) systBDTMuonUp  [1] = histo_VVV_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMuonScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMuonDown[1] = histo_VVV_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_BDTMuonScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMuonUp  [2] = histo_WZ_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_BDTMuonScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMuonDown[2] = histo_WZ_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMuonScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMuonUp  [3] = histo_ZZ_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMuonScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMuonDown[3] = histo_ZZ_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_BDTMuonScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMuonUp  [4] = histo_Higgs_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_BDTMuonScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMuonDown[4] = histo_Higgs_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
   
    double systBDTElectronUp  [5] = {1.0,1.0,1.0,1.0,1.0};
    double systBDTElectronDown[5] = {1.0,1.0,1.0,1.0,1.0};
    if(histo_Zg->GetBinContent(nb)	  > 0 && histo_Zg_CMS_BDTElectronScaleBoundingUp 	->GetBinContent(nb) > 0) systBDTElectronUp  [0] = histo_Zg_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    if(histo_Zg->GetBinContent(nb)	  > 0 && histo_Zg_CMS_BDTElectronScaleBoundingDown	->GetBinContent(nb) > 0) systBDTElectronDown[0] = histo_Zg_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTElectronScaleBoundingUp 	->GetBinContent(nb) > 0) systBDTElectronUp  [1] = histo_VVV_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTElectronScaleBoundingDown	->GetBinContent(nb) > 0) systBDTElectronDown[1] = histo_VVV_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_BDTElectronScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTElectronUp  [2] = histo_WZ_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_BDTElectronScaleBoundingDown	->GetBinContent(nb) > 0) systBDTElectronDown[2] = histo_WZ_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTElectronScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTElectronUp  [3] = histo_ZZ_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTElectronScaleBoundingDown	->GetBinContent(nb) > 0) systBDTElectronDown[3] = histo_ZZ_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_BDTElectronScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTElectronUp  [4] = histo_Higgs_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_BDTElectronScaleBoundingDown	->GetBinContent(nb) > 0) systBDTElectronDown[4] = histo_Higgs_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
   
    double systBDTMETUp  [5] = {1.0,1.0,1.0,1.0,1.0};
    double systBDTMETDown[5] = {1.0,1.0,1.0,1.0,1.0};
    if(histo_Zg->GetBinContent(nb)	  > 0 && histo_Zg_CMS_BDTMETScaleBoundingUp 	->GetBinContent(nb) > 0) systBDTMETUp  [0] = histo_Zg_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    if(histo_Zg->GetBinContent(nb)	  > 0 && histo_Zg_CMS_BDTMETScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMETDown[0] = histo_Zg_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMETScaleBoundingUp 	->GetBinContent(nb) > 0) systBDTMETUp  [1] = histo_VVV_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMETScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMETDown[1] = histo_VVV_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_BDTMETScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMETUp  [2] = histo_WZ_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_BDTMETScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMETDown[2] = histo_WZ_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMETScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMETUp  [3] = histo_ZZ_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMETScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMETDown[3] = histo_ZZ_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_BDTMETScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMETUp  [4] = histo_Higgs_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_BDTMETScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMETDown[4] = histo_Higgs_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
   
    double systPU[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_PUBoundingUp     ->GetBinContent(nb) > 0) systPU[0] = histo_Zg_CMS_PUBoundingUp->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    else if(histo_Zg->GetBinContent(nb)   > 0 && histo_Zg_CMS_PUBoundingDown   ->GetBinContent(nb) > 0) systPU[0] = histo_Zg->GetBinContent(nb)/histo_Zg_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)  > 0 && histo_VVV_CMS_PUBoundingUp    ->GetBinContent(nb) > 0) systPU[1] = histo_VVV_CMS_PUBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb  )> 0 && histo_VVV_CMS_PUBoundingDown  ->GetBinContent(nb) > 0) systPU[1] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_PUBoundingUp     ->GetBinContent(nb) > 0) systPU[2] = histo_WZ_CMS_PUBoundingUp->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb)   > 0 && histo_WZ_CMS_PUBoundingDown   ->GetBinContent(nb) > 0) systPU[2] = histo_WZ->GetBinContent(nb)/histo_WZ_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_PUBoundingUp     ->GetBinContent(nb) > 0) systPU[3] = histo_ZZ_CMS_PUBoundingUp->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb)   > 0 && histo_ZZ_CMS_PUBoundingDown   ->GetBinContent(nb) > 0) systPU[3] = histo_ZZ->GetBinContent(nb)/histo_ZZ_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_PUBoundingUp  ->GetBinContent(nb) > 0) systPU[4] = histo_Higgs_CMS_PUBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_PUBoundingDown->GetBinContent(nb) > 0) systPU[4] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_PUBoundingDown->GetBinContent(nb);
    for(int npu=0; npu<5; npu++) if(systPU[npu] > 1.02) systPU[npu] = 1.02;
    for(int npu=0; npu<5; npu++) if(systPU[npu] < 0.98) systPU[npu] = 0.98;

    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"MitZHAnalysis/datacards%s/histo_limits_wz3l%2s_%4s_bin%d.txt",subdirectory.c_str(), finalStateName,ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
    newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    newcardShape << Form("process Zg VVV WZ ZZ FakeM FakeE Higgs\n");
    if(histo_Higgs->GetSumOfWeights() > 0)
    newcardShape << Form("process 1 2 6 3 4 5 0\n");
    else if(isWZhinv)
    newcardShape << Form("process 1 2 6 3 4 5 7\n");
    else
    newcardShape << Form("process 1 2 0 3 4 5 6\n");
    newcardShape << Form("rate %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",TMath::Max(histo_Zg->GetBinContent(nb),0.0),TMath::Max(histo_VVV->GetBinContent(nb),0.0),TMath::Max(histo_WZ->GetBinContent(nb),0.0)*theWZSF,TMath::Max(histo_ZZ->GetBinContent(nb),0.0),TMath::Max(histo_FakeM->GetBinContent(nb),0.0),TMath::Max(histo_FakeE->GetBinContent(nb),0.0),TMath::Max(histo_Higgs->GetBinContent(nb),0.0));
    newcardShape << Form("lumi_%4s                               lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE); 		
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3],systLepEffM[4]);
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3],systLepEffE[4]);
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3],systLepResM[4]);
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3],systLepResE[4]);
    newcardShape << Form("CMS_pu2016                             lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",systPU[0],systPU[1],systPU[2],systPU[3],systPU[4]);
    newcardShape << Form("CMS_scale_met                          lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",systMet[0],systMet[1],systMet[2],systMet[3],systMet[4]);
    newcardShape << Form("CMS_scale_j                            lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",systJes[0],systJes[1],systJes[2],systJes[3],systJes[4]);
    newcardShape << Form("CMS_trigger2016                        lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",1.01,1.01,1.01,1.01,1.01);
    newcardShape << Form("CMS_eff_b_2016                         lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",systBtag[0],systBtag[1],systBtag[2],systBtag[3],systBtag[4]);
    newcardShape << Form("pdf_qqbar                              lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  %7.5f\n",TMath::Max(systPDF[0],1.01),TMath::Max(systPDF[1],1.01),TMath::Max(systPDF[2],1.01),TMath::Max(systPDF[3],1.01),TMath::Max(systPDF[4],1.01));
    newcardShape << Form("QCDscale_VVV		                 lnN    -     %7.5f   -     -	  -    -    -  \n",systQCDScale[1]); 	   
    newcardShape << Form("QCDscale_Zg		                 lnN  %7.5f     -     -     -     -    -    -  \n",systQCDScale[0]); 	   
    newcardShape << Form("QCDscale_WZ		                 lnN    -       -   %7.5f   -     -    -    -  \n",systQCDScale[2]); 	   
    newcardShape << Form("QCDscale_ZZ		                 lnN    -       -     -   %7.5f   -    -    -  \n",systQCDScale[3]); 	   
    newcardShape << Form("QCDscale_ggH		                 lnN    -       -     -     -     -    -  %7.5f\n",systQCDScale[4]); 	   
    newcardShape << Form("CMS_wz3l_FakeMSyst_%4s                 lnN    -       -     -	    -   %7.5f  -    -  \n",ECMsb.Data(),1.30);  	
    newcardShape << Form("CMS_wz3l_FakeESyst_%4s                 lnN    -       -     -	    -     -  %7.5f  -  \n",ECMsb.Data(),1.30);  	
    if(MVAVarType == 3 || MVAVarType==4) {
    newcardShape << Form("CMS_BDT_scale_m                        lnN	%7.5f/%7.5f  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f - -  %7.5f/%7.5f   \n", systBDTMuonUp[0], systBDTMuonDown[0], systBDTMuonUp[1], systBDTMuonDown[1], systBDTMuonUp[2], systBDTMuonDown[2], systBDTMuonUp[3], systBDTMuonDown[3], systBDTMuonUp[4], systBDTMuonDown[4]);	    
    newcardShape << Form("CMS_BDT_scale_e                        lnN	%7.5f/%7.5f  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f - -  %7.5f/%7.5f   \n", systBDTElectronUp[0], systBDTElectronDown[0], systBDTElectronUp[1], systBDTElectronDown[1], systBDTElectronUp[2], systBDTElectronDown[2], systBDTElectronUp[3], systBDTElectronDown[3], systBDTElectronUp[4], systBDTElectronDown[4]);	    
    newcardShape << Form("CMS_BDT_scale_MET                      lnN	%7.5f/%7.5f  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f - -  %7.5f/%7.5f   \n", systBDTMETUp[0], systBDTMETDown[0], systBDTMETUp[1], systBDTMETDown[1], systBDTMETUp[2], systBDTMETDown[2], systBDTMETUp[3], systBDTMETDown[3], systBDTMETUp[4], systBDTMETDown[4]);	    
    }
    if(histo_Zg->GetBinContent(nb)        > 0) newcardShape << Form("CMS_wz3l%s_MVAZgStatBounding_%s_Bin%d	  lnN    %7.5f   -    -    -	-    -    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_Zg    ->GetBinError(nb)/histo_Zg    ->GetBinContent(nb));
    if(histo_VVV->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%d       lnN      -  %7.5f   -    -	-    -    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_VVV   ->GetBinError(nb)/histo_VVV   ->GetBinContent(nb));
    if(histo_WZ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAWZStatBounding_%s_Bin%d	  lnN      -     -  %7.5f  -	-    -    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_WZ    ->GetBinError(nb)/histo_WZ    ->GetBinContent(nb));
    if(histo_ZZ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAZZStatBounding_%s_Bin%d	  lnN      -     -    -  %7.5f  -    -    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ZZ    ->GetBinError(nb)/histo_ZZ    ->GetBinContent(nb));
    if(histo_FakeM->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%d	  lnN      -     -    -    -  %7.5f  -    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_FakeM ->GetBinError(nb)/histo_FakeM ->GetBinContent(nb));
    if(histo_FakeE->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAFameEStatBounding_%s_Bin%d	  lnN      -     -    -    -	-  %7.5f  -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_FakeE ->GetBinError(nb)/histo_FakeE ->GetBinContent(nb));
    if(histo_Higgs->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%d     lnN      -     -    -    -	-    -  %7.5f\n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_Higgs ->GetBinError(nb)/histo_Higgs ->GetBinContent(nb));

    if(nb != 1 && isWZhinv){
      if(useZZWZEWKUnc){
      newcardShape << Form("CMS_hinv_vvnorm_bin%d rateParam  * WZ 1 [0.1,10]\n",nb-1);
      } else {
      newcardShape << Form("CMS_hinv_wznorm_bin%d rateParam  * WZ 1 [0.1,10]\n",nb-1);
      }
    }

    newcardShape.close();
  }
}
