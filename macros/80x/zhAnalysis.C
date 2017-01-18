#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom.h>
#include <iostream>
#include <fstream>
#include "TMVA/Reader.h"

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"

#include "MitAnalysisRunII/macros/80x/factors.h"
#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"
#include "MitZHAnalysis/macros/80x/zhMVA.h"

// 0 == sm, 7 == mh500, 24 = A_Mx-150_Mv-500, 55 == V_Mx-150_Mv-500

bool isMINIAOD = true;
int whichSkim = 4;
bool useZjetsTemplate = false;
bool usePureMC = true; 
bool useEMFromData = true;
bool useVVFromData = true;
bool useZZWZEWKUnc = true;
bool useCachedBDTSystematics=false;
double mcPrescale = 1.;
enum selType                     {ZSEL=0,  SIGSEL,   WWSEL,   WWLOOSESEL,   BTAGSEL,   WZSEL,   PRESEL,   CR1SEL,   CR2SEL,   CR12SEL,   TIGHTSEL,   DYSANESEL1,   DYSANESEL2,  nSelTypes};
TString selTypeName[nSelTypes]= {"ZSEL",  "SIGSEL", "WWSEL", "WWLOOSESEL", "BTAGSEL", "WZSEL", "PRESEL", "CR1SEL", "CR2SEL", "CR12SEL", "TIGHTSEL", "DYSANESEL1", "DYSANESEL2"};
enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};
const TString typeLepSel = "medium";

void zhAnalysis(
 unsigned int nJetsType = 1,
 bool isBlinded = false,
 Int_t typeSel = 3,
 Int_t plotModel = 0,
 bool verbose = true,
 string the_BDT_weights="",
 string subdirectory=""
 ){
  if(subdirectory!="" && subdirectory.c_str()[0]!='/') subdirectory = "/"+subdirectory;
  system(("mkdir -p MitZHAnalysis/datacards"+subdirectory).c_str());
  system(("mkdir -p MitZHAnalysis/plots"+subdirectory).c_str());
  bool makeMVAtrees=false;
  bool useBDT=false;
  if(makeMVAtrees) system("mkdir -p MitZHAnalysis/mva");
  Int_t period = 1;
  TString filesPathDA = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/met_";
  TString filesPathMC  = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/met_";
  TString filesPathDMMC  = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/";
  Double_t lumi = 36.5;
  TString processTag = "";

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev, signalName_;
  vector<Int_t> infilecatv, signalIndex_;  

  TString puPath = "";
  TString zjetsTemplatesPath = "";
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
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg.root",filesPathMC.Data()));					      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));infilecatv.push_back(1);

  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root",filesPathMC.Data()));	              infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                            infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data())); 	      infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));           infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infilecatv.push_back(1);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));        infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(2);
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
  infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root",filesPathMC.Data())); 		      infilecatv.push_back(5);

  for(int ifile=0; ifile<(int)infilenamev.size(); ifile++) {
    signalIndex_.push_back(-1); // Populate vector of signal indices with -1 for the non-MC-signal files
  }

  // Monte Carlo signals
  { // Model 0: standard model Higgs (125) with glu-glu
    int mH=125;
    signalName_.push_back("sm");
    infilenamev.push_back(Form("%sZH_ZToMM_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH)); infilecatv.push_back(6); signalIndex_.push_back(0);
    infilenamev.push_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH)); infilecatv.push_back(6); signalIndex_.push_back(0);
    infilenamev.push_back(Form("%sggZH_HToInv_ZToLL_M125_13TeV_powheg_pythia8.root",filesPathDMMC.Data()));       infilecatv.push_back(7); signalIndex_.push_back(0);
  }  // Models 1 thru 8: standard-model-like Higgs mass points without glu-glu (8 models)
  { int mH_[10]={110, 125, 150, 200, 300, 400, 500, 600, 800, 1000}; int iH=0; for(int i=1; i<=10; i++) { int mH = mH_[iH]; 
    signalName_.push_back(Form("mh%d", mH));
    infilenamev.push_back(Form("%sZH_ZToMM_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH)); 
    infilenamev.push_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH));
    infilecatv.push_back(6); signalIndex_.push_back(iH+1);
    infilecatv.push_back(6); signalIndex_.push_back(iH+1);
    iH++;
  }}

  { // dark matter models () ls -l /scratch5/ceballos/ntuples_weightsMC_80x/|grep DarkMatter_MonoZToLL|awk '{printf("    signalName\_.push_back(\"%s\"); infileName\_.push_back(Form(\"%s\", filesPathDMMC.Data())); infileCategory\_.push\_back(6); signalIndex\_.push\_back(i); i++;\n",$9,$9)}'
    int i=signalName_.size();
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-1000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-1995_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-1995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-100_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-20_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-50_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-200_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-500_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-100_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-2000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-200_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-20_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-300_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-500_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-2000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-500_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-995_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-200_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-300_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-50_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-95_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-95_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-1000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-1995_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-1995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-100_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-20_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-200_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-295_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-500_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-1000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-100_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-10_1-gDMgQ"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-10_1-gDMgQ_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-2000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-200_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-20_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-300_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-50_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-2000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-995_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-10_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-200_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-300_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-5000_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-50_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-95_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-95_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
  }
/*
  { // ls -l /scratch5/ceballos/ntuples_weightsMC_80x/|grep monoz_med|awk '{printf("    signalName\_.push_back(\"%s\"); infileName\_.push_back(Form(\"%s\", filesPathDMMC.Data())); infileCategory\_.push\_back(6); signalIndex\_.push\_back(i); i++;\n",$9,$9)}'
    int i=signalName_.size();
    signalName_.push_back("pseudoscalarmonoz_med-10_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-10_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-110_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-110_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-160_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-160_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-210_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-210_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-260_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-260_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-310_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-310_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-360_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-360_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-410_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-410_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-500_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-500_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-600_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-600_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-60_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-60_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("pseudoscalarmonoz_med-700_dm-50"); infilenamev.push_back(Form("%spseudoscalarmonoz_med-700_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-1000_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-1000_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-100_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-100_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-10_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-10_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-110_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-110_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-1300_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-1300_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-1500_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-1500_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-160_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-160_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-1800_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-1800_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-2000_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-2000_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-210_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-210_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-25_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-25_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-260_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-260_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-300_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-300_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-410_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-410_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-600_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-600_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-600_dm-50_gq-0.25"); infilenamev.push_back(Form("%sscalarmonoz_med-600_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-60_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-60_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("scalarmonoz_med-700_dm-50"); infilenamev.push_back(Form("%sscalarmonoz_med-700_dm-50.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-1000_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-1000_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-100_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-100_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-1300_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-1300_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-1500_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-1500_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-1800_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-1800_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-2000_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-2000_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-300_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-300_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
    signalName_.push_back("vectormonoz_med-600_dm-50_gq-0.25"); infilenamev.push_back(Form("%svectormonoz_med-600_dm-50_gq-0.25.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
  }
*/

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}
  
  //signalName_.clear();infilenamev.clear();infilecatv.clear();signalIndex_.clear();int i = 0;
  //signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1"); infilenamev.push_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data())); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
  //signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-2"); infilenamev.push_back(Form("/scratch5/ceballos/test/NeroNtuples_skim.root")); infilecatv.push_back(6); signalIndex_.push_back(i); i++;
  //infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8.root",filesPathMC.Data()));  				   infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathMC.Data(),600)); infilecatv.push_back(5);
  //infilenamev.push_back(Form("/home/ceballos/cms/hist/tt_all/t2mit/filefi/044/ZZTo2L2Nu_13TeV_powheg_pythia8.root")); infilecatv.push_back(5);
  //infilenamev.push_back(Form("/home/ceballos/cms/hist/tt_all/t2mit/filefi/044/ZZTo2L2Nu_13TeV_powheg_pythia8_JSON.root")); infilecatv.push_back(5);

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

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;

  TFile *fTrackElectronReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_13p0ifb.root"));
  TH2D *fhDeltrksf= (TH2D*)(fTrackElectronReco_SF->Get("scalefactors_Reco_Electron")); assert(fhDeltrksf); fhDeltrksf->SetDirectory(0);
  delete fTrackElectronReco_SF;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_13p0ifb.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_Medium_Electron"));
  TH2D *fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF->SetDirectory(0);
  delete fElSF;

  TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/trackMuReco_SF.root"));
  TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptg10")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
  //TH1D *fhDmutrksfptl10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptl10")); assert(fhDmutrksfptl10); fhDmutrksfptl10->SetDirectory(0);
  delete fTrackMuonReco_SF;

  //TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x.root"));
  //TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("scalefactors_Tight_Muon")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonID_Z_RunBCD_prompt80X_7p65.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  delete fMuSF;

  TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonIso_Z_RunBCD_prompt80X_7p65.root"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  delete fMuIsoSF;

  // Dilepton trigger efficiencies
  TFile *fEffDilepTrigs = TFile::Open(Form("MitAnalysisRunII/data/80x/dilepton_trigger_efficiencies_80x.root")); 
  TH2D *fhDEffDimuonSoup_pt0     = (TH2D*)(fEffDilepTrigs->Get("h_dimuon_eff_0"));
  TH2D *fhDEffDimuonSoup_pt1     = (TH2D*)(fEffDilepTrigs->Get("h_dimuon_eff_1"));
  TH2D *fhDEffDimuonSoup_pt2     = (TH2D*)(fEffDilepTrigs->Get("h_dimuon_eff_2"));
  TH2D *fhDEffDimuonSoup_pt3     = (TH2D*)(fEffDilepTrigs->Get("h_dimuon_eff_3"));
  TH2D *fhDEffDielectronSoup_pt0 = (TH2D*)(fEffDilepTrigs->Get("h_dielectron_eff_0"));
  TH2D *fhDEffDielectronSoup_pt1 = (TH2D*)(fEffDilepTrigs->Get("h_dielectron_eff_1"));
  TH2D *fhDEffDielectronSoup_pt2 = (TH2D*)(fEffDilepTrigs->Get("h_dielectron_eff_2"));
  TH2D *fhDEffDielectronSoup_pt3 = (TH2D*)(fEffDilepTrigs->Get("h_dielectron_eff_3"));
  assert(fhDEffDimuonSoup_pt0    );  
  assert(fhDEffDimuonSoup_pt1    ); 
  assert(fhDEffDimuonSoup_pt2    ); 
  assert(fhDEffDimuonSoup_pt3    ); 
  assert(fhDEffDielectronSoup_pt0); 
  assert(fhDEffDielectronSoup_pt1); 
  assert(fhDEffDielectronSoup_pt2); 
  assert(fhDEffDielectronSoup_pt3); 
  fhDEffDimuonSoup_pt0     ->SetDirectory(0); 
  fhDEffDimuonSoup_pt1     ->SetDirectory(0); 
  fhDEffDimuonSoup_pt2     ->SetDirectory(0); 
  fhDEffDimuonSoup_pt3     ->SetDirectory(0); 
  fhDEffDielectronSoup_pt0 ->SetDirectory(0); 
  fhDEffDielectronSoup_pt1 ->SetDirectory(0); 
  fhDEffDielectronSoup_pt2 ->SetDirectory(0); 
  fhDEffDielectronSoup_pt3 ->SetDirectory(0); 
  delete fEffDilepTrigs;

  TString ECMsb  = "13TeV2016";
  
  // MVA variable types:
  // 1: MET only
  // 2: MET x mll
  // 3: classifier only
  // 4: MET x classifier

  //const int MVAVarType = 0; const int nBinMVA = 8; Double_t xbins[nBinMVA+1] = {0, 50, 200, 250, 300, 400, 600, 800, 1000}; TString addChan = "";
  //const int MVAVarType = 0; const int nBinMVA = 14; Double_t xbins[nBinMVA+1] = {0, 50, 200, 225, 250, 275, 300, 350, 400, 500, 600, 700, 800, 900, 1000}; TString addChan = "";
  //const int MVAVarType = 1; const int nBinMVA = 8; Double_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 350}; TString addChan = "1";
  // const int MVAVarType = 1; const int nBinMVA = 12; Double_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600}; TString addChan = "1";
  //const int MVAVarType = 1; const int nBinMVA = 11; Double_t xbins[nBinMVA+1] = {0, 50, 100, 160, 240, 320, 400, 480, 560, 640, 800, 1200}; TString addChan = "1";
  //const int MVAVarType = 1; const int nBinMVA = 21; Double_t xbins[nBinMVA+1] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1200}; TString addChan = "1";
  //const int MVAVarType = 2; const int nBinMVA = 20; Double_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 350,
  //                                                                                         1125,1150,1175,1200,1250,1350,
  //											     2125,2150,2175,2200,2250,2350}; TString addChan = "2";
  const int MVAVarType = 3; const int nBinMVA = 15; Double_t xbins[nBinMVA+1] =  {-2, -1, 0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.4}; TString addChan = "3";
  //const int MVAVarType = 4; const int nBinMVA = 26; Double_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 350,
  //                                                                                         1125,1150,1175,1200,1250,1350,
  //                                                                                         2125,2150,2175,2200,2250,2350,
  //                                                                                         3125,3150,3175,3200,3250,3350}; TString addChan = "4";
  
  if (MVAVarType==3 || MVAVarType==4) useBDT=true;

  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  if     (MVAVarType == 0) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt50_mt.root";
  //if     (MVAVarType == 0) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt50_mt_13bins.root";
  else if(MVAVarType == 1) zjetsTemplatesPath = "MitZHAnalysis/data/80x/zjets_13TeV_25ns_metgt50_met.root";
  else useZjetsTemplate = false;

  TH1D *fhDZjets;
  TH1D *fhDZjetsSyst;
  if(useZjetsTemplate){
    TFile *fZjetsTemplatesFile = TFile::Open(Form("%s",zjetsTemplatesPath.Data()));

    if     (nJetsType == 0) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets0"));
    else if(nJetsType == 1) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
    else if(nJetsType == 2) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets2"));
    assert(fhDZjets);
    fhDZjets->SetDirectory(0);

    if     (nJetsType == 0) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
    else if(nJetsType == 1) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets0"));
    else if(nJetsType == 2) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
    assert(fhDZjetsSyst);
    fhDZjetsSyst->SetDirectory(0);

    delete fZjetsTemplatesFile;
  }

  TFile *fVJetsKfactorFile = TFile::Open(Form("MitAnalysisRunII/data/80x/kfactors_vjets.root"));
  TH1D *fhDVjetsNum = (TH1D*)(fVJetsKfactorFile->Get("EWKcorr/Z"));
  TH1D *fhDVjetsDen = (TH1D*)(fVJetsKfactorFile->Get("ZJets_LO/inv_pt"));
  assert(fhDVjetsNum);
  assert(fhDVjetsDen);
  fhDVjetsNum->SetDirectory(0);
  fhDVjetsDen->SetDirectory(0);
  delete fVJetsKfactorFile;

  const int numberCuts = 11;
  TH1D* histoZHSEL[4];
  histoZHSEL[0] = new TH1D("histoZHSEL_0", "histoZHSEL_0", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[1] = new TH1D("histoZHSEL_1", "histoZHSEL_1", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[2] = new TH1D("histoZHSEL_2", "histoZHSEL_2", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[3] = new TH1D("histoZHSEL_3", "histoZHSEL_3", numberCuts+1, -0.5, numberCuts+0.5);
  TString cutName[numberCuts+1] = {"ptl>20/20","3rd lepton veto","btag-veto","tauVeto","Njets","Z mass","ptll>60","MET>100","dPhi(Z-MET)>2.8","|ptll-MET|/ptll<0.4","dPhiJetMet>0.5","all"};

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 39;
  const int histBins = 8;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "...DY", "...WZ", "....ZZ", "...VVV", "....ZH", "..ggZH"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot ==  0) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =1000.0;}
    else if(thePlot ==  1) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot ==  2) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot ==  3) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot ==  4) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot ==  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot ==  6) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot ==  7) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot ==  8) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot == 10) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot == 11) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot == 12) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 13) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 14) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot == 15) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot == 16) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot == 17) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot == 18) {nBinPlot = 500; xminPlot =-0.5; xmaxPlot = 499.5;}
    else if(thePlot == 19) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot == 20) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   2.5;}
    else if(thePlot == 21) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =1000.0;}
    else if(thePlot == 22) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot == 23) {nBinPlot =  32; xminPlot =-0.1; xmaxPlot =   3.1;} // Delta phi jet met agreement
    else if(thePlot == 24) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =   2.0;} // Calo and PF met agreement
    else if(thePlot == 25) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;} // MET
    else if(thePlot == 26) {nBinPlot =  32; xminPlot =-0.1; xmaxPlot =   3.1;} // Delta phi jet met agreement
    else if(thePlot == 27) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =   2.0;} // Calo and PF met agreement
    else if(thePlot == 28) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;} // MET
    else if(thePlot == 29) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;} // Jet multiplicity
    else if(thePlot == 30) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;} // B-tagged jets
    else if(thePlot == 31) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;} // Lepton multiplicity
    else if(thePlot == 32) {nBinPlot =  32; xminPlot =-0.1; xmaxPlot =   3.1;} // Delta phi Z MET
    else if(thePlot == 33) {nBinPlot =  60; xminPlot =40.0; xmaxPlot = 100.0;}
    else if(thePlot == 34) {nBinPlot =  60; xminPlot = 0.0; xmaxPlot =   3.0;}
    else if(thePlot == 35) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot == 36) {nBinPlot = 100; xminPlot =-1.0; xmaxPlot =   1.0;}
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
  TH1D* histo_ZH_hinv_CMS_PDFBounding[nSigModels][102];
  TH1D* histo_VVV_CMS_PDFBounding[102];
  TH1D* histo_WZ_CMS_PDFBounding[102];
  TH1D* histo_ZZ_CMS_PDFBounding[102];
  TH1D* histo_ggZH_hinv_CMS_PDFBounding[102];
  for(int nb=0; nb<102; nb++) {
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

  TH1D* histo_VVV_CMS_MVAJESBoundingUp      	= new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown    	= new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_WZ_CMS_scale_jUp")  , Form("histo_WZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_WZ_CMS_scale_jDown"), Form("histo_WZ_CMS_scale_jDown"), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_ZZ_CMS_scale_jUp")  , Form("histo_ZZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_ZZ_CMS_scale_jDown"), Form("histo_ZZ_CMS_scale_jDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jUp")  , Form("histo_ggZH_hinv_CMS_scale_jUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jDown"), Form("histo_ggZH_hinv_CMS_scale_jDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_BDTMuonScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_muonUp")  , Form("histo_VVV_CMS_bdt_muonUp")  , nBinMVA, xbins);histo_VVV_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMuonScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_muonDown"), Form("histo_VVV_CMS_bdt_muonDown"), nBinMVA, xbins);histo_VVV_CMS_BDTMuonScaleBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_BDTMuonScaleBoundingUp    	= new TH1D( Form("histo_WZ_CMS_bdt_muonUp")  , Form("histo_WZ_CMS_bdt_muonUp")  , nBinMVA, xbins);  histo_WZ_CMS_BDTMuonScaleBoundingUp   ->Sumw2();
  TH1D* histo_WZ_CMS_BDTMuonScaleBoundingDown  	= new TH1D( Form("histo_WZ_CMS_bdt_muonDown"), Form("histo_WZ_CMS_bdt_muonDown"), nBinMVA, xbins);  histo_WZ_CMS_BDTMuonScaleBoundingDown ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_muonUp")  , Form("histo_ZZ_CMS_bdt_muonUp")  , nBinMVA, xbins);  histo_ZZ_CMS_BDTMuonScaleBoundingUp   ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_muonDown"), Form("histo_ZZ_CMS_bdt_muonDown"), nBinMVA, xbins);  histo_ZZ_CMS_BDTMuonScaleBoundingDown ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTMuonScaleBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_muonUp")  , Form("histo_ggZH_hinv_CMS_bdt_muonUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTMuonScaleBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_muonDown"), Form("histo_ggZH_hinv_CMS_bdt_muonDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMuonScaleBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_BDTMuonScaleBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_BDTMuonScaleBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_BDTElectronScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_electronUp")  , Form("histo_VVV_CMS_bdt_electronUp")  , nBinMVA, xbins);             histo_VVV_CMS_BDTElectronScaleBoundingUp	 ->Sumw2();
  TH1D* histo_VVV_CMS_BDTElectronScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_electronDown"), Form("histo_VVV_CMS_bdt_electronDown"), nBinMVA, xbins);             histo_VVV_CMS_BDTElectronScaleBoundingDown	 ->Sumw2();
  TH1D* histo_WZ_CMS_BDTElectronScaleBoundingUp    	= new TH1D( Form("histo_WZ_CMS_bdt_electronUp")  , Form("histo_WZ_CMS_bdt_electronUp")  , nBinMVA, xbins);               histo_WZ_CMS_BDTElectronScaleBoundingUp	 ->Sumw2();
  TH1D* histo_WZ_CMS_BDTElectronScaleBoundingDown  	= new TH1D( Form("histo_WZ_CMS_bdt_electronDown"), Form("histo_WZ_CMS_bdt_electronDown"), nBinMVA, xbins);               histo_WZ_CMS_BDTElectronScaleBoundingDown	 ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_electronUp")  , Form("histo_ZZ_CMS_bdt_electronUp")  , nBinMVA, xbins);               histo_ZZ_CMS_BDTElectronScaleBoundingUp	 ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_electronDown"), Form("histo_ZZ_CMS_bdt_electronDown"), nBinMVA, xbins);               histo_ZZ_CMS_BDTElectronScaleBoundingDown	 ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTElectronScaleBoundingUp  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_electronUp")  , Form("histo_ggZH_hinv_CMS_bdt_electronUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTElectronScaleBoundingDown= new TH1D( Form("histo_ggZH_hinv_CMS_bdt_electronDown"), Form("histo_ggZH_hinv_CMS_bdt_electronDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTElectronScaleBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_BDTElectronScaleBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_BDTElectronScaleBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_BDTMETScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_METUp")  , Form("histo_VVV_CMS_bdt_METUp")  , nBinMVA, xbins); histo_VVV_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMETScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_METDown"), Form("histo_VVV_CMS_bdt_METDown"), nBinMVA, xbins); histo_VVV_CMS_BDTMETScaleBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_BDTMETScaleBoundingUp    	= new TH1D( Form("histo_WZ_CMS_bdt_METUp")  , Form("histo_WZ_CMS_bdt_METUp")  , nBinMVA, xbins);   histo_WZ_CMS_BDTMETScaleBoundingUp	->Sumw2();
  TH1D* histo_WZ_CMS_BDTMETScaleBoundingDown  	= new TH1D( Form("histo_WZ_CMS_bdt_METDown"), Form("histo_WZ_CMS_bdt_METDown"), nBinMVA, xbins);   histo_WZ_CMS_BDTMETScaleBoundingDown ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMETScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_METUp")  , Form("histo_ZZ_CMS_bdt_METUp")  , nBinMVA, xbins);   histo_ZZ_CMS_BDTMETScaleBoundingUp	->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMETScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_METDown"), Form("histo_ZZ_CMS_bdt_METDown"), nBinMVA, xbins);   histo_ZZ_CMS_BDTMETScaleBoundingDown ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTMETScaleBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_METUp")  , Form("histo_ggZH_hinv_CMS_bdt_METUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTMETScaleBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_METDown"), Form("histo_ggZH_hinv_CMS_bdt_METDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMETScaleBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_BDTMETScaleBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_BDTMETScaleBoundingDown[nSigModels];

  TH1D* histo_VVV_CMS_BDTJetScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_JESUp")  , Form("histo_VVV_CMS_bdt_JESUp")  , nBinMVA, xbins); histo_VVV_CMS_BDTJetScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTJetScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_JESDown"), Form("histo_VVV_CMS_bdt_JESDown"), nBinMVA, xbins); histo_VVV_CMS_BDTJetScaleBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_BDTJetScaleBoundingUp    	= new TH1D( Form("histo_WZ_CMS_bdt_JESUp")  , Form("histo_WZ_CMS_bdt_JESUp")  , nBinMVA, xbins);   histo_WZ_CMS_BDTJetScaleBoundingUp	->Sumw2();
  TH1D* histo_WZ_CMS_BDTJetScaleBoundingDown  	= new TH1D( Form("histo_WZ_CMS_bdt_JESDown"), Form("histo_WZ_CMS_bdt_JESDown"), nBinMVA, xbins);   histo_WZ_CMS_BDTJetScaleBoundingDown ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTJetScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_JESUp")  , Form("histo_ZZ_CMS_bdt_JESUp")  , nBinMVA, xbins);   histo_ZZ_CMS_BDTJetScaleBoundingUp	->Sumw2();
  TH1D* histo_ZZ_CMS_BDTJetScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_JESDown"), Form("histo_ZZ_CMS_bdt_JESDown"), nBinMVA, xbins);   histo_ZZ_CMS_BDTJetScaleBoundingDown ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTJetScaleBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_JESUp")  , Form("histo_ggZH_hinv_CMS_bdt_JESUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTJetScaleBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_BDTJetScaleBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_JESDown"), Form("histo_ggZH_hinv_CMS_bdt_JESDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTJetScaleBoundingDown->Sumw2();
  TH1D* histo_ZH_hinv_CMS_BDTJetScaleBoundingUp[nSigModels];  
  TH1D* histo_ZH_hinv_CMS_BDTJetScaleBoundingDown[nSigModels];

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

  TH1D* histo_WZ_CMS_EWKCorrUp    	        = new TH1D( Form("histo_WZ_EWKCorrUp")  , Form("histo_WZ_EWKCorrUp")  , nBinMVA, xbins); histo_WZ_CMS_EWKCorrUp  ->Sumw2();
  TH1D* histo_WZ_CMS_EWKCorrDown                = new TH1D( Form("histo_WZ_EWKCorrDown"), Form("histo_WZ_EWKCorrDown"), nBinMVA, xbins); histo_WZ_CMS_EWKCorrDown->Sumw2();
  TH1D* histo_ZZ_CMS_EWKCorrUp                  = new TH1D( Form("histo_ZZ_EWKCorrUp")  , Form("histo_ZZ_EWKCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_EWKCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_EWKCorrDown                = new TH1D( Form("histo_ZZ_EWKCorrDown"), Form("histo_ZZ_EWKCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_EWKCorrDown->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrUp                   = new TH1D( Form("histo_ZZ_ggCorrUp")  , Form("histo_ZZ_ggCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_ggCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrDown                 = new TH1D( Form("histo_ZZ_ggCorrDown"), Form("histo_ZZ_ggCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_ggCorrDown->Sumw2();
  TH1D* histo_Zjets_CMS_ZjetsSystUp    	        = new TH1D( Form("histo_Zjets_ZjetsSystUp")  , Form("histo_Zjets_ZjetsSystUp")  , nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_ZjetsSystDown           = new TH1D( Form("histo_Zjets_ZjetsSystDown"), Form("histo_Zjets_ZjetsSystDown"), nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystDown->Sumw2();

  for(int nModel=0; nModel<nSigModels; nModel++) { 
    histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]     = new TH1D( Form("histo_ZH_hinv_%s_%sUp",   signalName_[nModel].Data(), effMName), Form("histo_ZH_hinv_%s_%sUp",  signalName_[nModel].Data(), effMName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]   = new TH1D( Form("histo_ZH_hinv_%s_%sDown", signalName_[nModel].Data(), effMName), Form("histo_ZH_hinv_%s_%sDown",signalName_[nModel].Data(), effMName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffMBoundingAvg [nModel]           = new TH1D( Form("histo_ZH_hinv_%s_%sAvg",             signalName_[nModel].Data(), effMName), Form("histo_ZH_hinv_%s_%sAvg" ,          signalName_[nModel].Data(), effMName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingUp [nModel]            = new TH1D( Form("histo_ZH_hinv_%s_%sUp",              signalName_[nModel].Data(), effEName), Form("histo_ZH_hinv_%s_%sUp"  ,          signalName_[nModel].Data(), effEName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingDown [nModel]          = new TH1D( Form("histo_ZH_hinv_%s_%sDown",            signalName_[nModel].Data(), effEName), Form("histo_ZH_hinv_%s_%sDown",          signalName_[nModel].Data(), effEName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingAvg [nModel]           = new TH1D( Form("histo_ZH_hinv_%s_%sAvg",             signalName_[nModel].Data(), effEName), Form("histo_ZH_hinv_%s_%sAvg" ,          signalName_[nModel].Data(), effEName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAMETBoundingUp [nModel]                = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_metUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_metUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAMETBoundingDown [nModel]              = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_metDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_metDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVAJESBoundingUp [nModel]                = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_jUp"    , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_jUp"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAJESBoundingDown [nModel]              = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_jDown"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_jDown"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTMuonScaleBoundingUp [nModel]                = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_muonUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_muonUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMuonScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTMuonScaleBoundingDown [nModel]              = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_muonDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_muonDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMuonScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTElectronScaleBoundingUp [nModel]                = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_electronUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_electronUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTElectronScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTElectronScaleBoundingDown [nModel]              = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_electronDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_electronDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTElectronScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTMETScaleBoundingUp [nModel]                = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_METUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_METUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMETScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTMETScaleBoundingDown [nModel]              = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_METDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_METDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMETScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTJetScaleBoundingUp [nModel]                = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_JESUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_JESUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTJetScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTJetScaleBoundingDown [nModel]              = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_JESDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_JESDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTJetScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_PUBoundingUp [nModel]                    = new TH1D( Form("histo_ZH_hinv_%s_CMS_puUp"         , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_puUp"         , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_PUBoundingDown [nModel]                  = new TH1D( Form("histo_ZH_hinv_%s_CMS_puDown"       , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_puDown"       , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingDown[nModel]->Sumw2();

  }

  double bgdDecay[nSigModels][nSelTypes*4][histBins],weiDecay[nSigModels][nSelTypes*4][histBins];
  for(int nModel=0; nModel<nSigModels; nModel++) { for(unsigned int i=0; i<nSelTypes*4; i++) { for(int j=0; j<histBins; j++) {       
    bgdDecay[nModel][i][j] = 0.0; weiDecay[nModel][i][j] = 0.0; 
  }}}
  TFile *mva_trees;
  TTree *data_mva_tree, *Zjets_mva_tree, *EM_mva_tree, *WZ_mva_tree, *ZZ_mva_tree, *VVV_mva_tree, *signal_mva_trees[nSigModels];
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
           mva_weight,
           aux_PDFscale[102],
           aux_PUscale,
           aux_QCDscale_r1f2,
           aux_QCDscale_r1f5,
           aux_QCDscale_r2f1,
           aux_QCDscale_r2f2,
           aux_QCDscale_r5f1,
           aux_QCDscale_r5f5,
           aux_MET_JESup,
           aux_MET_JESdown;
  UChar_t  mva_njets,
           mva_ntaus,
           aux_njets_JESup,
           aux_njets_JESdown;
  Bool_t   mva_btag_veto,
           mva_3lveto;
  if(makeMVAtrees) {
    mva_trees=new TFile("MitZHAnalysis/mva/mva_input_trees.root", "RECREATE");
    data_mva_tree = new TTree("bkg_mva_tree_Zjets", "MVA input tree with data events");
    data_mva_tree->Branch( "mva_balance"             , &mva_balance              , "mva_balance/F"           ); 
    data_mva_tree->Branch( "mva_cos_theta_star_l1"   , &mva_cos_theta_star_l1    , "mva_cos_theta_star_l1/F" ); 
    data_mva_tree->Branch( "mva_cos_theta_CS_l1"     , &mva_cos_theta_CS_l1      , "mva_cos_theta_CS_l1/F"   ); 
    data_mva_tree->Branch( "mva_delphi_ptll_MET"     , &mva_delphi_ptll_MET      , "mva_delphi_ptll_MET/F"   ); 
    data_mva_tree->Branch( "mva_delphi_ll"           , &mva_delphi_ll            , "mva_delphi_ll/F"         ); 
    data_mva_tree->Branch( "mva_delphi_jet_MET"      , &mva_delphi_jet_MET       , "mva_delphi_jet_MET/F"    ); 
    data_mva_tree->Branch( "mva_deltaR_ll"           , &mva_deltaR_ll            , "mva_deltaR_ll/F"         ); 
    data_mva_tree->Branch( "mva_etall"               , &mva_etall                , "mva_etall/F"             ); 
    data_mva_tree->Branch( "mva_etal1"               , &mva_etal1                , "mva_etal1/F"             ); 
    data_mva_tree->Branch( "mva_etal2"               , &mva_etal2                , "mva_etal2/F"             ); 
    data_mva_tree->Branch( "mva_MET"                 , &mva_MET                  , "mva_MET/F"               ); 
    data_mva_tree->Branch( "mva_mll_minus_mZ"        , &mva_mll_minus_mZ         , "mva_mll_minus_mZ/F"      ); 
    data_mva_tree->Branch( "mva_mTjetMET"            , &mva_mTjetMET             , "mva_mTjetMET/F"          ); 
    data_mva_tree->Branch( "mva_mTll"                , &mva_mTll                 , "mva_mTll/F"              ); 
    data_mva_tree->Branch( "mva_mTl1MET"             , &mva_mTl1MET              , "mva_mTl1MET/F"           ); 
    data_mva_tree->Branch( "mva_mTl2MET"             , &mva_mTl2MET              , "mva_mTl2MET/F"           ); 
    data_mva_tree->Branch( "mva_ptll"                , &mva_ptll                 , "mva_ptll/F"              ); 
    data_mva_tree->Branch( "mva_ptl1"                , &mva_ptl1                 , "mva_ptl1/F"              ); 
    data_mva_tree->Branch( "mva_ptl2"                , &mva_ptl2                 , "mva_ptl2/F"              ); 
    data_mva_tree->Branch( "mva_ptl1mptl2_over_ptll" , &mva_ptl1mptl2_over_ptll  , "mva_ptl1mptl2_over_ptll/F"); 
    data_mva_tree->Branch( "mva_weight"              , &mva_weight               , "mva_weight/F"            ); 
    data_mva_tree->Branch( "mva_njets"               , &mva_njets                , "mva_njets/b"             ); 
    data_mva_tree->Branch( "mva_ntaus"               , &mva_ntaus                , "mva_ntaus/b"             ); 
    data_mva_tree->Branch( "mva_3lveto"              , &mva_3lveto               , "mva_3lveto/O"            ); 
    data_mva_tree->Branch( "mva_btag_veto"           , &mva_btag_veto            , "mva_btag_veto/O"         ); 
    Zjets_mva_tree = (TTree*) data_mva_tree->CloneTree(); Zjets_mva_tree  ->SetName("bkg_mva_tree_Zjets" ); Zjets_mva_tree  ->SetTitle( "MVA input tree with Drell-Yan background events" );
    Zjets_mva_tree->Branch( "aux_PDFscale"            , &aux_PDFscale);
    Zjets_mva_tree->Branch( "aux_PUscale"             , &aux_PUscale              , "aux_PUscale/F"           );
    Zjets_mva_tree->Branch( "aux_QCDscale_r1f2"       , &aux_QCDscale_r1f2        , "aux_QCDscale_r1f2/F"     );
    Zjets_mva_tree->Branch( "aux_QCDscale_r1f5"       , &aux_QCDscale_r1f5        , "aux_QCDscale_r1f5/F"     );
    Zjets_mva_tree->Branch( "aux_QCDscale_r2f1"       , &aux_QCDscale_r2f1        , "aux_QCDscale_r2f1/F"     );
    Zjets_mva_tree->Branch( "aux_QCDscale_r2f2"       , &aux_QCDscale_r2f2        , "aux_QCDscale_r2f2/F"     );
    Zjets_mva_tree->Branch( "aux_QCDscale_r5f1"       , &aux_QCDscale_r5f1        , "aux_QCDscale_r5f1/F"     );
    Zjets_mva_tree->Branch( "aux_QCDscale_r5f5"       , &aux_QCDscale_r5f5        , "aux_QCDscale_r5f5/F"     );
    Zjets_mva_tree->Branch( "aux_njets_JESup"         , &aux_njets_JESup          , "aux_njets_JESup/b"       ); 
    Zjets_mva_tree->Branch( "aux_njets_JESdown"       , &aux_njets_JESdown        , "aux_njets_JESdown/b"     ); 
    Zjets_mva_tree->Branch( "aux_MET_JESup"           , &aux_MET_JESup            , "aux_MET_JESup/F"         );
    Zjets_mva_tree->Branch( "aux_MET_JESdown"         , &aux_MET_JESdown          , "aux_MET_JESdown/F"       );
    EM_mva_tree    = (TTree*)Zjets_mva_tree->CloneTree(); EM_mva_tree     ->SetName("bkg_mva_tree_EM" );    EM_mva_tree     ->SetTitle( "MVA input tree with WW/top background events" );
    WZ_mva_tree    = (TTree*)Zjets_mva_tree->CloneTree(); WZ_mva_tree     ->SetName("bkg_mva_tree_WZ" );    WZ_mva_tree     ->SetTitle( "MVA input tree with WZ background events"     );
    ZZ_mva_tree    = (TTree*)Zjets_mva_tree->CloneTree(); ZZ_mva_tree     ->SetName("bkg_mva_tree_ZZ" );    ZZ_mva_tree     ->SetTitle( "MVA input tree with ZZ background events"     );
    VVV_mva_tree   = (TTree*)Zjets_mva_tree->CloneTree(); VVV_mva_tree    ->SetName("bkg_mva_tree_VVV");    VVV_mva_tree    ->SetTitle( "MVA input tree with VVV background events"    );
    for(int nModel=0; nModel<nSigModels; nModel++) {
      signal_mva_trees[nModel] = (TTree*)Zjets_mva_tree->CloneTree(); 
      signal_mva_trees[nModel]->SetName( Form("signal_mva_tree_%s", signalName_[nModel].Data()));
      signal_mva_trees[nModel]->SetTitle(Form("MVA input tree with signal events (%s)", signalName_[nModel].Data()));
    }
  }
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
      reader->AddVariable( "ptl1mptl2_over_ptll"               , &mva_ptl1mptl2_over_ptll);
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
      reader->AddVariable( "ptl1mptl2_over_ptll"               , &mva_ptl1mptl2_over_ptll);
    }
    reader->BookMVA("BDT", the_BDT_weights);
  }
  
  const unsigned int num_bdt_toys=1000;
  float bdt_toy_scale[num_bdt_toys];
  TRandom toy_machine(3393); // Random seed is Dylan's birthday :-)
  // Pointers for TH2 objects to store the toy BDT shapes
  TH2F* histo_bdt_toys_electronScale_VVV, *histo_bdt_toys_electronScale_WZ, *histo_bdt_toys_electronScale_ZZ, *histo_bdt_toys_electronScale_ggZH_hinv, *histo_bdt_toys_electronScale_ZH_hinv[nSigModels];
  TH2F* histo_bdt_toys_muonScale_VVV, *histo_bdt_toys_muonScale_WZ, *histo_bdt_toys_muonScale_ZZ, *histo_bdt_toys_muonScale_ggZH_hinv, *histo_bdt_toys_muonScale_ZH_hinv[nSigModels];
  TH2F* histo_bdt_toys_METScale_VVV, *histo_bdt_toys_METScale_WZ, *histo_bdt_toys_METScale_ZZ, *histo_bdt_toys_METScale_ggZH_hinv, *histo_bdt_toys_METScale_ZH_hinv[nSigModels];
  // Pointers for TH1 arrays to store the bin yields
  TH1F *bdt_toy_binyields_electronScale_VVV[nBinMVA], *bdt_toy_binyields_electronScale_WZ[nBinMVA], *bdt_toy_binyields_electronScale_ZZ[nBinMVA], *bdt_toy_binyields_electronScale_ZH_hinv[nSigModels][nBinMVA], *bdt_toy_binyields_electronScale_ggZH_hinv[nBinMVA];
  TH1F *bdt_toy_binyields_muonScale_VVV[nBinMVA], *bdt_toy_binyields_muonScale_WZ[nBinMVA], *bdt_toy_binyields_muonScale_ZZ[nBinMVA], *bdt_toy_binyields_muonScale_ZH_hinv[nSigModels][nBinMVA], *bdt_toy_binyields_muonScale_ggZH_hinv[nBinMVA];
  TH1F *bdt_toy_binyields_METScale_VVV[nBinMVA], *bdt_toy_binyields_METScale_WZ[nBinMVA], *bdt_toy_binyields_METScale_ZZ[nBinMVA], *bdt_toy_binyields_METScale_ZH_hinv[nSigModels][nBinMVA], *bdt_toy_binyields_METScale_ggZH_hinv[nBinMVA];
  // Pointers for TH1's to store the relative up/down systematics (this is what will be cached)
  TH1F *bdt_syst_electronScaleUp_VVV, *bdt_syst_electronScaleUp_WZ, *bdt_syst_electronScaleUp_ZZ, *bdt_syst_electronScaleUp_ZH_hinv[nSigModels], *bdt_syst_electronScaleUp_ggZH_hinv;
  TH1F *bdt_syst_muonScaleUp_VVV, *bdt_syst_muonScaleUp_WZ, *bdt_syst_muonScaleUp_ZZ, *bdt_syst_muonScaleUp_ZH_hinv[nSigModels], *bdt_syst_muonScaleUp_ggZH_hinv;
  TH1F *bdt_syst_METScaleUp_VVV, *bdt_syst_METScaleUp_WZ, *bdt_syst_METScaleUp_ZZ, *bdt_syst_METScaleUp_ZH_hinv[nSigModels], *bdt_syst_METScaleUp_ggZH_hinv;
  TH1F *bdt_syst_electronScaleDown_VVV, *bdt_syst_electronScaleDown_WZ, *bdt_syst_electronScaleDown_ZZ, *bdt_syst_electronScaleDown_ZH_hinv[nSigModels], *bdt_syst_electronScaleDown_ggZH_hinv;
  TH1F *bdt_syst_muonScaleDown_VVV, *bdt_syst_muonScaleDown_WZ, *bdt_syst_muonScaleDown_ZZ, *bdt_syst_muonScaleDown_ZH_hinv[nSigModels], *bdt_syst_muonScaleDown_ggZH_hinv;
  TH1F *bdt_syst_METScaleDown_VVV, *bdt_syst_METScaleDown_WZ, *bdt_syst_METScaleDown_ZZ, *bdt_syst_METScaleDown_ZH_hinv[nSigModels], *bdt_syst_METScaleDown_ggZH_hinv;
    // Pointers for TH2's which hold the 2D BDT toy envelope maps
  TH2F *bdt_toy_envelope_electronScale_VVV, *bdt_toy_envelope_electronScale_WZ, *bdt_toy_envelope_electronScale_ZZ, *bdt_toy_envelope_electronScale_ZH_hinv[nSigModels], *bdt_toy_envelope_electronScale_ggZH_hinv;
  TH2F *bdt_toy_envelope_muonScale_VVV, *bdt_toy_envelope_muonScale_WZ, *bdt_toy_envelope_muonScale_ZZ, *bdt_toy_envelope_muonScale_ZH_hinv[nSigModels], *bdt_toy_envelope_muonScale_ggZH_hinv;
  TH2F *bdt_toy_envelope_METScale_VVV, *bdt_toy_envelope_METScale_WZ, *bdt_toy_envelope_METScale_ZZ, *bdt_toy_envelope_METScale_ZH_hinv[nSigModels], *bdt_toy_envelope_METScale_ggZH_hinv;
  TFile* cached_BDT_systematics;
  if (useBDT) { 
    char inputCachedBDTSysts[200]; // Filename for the cached bdt systematics, we will either read them from here or write fresh ones to here
    sprintf(inputCachedBDTSysts,"MitZHAnalysis/plots%s/zll%shinv%s_%s_BDTsyst_%s.root", subdirectory.c_str(), addChan.Data(), finalStateName, signalName_[plotModel].Data(), ECMsb.Data());
    if(useCachedBDTSystematics) {
      printf("Opening file with cached BDT systematics: \"%s\"\n", inputCachedBDTSysts);
      cached_BDT_systematics = new TFile(inputCachedBDTSysts,"read"); assert(cached_BDT_systematics->IsOpen());
      bdt_syst_electronScaleUp_VVV       = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleUp_VVV"      );
      bdt_syst_electronScaleUp_WZ        = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleUp_WZ"       );
      bdt_syst_electronScaleUp_ZZ        = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleUp_ZZ"       );
      bdt_syst_electronScaleUp_ggZH_hinv = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleUp_ggZH_hinv");
      bdt_syst_electronScaleUp_ZH_hinv[plotModel] = (TH1F*) cached_BDT_systematics->Get(Form("bdt_syst_electronScaleUp_ZH_hinv_%s",signalName_[plotModel].Data()));
      bdt_syst_muonScaleUp_VVV           = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleUp_VVV"          );
      bdt_syst_muonScaleUp_WZ            = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleUp_WZ"           );
      bdt_syst_muonScaleUp_ZZ            = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleUp_ZZ"           );
      bdt_syst_muonScaleUp_ggZH_hinv     = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleUp_ggZH_hinv"    );
      bdt_syst_muonScaleUp_ZH_hinv[plotModel] = (TH1F*) cached_BDT_systematics->Get(Form("bdt_syst_muonScaleUp_ZH_hinv_%s",signalName_[plotModel].Data()));
      bdt_syst_METScaleUp_VVV            = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleUp_VVV"           );
      bdt_syst_METScaleUp_WZ             = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleUp_WZ"            );
      bdt_syst_METScaleUp_ZZ             = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleUp_ZZ"            );
      bdt_syst_METScaleUp_ggZH_hinv      = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleUp_ggZH_hinv"     );
      bdt_syst_METScaleUp_ZH_hinv[plotModel] = (TH1F*) cached_BDT_systematics->Get(Form("bdt_syst_METScaleUp_ZH_hinv_%s",signalName_[plotModel].Data()));
      bdt_syst_electronScaleDown_VVV       = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleDown_VVV"      );
      bdt_syst_electronScaleDown_WZ        = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleDown_WZ"       );
      bdt_syst_electronScaleDown_ZZ        = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleDown_ZZ"       );
      bdt_syst_electronScaleDown_ggZH_hinv = (TH1F*) cached_BDT_systematics->Get("bdt_syst_electronScaleDown_ggZH_hinv");
      bdt_syst_electronScaleDown_ZH_hinv[plotModel] = (TH1F*) cached_BDT_systematics->Get(Form("bdt_syst_electronScaleDown_ZH_hinv_%s",signalName_[plotModel].Data()));
      bdt_syst_muonScaleDown_VVV           = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleDown_VVV"          );
      bdt_syst_muonScaleDown_WZ            = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleDown_WZ"           );
      bdt_syst_muonScaleDown_ZZ            = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleDown_ZZ"           );
      bdt_syst_muonScaleDown_ggZH_hinv     = (TH1F*) cached_BDT_systematics->Get("bdt_syst_muonScaleDown_ggZH_hinv"    );
      bdt_syst_muonScaleDown_ZH_hinv[plotModel] = (TH1F*) cached_BDT_systematics->Get(Form("bdt_syst_muonScaleDown_ZH_hinv_%s",signalName_[plotModel].Data()));
      bdt_syst_METScaleDown_VVV            = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleDown_VVV"           );
      bdt_syst_METScaleDown_WZ             = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleDown_WZ"            );
      bdt_syst_METScaleDown_ZZ             = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleDown_ZZ"            );
      bdt_syst_METScaleDown_ggZH_hinv      = (TH1F*) cached_BDT_systematics->Get("bdt_syst_METScaleDown_ggZH_hinv"     );
      bdt_syst_METScaleDown_ZH_hinv[plotModel] = (TH1F*) cached_BDT_systematics->Get(Form("bdt_syst_METScaleDown_ZH_hinv_%s",signalName_[plotModel].Data()));
      bdt_syst_electronScaleUp_VVV       ->SetDirectory(0);
      bdt_syst_electronScaleUp_WZ        ->SetDirectory(0);
      bdt_syst_electronScaleUp_ZZ        ->SetDirectory(0);
      bdt_syst_electronScaleUp_ggZH_hinv ->SetDirectory(0);
      bdt_syst_electronScaleUp_ZH_hinv[plotModel] ->SetDirectory(0);
      bdt_syst_muonScaleUp_VVV           ->SetDirectory(0);
      bdt_syst_muonScaleUp_WZ            ->SetDirectory(0);
      bdt_syst_muonScaleUp_ZZ            ->SetDirectory(0);
      bdt_syst_muonScaleUp_ggZH_hinv     ->SetDirectory(0);
      bdt_syst_muonScaleUp_ZH_hinv[plotModel] ->SetDirectory(0);
      bdt_syst_METScaleUp_VVV            ->SetDirectory(0);
      bdt_syst_METScaleUp_WZ             ->SetDirectory(0);
      bdt_syst_METScaleUp_ZZ             ->SetDirectory(0);
      bdt_syst_METScaleUp_ggZH_hinv      ->SetDirectory(0);
      bdt_syst_METScaleUp_ZH_hinv[plotModel] ->SetDirectory(0);
      bdt_syst_electronScaleDown_VVV       ->SetDirectory(0);
      bdt_syst_electronScaleDown_WZ        ->SetDirectory(0);
      bdt_syst_electronScaleDown_ZZ        ->SetDirectory(0);
      bdt_syst_electronScaleDown_ggZH_hinv ->SetDirectory(0);
      bdt_syst_electronScaleDown_ZH_hinv[plotModel] ->SetDirectory(0);
      bdt_syst_muonScaleDown_VVV           ->SetDirectory(0);
      bdt_syst_muonScaleDown_WZ            ->SetDirectory(0);
      bdt_syst_muonScaleDown_ZZ            ->SetDirectory(0);
      bdt_syst_muonScaleDown_ggZH_hinv     ->SetDirectory(0);
      bdt_syst_muonScaleDown_ZH_hinv[plotModel] ->SetDirectory(0);
      bdt_syst_METScaleDown_VVV            ->SetDirectory(0);
      bdt_syst_METScaleDown_WZ             ->SetDirectory(0);
      bdt_syst_METScaleDown_ZZ             ->SetDirectory(0);
      bdt_syst_METScaleDown_ggZH_hinv      ->SetDirectory(0);
      bdt_syst_METScaleDown_ZH_hinv[plotModel] ->SetDirectory(0);
      cached_BDT_systematics->Close();
    } else { // If not using cached BDT systematics, construct all of the histos we are going to use to compute them
      cached_BDT_systematics = new TFile(inputCachedBDTSysts,"recreate");
      assert(cached_BDT_systematics->IsOpen());
      bdt_toy_envelope_electronScale_VVV        = new TH2F("bdt_toy_envelope_electronScale_VVV"       , "bdt_toy_envelope_electronScale_VVV"       , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_electronScale_WZ         = new TH2F("bdt_toy_envelope_electronScale_WZ"        , "bdt_toy_envelope_electronScale_WZ"        , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_electronScale_ZZ         = new TH2F("bdt_toy_envelope_electronScale_ZZ"        , "bdt_toy_envelope_electronScale_ZZ"        , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_electronScale_ggZH_hinv  = new TH2F("bdt_toy_envelope_electronScale_ggZH_hinv" , "bdt_toy_envelope_electronScale_ggZH_hinv" , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_muonScale_VVV        = new TH2F("bdt_toy_envelope_muonScale_VVV"       , "bdt_toy_envelope_muonScale_VVV"       , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_muonScale_WZ         = new TH2F("bdt_toy_envelope_muonScale_WZ"        , "bdt_toy_envelope_muonScale_WZ"        , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_muonScale_ZZ         = new TH2F("bdt_toy_envelope_muonScale_ZZ"        , "bdt_toy_envelope_muonScale_ZZ"        , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_muonScale_ggZH_hinv  = new TH2F("bdt_toy_envelope_muonScale_ggZH_hinv" , "bdt_toy_envelope_muonScale_ggZH_hinv" , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_METScale_VVV        = new TH2F("bdt_toy_envelope_METScale_VVV"       , "bdt_toy_envelope_METScale_VVV"       , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_METScale_WZ         = new TH2F("bdt_toy_envelope_METScale_WZ"        , "bdt_toy_envelope_METScale_WZ"        , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_METScale_ZZ         = new TH2F("bdt_toy_envelope_METScale_ZZ"        , "bdt_toy_envelope_METScale_ZZ"        , nBinMVA, xbins, 1000 , 0, 4); 
      bdt_toy_envelope_METScale_ggZH_hinv  = new TH2F("bdt_toy_envelope_METScale_ggZH_hinv" , "bdt_toy_envelope_METScale_ggZH_hinv" , nBinMVA, xbins, 1000 , 0, 4); 
      histo_bdt_toys_electronScale_VVV       = new TH2F("histo_bdt_toys_electronScale_VVV"      , "histo_bdt_toys_electronScale_VVV"       , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_VVV->Sumw2();
      histo_bdt_toys_electronScale_WZ        = new TH2F("histo_bdt_toys_electronScale_WZ"       , "histo_bdt_toys_electronScale_WZ"        , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_WZ->Sumw2();
      histo_bdt_toys_electronScale_ZZ        = new TH2F("histo_bdt_toys_electronScale_ZZ"       , "histo_bdt_toys_electronScale_ZZ"        , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_ZZ->Sumw2();
      histo_bdt_toys_electronScale_ggZH_hinv = new TH2F("histo_bdt_toys_electronScale_ggZH_hinv", "histo_bdt_toys_electronScale_ggZH_hinv" , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_ggZH_hinv->Sumw2();
      histo_bdt_toys_muonScale_VVV           = new TH2F("histo_bdt_toys_muonScale_VVV"          , "histo_bdt_toys_muonScale_VVV"           , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_VVV->Sumw2();
      histo_bdt_toys_muonScale_WZ            = new TH2F("histo_bdt_toys_muonScale_WZ"           , "histo_bdt_toys_muonScale_WZ"            , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_WZ->Sumw2();
      histo_bdt_toys_muonScale_ZZ            = new TH2F("histo_bdt_toys_muonScale_ZZ"           , "histo_bdt_toys_muonScale_ZZ"            , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_ZZ->Sumw2();
      histo_bdt_toys_muonScale_ggZH_hinv     = new TH2F("histo_bdt_toys_muonScale_ggZH_hinv"    , "histo_bdt_toys_muonScale_ggZH_hinv"     , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_ggZH_hinv->Sumw2();
      histo_bdt_toys_METScale_VVV            = new TH2F("histo_bdt_toys_METScale_VVV"           , "histo_bdt_toys_METScale_VVV"            , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_VVV->Sumw2();
      histo_bdt_toys_METScale_WZ             = new TH2F("histo_bdt_toys_METScale_WZ"            , "histo_bdt_toys_METScale_WZ"             , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_WZ->Sumw2();
      histo_bdt_toys_METScale_ZZ             = new TH2F("histo_bdt_toys_METScale_ZZ"            , "histo_bdt_toys_METScale_ZZ"             , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_ZZ->Sumw2();
      histo_bdt_toys_METScale_ggZH_hinv      = new TH2F("histo_bdt_toys_METScale_ggZH_hinv"     , "histo_bdt_toys_METScale_ggZH_hinv"      , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_ggZH_hinv->Sumw2();
      bdt_syst_electronScaleUp_VVV       = new TH1F("bdt_syst_electronScaleUp_VVV", "bdt_syst_electronScaleUp_VVV", nBinMVA, xbins);
      bdt_syst_electronScaleUp_WZ        = new TH1F("bdt_syst_electronScaleUp_WZ", "bdt_syst_electronScaleUp_WZ", nBinMVA, xbins);
      bdt_syst_electronScaleUp_ZZ        = new TH1F("bdt_syst_electronScaleUp_ZZ", "bdt_syst_electronScaleUp_ZZ", nBinMVA, xbins);
      bdt_syst_electronScaleUp_ggZH_hinv = new TH1F("bdt_syst_electronScaleUp_ggZH_hinv", "bdt_syst_electronScaleUp_ggZH_hinv", nBinMVA, xbins);
      bdt_syst_muonScaleUp_VVV           = new TH1F("bdt_syst_muonScaleUp_VVV", "bdt_syst_muonScaleUp_VVV", nBinMVA, xbins);
      bdt_syst_muonScaleUp_WZ            = new TH1F("bdt_syst_muonScaleUp_WZ", "bdt_syst_muonScaleUp_WZ", nBinMVA, xbins);
      bdt_syst_muonScaleUp_ZZ            = new TH1F("bdt_syst_muonScaleUp_ZZ", "bdt_syst_muonScaleUp_ZZ", nBinMVA, xbins);
      bdt_syst_muonScaleUp_ggZH_hinv     = new TH1F("bdt_syst_muonScaleUp_ggZH_hinv", "bdt_syst_muonScaleUp_ggZH_hinv", nBinMVA, xbins);
      bdt_syst_METScaleUp_VVV            = new TH1F("bdt_syst_METScaleUp_VVV", "bdt_syst_METScaleUp_VVV", nBinMVA, xbins);
      bdt_syst_METScaleUp_WZ             = new TH1F("bdt_syst_METScaleUp_WZ", "bdt_syst_METScaleUp_WZ", nBinMVA, xbins);
      bdt_syst_METScaleUp_ZZ             = new TH1F("bdt_syst_METScaleUp_ZZ", "bdt_syst_METScaleUp_ZZ", nBinMVA, xbins);
      bdt_syst_METScaleUp_ggZH_hinv      = new TH1F("bdt_syst_METScaleUp_ggZH_hinv", "bdt_syst_METScaleUp_ggZH_hinv", nBinMVA, xbins);
      bdt_syst_electronScaleDown_VVV       = new TH1F("bdt_syst_electronScaleDown_VVV", "bdt_syst_electronScaleDown_VVV", nBinMVA, xbins);
      bdt_syst_electronScaleDown_WZ        = new TH1F("bdt_syst_electronScaleDown_WZ", "bdt_syst_electronScaleDown_WZ", nBinMVA, xbins);
      bdt_syst_electronScaleDown_ZZ        = new TH1F("bdt_syst_electronScaleDown_ZZ", "bdt_syst_electronScaleDown_ZZ", nBinMVA, xbins);
      bdt_syst_electronScaleDown_ggZH_hinv = new TH1F("bdt_syst_electronScaleDown_ggZH_hinv", "bdt_syst_electronScaleDown_ggZH_hinv", nBinMVA, xbins);
      bdt_syst_muonScaleDown_VVV           = new TH1F("bdt_syst_muonScaleDown_VVV", "bdt_syst_muonScaleDown_VVV", nBinMVA, xbins);
      bdt_syst_muonScaleDown_WZ            = new TH1F("bdt_syst_muonScaleDown_WZ", "bdt_syst_muonScaleDown_WZ", nBinMVA, xbins);
      bdt_syst_muonScaleDown_ZZ            = new TH1F("bdt_syst_muonScaleDown_ZZ", "bdt_syst_muonScaleDown_ZZ", nBinMVA, xbins);
      bdt_syst_muonScaleDown_ggZH_hinv     = new TH1F("bdt_syst_muonScaleDown_ggZH_hinv", "bdt_syst_muonScaleDown_ggZH_hinv", nBinMVA, xbins);
      bdt_syst_METScaleDown_VVV            = new TH1F("bdt_syst_METScaleDown_VVV", "bdt_syst_METScaleDown_VVV", nBinMVA, xbins);
      bdt_syst_METScaleDown_WZ             = new TH1F("bdt_syst_METScaleDown_WZ", "bdt_syst_METScaleDown_WZ", nBinMVA, xbins);
      bdt_syst_METScaleDown_ZZ             = new TH1F("bdt_syst_METScaleDown_ZZ", "bdt_syst_METScaleDown_ZZ", nBinMVA, xbins);
      bdt_syst_METScaleDown_ggZH_hinv      = new TH1F("bdt_syst_METScaleDown_ggZH_hinv", "bdt_syst_METScaleDown_ggZH_hinv", nBinMVA, xbins);
      bdt_syst_electronScaleUp_VVV       ->SetDirectory(0);
      bdt_syst_electronScaleUp_WZ        ->SetDirectory(0);
      bdt_syst_electronScaleUp_ZZ        ->SetDirectory(0);
      bdt_syst_electronScaleUp_ggZH_hinv ->SetDirectory(0);
      bdt_syst_muonScaleUp_VVV           ->SetDirectory(0);
      bdt_syst_muonScaleUp_WZ            ->SetDirectory(0);
      bdt_syst_muonScaleUp_ZZ            ->SetDirectory(0);
      bdt_syst_muonScaleUp_ggZH_hinv     ->SetDirectory(0);
      bdt_syst_METScaleUp_VVV            ->SetDirectory(0);
      bdt_syst_METScaleUp_WZ             ->SetDirectory(0);
      bdt_syst_METScaleUp_ZZ             ->SetDirectory(0);
      bdt_syst_METScaleUp_ggZH_hinv      ->SetDirectory(0);
      bdt_syst_electronScaleDown_VVV       ->SetDirectory(0);
      bdt_syst_electronScaleDown_WZ        ->SetDirectory(0);
      bdt_syst_electronScaleDown_ZZ        ->SetDirectory(0);
      bdt_syst_electronScaleDown_ggZH_hinv ->SetDirectory(0);
      bdt_syst_muonScaleDown_VVV           ->SetDirectory(0);
      bdt_syst_muonScaleDown_WZ            ->SetDirectory(0);
      bdt_syst_muonScaleDown_ZZ            ->SetDirectory(0);
      bdt_syst_muonScaleDown_ggZH_hinv     ->SetDirectory(0);
      bdt_syst_METScaleDown_VVV            ->SetDirectory(0);
      bdt_syst_METScaleDown_WZ             ->SetDirectory(0);
      bdt_syst_METScaleDown_ZZ             ->SetDirectory(0);
      bdt_syst_METScaleDown_ggZH_hinv      ->SetDirectory(0);
      for(int nModel=0; nModel<nSigModels; nModel++) { 
        bdt_toy_envelope_electronScale_ZH_hinv[nModel] = new TH2F(Form("bdt_toy_envelope_electronScale_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_toy_envelope_electronScale_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins, 1000 , 0, 4); 
        bdt_toy_envelope_muonScale_ZH_hinv[nModel]     = new TH2F(Form("bdt_toy_envelope_muonScale_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_toy_envelope_muonScale_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins, 1000 , 0, 4); 
        bdt_toy_envelope_METScale_ZH_hinv[nModel]      = new TH2F(Form("bdt_toy_envelope_METScale_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_toy_envelope_METScale_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins, 1000 , 0, 4); 
        histo_bdt_toys_electronScale_ZH_hinv[nModel]   = new TH2F(Form("histo_bdt_toys_electronScale_ZH_hinv_%s",signalName_[nModel].Data()),   Form("histo_bdt_toys_electronScale_ZH_hinv_%s",signalName_[nModel].Data()),   nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_ZH_hinv[nModel]->Sumw2();
        histo_bdt_toys_muonScale_ZH_hinv[nModel]   = new TH2F(Form("histo_bdt_toys_muonScale_ZH_hinv_%s",signalName_[nModel].Data()),   Form("histo_bdt_toys_muonScale_ZH_hinv_%s",signalName_[nModel].Data()),   nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_ZH_hinv[nModel]->Sumw2();
        histo_bdt_toys_METScale_ZH_hinv[nModel]   = new TH2F(Form("histo_bdt_toys_METScale_ZH_hinv_%s",signalName_[nModel].Data()),   Form("histo_bdt_toys_METScale_ZH_hinv_%s",signalName_[nModel].Data()),   nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_ZH_hinv[nModel]->Sumw2();
        bdt_syst_electronScaleUp_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_electronScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_electronScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
        bdt_syst_muonScaleUp_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_muonScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_muonScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
        bdt_syst_METScaleUp_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_METScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_METScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
        bdt_syst_electronScaleDown_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_electronScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_electronScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
        bdt_syst_muonScaleDown_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_muonScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_muonScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
        bdt_syst_METScaleDown_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_METScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_METScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
        bdt_syst_electronScaleUp_ZH_hinv[nModel] ->SetDirectory(0);
        bdt_syst_muonScaleUp_ZH_hinv[nModel] ->SetDirectory(0);
        bdt_syst_METScaleUp_ZH_hinv[nModel] ->SetDirectory(0);
        bdt_syst_electronScaleDown_ZH_hinv[nModel] ->SetDirectory(0);
        bdt_syst_muonScaleDown_ZH_hinv[nModel] ->SetDirectory(0);
        bdt_syst_METScaleDown_ZH_hinv[nModel] ->SetDirectory(0);
      }
      for(int nb=1; nb<nBinMVA; nb++) {  
        bdt_syst_electronScaleUp_VVV       ->SetBinContent(nb,1.);
        bdt_syst_electronScaleUp_WZ        ->SetBinContent(nb,1.);
        bdt_syst_electronScaleUp_ZZ        ->SetBinContent(nb,1.);
        bdt_syst_electronScaleUp_ggZH_hinv ->SetBinContent(nb,1.);
        bdt_syst_muonScaleUp_VVV           ->SetBinContent(nb,1.);
        bdt_syst_muonScaleUp_WZ            ->SetBinContent(nb,1.);
        bdt_syst_muonScaleUp_ZZ            ->SetBinContent(nb,1.);
        bdt_syst_muonScaleUp_ggZH_hinv     ->SetBinContent(nb,1.);
        bdt_syst_METScaleUp_VVV            ->SetBinContent(nb,1.);
        bdt_syst_METScaleUp_WZ             ->SetBinContent(nb,1.);
        bdt_syst_METScaleUp_ZZ             ->SetBinContent(nb,1.);
        bdt_syst_METScaleUp_ggZH_hinv      ->SetBinContent(nb,1.);
        bdt_syst_electronScaleDown_VVV       ->SetBinContent(nb,1.);
        bdt_syst_electronScaleDown_WZ        ->SetBinContent(nb,1.);
        bdt_syst_electronScaleDown_ZZ        ->SetBinContent(nb,1.);
        bdt_syst_electronScaleDown_ggZH_hinv ->SetBinContent(nb,1.);
        bdt_syst_muonScaleDown_VVV           ->SetBinContent(nb,1.);
        bdt_syst_muonScaleDown_WZ            ->SetBinContent(nb,1.);
        bdt_syst_muonScaleDown_ZZ            ->SetBinContent(nb,1.);
        bdt_syst_muonScaleDown_ggZH_hinv     ->SetBinContent(nb,1.);
        bdt_syst_METScaleDown_VVV            ->SetBinContent(nb,1.);
        bdt_syst_METScaleDown_WZ             ->SetBinContent(nb,1.);
        bdt_syst_METScaleDown_ZZ             ->SetBinContent(nb,1.);
        bdt_syst_METScaleDown_ggZH_hinv      ->SetBinContent(nb,1.);
        for(int nModel=0; nModel<nSigModels; nModel++) {
          bdt_syst_electronScaleUp_ZH_hinv[nModel] ->SetBinContent(nb,1.);
          bdt_syst_muonScaleUp_ZH_hinv[nModel] ->SetBinContent(nb,1.);
          bdt_syst_METScaleUp_ZH_hinv[nModel] ->SetBinContent(nb,1.);
          bdt_syst_electronScaleDown_ZH_hinv[nModel] ->SetBinContent(nb,1.);
          bdt_syst_muonScaleDown_ZH_hinv[nModel] ->SetBinContent(nb,1.);
          bdt_syst_METScaleDown_ZH_hinv[nModel] ->SetBinContent(nb,1.);
        }
      }
      // Store a value from Gaussian distribution with mean=0, sigma=1
      // These values will get multiplied by the individual nuisance sizes
      for(unsigned int i_toy=0; i_toy<num_bdt_toys; i_toy++) {
        bdt_toy_scale[i_toy]=toy_machine.Gaus(0,1);
      }
    }
  }

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
    if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
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
    if(infilecatv[ifile] == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }
    else {
    }
    
    char weightDef[256];
    int initPDFTag = -1;
    if(the_PDF_tree) {
      the_PDF_tree->SetBranchAddress("weightDef", &weightDef);
      for (int i=0; i<the_PDF_tree->GetEntries(); ++i) {
        the_PDF_tree->GetEntry(i);
        char **tokensPDF;
        size_t numtokensPDF;
        tokensPDF = strsplit(weightDef, " = ", &numtokensPDF);
        for (int k = 0; k < (int)numtokensPDF; k++) if(strcmp(tokensPDF[k],"292201") == 0||
	                                               strcmp(tokensPDF[k],"292001") == 0||
						       strcmp(tokensPDF[k],"260001") == 0) {initPDFTag = i; break;}
	if(initPDFTag != -1) break;
      }
    }
    if(infilecatv[ifile] != 0 && initPDFTag == -1) {
      printf("PDFTAG PROBLEM\n");
      if(the_PDF_tree) {
        printf("PDFTree Entries: %d\n",(int)the_PDF_tree->GetEntries());
        for (int i=0; i<the_PDF_tree->GetEntries(); ++i) {
          the_PDF_tree->GetEntry(i);
	  //printf("PDF(%d): %s\n",i,weightDef);
        }
      }
      else {
        printf("PDFTree not available\n");
      }
      //return;
    }

    bool errorMsgQCDscale = false;
    unsigned int selBit_= 0;
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    histoZHSEL[0]->Scale(0.0);
    histoZHSEL[1]->Scale(0.0);
    histoZHSEL[2]->Scale(0.0);
    histoZHSEL[3]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;
      the_input_tree->GetEntry(i);

      Bool_t passFilter[4] = {kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      if(infilecatv[ifile] == 0) {
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

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 
      TLorentzVector dilepMET(dilep + (*((TLorentzVector*)(*eventMet.p4)[0]))); 

      vector<int> idJet,idJetUp,idJetDown,idBJet;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double mTJetMET = -1;
      double dPhiJetDiLep = -1.0;
      TLorentzVector dilepJet = dilep;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 15) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.4) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1 && ((TLorentzVector*)(*eventJets.p4)[nj])->Pt()> 30) {
          dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
          mTJetMET = TMath::Sqrt(2.0*((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])))))); 
        }

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20) { 
	   if ((float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];
           if ((float)(*eventJets.bDiscr)[nj] > 0.8) idBJet.push_back(nj);
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
           (double)(*eventTaus.iso)[ntau] < 5.0){
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
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typePair = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typePair = 2;
      

      // Calculate a lot of physics quantities used for the rectangular selection

      double dPhiDiLepMET = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))); // TMath::Abs((*(TLorentzVector*)(*eventMet.p4)[0]).DeltaPhi(*eventMet.trackMet));
      double ptFrac = TMath::Abs(dilep.Pt()-((TLorentzVector*)(*eventMet.p4)[0])->Pt())/dilep.Pt(); // TMath::Abs(dilepJet.Pt()-((TLorentzVector*)(*eventMet.p4)[0])->Pt())/dilepJet.Pt();
      double deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));

      double caloMinusPFMETRel = TMath::Abs( (double)eventMet.CaloMet->Pt() - ((TLorentzVector*)(*eventMet.p4)[0])->Pt() ) / ((TLorentzVector*)(*eventMet.p4)[0])->Pt();
      
      TVector2 metv(((TLorentzVector*)(*eventMet.p4)[0])->Px(), ((TLorentzVector*)(*eventMet.p4)[0])->Py());
      TVector2 dilv(dilep.Px(), dilep.Py());
      TVector2 utv = -1.*(metv+dilv);
      double phiv = utv.DeltaPhi(dilv);
      double the_upara = TMath::Abs(utv.Mod()*TMath::Cos(phiv))/dilep.Pt();
      
      // Helicity angle calculation
      double cos_theta_star_l1 = cos_theta_star( *(TLorentzVector*)(*eventLeptons.p4)[idLep[0]], *(TLorentzVector*)(*eventLeptons.p4)[idLep[1]], dilepMET);
      
      bool passZMass     = dilep.M() > 76.1876 && dilep.M() < 106.1876;
      bool passNjets     = idJet.size() <= nJetsType;

      double metMIN = 100; double mtMIN = 200; double metTIGHT = 100;
      
      bool passMETTight  = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > metTIGHT;

      if(MVAVarType == 0)                                            {metMIN = 50; mtMIN = 200;}
      else if(MVAVarType == 1 || MVAVarType == 2 || MVAVarType == 4) {metMIN = 50; mtMIN = 0;}
      else if(MVAVarType ==3) { metMIN = 80; mtMIN = 0;}
      bool passMET = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > metMIN;
      bool passMT = mtW > mtMIN || !passMETTight;
      if(infilecatv[ifile] == 0 && isBlinded) passMET = passMET && ((TLorentzVector*)(*eventMet.p4)[0])->Pt() < 100;

      bool passPTFrac    = ptFrac < 0.4;
      bool passDPhiZMET  = dPhiDiLepMET > 2.6;
      //bool passBtagVeto  = bDiscrMax < 0.800 && idSoft.size() == 0;
      bool passBtagVeto  = bDiscrMax < 0.800;
      bool passPTLL      = dilep.Pt() > 60;
      bool pass3rdLVeto  = idLep.size() == numberOfLeptons && TMath::Abs(signQ) == 0;
      double dphill = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]));
      double detall = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()-((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta());
      double drll = sqrt(dphill*dphill+detall*detall);
      bool passDelphiLL  = drll < 1.8;//dphill < TMath::Pi()/2.;

      bool passZMassLarge = TMath::Abs(dilep.M()-91.1876) < 30.0;
      bool passZMassSB    = (dilep.M() > 110.0 && dilep.M() < 200.0);

      bool passDPhiJetMET = dPhiJetMET == -1 || dPhiJetMET >= 0.5;
      bool passTauVeto    = numberGoodTaus == 0;

      bool passNMinusOne[11] = {
                  passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT &&	       passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass &&  	    passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
                  passZMass && passNjets &&	       passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
        passMT && passZMass && passNjets && passMET &&  	     passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac		  && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET		  && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto	      &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
	passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto  	       && passDPhiJetMET && passTauVeto && passMETTight,
	passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL &&		    passTauVeto && passMETTight,
	passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET		&& passMETTight
                               };


      bool passAllCuts[nSelTypes] = {                   
            		   passZMass && passNjets										  &&  pass3rdLVeto						   ,	 // ZSEL
        		   passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,	 // SIGSEL
        passZMassLarge && !passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,	 // WWSEL
        passZMassSB    && !passZMass && passNjets && passMET    				     && !passBtagVeto		  &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,	 // WWLOOSESEL
        		   passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET && !passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,	 // BTAGSEL
        		   passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL && !pass3rdLVeto,							 // WZSEL
	   		   passZMass && passNjets && passMET    		     && passDPhiZMET		      && passPTLL &&  pass3rdLVeto && ((TLorentzVector*)(*eventMet.p4)[0])->Pt() < 100., // PRESEL
	   		   passZMass && passNjets && passMT && passMET &&!passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,	 // CR1SEL
	   		   passZMass && passNjets && passMT && passMET && passPTFrac &&!passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,	 // CR2SEL
	   		   passZMass && passNjets && passMT && passMET &&!passPTFrac &&!passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,	 // CR12SEL
        		   passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight, // TIGHTSEL
	   		   passZMass && passNjets &&           passMET &&!passPTFrac 		    		      && passPTLL &&  pass3rdLVeto                                                  && passMETTight, // DYSANESEL1
	   		   passZMass && passNjets &&           passMET               &&!passDPhiZMET   		      && passPTLL &&  pass3rdLVeto                                                  && passMETTight  // DYSANESEL2
                                    };
     // Evaluate nominal BDT value
     double bdt_value=-1;
     if(useBDT) {
       TLorentzVector lepton1 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[0]]),
                      lepton2 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[1]]),
                      MET     = *((TLorentzVector*)(*eventMet.p4)[0]),
                      jet1    = idJet.size() > 0 ? *((TLorentzVector*)(*eventJets.p4)[idJet[0]]) : TLorentzVector(0,0,0,0);
        bdt_value = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0,0,0,0);
     }
     if(MVAVarType==3) {
       passAllCuts[WZSEL]  = passZMassLarge && passNjets && passMET && passBtagVeto  && !pass3rdLVeto;
       passAllCuts[PRESEL] = passZMassLarge && passNjets && passMET && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET;
       passAllCuts[SIGSEL] = (passAllCuts[PRESEL] && bdt_value>0 && passMETTight) || (passZMass && passNjets && passMT && passMET && !passMETTight && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto);
       passAllCuts[TIGHTSEL] = passAllCuts[PRESEL] && bdt_value>0 && passMETTight;
     }
     bool passEvolFilter[numberCuts] = {pass3rdLVeto,passBtagVeto,passTauVeto,passNjets,passZMass,passPTLL,passMETTight,passDPhiZMET,passPTFrac,passDPhiJetMET,passDelphiLL&&passMT};
     //bool passEvolFilter[numberCuts] = {pass3rdLVeto,passBtagVeto,passTauVeto,passNjets,passZMass,passPTLL,true,true,true,true,true};

     //if(typePair!=0)printf("LLL %d %d %llu\n",eventEvent.runNum,eventEvent.lumiNum,eventEvent.eventNum);
     int sumEvol = 0;
     bool totalSel = kTRUE;
     for(int isel=0; isel<numberCuts; isel++) {
       totalSel = totalSel && passEvolFilter[isel];
       if(totalSel == kTRUE) sumEvol++;
       //if(totalSel == kTRUE && isel == 2&&typePair!=0) printf("TTT %d %d %llu\n",eventEvent.runNum,eventEvent.lumiNum,eventEvent.eventNum);
       //if(totalSel == kTRUE && isel == 7&&typePair!=0) printf("JJJ %d %d %llu\n",eventEvent.runNum,eventEvent.lumiNum,eventEvent.eventNum);
     }
     double mtWSyst[2] = {TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt()*(1.0 - cos(deltaPhiDileptonMet))),
                          TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt()*(1.0 - cos(deltaPhiDileptonMet)))};
     // Syst cuts for MVA
     bool passSystCuts[nSystTypes] = {
          passZMass && idJetUp.size() <= nJetsType  && passMET && passMT                                                                        								 && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && idJetDown.size()<= nJetsType && passMET && passMT                                                                        								 && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && passNjets && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() > metMIN && (mtWSyst[0] > mtMIN || ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() < metTIGHT) && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && passNjets && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > metMIN && (mtWSyst[1] > mtMIN || ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() < metTIGHT) && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto
     };
     if(MVAVarType==3) {
         // (Cuts in common between DY bin and signal region) && ((Exclusive BDT signal region cuts || Exclusive DY bin cuts))  
         passSystCuts[0] = (idJetUp.size()   <= nJetsType && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && passMETTight && bdt_value>=0) || (passZMass && passMET && !passMETTight && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
         passSystCuts[1] = (idJetDown.size() <= nJetsType && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && passMETTight && bdt_value>=0) || (passZMass && passMET && !passMETTight && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
         passSystCuts[2] = (passNjets && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() > metTIGHT && bdt_value>=0) || (passZMass && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() > metMIN && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() <= metTIGHT && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
         passSystCuts[3] = (passNjets && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > metTIGHT && bdt_value>=0) || (passZMass && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > metMIN && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() <= metTIGHT && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
     }
      
      // begin event weighting
      vector<int>zzBoson;
      vector<bool> isGenDupl;double bosonPtMin = 1000000000; bool isBosonFound = false;vector<bool> isNeuDupl;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23) zzBoson.push_back(ngen0);
        if((TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23||TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) &&
	   ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() < bosonPtMin) {bosonPtMin = ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt(); isBosonFound = true;}
        // begin neutrinos
        isNeuDupl.push_back(0);
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 12 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 14 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 16) isNeuDupl[ngen0] = 1;
	else {
          for(int ngen1=ngen0+1; ngen1<eventMonteCarlo.p4->GetEntriesFast(); ngen1++) {
	    if((int)(*eventMonteCarlo.pdgId)[ngen0] != (int)(*eventMonteCarlo.pdgId)[ngen1]) continue;
            if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen1])) < 0.02) {
	      isNeuDupl[ngen0] = 1;
	      break;
	    }
          }
        }
	// begin leptons	
        isGenDupl.push_back(0);
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) isGenDupl[ngen0] = 1;
	else {
          for(int ngen1=ngen0+1; ngen1<eventMonteCarlo.p4->GetEntriesFast(); ngen1++) {
	    if((int)(*eventMonteCarlo.pdgId)[ngen0] != (int)(*eventMonteCarlo.pdgId)[ngen1]) continue;
            if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen1])) < 0.02) {
	      isGenDupl[ngen0] = 1;
	      break;
	    }
          }
	}
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
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.1) {
	    isGenLepton = true;
	    break;
	  }
	}
	if(isGenLepton == true) {isGenLep.push_back(1); goodIsGenLep++;}
	else                    {isGenLep.push_back(0);}
      }

      // trigger efficiency
      double trigEff = 1.0;
      if(infilecatv[ifile] != 0) { 
        if(typePair==1) {
          int nbin = fhDEffDimuonSoup_pt0->FindBin( TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) , TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) );
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  40) trigEff=fhDEffDimuonSoup_pt0->GetBinContent(nbin);
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() >= 40) trigEff=fhDEffDimuonSoup_pt1->GetBinContent(nbin);
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() >= 40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  40) trigEff=fhDEffDimuonSoup_pt2->GetBinContent(nbin);
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() >= 40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() >= 40) trigEff=fhDEffDimuonSoup_pt3->GetBinContent(nbin);
        } else if(typePair==2) {
          int nbin = fhDEffDielectronSoup_pt0->FindBin( TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) , TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) );
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  40) trigEff=fhDEffDielectronSoup_pt0->GetBinContent(nbin);
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() >= 40) trigEff=fhDEffDielectronSoup_pt1->GetBinContent(nbin);
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() >= 40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  40) trigEff=fhDEffDielectronSoup_pt2->GetBinContent(nbin);
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() >= 40 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() >= 40) trigEff=fhDEffDielectronSoup_pt3->GetBinContent(nbin);
        }
      }
      //if(infilecatv[ifile] != 0) {
      //  trigEff = trigLookup.GetExpectedTriggerEfficiency(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),
      //  						  ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),
      //  						 TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]));
      //}
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
		typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,false);
        }
      }

      // fake rate
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        printf("NEED TO WORK ON IT IF WE WANT TO USE IT\n");return;
        if     (theCategory == 5){ // remove W+jets from MC
          fakeSF = 0.0;
        }
        else if(theCategory == 2 && goodIsTight != idTight.size()){ // remove Z+jets from MC as fakeable objects
          fakeSF = 0.0;
        }
        else if((infilecatv[ifile] == 0 || infilecatv[ifile] == 6 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add W+jets from data
          for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    effSF = effSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    theCategory = 5;
          }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) effSF =  1.0 * effSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) effSF = -1.0 * effSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) effSF = -1.0 * effSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) effSF =  1.0 * effSF; // single fake, data
        }
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 6 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 6){ // data or W+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infilecatv[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff;
      //printf("totalWeight: %f * %f * %f * %f * %f * %f * %f = %f\n",mcWeight,theLumi,puWeight,effSF,fakeSF,theMCPrescale,trigEff,totalWeight);

      //if(totalWeight == 0) continue;

      // DY
      if(theCategory == 2 && zzBoson.size() >= 1) {
        double myptz = ((TLorentzVector*)(*eventMonteCarlo.p4)[zzBoson[0]])->Pt();
        Int_t binpt = fhDVjetsNum->GetXaxis()->FindBin(myptz);
	if     (binpt <= 0) binpt = 1;
	else if(binpt > fhDVjetsNum->GetNbinsX()) binpt = fhDVjetsNum->GetNbinsX();
        totalWeight = totalWeight * fhDVjetsNum->GetBinContent(binpt)/fhDVjetsDen->GetBinContent(binpt) ;
      }

      // ZZ
      double the_rho = 0.0; if(the_rhoP4.P() > 0) the_rho = the_rhoP4.Pt()/the_rhoP4.P();
      double theZZCorr[2] {1,1};
      if(theCategory == 4 && infilenamev[ifile].Contains("GluGlu") == kFALSE) {
	theZZCorr[0] = weightEWKCorr(bosonPtMin,1);

        //float GENdPhiZZ = 5;
	//if(zzBoson.size() >= 2) GENdPhiZZ = TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[zzBoson[0]])->DeltaPhi(*((TLorentzVector*)(*eventMonteCarlo.p4)[zzBoson[1]])));
	//theZZCorr[1] = kfactor_qqZZ_qcd_dPhi(GENdPhiZZ);
        float GENmZZ = 0.0;
	if(zzBoson.size() >= 2) GENmZZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[1])) ) ).M();
	theZZCorr[1] = kfactor_qqZZ_qcd_M(GENmZZ);
        //float GENptZZ = 0.0;
	//if(zzBoson.size() >= 2) GENptZZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[1])) ) ).Pt();
	//theZZCorr[1] = kfactor_qqZZ_qcd_M(GENptZZ);

        totalWeight = totalWeight * (theZZCorr[0]*theZZCorr[1]);
      }
      // end event weighting
      //totalWeight = 1;

      if(makeMVAtrees) {
        // Save values for MVA trees
        mva_balance             = ptFrac;
        mva_delphi_ptll_MET     = dPhiDiLepMET; 
        mva_cos_theta_star_l1   = cos_theta_star_l1;
        mva_cos_theta_CS_l1     = cos_theta_collins_soper(*(TLorentzVector*)(*eventLeptons.p4)[idLep[0]],*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]);
        mva_deltaR_ll           = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaR(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]])); 
        mva_delphi_ll           = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]])); 
        mva_delphi_jet_MET      = dPhiJetMET;
        mva_etal1               = ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta();
        mva_etal2               = ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta();
        mva_etall               = dilep.Eta(); 
        mva_MET                 = ((TLorentzVector*)(*eventMet.p4)[0])->Pt(); 
        mva_mll_minus_mZ        = TMath::Abs(dilep.M() - 91.1876); 
        mva_mTjetMET            = mTJetMET;
        mva_mTll                = mtW; 
        mva_mTl1MET             = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])))))); 
        mva_mTl2MET             = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])))))); 
        mva_ptll                = dilep.Pt(); 
        mva_ptl1                = ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(); 
        mva_ptl2                = ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(); 
        mva_ptl1mptl2_over_ptll = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() - ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt()) / dilep.Pt();
        mva_weight              = totalWeight; 
        mva_njets               = idJet.size(); 
        mva_ntaus               = (unsigned char) numberGoodTaus; 
        mva_btag_veto           = passBtagVeto; 
        mva_3lveto              = pass3rdLVeto;
      }
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[SIGSEL]) sumEventsProcess[ifile] += totalWeight;

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
	    if     (thePlot ==  0 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(mtW,999.999);}
	    else if(thePlot ==  1 && passNMinusOne[1])       {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.M()-91.1876),99.999);}
	    else if(thePlot ==  2 && passNMinusOne[2])       {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	    else if(thePlot ==  3 && passNMinusOne[3])       {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),399.999);}
	    else if(thePlot ==  4 && passNMinusOne[4])       {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
	    else if(thePlot ==  5 && passNMinusOne[5])       {makePlot = true;theVar = dPhiDiLepMET;}
	    else if(thePlot ==  6 && passNMinusOne[6])       {makePlot = true;theVar = TMath::Max(TMath::Min(bDiscrMax,0.999),0.001);}
	    else if(thePlot ==  7 && passNMinusOne[7])       {makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
	    else if(thePlot ==  8 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
	    else if(thePlot ==  9 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
	    else if(thePlot == 10 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min((double)eventEvent.rho,39.999);}
	    else if(thePlot == 11 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	    else if(thePlot == 12 && passNMinusOne[9])       {makePlot = true;theVar = dPhiJetMET;}
	    else if(thePlot == 13 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = dPhiLepMETMin;}
	    else if(thePlot == 14 && passNMinusOne[8])       {makePlot = true;theVar = dphill;}
	    else if(thePlot == 15 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	    else if(thePlot == 16 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
	    else if(thePlot == 17 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
	    else if(thePlot == 18 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = (double)(numberGoodGenLep[0]+10*numberGoodGenLep[1]+100*numberGoodGenLep[2]);}
	    else if(thePlot == 19 && passNMinusOne[10])      {makePlot = true;theVar = TMath::Min((double)numberGoodTaus,3.499);}
	    else if(thePlot == 20 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.Eta()),2.499);}
	    else if(thePlot == 21 && passNMinusOne[0])       {makePlot = true;theVar = TMath::Min(mtW,999.999);}
	    else if(thePlot == 22 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(the_rho,0.999);}
	    else if(thePlot == 23 && passAllCuts[DYSANESEL1]){makePlot = true;theVar = TMath::Min(TMath::Max(dPhiJetMET,-0.05),3.099);}
	    else if(thePlot == 24 && passAllCuts[DYSANESEL1]){makePlot = true;theVar = TMath::Min(caloMinusPFMETRel,1.999);}
	    else if(thePlot == 25 && passAllCuts[DYSANESEL1]){makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),499.999);}
	    else if(thePlot == 26 && passAllCuts[DYSANESEL2]){makePlot = true;theVar = TMath::Min(TMath::Max(dPhiJetMET,-0.05),3.099);}
	    else if(thePlot == 27 && passAllCuts[DYSANESEL2]){makePlot = true;theVar = TMath::Min(caloMinusPFMETRel,1.999);}
	    else if(thePlot == 28 && passAllCuts[DYSANESEL2]){makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),499.999);}
	    else if(thePlot == 29 && passZMass)              {makePlot = true;theVar = idJet.size();}
	    else if(thePlot == 30 && passZMass)              {makePlot = true;theVar = idBJet.size();}
	    else if(thePlot == 31 && passZMass)              {makePlot = true;theVar = idLep.size();;}
	    else if(thePlot == 32 && passZMass)              {makePlot = true;theVar = TMath::Min(TMath::Max(dPhiDiLepMET,-0.05),3.099);}
	    else if(thePlot == 33 && passNMinusOne[3])       {makePlot = true;theVar = (double)((TLorentzVector*)(*eventMet.p4)[0])->Pt();}
	    else if(thePlot == 34 && passNMinusOne[8])       {makePlot = true;theVar = TMath::Min(drll,2.999);}
	    else if(thePlot == 35 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(dilep.Pt()/mtW,0.999);}
            else if(thePlot == 36 && passAllCuts[TIGHTSEL] && ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 150) {makePlot=true;theVar = TMath::Min(1., TMath::Max(-1.,bdt_value));}
	    if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
	  }
        }
      }

      if(typeSel == typePair || typeSel == 3) {
	double MVAVar = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, ((TLorentzVector*)(*eventMet.p4)[0])->Pt(), mtW, dilep.M(), bdt_value, xbins[nBinMVA]);
	double MVAVarMETSyst[2] = {MVAVar, MVAVar};	
        if(MVAVarType==1) {
          MVAVarMETSyst[0] = TMath::Min(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt(),xbins[nBinMVA]-0.001);
          MVAVarMETSyst[1] = TMath::Min(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt(),xbins[nBinMVA]-0.001);
        }

	// Avoid QCD scale weights that are anomalous high
	double maxQCDscale = (TMath::Abs((double)eventMonteCarlo.r1f2)+TMath::Abs((double)eventMonteCarlo.r1f5)+TMath::Abs((double)eventMonteCarlo.r2f1)+
                              TMath::Abs((double)eventMonteCarlo.r2f2)+TMath::Abs((double)eventMonteCarlo.r5f1)+TMath::Abs((double)eventMonteCarlo.r5f5))/6.0;
        // Throw O(1000) toys for the BDT nuisance evaluations
        if(useBDT && passAllCuts[SIGSEL] && !useCachedBDTSystematics) {
          double bdt_toy_value_muonScale, bdt_toy_value_electronScale, bdt_toy_value_METScale, MVAVar_toy_muonScale, MVAVar_toy_electronScale, MVAVar_toy_METScale;
          TLorentzVector lepton1 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[0]]),
                         lepton2 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[1]]),
                         MET     = *((TLorentzVector*)(*eventMet.p4)[0]),
                         jet1    = idJet.size() > 0 ? *((TLorentzVector*)(*eventJets.p4)[idJet[0]]) : TLorentzVector(0,0,0,0);
          double lepton1_scale_variation, lepton2_scale_variation;
          for(unsigned int i_toy=0; i_toy<num_bdt_toys; i_toy++) {
            // BDT variation with the muon scale variation (flat 1%)
            bdt_toy_value_muonScale = bdt_value; MVAVar_toy_muonScale = MVAVar;
            if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13 || TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) { // don't perform expensive BDT evaluation unless there is a muon
              lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13 ? 0.01 * bdt_toy_scale[i_toy] : 0;
              lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13 ? 0.01 * bdt_toy_scale[i_toy] : 0;
              bdt_toy_value_muonScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation, 0, 0);
  	      MVAVar_toy_muonScale = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, ((TLorentzVector*)(*eventMet.p4)[0])->Pt(), mtW, dilep.M(), bdt_toy_value_muonScale, xbins[nBinMVA]);
            }
            // BDT variation with the electron scale variation (flat 1%)
            bdt_toy_value_electronScale = bdt_value; MVAVar_toy_electronScale = MVAVar;
            if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11 || TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) { // don't perform expensive BDT evaluation unless there is an electron
              lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11 ? 0.01 * bdt_toy_scale[i_toy] : 0;
              lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11 ? 0.01 * bdt_toy_scale[i_toy] : 0;
              bdt_toy_value_electronScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation, 0, 0);
    	      MVAVar_toy_electronScale = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, ((TLorentzVector*)(*eventMet.p4)[0])->Pt(), mtW, dilep.M(), bdt_toy_value_electronScale, xbins[nBinMVA]);
            }
            // BDT variation with the MET scale variation (from Nero)
            if(bdt_toy_scale[i_toy] >= 0) 
              bdt_toy_value_METScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, bdt_toy_scale[i_toy] * (((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt() / MET.Pt() - 1.), 0);
            else
              bdt_toy_value_METScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, -bdt_toy_scale[i_toy] * (((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])  ->Pt() / MET.Pt() - 1.), 0);
  	        MVAVar_toy_METScale = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, ((TLorentzVector*)(*eventMet.p4)[0])->Pt(), mtW, dilep.M(), bdt_toy_value_METScale, xbins[nBinMVA]);
            // debug
            //if(passAllCuts[SIGSEL] && theCategory == 4 && typePair!=3) printf("ZZ event METscale toy %d changes BDT value / MVA var (%f, %f) => (%f, %f)\n", i_toy, bdt_value, MVAVar, bdt_toy_value_METScale, MVAVar_toy_METScale);



            if(theCategory == 3)      { histo_bdt_toys_muonScale_WZ       ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_WZ       ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_WZ       ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 4) { histo_bdt_toys_muonScale_ZZ       ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_ZZ       ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_ZZ       ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 5) { histo_bdt_toys_muonScale_VVV      ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_VVV      ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_VVV      ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 6) { histo_bdt_toys_muonScale_ZH_hinv[nModel]  ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_ZH_hinv[nModel]  ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_ZH_hinv[nModel]  ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 7) { histo_bdt_toys_muonScale_ggZH_hinv->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_ggZH_hinv->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_ggZH_hinv->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
          }
        }

        if     (theCategory == 0){
	  if(passAllCuts[SIGSEL]) histo_Data->Fill(MVAVar,totalWeight);
	  if(passAllCuts[TIGHTSEL] && verbose) {
            if      (typePair==1) printf("mm data event: ");
            else if (typePair==2) printf("ee data event: " );
            else if (typePair==3) printf("em data event: " );
            printf("runnumber %d lumisection %d eventnumber %lld ptll %f MET %f njets %d l1_pt %f l1_eta %f l2_pt %f l2_eta %f balance %f\n",
              eventEvent.runNum,
              eventEvent.lumiNum,
              eventEvent.eventNum,
              dilep.Pt(),
              ((TLorentzVector*)(*eventMet.p4)[0])->Pt(),
              (int)idJet.size(),
              ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(), ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),
              ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(), ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),
              ptFrac
            );
          }
          if(makeMVAtrees) data_mva_tree->Fill();
        }
        else if(theCategory == 1){
	  if(passAllCuts[SIGSEL]) histo_EM   ->Fill(MVAVar,totalWeight);
	  if(passAllCuts[SIGSEL]) histo_EMNoW->Fill(MVAVar,1.);
          if(makeMVAtrees) EM_mva_tree->Fill();
        }
        else if(theCategory == 2){
	  if(passAllCuts[SIGSEL]) histo_Zjets   ->Fill(MVAVar,totalWeight);
	  if(passAllCuts[SIGSEL]) histo_ZjetsNoW->Fill(MVAVar,1.);
          if(makeMVAtrees) Zjets_mva_tree->Fill();
        }
        else if(theCategory == 3){
	  if(passAllCuts[SIGSEL]) {
	     histo_WZ              ->Fill(MVAVar,totalWeight);
	     histo_WZNoW           ->Fill(MVAVar,1.);
	     histo_WZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*1.10);
	     histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_WZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_WZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
          if(makeMVAtrees) WZ_mva_tree->Fill();
	}
        else if(theCategory == 4){
	  if(passAllCuts[SIGSEL]) {
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
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_ZZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
          }
          if(passSystCuts[JESUP])  histo_ZZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
          if(makeMVAtrees) ZZ_mva_tree->Fill();
        }
        else if(theCategory == 5){
	  if(passAllCuts[SIGSEL]) {
	     histo_VVV   ->Fill(MVAVar,totalWeight);
	     histo_VVVNoW->Fill(MVAVar,1.);
	     histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
          }
          if(passSystCuts[JESUP])  histo_VVV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VVV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
          if(makeMVAtrees) VVV_mva_tree->Fill();
        }
        else if(theCategory == 6){
	  if(passAllCuts[SIGSEL]) {
	     histo_ZH_hinv[nModel]   ->Fill(MVAVar,totalWeight);
	     histo_ZH_hinvNoW[nModel]->Fill(MVAVar,1.);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ZH_hinv_CMS_QCDScaleBounding[nModel][5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_ZH_hinv_CMS_PDFBounding[nModel][npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_ZH_hinv_CMS_PDFBounding[nModel][npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingAvg [nModel]->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingAvg [nModel]->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingUp  [nModel]->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingUp  [nModel]->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Fill(MVAVar,totalWeight*0.98);
             histo_ZH_hinv_CMS_PUBoundingUp  [nModel]->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZH_hinv_CMS_PUBoundingDown[nModel]->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_ZH_hinv_CMS_MVAJESBoundingUp  [nModel]->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZH_hinv_CMS_MVAMETBoundingUp  [nModel]->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->Fill(MVAVarMETSyst[1],totalWeight);
          if(makeMVAtrees) signal_mva_trees[nModel]->Fill();
        }
        else if(theCategory == 7){
	  if(passAllCuts[SIGSEL]) {
	     histo_ggZH_hinv   ->Fill(MVAVar,totalWeight);
	     histo_ggZH_hinvNoW->Fill(MVAVar,1.0);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_ggZH_hinv_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_ggZH_hinv_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_ggZH_hinv_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ggZH_hinv_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ggZH_hinv_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ggZH_hinv_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
          if(makeMVAtrees) signal_mva_trees[nModel]->Fill();
        }
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      }
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);
    for(int nc=0; nc<numberCuts+1; nc++){
      printf("(%20s): %10.2f %10.2f %10.2f %10.2f\n",cutName[nc].Data(),histoZHSEL[0]->GetBinContent(nc+1),histoZHSEL[1]->GetBinContent(nc+1),histoZHSEL[2]->GetBinContent(nc+1),histoZHSEL[3]->GetBinContent(nc+1));
    }
    the_input_file->Close();
  } // end of chain
  if(makeMVAtrees) mva_trees->Close();

  // "-1" to remove the Higgs contribution
  double sumEvents = 0;
  for(int np=1; np<histBins-1; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  //printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  //printf(" = %.3f\n",sumEvents);

  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Printing yields and statistical uncertainties for all signal models\n\n");
  for(int nModel=0; nModel<nSigModels; nModel++) {
    if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
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
  double kEoverM  = sqrt(bgdDecay[0][ZSEL+nSelTypes*2][0]/bgdDecay[0][ZSEL+nSelTypes*1][0]);
  double NemFact[3] = {(1.0/kEoverM)*0.5, (kEoverM)*0.5, (kEoverM+1.0/kEoverM)*0.5};

  // This is the default method
  double NemFact_FromMLLSB[3] = {NemFact[0], NemFact[1], NemFact[2]}; double NemFact_FromMLLSBE[3] = {1.0, 1.0, 1.0};
  if(bgdDecay[0][WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWLOOSESEL+nSelTypes*1][0]) > 0) {
    NemFact_FromMLLSB[0]  = (bgdDecay[0][WWLOOSESEL+nSelTypes*1][0])/bgdDecay[0][WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[0] = sqrt(1./bgdDecay[0][WWLOOSESEL+nSelTypes*1][0]+1./bgdDecay[0][WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[0];
  }
  if(bgdDecay[0][WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWLOOSESEL+nSelTypes*2][0]) > 0) {
    NemFact_FromMLLSB[1]  = (bgdDecay[0][WWLOOSESEL+nSelTypes*2][0])/bgdDecay[0][WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[1] = sqrt(1./bgdDecay[0][WWLOOSESEL+nSelTypes*2][0]+1./bgdDecay[0][WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[1];
  }
  if(bgdDecay[0][WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWLOOSESEL+nSelTypes*3][0]) > 0) {
    NemFact_FromMLLSB[2]  = (bgdDecay[0][WWLOOSESEL+nSelTypes*3][0])/bgdDecay[0][WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[2] = sqrt(1./bgdDecay[0][WWLOOSESEL+nSelTypes*3][0]+1./bgdDecay[0][WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[2];
  }
  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Computing the flavor k-factors used for the EM background\n\n");

  printf("(mm) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[0],NemFact_FromMLLSB[0],NemFact_FromMLLSBE[0]);
  printf("(ee) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[1],NemFact_FromMLLSB[1],NemFact_FromMLLSBE[1]);
  printf("(ll) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[2],NemFact_FromMLLSB[2],NemFact_FromMLLSBE[2]);

  // There uncertainties: closure test, different between default and alternative method, and data statistics
  // Only the second item is used as systematic uncertainty
  double EMSystTotal[3] = {1.0,1.0,1.0}; double EMSyst[2][3] = {0.0,0.0,0.0,0.0,0.0,0.0};

  EMSyst[0][0] = bgdDecay[0][SIGSEL+nSelTypes*1][1]/bgdDecay[0][SIGSEL+nSelTypes*0][1]/NemFact[0];
  if(EMSyst[0][0] < 1.0) EMSyst[0][0] = 1/EMSyst[0][0]; EMSyst[0][0] = EMSyst[0][0] - 1.0;

  EMSyst[0][1] = bgdDecay[0][SIGSEL+nSelTypes*2][1]/bgdDecay[0][SIGSEL+nSelTypes*0][1]/NemFact[1];
  if(EMSyst[0][1] < 1.0) EMSyst[0][1] = 1/EMSyst[0][1]; EMSyst[0][1] = EMSyst[0][1] - 1.0;

  EMSyst[0][2] = bgdDecay[0][SIGSEL+nSelTypes*3][1]/bgdDecay[0][SIGSEL+nSelTypes*0][1]/NemFact[2];
  if(EMSyst[0][2] < 1.0) EMSyst[0][2] = 1/EMSyst[0][2]; EMSyst[0][2] = EMSyst[0][2] - 1.0;

  EMSyst[1][0] = NemFact_FromMLLSB[0]/NemFact[0];
  if(EMSyst[1][0] < 1.0) EMSyst[1][0] = 1/EMSyst[1][0]; EMSyst[1][0] = EMSyst[1][0] - 1.0;

  EMSyst[1][1] = NemFact_FromMLLSB[1]/NemFact[1];
  if(EMSyst[1][1] < 1.0) EMSyst[1][1] = 1/EMSyst[1][1]; EMSyst[1][1] = EMSyst[1][1] - 1.0;

  EMSyst[1][2] = NemFact_FromMLLSB[2]/NemFact[2];
  if(EMSyst[1][2] < 1.0) EMSyst[1][2] = 1/EMSyst[1][2]; EMSyst[1][2] = EMSyst[1][2] - 1.0;

  if(bgdDecay[0][SIGSEL][0] > 0) EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1/bgdDecay[0][SIGSEL][0]);
  else                        EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1.0);

  if(bgdDecay[0][SIGSEL][0] > 0) EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1/bgdDecay[0][SIGSEL][0]);
  else                        EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1.0);

  if(bgdDecay[0][SIGSEL][0] > 0) EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1/bgdDecay[0][SIGSEL][0]);
  else                        EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1.0);

  double EMbkg = bgdDecay[0][SIGSEL][2]+bgdDecay[0][SIGSEL][3]+bgdDecay[0][SIGSEL][4]+bgdDecay[0][SIGSEL][5];

  printf("(mm) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][SIGSEL+nSelTypes*1][1],sqrt(weiDecay[0][SIGSEL+nSelTypes*1][4]),
	 bgdDecay[0][SIGSEL][1]*NemFact[0]  ,sqrt(weiDecay[0][SIGSEL][4])*NemFact[0],bgdDecay[0][SIGSEL][0],EMbkg,EMSystTotal[0],EMSyst[0][0],EMSyst[1][0],sqrt(EMSystTotal[0]*EMSystTotal[0] - EMSyst[0][0]*EMSyst[0][0] - EMSyst[1][0]*EMSyst[1][0]));

  printf("(ee) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][SIGSEL+nSelTypes*2][1],sqrt(weiDecay[0][SIGSEL+nSelTypes*2][4]),
	 bgdDecay[0][SIGSEL][1]*NemFact[1]  ,sqrt(weiDecay[0][SIGSEL][4])*NemFact[1],bgdDecay[0][SIGSEL][0],EMbkg,EMSystTotal[1],EMSyst[0][1],EMSyst[1][1],sqrt(EMSystTotal[1]*EMSystTotal[1] - EMSyst[0][1]*EMSyst[0][1] - EMSyst[1][1]*EMSyst[1][1]));

  printf("(ll) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][SIGSEL+nSelTypes*3][1],sqrt(weiDecay[0][SIGSEL+nSelTypes*3][4]),
	 bgdDecay[0][SIGSEL][1]*NemFact[2]  ,sqrt(weiDecay[0][SIGSEL][4])*NemFact[2],bgdDecay[0][SIGSEL][0],EMbkg,EMSystTotal[2],EMSyst[0][2],EMSyst[1][2],sqrt(EMSystTotal[2]*EMSystTotal[2] - EMSyst[0][2]*EMSyst[0][2] - EMSyst[1][2]*EMSyst[1][2]));

  double systEM[2] = {1.0, 1.0};
  if(histo_EM->GetSumOfWeights() > 1) systEM[0] = 1. + 1./sqrt(histo_EM->GetSumOfWeights());
  else  			      systEM[0] = 2.;
  if(useEMFromData == true){
    EMbkg = bgdDecay[0][TIGHTSEL][2]+bgdDecay[0][TIGHTSEL][3]+bgdDecay[0][TIGHTSEL][4]+bgdDecay[0][TIGHTSEL][5];
    double EMNormFact[4] = {((bgdDecay[0][TIGHTSEL][0]-EMbkg)*NemFact[0])/bgdDecay[0][TIGHTSEL+nSelTypes*(1)][1],
                            ((bgdDecay[0][TIGHTSEL][0]-EMbkg)*NemFact[1])/bgdDecay[0][TIGHTSEL+nSelTypes*(2)][1],
                            ((bgdDecay[0][TIGHTSEL][0]-EMbkg)*NemFact[2])/bgdDecay[0][TIGHTSEL+nSelTypes*(3)][1],
                            ((bgdDecay[0][TIGHTSEL][0]-EMbkg)*1.00000000)/bgdDecay[0][TIGHTSEL+nSelTypes*(0)][1]};

    printf("EM(1): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*1][1],bgdDecay[0][TIGHTSEL][0],EMbkg,NemFact[0],bgdDecay[0][TIGHTSEL+nSelTypes*(1)][1],bgdDecay[0][TIGHTSEL+nSelTypes*1][1]*EMNormFact[0],bgdDecay[0][TIGHTSEL+nSelTypes*1][1]*EMNormFact[0]*EMSystTotal[0]);
    printf("EM(2): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*2][1],bgdDecay[0][TIGHTSEL][0],EMbkg,NemFact[1],bgdDecay[0][TIGHTSEL+nSelTypes*(2)][1],bgdDecay[0][TIGHTSEL+nSelTypes*2][1]*EMNormFact[1],bgdDecay[0][TIGHTSEL+nSelTypes*2][1]*EMNormFact[1]*EMSystTotal[1]);
    printf("EM(3): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*3][1],bgdDecay[0][TIGHTSEL][0],EMbkg,NemFact[2],bgdDecay[0][TIGHTSEL+nSelTypes*(3)][1],bgdDecay[0][TIGHTSEL+nSelTypes*3][1]*EMNormFact[2],bgdDecay[0][TIGHTSEL+nSelTypes*3][1]*EMNormFact[2]*EMSystTotal[2]);
    printf("EM(4): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*0][1],bgdDecay[0][TIGHTSEL][0],EMbkg,1.00000000,bgdDecay[0][TIGHTSEL+nSelTypes*(0)][1],bgdDecay[0][TIGHTSEL+nSelTypes*0][1]*EMNormFact[3],bgdDecay[0][TIGHTSEL+nSelTypes*0][1]*EMNormFact[3]*0.000000000000);

    //systEM[0] = 1. + EMSystTotal[typeSel-1];
    //systEM[0] = 1. + EMSyst[1][typeSel-1];
    systEM[0] = 1. + 0.20;
    systEM[1] = TMath::Max(bgdDecay[0][TIGHTSEL][0],1.0);
    if(bgdDecay[0][TIGHTSEL][0] > 0) {
      histo_EM->Scale(EMNormFact[typeSel-1]);
      histo_EM->SetBinContent(1,histo_EM->GetBinContent(1)*EMNormFact[3]/EMNormFact[typeSel-1]);
    }
  }

  // DY background estimation
  if(useZjetsTemplate){
    histo_Zjets_CMS_ZjetsSystUp->Add(fhDZjets);
    histo_Zjets->Scale(0.0);
    histo_Zjets->Add(fhDZjets);
    histo_ZjetsNoW->Scale(0.0);
    histo_ZjetsNoW->Add(fhDZjets);
    double ZJetsNorm[3] = {0.56, 0.52, 1.};
    if     (MVAVarType == 0 || MVAVarType == 1 ||MVAVarType == 2 || MVAVarType == 3) {
      ZJetsNorm[0] = 1.; ZJetsNorm[1] = 1.; ZJetsNorm[2] = 1.; 
    }
    if     (nJetsType == 0) histo_Zjets->Scale(TMath::Abs(ZJetsNorm[0]/histo_Zjets->GetSumOfWeights()));
    else if(nJetsType == 1) histo_Zjets->Scale(TMath::Abs(ZJetsNorm[1]/histo_Zjets->GetSumOfWeights()));
    else                    histo_Zjets->Scale(TMath::Abs(ZJetsNorm[2]/histo_Zjets->GetSumOfWeights()));
    
    histo_Zjets_CMS_ZjetsSystUp->Scale(0.0);
    histo_Zjets_CMS_ZjetsSystUp->Add(fhDZjetsSyst);
    if     (nJetsType == 0) histo_Zjets_CMS_ZjetsSystUp->Scale(TMath::Abs(ZJetsNorm[0]/histo_Zjets_CMS_ZjetsSystUp->GetSumOfWeights()));
    else if(nJetsType == 1) histo_Zjets_CMS_ZjetsSystUp->Scale(TMath::Abs(ZJetsNorm[1]/histo_Zjets_CMS_ZjetsSystUp->GetSumOfWeights()));
    else                    histo_Zjets_CMS_ZjetsSystUp->Scale(TMath::Abs(ZJetsNorm[2]/histo_Zjets_CMS_ZjetsSystUp->GetSumOfWeights()));
  }

  // computing DY scale factor using the first bin
  if(MVAVarType == 0 || MVAVarType == 1 || MVAVarType == 2 || MVAVarType==3 || MVAVarType==4){
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("Computing the Drell-Yan data/MC scale factor using the first shape bin\n\n");
    int theBin = 2;
    double the_bck = histo_EM->GetBinContent(theBin) + histo_WZ->GetBinContent(theBin) + histo_ZZ->GetBinContent(theBin) + histo_VVV->GetBinContent(theBin);
    double the_data = histo_Data->GetBinContent(theBin);
    double the_sf;
    //if(MVAVarType!=3) { 
    the_sf = (the_data-the_bck)/histo_Zjets->GetBinContent(theBin);
    printf("DY SF: data/bck/DY = %f %f %f ==> %f\n",the_data,the_bck,histo_Zjets->GetBinContent(theBin),the_sf);
    //} else {
    //  the_sf=1.068889;
    //}
    histo_Zjets->Scale(the_sf);
    histo_Zjets_CMS_ZjetsSystUp->Scale(the_sf);
  }

  double nZjetsNorm = histo_Zjets->GetSumOfWeights();
  if(nZjetsNorm > 0){
    for(int i=1; i<=histo_Zjets->GetNbinsX(); i++){
      histo_Zjets->SetBinContent(i,TMath::Max(histo_Zjets->GetBinContent(i),0.000001));
    }
    histo_Zjets->Scale(nZjetsNorm/histo_Zjets->GetSumOfWeights());
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
  
  printf("EWK Corr: WZ(%f/%f/%f) ZZ(%f/%f/%f) ggZZ(%f/%f/%f)\n",
                     histo_WZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_WZ->GetSumOfWeights(),histo_WZ_CMS_EWKCorrDown->GetSumOfWeights(),
                     histo_ZZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_EWKCorrDown->GetSumOfWeights(),
		     histo_ZZ_CMS_ggCorrUp ->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_ggCorrDown ->GetSumOfWeights());

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
    sprintf(output,"MitZHAnalysis/plots%s/histo%szh%s_nice_%s_%d.root",subdirectory.c_str(),addChan.Data(),finalStateName, signalName_[plotModel].Data(),thePlot);	  
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
      if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
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
    if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
    sprintf(outputLimits,"MitZHAnalysis/plots%s/zll%shinv%s_%s_input_%s.root", subdirectory.c_str(), addChan.Data(), finalStateName, signalName_[nModel].Data(), ECMsb.Data());
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
      printf("uncertainties Zjets\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystUp  ->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystDown->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
    } 

    if(useBDT && !useCachedBDTSystematics) {
      for(int nb=1; nb<=nBinMVA; nb++) {
        for(unsigned int i_toy=0; i_toy<num_bdt_toys; i_toy++) {
          int nbin2d = histo_bdt_toys_electronScale_VVV->GetBin(nb, i_toy);
          if(histo_VVV            ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_VVV      ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_VVV            ->GetBinContent(nbin2d) / histo_VVV            ->GetBinContent(nb));
          if(histo_WZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_WZ       ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_WZ             ->GetBinContent(nbin2d) / histo_WZ             ->GetBinContent(nb));
          if(histo_ZZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_ZZ       ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_ZZ             ->GetBinContent(nbin2d) / histo_ZZ             ->GetBinContent(nb));
          if(histo_ZH_hinv[nModel]->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_ZH_hinv[nModel]  ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_ZH_hinv[nModel]->GetBinContent(nbin2d) / histo_ZH_hinv[nModel]->GetBinContent(nb));
          if(histo_ggZH_hinv      ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_ggZH_hinv->Fill( xbins[nb-1], histo_bdt_toys_electronScale_ggZH_hinv      ->GetBinContent(nbin2d) / histo_ggZH_hinv      ->GetBinContent(nb));
          if(histo_VVV            ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_VVV      ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_VVV            ->GetBinContent(nbin2d) / histo_VVV            ->GetBinContent(nb));
          if(histo_WZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_WZ       ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_WZ             ->GetBinContent(nbin2d) / histo_WZ             ->GetBinContent(nb));
          if(histo_ZZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_ZZ       ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_ZZ             ->GetBinContent(nbin2d) / histo_ZZ             ->GetBinContent(nb));
          if(histo_ZH_hinv[nModel]->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_ZH_hinv[nModel]  ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_ZH_hinv[nModel]->GetBinContent(nbin2d) / histo_ZH_hinv[nModel]->GetBinContent(nb));
          if(histo_ggZH_hinv      ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_ggZH_hinv->Fill( xbins[nb-1], histo_bdt_toys_muonScale_ggZH_hinv      ->GetBinContent(nbin2d) / histo_ggZH_hinv      ->GetBinContent(nb));
          if(histo_VVV            ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_VVV      ->Fill( xbins[nb-1], histo_bdt_toys_METScale_VVV            ->GetBinContent(nbin2d) / histo_VVV            ->GetBinContent(nb));
          if(histo_WZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_WZ       ->Fill( xbins[nb-1], histo_bdt_toys_METScale_WZ             ->GetBinContent(nbin2d) / histo_WZ             ->GetBinContent(nb));
          if(histo_ZZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_ZZ       ->Fill( xbins[nb-1], histo_bdt_toys_METScale_ZZ             ->GetBinContent(nbin2d) / histo_ZZ             ->GetBinContent(nb));
          if(histo_ZH_hinv[nModel]->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_ZH_hinv[nModel]  ->Fill( xbins[nb-1], histo_bdt_toys_METScale_ZH_hinv[nModel]->GetBinContent(nbin2d) / histo_ZH_hinv[nModel]->GetBinContent(nb));
          if(histo_ggZH_hinv      ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_ggZH_hinv->Fill( xbins[nb-1], histo_bdt_toys_METScale_ggZH_hinv      ->GetBinContent(nbin2d) / histo_ggZH_hinv      ->GetBinContent(nb));
        }
        bdt_toy_binyields_electronScale_VVV[nb-1]             = (TH1F*) bdt_toy_envelope_electronScale_VVV            ->ProjectionY(Form("bdt_toy_binyields_electronScale_VVV_bin%d", nb)      , nb, nb);
        bdt_toy_binyields_electronScale_WZ[nb-1]              = (TH1F*) bdt_toy_envelope_electronScale_WZ             ->ProjectionY(Form("bdt_toy_binyields_electronScale_WZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_electronScale_ZZ[nb-1]              = (TH1F*) bdt_toy_envelope_electronScale_ZZ             ->ProjectionY(Form("bdt_toy_binyields_electronScale_ZZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_electronScale_ZH_hinv[nModel][nb-1] = (TH1F*) bdt_toy_envelope_electronScale_ZH_hinv[nModel]->ProjectionY(Form("bdt_toy_binyields_electronScale_ZH_hinv_%s_bin%d",signalName_[nModel].Data(),nb), nb, nb);
        bdt_toy_binyields_electronScale_ggZH_hinv[nb-1]       = (TH1F*) bdt_toy_envelope_electronScale_ggZH_hinv      ->ProjectionY(Form("bdt_toy_binyields_electronScale_ggZH_hinv_bin%d", nb), nb, nb);
        bdt_toy_binyields_muonScale_VVV[nb-1]                 = (TH1F*) bdt_toy_envelope_muonScale_VVV                ->ProjectionY(Form("bdt_toy_binyields_muonScale_VVV_bin%d", nb)      , nb, nb);
        bdt_toy_binyields_muonScale_WZ[nb-1]                  = (TH1F*) bdt_toy_envelope_muonScale_WZ                 ->ProjectionY(Form("bdt_toy_binyields_muonScale_WZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_muonScale_ZZ[nb-1]                  = (TH1F*) bdt_toy_envelope_muonScale_ZZ                 ->ProjectionY(Form("bdt_toy_binyields_muonScale_ZZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_muonScale_ZH_hinv[nModel][nb-1]     = (TH1F*) bdt_toy_envelope_muonScale_ZH_hinv[nModel]    ->ProjectionY(Form("bdt_toy_binyields_muonScale_ZH_hinv_%s_bin%d", signalName_[nModel].Data(),nb), nb, nb);
        bdt_toy_binyields_muonScale_ggZH_hinv[nb-1]           = (TH1F*) bdt_toy_envelope_muonScale_ggZH_hinv          ->ProjectionY(Form("bdt_toy_binyields_muonScale_ggZH_hinv_bin%d", nb), nb, nb);
        bdt_toy_binyields_METScale_VVV[nb-1]                  = (TH1F*) bdt_toy_envelope_METScale_VVV                 ->ProjectionY(Form("bdt_toy_binyields_METScale_VVV_bin%d", nb)      , nb, nb);
        bdt_toy_binyields_METScale_WZ[nb-1]                   = (TH1F*) bdt_toy_envelope_METScale_WZ                  ->ProjectionY(Form("bdt_toy_binyields_METScale_WZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_METScale_ZZ[nb-1]                   = (TH1F*) bdt_toy_envelope_METScale_ZZ                  ->ProjectionY(Form("bdt_toy_binyields_METScale_ZZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_METScale_ZH_hinv[nModel][nb-1]      = (TH1F*) bdt_toy_envelope_METScale_ZH_hinv[nModel]     ->ProjectionY(Form("bdt_toy_binyields_METScale_ZH_hinv_%s_bin%d", signalName_[nModel].Data(),nb), nb, nb);
        bdt_toy_binyields_METScale_ggZH_hinv[nb-1]            = (TH1F*) bdt_toy_envelope_METScale_ggZH_hinv           ->ProjectionY(Form("bdt_toy_binyields_METScale_ggZH_hinv_bin%d", nb), nb, nb);
        
        // Since not using cached results, finally calculate the bin-by-bin % systematics (up and down) for the BDT 
        double toy_mean_electronScale, toy_mean_muonScale, toy_mean_METScale;
        double toy_rms_electronScale, toy_rms_muonScale, toy_rms_METScale;
        if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 ) { 
          toy_mean_electronScale = bdt_toy_binyields_electronScale_ZH_hinv[nModel][nb-1]->GetMean();
          toy_mean_muonScale = bdt_toy_binyields_muonScale_ZH_hinv[nModel][nb-1]->GetMean();
          toy_mean_METScale = bdt_toy_binyields_METScale_ZH_hinv[nModel][nb-1]->GetMean();
          toy_rms_electronScale = bdt_toy_binyields_electronScale_ZH_hinv[nModel][nb-1]->GetRMS();
          toy_rms_muonScale = bdt_toy_binyields_muonScale_ZH_hinv[nModel][nb-1]->GetRMS();
          toy_rms_METScale = bdt_toy_binyields_METScale_ZH_hinv[nModel][nb-1]->GetRMS();
          bdt_syst_electronScaleUp_ZH_hinv[nModel]   ->SetBinContent(nb, 1. + 0.5 * (toy_mean_electronScale - 1. + toy_rms_electronScale));
          bdt_syst_electronScaleDown_ZH_hinv[nModel] ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_electronScale + toy_rms_electronScale));
          bdt_syst_muonScaleUp_ZH_hinv[nModel]       ->SetBinContent(nb, 1. + 0.5 * (toy_mean_muonScale - 1. + toy_rms_muonScale));
          bdt_syst_muonScaleDown_ZH_hinv[nModel]     ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_muonScale + toy_rms_muonScale));
          bdt_syst_METScaleUp_ZH_hinv[nModel]        ->SetBinContent(nb, 1. + 0.5 * (toy_mean_METScale - 1. + toy_rms_METScale));
          bdt_syst_METScaleDown_ZH_hinv[nModel]      ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_METScale + toy_rms_METScale));
        }
        if(histo_VVV->GetBinContent(nb) > 0 ) {
          toy_mean_electronScale = bdt_toy_binyields_electronScale_VVV[nb-1]->GetMean();
          toy_mean_muonScale = bdt_toy_binyields_muonScale_VVV[nb-1]->GetMean();
          toy_mean_METScale = bdt_toy_binyields_METScale_VVV[nb-1]->GetMean();
          toy_rms_electronScale = bdt_toy_binyields_electronScale_VVV[nb-1]->GetRMS();
          toy_rms_muonScale = bdt_toy_binyields_muonScale_VVV[nb-1]->GetRMS();
          toy_rms_METScale = bdt_toy_binyields_METScale_VVV[nb-1]->GetRMS();
          bdt_syst_electronScaleUp_VVV   ->SetBinContent(nb, 1. + 0.5 * (toy_mean_electronScale - 1. + toy_rms_electronScale));
          bdt_syst_electronScaleDown_VVV ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_electronScale + toy_rms_electronScale));
          bdt_syst_muonScaleUp_VVV       ->SetBinContent(nb, 1. + 0.5 * (toy_mean_muonScale - 1. + toy_rms_muonScale));
          bdt_syst_muonScaleDown_VVV     ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_muonScale + toy_rms_muonScale));
          bdt_syst_METScaleUp_VVV        ->SetBinContent(nb, 1. + 0.5 * (toy_mean_METScale - 1. + toy_rms_METScale));
          bdt_syst_METScaleDown_VVV      ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_METScale + toy_rms_METScale));
        }
        if(histo_WZ->GetBinContent(nb) > 0 ) {
          toy_mean_electronScale = bdt_toy_binyields_electronScale_WZ[nb-1]->GetMean();
          toy_mean_muonScale = bdt_toy_binyields_muonScale_WZ[nb-1]->GetMean();
          toy_mean_METScale = bdt_toy_binyields_METScale_WZ[nb-1]->GetMean();
          toy_rms_electronScale = bdt_toy_binyields_electronScale_WZ[nb-1]->GetRMS();
          toy_rms_muonScale = bdt_toy_binyields_muonScale_WZ[nb-1]->GetRMS();
          toy_rms_METScale = bdt_toy_binyields_METScale_WZ[nb-1]->GetRMS();
          bdt_syst_electronScaleUp_WZ   ->SetBinContent(nb, 1. + 0.5 * (toy_mean_electronScale - 1. + toy_rms_electronScale));
          bdt_syst_electronScaleDown_WZ ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_electronScale + toy_rms_electronScale));
          bdt_syst_muonScaleUp_WZ       ->SetBinContent(nb, 1. + 0.5 * (toy_mean_muonScale - 1. + toy_rms_muonScale));
          bdt_syst_muonScaleDown_WZ     ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_muonScale + toy_rms_muonScale));
          bdt_syst_METScaleUp_WZ        ->SetBinContent(nb, 1. + 0.5 * (toy_mean_METScale - 1. + toy_rms_METScale));
          bdt_syst_METScaleDown_WZ      ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_METScale + toy_rms_METScale));
        }
        if(histo_ZZ->GetBinContent(nb) > 0 ) {
          toy_mean_electronScale = bdt_toy_binyields_electronScale_ZZ[nb-1]->GetMean();
          toy_mean_muonScale = bdt_toy_binyields_muonScale_ZZ[nb-1]->GetMean();
          toy_mean_METScale = bdt_toy_binyields_METScale_ZZ[nb-1]->GetMean();
          toy_rms_electronScale = bdt_toy_binyields_electronScale_ZZ[nb-1]->GetRMS();
          toy_rms_muonScale = bdt_toy_binyields_muonScale_ZZ[nb-1]->GetRMS();
          toy_rms_METScale = bdt_toy_binyields_METScale_ZZ[nb-1]->GetRMS();
          bdt_syst_electronScaleUp_ZZ   ->SetBinContent(nb, 1. + 0.5 * (toy_mean_electronScale - 1. + toy_rms_electronScale));
          bdt_syst_electronScaleDown_ZZ ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_electronScale + toy_rms_electronScale));
          bdt_syst_muonScaleUp_ZZ       ->SetBinContent(nb, 1. + 0.5 * (toy_mean_muonScale - 1. + toy_rms_muonScale));
          bdt_syst_muonScaleDown_ZZ     ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_muonScale + toy_rms_muonScale));
          bdt_syst_METScaleUp_ZZ        ->SetBinContent(nb, 1. + 0.5 * (toy_mean_METScale - 1. + toy_rms_METScale));
          bdt_syst_METScaleDown_ZZ      ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_METScale + toy_rms_METScale));
        }
        if(histo_ggZH_hinv->GetBinContent(nb) > 0 ) {
          toy_mean_electronScale = bdt_toy_binyields_electronScale_ggZH_hinv[nb-1]->GetMean();
          toy_mean_muonScale = bdt_toy_binyields_muonScale_ggZH_hinv[nb-1]->GetMean();
          toy_mean_METScale = bdt_toy_binyields_METScale_ggZH_hinv[nb-1]->GetMean();
          toy_rms_electronScale = bdt_toy_binyields_electronScale_ggZH_hinv[nb-1]->GetRMS();
          toy_rms_muonScale = bdt_toy_binyields_muonScale_ggZH_hinv[nb-1]->GetRMS();
          toy_rms_METScale = bdt_toy_binyields_METScale_ggZH_hinv[nb-1]->GetRMS();
          bdt_syst_electronScaleUp_ggZH_hinv   ->SetBinContent(nb, 1. + 0.5 * (toy_mean_electronScale - 1. + toy_rms_electronScale));
          bdt_syst_electronScaleDown_ggZH_hinv ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_electronScale + toy_rms_electronScale));
          bdt_syst_muonScaleUp_ggZH_hinv       ->SetBinContent(nb, 1. + 0.5 * (toy_mean_muonScale - 1. + toy_rms_muonScale));
          bdt_syst_muonScaleDown_ggZH_hinv     ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_muonScale + toy_rms_muonScale));
          bdt_syst_METScaleUp_ggZH_hinv        ->SetBinContent(nb, 1. + 0.5 * (toy_mean_METScale - 1. + toy_rms_METScale));
          bdt_syst_METScaleDown_ggZH_hinv      ->SetBinContent(nb, 1. - 0.5 * (1. - toy_mean_METScale + toy_rms_METScale));
        }
      }
    }

    histo_Data             ->Write();
    histo_ZH_hinv[nModel]  ->Write();
    histo_Zjets            ->Write();
    histo_VVV              ->Write();
    histo_WZ               ->Write();
    histo_ZZ               ->Write();
    histo_EM               ->Write();
    histo_ggZH_hinv        ->Write();
    
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]	        ->Write();
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]             ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingUp	                        ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingDown                        ->Write();
    histo_WZ_CMS_MVAWZStatBoundingUp	                        ->Write();
    histo_WZ_CMS_MVAWZStatBoundingDown	                        ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingUp	                        ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingDown	                        ->Write();
    histo_EM_CMS_MVAEMStatBoundingUp	                        ->Write();
    histo_EM_CMS_MVAEMStatBoundingDown	                        ->Write();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp                   ->Write();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown                 ->Write();
    
    histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]              ->Write();
    histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]            ->Write();
    histo_VVV_CMS_MVALepEffMBoundingUp	                        ->Write();
    histo_VVV_CMS_MVALepEffMBoundingDown	                ->Write();
    histo_WZ_CMS_MVALepEffMBoundingUp	                        ->Write();
    histo_WZ_CMS_MVALepEffMBoundingDown	                        ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingUp	                        ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingDown	                        ->Write();
    histo_ggZH_hinv_CMS_MVALepEffMBoundingUp                    ->Write();
    histo_ggZH_hinv_CMS_MVALepEffMBoundingDown                  ->Write();
    
    histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]              ->Write();
    histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]            ->Write();
    histo_VVV_CMS_MVALepEffEBoundingUp	                        ->Write();
    histo_VVV_CMS_MVALepEffEBoundingDown	                ->Write();
    histo_WZ_CMS_MVALepEffEBoundingUp	                        ->Write();
    histo_WZ_CMS_MVALepEffEBoundingDown	                        ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingUp	                        ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingDown	                        ->Write();
    histo_ggZH_hinv_CMS_MVALepEffEBoundingUp                    ->Write();
    histo_ggZH_hinv_CMS_MVALepEffEBoundingDown                  ->Write();
    
    histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]                  ->Write();
    histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]                ->Write();
    histo_VVV_CMS_MVAMETBoundingUp	                        ->Write();
    histo_VVV_CMS_MVAMETBoundingDown	                        ->Write();
    histo_WZ_CMS_MVAMETBoundingUp	                        ->Write();
    histo_WZ_CMS_MVAMETBoundingDown	                        ->Write();
    histo_ZZ_CMS_MVAMETBoundingUp	                        ->Write();
    histo_ZZ_CMS_MVAMETBoundingDown	                        ->Write();
    histo_ggZH_hinv_CMS_MVAMETBoundingUp                        ->Write();
    histo_ggZH_hinv_CMS_MVAMETBoundingDown                      ->Write();
    
    histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]	                ->Write();
    histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]                ->Write(); 
    histo_VVV_CMS_MVAJESBoundingUp	                        ->Write();
    histo_VVV_CMS_MVAJESBoundingDown	                        ->Write();
    histo_WZ_CMS_MVAJESBoundingUp 	                        ->Write();
    histo_WZ_CMS_MVAJESBoundingDown	                        ->Write();
    histo_ZZ_CMS_MVAJESBoundingUp 	                        ->Write();
    histo_ZZ_CMS_MVAJESBoundingDown	                        ->Write();
    histo_ggZH_hinv_CMS_MVAJESBoundingUp	                ->Write();
    histo_ggZH_hinv_CMS_MVAJESBoundingDown                      ->Write();
    
    histo_ZH_hinv_CMS_PUBoundingUp[nModel]                      ->Write();
    histo_ZH_hinv_CMS_PUBoundingDown[nModel]                    ->Write();
    histo_VVV_CMS_PUBoundingUp	                                ->Write();
    histo_VVV_CMS_PUBoundingDown	                        ->Write();
    histo_WZ_CMS_PUBoundingUp	                                ->Write();
    histo_WZ_CMS_PUBoundingDown	                                ->Write();
    histo_ZZ_CMS_PUBoundingUp	                                ->Write();
    histo_ZZ_CMS_PUBoundingDown	                                ->Write();
    histo_ggZH_hinv_CMS_PUBoundingUp                            ->Write();
    histo_ggZH_hinv_CMS_PUBoundingDown                          ->Write();
    
    histo_WZ_CMS_EWKCorrUp	                                ->Write();
    histo_WZ_CMS_EWKCorrDown	                                ->Write();
    histo_ZZ_CMS_EWKCorrUp	                                ->Write();
    histo_ZZ_CMS_EWKCorrDown	                                ->Write();
    histo_ZZ_CMS_ggCorrUp	                                ->Write();
    histo_ZZ_CMS_ggCorrDown	                                ->Write();
    
    histo_Zjets_CMS_ZjetsSystUp	                                ->Write();
    histo_Zjets_CMS_ZjetsSystDown	                        ->Write();
    if(useBDT) {
      bdt_syst_electronScaleUp_VVV       ->Write();
      bdt_syst_electronScaleUp_WZ        ->Write();
      bdt_syst_electronScaleUp_ZZ        ->Write();
      bdt_syst_electronScaleUp_ggZH_hinv ->Write();
      bdt_syst_electronScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_muonScaleUp_VVV           ->Write();
      bdt_syst_muonScaleUp_WZ            ->Write();
      bdt_syst_muonScaleUp_ZZ            ->Write();
      bdt_syst_muonScaleUp_ggZH_hinv     ->Write();
      bdt_syst_muonScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_METScaleUp_VVV            ->Write();
      bdt_syst_METScaleUp_WZ             ->Write();
      bdt_syst_METScaleUp_ZZ             ->Write();
      bdt_syst_METScaleUp_ggZH_hinv      ->Write();
      bdt_syst_METScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_electronScaleDown_VVV       ->Write();
      bdt_syst_electronScaleDown_WZ        ->Write();
      bdt_syst_electronScaleDown_ZZ        ->Write();
      bdt_syst_electronScaleDown_ggZH_hinv ->Write();
      bdt_syst_electronScaleDown_ZH_hinv[nModel] ->Write();
      bdt_syst_muonScaleDown_VVV           ->Write();
      bdt_syst_muonScaleDown_WZ            ->Write();
      bdt_syst_muonScaleDown_ZZ            ->Write();
      bdt_syst_muonScaleDown_ggZH_hinv     ->Write();
      bdt_syst_muonScaleDown_ZH_hinv[nModel] ->Write();
      bdt_syst_METScaleDown_VVV            ->Write();
      bdt_syst_METScaleDown_WZ             ->Write();
      bdt_syst_METScaleDown_ZZ             ->Write();
      bdt_syst_METScaleDown_ggZH_hinv      ->Write();
      bdt_syst_METScaleDown_ZH_hinv[nModel] ->Write();   
    }
    outFileLimits->Close();
    if(useBDT && !useCachedBDTSystematics) { // Record the BDT shape systs in a separate file for the caching
      cached_BDT_systematics->cd();
      histo_bdt_toys_electronScale_VVV       ->Write();
      histo_bdt_toys_electronScale_WZ        ->Write();
      histo_bdt_toys_electronScale_ZZ        ->Write();
      histo_bdt_toys_electronScale_ggZH_hinv ->Write();
      histo_bdt_toys_electronScale_ZH_hinv[nModel]->Write();
      histo_bdt_toys_muonScale_VVV           ->Write();
      histo_bdt_toys_muonScale_WZ            ->Write();
      histo_bdt_toys_muonScale_ZZ            ->Write();
      histo_bdt_toys_muonScale_ggZH_hinv     ->Write();
      histo_bdt_toys_muonScale_ZH_hinv[nModel]    ->Write();
      histo_bdt_toys_METScale_VVV            ->Write();
      histo_bdt_toys_METScale_WZ             ->Write();
      histo_bdt_toys_METScale_ZZ             ->Write();
      histo_bdt_toys_METScale_ggZH_hinv      ->Write();
      histo_bdt_toys_METScale_ZH_hinv[nModel]     ->Write();
      bdt_toy_envelope_electronScale_VVV->Write();
      bdt_toy_envelope_electronScale_WZ->Write();
      bdt_toy_envelope_electronScale_ZZ->Write();
      bdt_toy_envelope_electronScale_ZH_hinv[nModel]->Write();
      bdt_toy_envelope_electronScale_ggZH_hinv->Write();
      bdt_toy_envelope_muonScale_VVV->Write();
      bdt_toy_envelope_muonScale_WZ->Write();
      bdt_toy_envelope_muonScale_ZZ->Write();
      bdt_toy_envelope_muonScale_ZH_hinv[nModel]->Write();
      bdt_toy_envelope_muonScale_ggZH_hinv->Write();
      bdt_toy_envelope_METScale_VVV->Write();
      bdt_toy_envelope_METScale_WZ->Write();
      bdt_toy_envelope_METScale_ZZ->Write();
      bdt_toy_envelope_METScale_ZH_hinv[nModel]->Write();
      bdt_toy_envelope_METScale_ggZH_hinv->Write();
      bdt_syst_electronScaleUp_VVV       ->Write();
      bdt_syst_electronScaleUp_WZ        ->Write();
      bdt_syst_electronScaleUp_ZZ        ->Write();
      bdt_syst_electronScaleUp_ggZH_hinv ->Write();
      bdt_syst_electronScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_muonScaleUp_VVV           ->Write();
      bdt_syst_muonScaleUp_WZ            ->Write();
      bdt_syst_muonScaleUp_ZZ            ->Write();
      bdt_syst_muonScaleUp_ggZH_hinv     ->Write();
      bdt_syst_muonScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_METScaleUp_VVV            ->Write();
      bdt_syst_METScaleUp_WZ             ->Write();
      bdt_syst_METScaleUp_ZZ             ->Write();
      bdt_syst_METScaleUp_ggZH_hinv      ->Write();
      bdt_syst_METScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_electronScaleDown_VVV       ->Write();
      bdt_syst_electronScaleDown_WZ        ->Write();
      bdt_syst_electronScaleDown_ZZ        ->Write();
      bdt_syst_electronScaleDown_ggZH_hinv ->Write();
      bdt_syst_electronScaleDown_ZH_hinv[nModel] ->Write();
      bdt_syst_muonScaleDown_VVV           ->Write();
      bdt_syst_muonScaleDown_WZ            ->Write();
      bdt_syst_muonScaleDown_ZZ            ->Write();
      bdt_syst_muonScaleDown_ggZH_hinv     ->Write();
      bdt_syst_muonScaleDown_ZH_hinv[nModel] ->Write();
      bdt_syst_METScaleDown_VVV            ->Write();
      bdt_syst_METScaleDown_WZ             ->Write();
      bdt_syst_METScaleDown_ZZ             ->Write();
      bdt_syst_METScaleDown_ggZH_hinv      ->Write();
      bdt_syst_METScaleDown_ZH_hinv[nModel] ->Write();
      cached_BDT_systematics->Close();
    }
    double lumiE = 1.030;
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
      for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_ZH_hinv_CMS_PDFBounding[nModel][npdf]->GetBinContent(nb)-histo_ZH_hinv[nModel]->GetBinContent(nb))/histo_ZH_hinv[nModel]->GetBinContent(nb));
      systPDF[0] = 1.0+sqrt(TMath::Max(histo_Diff->GetRMS()*histo_Diff->GetRMS()-pdfTotal[0]*pdfTotal[0],0.0));
      histo_Diff->Reset();
      for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
      systPDF[1] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_WZ_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WZ->GetBinContent(nb))/histo_WZ->GetBinContent(nb));
      systPDF[2] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ZZ->GetBinContent(nb))/histo_ZZ->GetBinContent(nb));
      systPDF[3] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=1; npdf<102; npdf++) {
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

      // BDT systematic study
      double systBDTMuonUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTMuonDown[5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTElectronUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTElectronDown[5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTMETUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTMETDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(useBDT) { 
        systBDTElectronUp  [0] = bdt_syst_electronScaleUp_ZH_hinv[nModel]   ->GetBinContent(nb);
        systBDTElectronDown[0] = bdt_syst_electronScaleDown_ZH_hinv[nModel] ->GetBinContent(nb);
        systBDTMuonUp  [0]     = bdt_syst_muonScaleUp_ZH_hinv[nModel]       ->GetBinContent(nb);
        systBDTMuonDown[0]     = bdt_syst_muonScaleDown_ZH_hinv[nModel]     ->GetBinContent(nb);
        systBDTMETUp  [0]      = bdt_syst_METScaleUp_ZH_hinv[nModel]        ->GetBinContent(nb);
        systBDTMETDown[0]      = bdt_syst_METScaleDown_ZH_hinv[nModel]      ->GetBinContent(nb);
        systBDTElectronUp  [1] = bdt_syst_electronScaleUp_VVV   ->GetBinContent(nb);
        systBDTElectronDown[1] = bdt_syst_electronScaleDown_VVV ->GetBinContent(nb);
        systBDTMuonUp  [1]     = bdt_syst_muonScaleUp_VVV       ->GetBinContent(nb);
        systBDTMuonDown[1]     = bdt_syst_muonScaleDown_VVV     ->GetBinContent(nb);
        systBDTMETUp  [1]      = bdt_syst_METScaleUp_VVV        ->GetBinContent(nb);
        systBDTMETDown[1]      = bdt_syst_METScaleDown_VVV      ->GetBinContent(nb);
        systBDTElectronUp  [2] = bdt_syst_electronScaleUp_WZ   ->GetBinContent(nb);
        systBDTElectronDown[2] = bdt_syst_electronScaleDown_WZ ->GetBinContent(nb);
        systBDTMuonUp  [2]     = bdt_syst_muonScaleUp_WZ       ->GetBinContent(nb);
        systBDTMuonDown[2]     = bdt_syst_muonScaleDown_WZ     ->GetBinContent(nb);
        systBDTMETUp  [2]      = bdt_syst_METScaleUp_WZ        ->GetBinContent(nb);
        systBDTMETDown[2]      = bdt_syst_METScaleDown_WZ      ->GetBinContent(nb);
        systBDTElectronUp  [3] = bdt_syst_electronScaleUp_ZZ   ->GetBinContent(nb);
        systBDTElectronDown[3] = bdt_syst_electronScaleDown_ZZ ->GetBinContent(nb);
        systBDTMuonUp  [3]     = bdt_syst_muonScaleUp_ZZ       ->GetBinContent(nb);
        systBDTMuonDown[3]     = bdt_syst_muonScaleDown_ZZ     ->GetBinContent(nb);
        systBDTMETUp  [3]      = bdt_syst_METScaleUp_ZZ        ->GetBinContent(nb);
        systBDTMETDown[3]      = bdt_syst_METScaleDown_ZZ      ->GetBinContent(nb);
        systBDTElectronUp  [4] = bdt_syst_electronScaleUp_ggZH_hinv   ->GetBinContent(nb);
        systBDTElectronDown[4] = bdt_syst_electronScaleDown_ggZH_hinv ->GetBinContent(nb);
        systBDTMuonUp  [4]     = bdt_syst_muonScaleUp_ggZH_hinv       ->GetBinContent(nb);
        systBDTMuonDown[4]     = bdt_syst_muonScaleDown_ggZH_hinv     ->GetBinContent(nb);
        systBDTMETUp  [4]      = bdt_syst_METScaleUp_ggZH_hinv        ->GetBinContent(nb);
        systBDTMETDown[4]      = bdt_syst_METScaleDown_ggZH_hinv      ->GetBinContent(nb);
        
      }
   
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
  
      double syst_EWKCorrUp[3]   = {1.0,1.0,1.0}; // WZ, ZZ, ggZZ
      double syst_EWKCorrDown[3] = {1.0,1.0,1.0};
      if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [0] = histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[0] = histo_WZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [1] = histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[1] = histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb) > 0) syst_EWKCorrUp  [2] = histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb) > 0) syst_EWKCorrDown[2] = histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);

      // If VV normalization is from data
      double systVV[2] = {0,0};
      if(useVVFromData){
        systVV[0] = sqrt(TMath::Power(syst_EWKCorrUp[0]-1,2)+TMath::Power(systQCDScale[2]-1,2)+TMath::Power(systPDF[2]-1,2));
	systVV[1] = sqrt(TMath::Power(syst_EWKCorrUp[1]-1,2)+TMath::Power(systQCDScale[3]-1,2)+TMath::Power(systPDF[3]-1,2));
	//printf("systVV_theo = %f(%f/%f/%f) %f(%f/%f/%f)\n",systVV[0],syst_EWKCorrUp[0]-1,systQCDScale[2]-1,systPDF[2]-1,
	//                                                   systVV[1],syst_EWKCorrUp[1]-1,systQCDScale[3]-1,systPDF[3]-1);
      }

      char outputLimitsShape[200];                                            
      sprintf(outputLimitsShape,"MitZHAnalysis/datacards%s/histo_limits_zll%shinv%s_%s_shape_%s_bin%d.txt",subdirectory.c_str(), addChan.Data(), finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb-1);
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
      newcardShape << Form("bin hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process ZH_hinv Zjets VVV WZ ZZ EM ggZH_hinv\n");
      newcardShape << Form("process 0 1 2 3 4 5 -1\n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",histo_ZH_hinv[nModel]->GetBinContent(nb),histo_Zjets->GetBinContent(nb),TMath::Max(histo_VVV->GetBinContent(nb),0.0),histo_WZ->GetBinContent(nb),histo_ZZ->GetBinContent(nb),histo_EM->GetBinContent(nb),nggZHEvt);
      newcardShape << Form("lumi_%4s                               lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE);		     
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3],systLepEffM[4]);
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3],systLepEffE[4]);
      if(MVAVarType != 3) {
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3],systLepResM[4]);
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3],systLepResE[4]);
      }
      newcardShape << Form("CMS_pu2016                             lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systPUUp[0],systPUDown[0],systPUUp[1],systPUDown[1],systPUUp[2],systPUDown[2],systPUUp[3],systPUDown[3],systPUUp[0],systPUDown[0]); // 0 --> 4
      newcardShape << Form("CMS_scale_met                          lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systMetUp[0],systMetDown[0],systMetUp[1],systMetDown[1],systMetUp[2],systMetDown[2],systMetUp[3],systMetDown[3],systMetUp[0],systMetDown[0]); // 0 --> 4
      newcardShape << Form("CMS_scale_j                            lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systJesUp[0],systJesDown[0],systJesUp[1],systJesDown[1],systJesUp[2],systJesDown[2],systJesUp[3],systJesDown[3],systJesUp[0],systJesDown[0]); // 0 --> 4		 
      newcardShape << Form("CMS_trigger2016                        lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",1.02,1.02,1.02,1.02,1.02);
      newcardShape << Form("UEPS			           lnN  1.030   -     -     -     -     -   1.030\n");
      newcardShape << Form("CMS_eff_b_mistag2016                   lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",syst_btag,syst_btag,syst_btag,syst_btag,syst_btag);
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
      newcardShape << Form("QCDscale_VV		                   lnN    -     -     -   %7.5f       %7.5f    -      -  \n",systQCDScale[2],systQCDScale[3]);		

      if(useVVFromData && nb != 1){
      if(useZZWZEWKUnc){
      newcardShape << Form("CMS_hinv_vvnorm_bin%d rateParam  * WZ 1 [0.1,10]\n",nb-1);		
      newcardShape << Form("CMS_hinv_vvnorm_bin%d rateParam  * ZZ 1 [0.1,10]\n",nb-1);	
      newcardShape << Form("CMS_zllhinv_ZZWZ_EWKCorr               lnN    -     -     -   -      %7.5f    -      -\n",1.+sqrt(0.1*0.1+(syst_EWKCorrUp[1]-1.0)*(syst_EWKCorrUp[1]-1.0)));		
      } else {
      newcardShape << Form("CMS_hinv_wznorm_bin%d rateParam  * WZ 1 [0.1,10]\n",nb-1);		
      newcardShape << Form("CMS_hinv_wznorm_bin%d param 1 %5.3f\n",nb-1, systVV[0]);		
      newcardShape << Form("CMS_hinv_zznorm_bin%d rateParam  * ZZ 1 [0.1,10]\n",nb-1);	
      newcardShape << Form("CMS_hinv_zznorm_bin%d param 1 %5.3f\n",nb-1, systVV[1]);	
      }
      }
      else if(!useVVFromData){
      newcardShape << Form("CMS_zllhinv_WZ_EWKCorr                 lnN    -     -     -   %7.5f/%7.5f   -      -      -  \n",syst_EWKCorrUp[0],syst_EWKCorrDown[0]);		
      newcardShape << Form("CMS_zllhinv_ZZ_EWKCorr                 lnN    -     -     -     -    %7.5f/%7.5f   -      -  \n",syst_EWKCorrUp[1],syst_EWKCorrDown[1]);		
      }

      newcardShape << Form("CMS_zllhinv_ggZZCorr                   lnN    -     -     -     -   %7.5f/%7.5f   -     -  \n",syst_EWKCorrUp[2],syst_EWKCorrDown[2]);		

      newcardShape << Form("CMS_zllhinv_WZ_lep2016                 lnN     -	 -     -   %7.5f   -	  -    -  \n",syst_WZl[0]);	    
      newcardShape << Form("CMS_zllhinv_WZ_tau2016                 lnN     -	 -     -   %7.5f   -	  -    -  \n",syst_WZl[1]);	    
      if(MVAVarType == 3 || MVAVarType==4) {
      newcardShape << Form("CMS_BDT_scale_m                   lnN	%7.5f/%7.5f -  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f -  %7.5f/%7.5f   \n", systBDTMuonUp[0], systBDTMuonDown[0], systBDTMuonUp[1], systBDTMuonDown[1], systBDTMuonUp[2], systBDTMuonDown[2], systBDTMuonUp[3], systBDTMuonDown[3], systBDTMuonUp[4], systBDTMuonDown[4]);	    
      newcardShape << Form("CMS_BDT_scale_e                   lnN	%7.5f/%7.5f -  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f -  %7.5f/%7.5f   \n", systBDTElectronUp[0], systBDTElectronDown[0], systBDTElectronUp[1], systBDTElectronDown[1], systBDTElectronUp[2], systBDTElectronDown[2], systBDTElectronUp[3], systBDTElectronDown[3], systBDTElectronUp[4], systBDTElectronDown[4]);	    
      newcardShape << Form("CMS_BDT_scale_MET                 lnN	%7.5f/%7.5f -  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f -  %7.5f/%7.5f   \n", systBDTMETUp[0], systBDTMETDown[0], systBDTMETUp[1], systBDTMETDown[1], systBDTMETUp[2], systBDTMETDown[2], systBDTMETUp[3], systBDTMETDown[3], systBDTMETUp[4], systBDTMETDown[4]);	    
      }
      if(nb>=3){
      newcardShape << Form("CMS_zllhinv_ZLLNorm2016_%s_%s              lnN	-   %7.5f   -	  -     -     -     -  \n",finalStateName,ECMsb.Data(),2.0);	    
      //newcardShape << Form("CMS_zllhinv_ZLLShape2016_%s_%s             lnN	-   %7.5f/%7.5f   -	-     -     -     -  \n",finalStateName,ECMsb.Data(),systZjetsUp[0],systZjetsDown[0]);	    
      }
      if(MVAVarType == 1 || MVAVarType == 2 || MVAVarType == 4){
      newcardShape << Form("CMS_zllhinv_ZLLFit2016_%s_%s               lnU	-   %7.5f   -	  -     -     -     -  \n",finalStateName,ECMsb.Data(),2.0);	    
      }
      if(nb>=2)
      newcardShape << Form("CMS_zllhinv_EMSyst2016_%s_%s               lnN	-     -     -	  -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),systEM[0]);	       

      newcardShape << Form("CMS_zllhinv_EMNorm2016_%s_%s               lnU	-     -     -	  -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),2.0);      
  
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding2016_%s_Bin%d    lnN    %7.5f -      -	 -    -    -	 -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ZH_hinv[nModel]->GetBinError(nb)  /histo_ZH_hinv[nModel]->GetBinContent(nb)  );
  
      if(histo_Zjets->GetBinContent(nb)           > 0) newcardShape << Form("CMS_zllhinv%s_MVAZjetsStatBounding2016_%s_Bin%d  lnN	-  %7.5f   -	-    -    -	-  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_Zjets->GetBinError(nb)    /histo_Zjets->GetBinContent(nb)    );
  
      if(histo_VVV->GetBinContent(nb)             > 0) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding2016_%s_Bin%d    lnN	-   -  %7.5f   -    -	 -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_VVV->GetBinError(nb)	   /histo_VVV->GetBinContent(nb)      );
  
      if(histo_WZ->GetBinContent(nb)              > 0) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding2016_%s_Bin%d     lnN      -  -     -  %7.5f  -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_WZ->GetBinError(nb)	/histo_WZ->GetBinContent(nb)	   );
  
      if(histo_ZZ->GetBinContent(nb)	          > 0) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding2016_%s_Bin%d     lnN	  -  -     -	-  %7.5f  -	-  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ZZ->GetBinError(nb)	    /histo_ZZ->GetBinContent(nb)       );
  
      if(histo_EM->GetBinContent(nb)              > 0) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding2016_%s_Bin%d     lnN	  -  -     -	-    -  %7.5f	-  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_EM->GetBinError(nb)	    /histo_EM->GetBinContent(nb)       );
  
      if(histo_ggZH_hinv->GetBinContent(nb)       > 0) newcardShape << Form("CMS_zllhinv%s_MVAggZHStatBounding2016_%s_Bin%d   lnN      -    -     -    -    -    -   %7.5f\n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ggZH_hinv->GetBinError(nb)/histo_ggZH_hinv->GetBinContent(nb));
    }

  } // end loop over models
  
}
