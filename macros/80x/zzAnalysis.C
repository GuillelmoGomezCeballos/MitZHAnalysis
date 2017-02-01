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
#include "MitZHAnalysis/macros/80x/zhMVA.h"

bool isMINIAOD = true;
int whichSkim = 2;
bool usePureMC = true; 
double mcPrescale = 1.0;
bool useZZWZEWKUnc = true;
enum selType                     {ZZSEL=0,  ZHWWSEL,   ZHSEL, nSelTypes};
TString selTypeName[nSelTypes]= {"ZZSEL",  "ZHWWSEL", "ZHSEL",};
enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};
const TString typeLepSel = "medium";

void zzAnalysis(
 string the_BDT_weights="",
 string subdirectory=""
 ){
  if(subdirectory!="" && subdirectory.c_str()[0]!='/') subdirectory = "/"+subdirectory;
  system(("mkdir -p MitZHAnalysis/datacards"+subdirectory).c_str());
  system(("mkdir -p MitZHAnalysis/plots"+subdirectory).c_str());
  bool printMCEventList=true;
  bool useBDT=false;

  Int_t period = 1;
  TString filesPathDA = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/met_";
  TString filesPathMC  = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/met_";
  Double_t lumi = 36.8;
  bool verbose = true;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infileName_;  
  vector<Int_t> infileCategory_;  

  TString puPath = "";
  TString triggerSuffix = "*";
  if(isMINIAOD) triggerSuffix = "";

  puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";

  // Data files
  if(isMINIAOD) {
    infileName_.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data())); infileCategory_.push_back(0);
    infileName_.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data())); infileCategory_.push_back(0);
    infileName_.push_back(Form("%sdata_Run2016D.root",filesPathDA.Data())); infileCategory_.push_back(0);
    infileName_.push_back(Form("%sdata_Run2016E.root",filesPathDA.Data())); infileCategory_.push_back(0);
    infileName_.push_back(Form("%sdata_Run2016F.root",filesPathDA.Data())); infileCategory_.push_back(0);
    infileName_.push_back(Form("%sdata_Run2016G.root",filesPathDA.Data())); infileCategory_.push_back(0);
    infileName_.push_back(Form("%sdata_Run2016H.root",filesPathDA.Data())); infileCategory_.push_back(0);
  } else {
  }

  if(usePureMC == true){
  infileName_.push_back(Form("%sWWTo2L2Nu_13TeV-powheg.root",filesPathMC.Data()));                                            infileCategory_.push_back(1);
  infileName_.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV.root",filesPathMC.Data()));					      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));        infileCategory_.push_back(1);
  infileName_.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));            infileCategory_.push_back(1);
  infileName_.push_back(Form("%sTTTo2L2Nu_13TeV-powheg.root",filesPathMC.Data()));                                            infileCategory_.push_back(1);
  infileName_.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));    infileCategory_.push_back(1);
  infileName_.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));infileCategory_.push_back(1);
  infileName_.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));		      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                  infileCategory_.push_back(1);
  infileName_.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root",filesPathMC.Data()));               infileCategory_.push_back(1);
  infileName_.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8.root",filesPathMC.Data()));              infileCategory_.push_back(1);
  infileName_.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));		              infileCategory_.push_back(1);
  infileName_.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                            infileCategory_.push_back(1);
  infileName_.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8.root",filesPathMC.Data()));  				      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                   infileCategory_.push_back(1);
  infileName_.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));			      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infileCategory_.push_back(1);
  infileName_.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));           infileCategory_.push_back(1);
  infileName_.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 		                      infileCategory_.push_back(1);
  }

  infileName_.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8.root",filesPathMC.Data()));					      infileCategory_.push_back(4);
  infileName_.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infileCategory_.push_back(4);
  infileName_.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	 	      infileCategory_.push_back(4);
  infileName_.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	 	      infileCategory_.push_back(4);
  infileName_.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infileCategory_.push_back(4);
  infileName_.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infileCategory_.push_back(4);
  infileName_.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infileCategory_.push_back(4);

  infileName_.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infileCategory_.push_back(5);
  infileName_.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infileCategory_.push_back(5);
  infileName_.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infileCategory_.push_back(5);
  infileName_.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));	      infileCategory_.push_back(5);

  infileName_.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data())); 		      infileCategory_.push_back(6);
  infileName_.push_back(Form("%sttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.root",filesPathMC.Data()));     infileCategory_.push_back(6);
  infileName_.push_back(Form("%sGluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root",filesPathMC.Data()));               infileCategory_.push_back(6);
  infileName_.push_back(Form("%sVBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root",filesPathMC.Data()));                 infileCategory_.push_back(6);

  //infileName_.clear();infileCategory_.clear();
  //infileName_.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infileCategory_.push_back(1);
  //infileName_.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	 infileCategory_.push_back(2);
  //infileName_.push_back(Form("%sTT_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));     infileCategory_.push_back(3);
  //infileName_.push_back(Form("%sTTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	 infileCategory_.push_back(4);

  if(infileName_.size() != infileCategory_.size()) {assert(0); return;}

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPU);    fhDPU    ->SetDirectory(0);
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

  TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/trackMuReco_SF.root"));
  TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptg10")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
  //TH1D *fhDmutrksfptl10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptl10")); assert(fhDmutrksfptl10); fhDmutrksfptl10->SetDirectory(0);
  delete fTrackMuonReco_SF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("scalefactors_Id_Muon")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  //TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonID_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  delete fMuSF;

  TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("scalefactors_Iso_Muon")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  //TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonIso_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  delete fMuIsoSF;

  //const int MVAVarType = 1; const int nBinMVA = 8; Float_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 350}; TString addChan = "1";
  const int MVAVarType = 1; const int nBinMVA = 12; Float_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600}; TString addChan = "1";
  //const int MVAVarType = 1; const int nBinMVA = 11; Float_t xbins[nBinMVA+1] = {0, 50, 100, 160, 240, 320, 400, 480, 560, 640, 800, 1200}; TString addChan = "1";
  //const int MVAVarType = 1; const int nBinMVA = 21; Float_t xbins[nBinMVA+1] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1200}; TString addChan = "1";
  //const int MVAVarType = 3; const int nBinMVA = 18; Float_t xbins[nBinMVA+1] =  {-2, -1, 0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4}; TString addChan = "3";
  //const int MVAVarType = 4; const int nBinMVA = 26; Float_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 350,
  if (MVAVarType==3 || MVAVarType==4) useBDT=true;
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

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 35;
  const int histBins = 7;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "Zgamma", "....WZ", "...ZZ", "...VVV", ".Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >=  1 && thePlot <=  1) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=  3 && thePlot <=  4) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   4.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  8 && thePlot <=  9) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >= 11 && thePlot <= 12) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >= 13 && thePlot <= 16) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 17 && thePlot <= 17) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >= 18 && thePlot <= 18) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 20 && thePlot <= 21) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 22 && thePlot <= 24) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 25 && thePlot <= 25) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 26 && thePlot <= 26) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot =   2.0;}
    else if(thePlot >= 29 && thePlot <= 29) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 30 && thePlot <= 30) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 31 && thePlot <= 33) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    TH1D* histos;
    if(thePlot != allPlots-1 && thePlot != 27 && thePlot != 28) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    else                                                        histos = new TH1D("histos", "histos", nBinMVA, xbins);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Clear();
  }

  TH1D *histo_Data     = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_Fake     = (TH1D*) histoMVA->Clone("histo_Fake");	 
  TH1D *histo_ZZ       = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_VVV      = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_Higgs    = (TH1D*) histoMVA->Clone("histo_Higgs");	 

  char finalStateName[6],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff2016_m");sprintf(momMName,"CMS_scale2016_m");
  sprintf(effEName,"CMS_eff2016_e");sprintf(momEName,"CMS_scale2016_e");
  sprintf(finalStateName,"llll1j");

  const int allStates = 4;
  double bgdDecay[nSelTypes*allStates][histBins],weiDecay[nSelTypes*allStates][histBins];
  for(unsigned int i=0; i<nSelTypes*allStates; i++) {
    for(int j=0; j<histBins; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }

  TString ECMsb  = "13TeV2016";

  TH1D* histo_Fake_CMS_MVAFakeStatBoundingUp       = new TH1D( Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Fake_CMS_MVAFakeStatBoundingUp  ->Sumw2();
  TH1D* histo_Fake_CMS_MVAFakeStatBoundingDown     = new TH1D( Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Fake_CMS_MVAFakeStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingUp     = new TH1D( Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingDown   = new TH1D( Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingDown->Sumw2();

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBounding[6];
  for(int nb=0; nb<6; nb++){
    histo_ZZ_CMS_QCDScaleBounding[nb]	     = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]       = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_Higgs_CMS_QCDScaleBounding[nb]     = new TH1D(Form("histo_Higgs_QCDScale_f%d",nb),      Form("histo_Higgs_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Higgs_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_ZZ_CMS_PDFBounding[102];
  TH1D* histo_VVV_CMS_PDFBounding[102];
  TH1D* histo_Higgs_CMS_PDFBounding[102];
  for(int nb=0; nb<102; nb++) {
    histo_ZZ_CMS_PDFBounding[nb]        = new TH1D(Form("histo_ZZ_PDF_f%d",nb),      Form("histo_ZZ_PDF_f%d",nb),     nBinMVA, xbins); histo_ZZ_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]       = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_Higgs_CMS_PDFBounding[nb]     = new TH1D(Form("histo_Higgs_PDF_f%d",nb),      Form("histo_Higgs_PDF_f%d",nb),     nBinMVA, xbins); histo_Higgs_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_Fake_CMS_MVAFakeStatBoundingBinUp[nBinMVA];
  TH1D* histo_Fake_CMS_MVAFakeStatBoundingBinDown[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++) {
    histo_Fake_CMS_MVAFakeStatBoundingBinUp[nb]    = new TH1D(Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%s_Bin%dUp"	      ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_Fake_CMS_MVAFakeStatBoundingBinDown[nb]  = new TH1D(Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_Fake_CMS_zz4l%s_MVAFakeStatBounding_%s_Bin%dDown"	    ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]        = new TH1D(Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      = new TH1D(Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_zz4l%s_MVAZZStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      = new TH1D(Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%s_Bin%dUp" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%s_Bin%dUp"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    = new TH1D(Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%s_Bin%dDown"	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_zz4l%s_MVAVVVStatBounding_%s_Bin%dDown"	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]  = new TH1D(Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%s_Bin%dUp"	       ,finalStateName,  ECMsb.Data(),nb), Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%s_Bin%dUp"	    ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]= new TH1D(Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%s_Bin%dDown"	     ,finalStateName,  ECMsb.Data(),nb), Form("histo_Higgs_CMS_zz4l%s_MVAHiggsStatBounding_%s_Bin%dDown"	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_Fake_CMS_MVAFakeStatBoundingBinUp[nb]    ->Sumw2();
    histo_Fake_CMS_MVAFakeStatBoundingBinDown[nb]  ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	   ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	   ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	   ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]  ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]->Sumw2();
  }

  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp    	   = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown  	   = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	   = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	   = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingUp       = new TH1D( Form("histo_Higgs_%sUp",effMName)  , Form("histo_Higgs_%sUp",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingDown     = new TH1D( Form("histo_Higgs_%sDown",effMName), Form("histo_Higgs_%sDown",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    	   = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	   = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingAvg      = new TH1D( Form("histo_Higgs_%sAvg",effMName)  , Form("histo_Higgs_%sAvg",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp    	   = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown  	   = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	   = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	   = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingUp       = new TH1D( Form("histo_Higgs_%sUp",effEName)  , Form("histo_Higgs_%sUp",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingDown     = new TH1D( Form("histo_Higgs_%sDown",effEName), Form("histo_Higgs_%sDown",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingDown->Sumw2();

  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    	   = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	   = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingAvg      = new TH1D( Form("histo_Higgs_%sAvg",effEName)  , Form("histo_Higgs_%sAvg",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_ZZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_Higgs_CMS_scale_metUp")  , Form("histo_Higgs_CMS_scale_metUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_Higgs_CMS_scale_metDown"), Form("histo_Higgs_CMS_scale_metDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_ZZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_ZZ_CMS_scale_jUp")  , Form("histo_ZZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_ZZ_CMS_scale_jDown"), Form("histo_ZZ_CMS_scale_jDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp      	= new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown    	= new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_Higgs_CMS_scale_jUp")  , Form("histo_Higgs_CMS_scale_jUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_Higgs_CMS_scale_jDown"), Form("histo_Higgs_CMS_scale_jDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_muonUp")  , Form("histo_ZZ_CMS_bdt_muonUp")  , nBinMVA, xbins);  histo_ZZ_CMS_BDTMuonScaleBoundingUp   ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_muonDown"), Form("histo_ZZ_CMS_bdt_muonDown"), nBinMVA, xbins);  histo_ZZ_CMS_BDTMuonScaleBoundingDown ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMuonScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_muonUp")  , Form("histo_VVV_CMS_bdt_muonUp")  , nBinMVA, xbins);histo_VVV_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMuonScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_muonDown"), Form("histo_VVV_CMS_bdt_muonDown"), nBinMVA, xbins);histo_VVV_CMS_BDTMuonScaleBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMuonScaleBoundingUp    	= new TH1D( Form("histo_Higgs_CMS_bdt_muonUp")  , Form("histo_Higgs_CMS_bdt_muonUp")  , nBinMVA, xbins);  histo_Higgs_CMS_BDTMuonScaleBoundingUp   ->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMuonScaleBoundingDown  	= new TH1D( Form("histo_Higgs_CMS_bdt_muonDown"), Form("histo_Higgs_CMS_bdt_muonDown"), nBinMVA, xbins);  histo_Higgs_CMS_BDTMuonScaleBoundingDown ->Sumw2();

  TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_electronUp")  , Form("histo_ZZ_CMS_bdt_electronUp")  , nBinMVA, xbins);  histo_ZZ_CMS_BDTElectronScaleBoundingUp   ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_electronDown"), Form("histo_ZZ_CMS_bdt_electronDown"), nBinMVA, xbins);  histo_ZZ_CMS_BDTElectronScaleBoundingDown ->Sumw2();
  TH1D* histo_VVV_CMS_BDTElectronScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_electronUp")  , Form("histo_VVV_CMS_bdt_electronUp")  , nBinMVA, xbins);histo_VVV_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTElectronScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_electronDown"), Form("histo_VVV_CMS_bdt_electronDown"), nBinMVA, xbins);histo_VVV_CMS_BDTElectronScaleBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_BDTElectronScaleBoundingUp    	= new TH1D( Form("histo_Higgs_CMS_bdt_electronUp")  , Form("histo_Higgs_CMS_bdt_electronUp")  , nBinMVA, xbins);  histo_Higgs_CMS_BDTElectronScaleBoundingUp   ->Sumw2();
  TH1D* histo_Higgs_CMS_BDTElectronScaleBoundingDown  	= new TH1D( Form("histo_Higgs_CMS_bdt_electronDown"), Form("histo_Higgs_CMS_bdt_electronDown"), nBinMVA, xbins);  histo_Higgs_CMS_BDTElectronScaleBoundingDown ->Sumw2();

  TH1D* histo_ZZ_CMS_BDTMETScaleBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_bdt_METUp")  , Form("histo_ZZ_CMS_bdt_METUp")  , nBinMVA, xbins);  histo_ZZ_CMS_BDTMETScaleBoundingUp   ->Sumw2();
  TH1D* histo_ZZ_CMS_BDTMETScaleBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_bdt_METDown"), Form("histo_ZZ_CMS_bdt_METDown"), nBinMVA, xbins);  histo_ZZ_CMS_BDTMETScaleBoundingDown ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMETScaleBoundingUp   	= new TH1D( Form("histo_VVV_CMS_bdt_METUp")  , Form("histo_VVV_CMS_bdt_METUp")  , nBinMVA, xbins);histo_VVV_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_BDTMETScaleBoundingDown 	= new TH1D( Form("histo_VVV_CMS_bdt_METDown"), Form("histo_VVV_CMS_bdt_METDown"), nBinMVA, xbins);histo_VVV_CMS_BDTMETScaleBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMETScaleBoundingUp    	= new TH1D( Form("histo_Higgs_CMS_bdt_METUp")  , Form("histo_Higgs_CMS_bdt_METUp")  , nBinMVA, xbins);  histo_Higgs_CMS_BDTMETScaleBoundingUp   ->Sumw2();
  TH1D* histo_Higgs_CMS_BDTMETScaleBoundingDown  	= new TH1D( Form("histo_Higgs_CMS_bdt_METDown"), Form("histo_Higgs_CMS_bdt_METDown"), nBinMVA, xbins);  histo_Higgs_CMS_BDTMETScaleBoundingDown ->Sumw2();

  TH1D* histo_ZZ_CMS_PUBoundingUp            	= new TH1D( Form("histo_ZZ_CMS_puUp")  , Form("histo_ZZ_CMS_puUp")  , nBinMVA, xbins); histo_ZZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingDown  	        = new TH1D( Form("histo_ZZ_CMS_puDown"), Form("histo_ZZ_CMS_puDown"), nBinMVA, xbins); histo_ZZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingUp           	= new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown         	= new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingUp            = new TH1D( Form("histo_Higgs_CMS_puUp")  , Form("histo_Higgs_CMS_puUp")  , nBinMVA, xbins); histo_Higgs_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingDown          = new TH1D( Form("histo_Higgs_CMS_puDown"), Form("histo_Higgs_CMS_puDown"), nBinMVA, xbins); histo_Higgs_CMS_PUBoundingDown->Sumw2();

  TH1D* histo_ZZ_CMS_EWKCorrUp                  = new TH1D( Form("histo_ZZ_EWKCorrUp")  , Form("histo_ZZ_EWKCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_EWKCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_EWKCorrDown                = new TH1D( Form("histo_ZZ_EWKCorrDown"), Form("histo_ZZ_EWKCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_EWKCorrDown->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrUp                   = new TH1D( Form("histo_ZZ_ggCorrUp")  , Form("histo_ZZ_ggCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_ggCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrDown                 = new TH1D( Form("histo_ZZ_ggCorrDown"), Form("histo_ZZ_ggCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_ggCorrDown->Sumw2();
  
  char outputEventList[200];
  sprintf(outputEventList,"MitZHAnalysis/datacards%s/eventlist_zz_MC.csv",subdirectory.c_str());
  ofstream eventList;
  if(printMCEventList) {
    eventList.open(outputEventList);
    eventList << "eventNumber,njets,realMET,fakeMET,balance,dPhiZMET,ZpT,dPhiJetMET,pdgl1,pdgl2,pdgl3,pdgl4\n";
  }

  unsigned int numberOfLeptons = 4;

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infileName_.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infileName_.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infileName_[ifile].Data());

    TFile *the_input_file = TFile::Open(infileName_[ifile].Data());
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

    //BareTaus eventTaus;
    //eventTaus.setBranchAddresses(the_input_tree);

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
    if(infileCategory_[ifile] == 0){
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
    if(infileCategory_[ifile] != 0 && initPDFTag == -1) {
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
    double theMCPrescale = mcPrescale;
    if(infileCategory_[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;

      the_input_tree->GetEntry(i);

      Bool_t passFilter[14] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,
                               kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      if(infileCategory_[ifile] == 0) {
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
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() <= 10 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[3]])->Pt() <= 10) continue;

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

      passFilter[4] = TMath::Abs(signQ) == 0;
      if(passFilter[4] == kFALSE) continue;

      double minMassll = 999.0;
      double minMassZ[2] = {99999.0, 99999.0};
      double deltaRllMin = 999.0;
      int tagZ[4] = {-1,-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
          double deltaRllAux = ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl1]]));
          if(deltaRllAux < deltaRllMin) deltaRllMin = deltaRllAux;

	  if((int)(*eventLeptons.pdgId)[idLep[nl0]] * (int)(*eventLeptons.pdgId)[idLep[nl1]] > 0) continue;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl1])) ) ));

	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl1]])){
	    if(TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ[0]-91.1876)) {
	      minMassZ[0] = dilepAux.M();tagZ[0]=nl0;tagZ[1]=nl1;
	    }
	  }

	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }
      
      // Removing events with no Z candidate at all
      if(tagZ[0] == -1) continue;

      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(tagZ[0] == (int)nl || tagZ[1] == (int)nl) continue;
	if(tagZ[2] == -1) tagZ[2] = nl;
	else              tagZ[3] = nl;
      }

      if(tagZ[0] == -1 || tagZ[1] == -1 || tagZ[2] == -1 || tagZ[3] == -1) printf("TAGZ PROBLEM!: %d %d %d %d\n",tagZ[0],tagZ[1],tagZ[2],tagZ[3]);

      TLorentzVector fourLep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
        		       ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
        		       ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) + 
        		       ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[3])) ) ));
      double mass4l = fourLep.M();
      int typeL[2] = {0,0};
      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 11) typeL[0]++;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 13) typeL[1]++;
        else {printf("ZZ4lPROBLEM!\n");assert(0);return;}
      }

      int type4l = 0;
      if     (typeL[0] == 4) type4l = 0;
      else if(typeL[1] == 4) type4l = 1;
      else if(typeL[0] == 2) type4l = 2;
      else		     type4l = 3;
      // Determine flavor of the pair (0 means e-mu pair, 1 means mu-mu, 2 means e-e)
      int typePair = 0;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==13) typePair = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==11) typePair = 2;

      TLorentzVector dilepZ (( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) ) ));
      TLorentzVector dilepLL(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[2]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[3]])) ) ));
      minMassZ[1] = dilepLL.M();
      double deltaPhiDileptonMet = TMath::Abs(dilepLL.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtHWW = TMath::Sqrt(2.0*dilepLL.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));

      TLorentzVector dilepWW;
      dilepWW.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px()+dilepLL.Px());
      dilepWW.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py()+dilepLL.Py());
      dilepWW.SetPz(0.0);

      TLorentzVector theFakeMETUp, theFakeMETDown, theZZllnnMET, theOtherZZllnnMET, dilepZll, dilepZnn;
      double dPhiDiLepMET = 0, ptFrac = 0;
      { // MET emulation l -> nu
        TLorentzVector theFakeMET[2];
        theFakeMET[0].SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Px());
        theFakeMET[0].SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Py());
        theFakeMET[1].SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
        theFakeMET[1].SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
        if(TMath::Abs(minMassZ[0]-91.1876) < TMath::Abs(minMassZ[1]-91.1876)) { // (1) => ll, (2) => nn
	  theZZllnnMET = theFakeMET[1];
          dilepZll = ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[0]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[1]])) );
          dilepZnn = ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[2]])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[tagZ[3]])) );
          theFakeMETUp  .SetPx(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
          theFakeMETUp  .SetPy(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
          theFakeMETDown.SetPx(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
          theFakeMETDown.SetPy(((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py()+((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
	  theOtherZZllnnMET = theFakeMET[0];
	} else if(type4l != 3) { // (1) => nn, (2) ==> ll
	  printf("IMPOSSIBLE!!!: %f %f %d %d %d %d %d\n",minMassZ[0],minMassZ[1],tagZ[0],tagZ[1],tagZ[2],tagZ[3],type4l);
	}
        dPhiDiLepMET = TMath::Abs(dilepZll.DeltaPhi(theZZllnnMET));
        ptFrac = TMath::Abs(dilepZll.Pt()-theZZllnnMET.Pt())/dilepZll.Pt();
      }
      
      vector<int> idJet,idJetUp,idJetDown;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1 && ((TLorentzVector*)(*eventJets.p4)[nj])->Pt()> 30) {
          dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(theZZllnnMET));
        }

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()      > 30) idJet.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*1.03 > 30) idJetUp.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*0.97 > 30) idJetDown.push_back(nj);
      }

      passFilter[5]  = typeL[0] == 4 || typeL[1] == 4 || (typeL[0] == 2 && typeL[1] == 2);
      passFilter[6]  = minMassll > 4.0;
      passFilter[7]  = minMassZ[0] > 76.1876 && minMassZ[0] < 106.1876 && minMassZ[1] > 76.1876 && minMassZ[1] < 106.1876;
      passFilter[8]  = minMassZ[0] > 76.1876 && minMassZ[0] < 106.1876;
      passFilter[9]  = minMassZ[1] > 10 && minMassZ[1] < 65;
      passFilter[10] = mass4l > 140.0 || type4l == 3;
      passFilter[11] = (double)((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 20.;
      passFilter[12] = bDiscrMax < 0.560 && idJet.size() <= 1;
      passFilter[13] = fourLep.Pt() > 30 || type4l == 3;

      bool passNMinusOne = passFilter[5];
      bool passZZLoose = passFilter[5] && passFilter[6];
      bool passZZSel = passFilter[5] && passFilter[6] && passFilter[7];
      bool passZHWWSelNMinusOne[6] = {passFilter[6] &&                  passFilter[9] && passFilter[10] && passFilter[11] && passFilter[12] && passFilter[13],
                                      passFilter[6] && passFilter[8] &&                  passFilter[10] && passFilter[11] && passFilter[12] && passFilter[13],
				      passFilter[6] && passFilter[8] && passFilter[9]                   && passFilter[11] && passFilter[12] && passFilter[13],
				      passFilter[6] && passFilter[8] && passFilter[9] && passFilter[10]                   && passFilter[12] && passFilter[13],
				      passFilter[6] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11]                   && passFilter[13],
				      passFilter[6] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11] && passFilter[12]};
      bool passZHWWSel = passFilter[6] && passFilter[8] && passFilter[9] && passFilter[10] && passFilter[11] && passFilter[12] && passFilter[13];

      double dphill = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]));
      double detall = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Eta()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Eta());
      double drll = sqrt(dphill*dphill+detall*detall);
      bool useZHcuts = true;
      bool passDelphiLL   = useZHcuts ? (drll < 1.8) : (true);
      bool passNjets      = idJet.size() <= 1;
      bool passMET        = theZZllnnMET.Pt() > 50;
      bool passPTFrac     = ptFrac < (useZHcuts ? 0.4 : 1.0);
      bool passDPhiZMET   = dPhiDiLepMET > (useZHcuts ? 2.6 : 2.0);
      bool passPTLL       = dilepZll.Pt() > 60;
      bool passDPhiJetMET = useZHcuts ? (dPhiJetMET == -1 || dPhiJetMET >= 0.5) : (true);
      bool passZZhinvSel = {passZZSel && passNjets && passMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL};
      if(printMCEventList && infileName_[ifile].Contains("GluGlu") == kFALSE && infileCategory_[ifile] == 4 && theZZllnnMET.Pt() > 100 && passZZhinvSel)
        eventList << Form("%lld,%d,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d\n", eventEvent.eventNum, (int)idJet.size(), ((TLorentzVector*)(*eventMet.p4)[0])->Pt(), theZZllnnMET.Pt(), ptFrac, dPhiDiLepMET, dilepZll.Pt(), dPhiJetMET, (int)(*eventLeptons.pdgId)[idLep[tagZ[0]]], (int)(*eventLeptons.pdgId)[idLep[tagZ[1]]], (int)(*eventLeptons.pdgId)[idLep[tagZ[2]]], (int)(*eventLeptons.pdgId)[idLep[tagZ[3]]]);
      
      // Evaluate nominal BDT value
      double bdt_value=-1;
      if(useBDT && passZZhinvSel) {
        TLorentzVector lepton1 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]),
                       lepton2 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]),
                       MET     = theZZllnnMET,
                       jet1    = idJet.size() > 0 ? *((TLorentzVector*)(*eventJets.p4)[idJet[0]]) : TLorentzVector(0,0,0,0);
         bdt_value = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0,0,0,0);
      }
      
      bool passSystCuts[nSystTypes] = {
          passZZSel && idJetUp.size() <= 1   && passMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL,
	  passZZSel && idJetDown.size() <= 1 && passMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL,
          passZZSel && passNjets             && passMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL,
	  passZZSel && passNjets             && passMET && passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passDelphiLL
      };
      
      bool passAllCuts[nSelTypes] = {passZZSel,passZHWWSel,passZZhinvSel};

      // Do the BDT nuisance evaluations
      double bdt_muonScaleDown=-1, bdt_muonScaleUp=-1, bdt_electronScaleDown=-1, bdt_electronScaleUp=-1, bdt_METScaleDown=-1, bdt_METScaleUp=-1, bdt_jetScaleDown=-1, bdt_jetScaleUp=-1;
      double MVAVar_muonScaleDown=-9999, MVAVar_muonScaleUp=-9999, MVAVar_electronScaleDown=-9999, MVAVar_electronScaleUp=-9999, MVAVar_METScaleDown=-9999, MVAVar_METScaleUp=-9999, MVAVar_jetScaleDown=-9999, MVAVar_jetScaleUp=-9999;
      if(useBDT && passZZhinvSel) {
        TLorentzVector lepton1 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]),
                       lepton2 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]),
                       MET     = theZZllnnMET,
                       jet1    = idJet.size() > 0 ? *((TLorentzVector*)(*eventJets.p4)[idJet[0]]) : TLorentzVector(0,0,0,0);
        double lepton1_scale_variation, lepton2_scale_variation, neutrino1_scale_variation, neutrino2_scale_variation;
        // BDT variation with the muon scale variation (flat 1%)
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==13 ? 0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==13 ? 0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==13 ? 0.01 : 0;
        neutrino2_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[3]]])==13 ? 0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
         bdt_muonScaleUp = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_muonScaleUp = getMVAVar(MVAVarType, passZZhinvSel, typePair, MET.Pt(), 0, dilepZll.M(), bdt_muonScaleUp, xbins[nBinMVA]);
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==13 ? -0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==13 ? -0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==13 ? -0.01 : 0;
        neutrino2_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[3]]])==13 ? -0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
        bdt_muonScaleDown = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_muonScaleDown = getMVAVar(MVAVarType, passZZhinvSel, typePair, MET.Pt(), 0, dilepZll.M(), bdt_muonScaleUp, xbins[nBinMVA]);
        // BDT variation with the electron scale variation (flat 1%)
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==11 ? 0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==11 ? 0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==11 ? 0.01 : 0;
        neutrino2_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[3]]])==11 ? 0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
        bdt_electronScaleUp = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_electronScaleUp = getMVAVar(MVAVarType, passZZhinvSel, typePair, MET.Pt(), 0, dilepZll.M(), bdt_electronScaleUp, xbins[nBinMVA]);
        lepton1_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[0]]])==11 ? -0.01 : 0;
        lepton2_scale_variation = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[1]]])==11 ? -0.01 : 0;
        neutrino1_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]])==11 ? -0.01 : 0;
        neutrino2_scale_variation  = TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[3]]])==11 ? -0.01 : 0;
        MET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
        MET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py() + (1.+neutrino1_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py() + (1.+neutrino2_scale_variation)*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
        bdt_electronScaleDown = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation,0,0);
        MVAVar_electronScaleDown = getMVAVar(MVAVarType, passZZhinvSel, typePair, MET.Pt(), 0, dilepZll.M(), bdt_electronScaleUp, xbins[nBinMVA]);
        // BDT variation with the MET scale variation (from Nero)
        MET = theFakeMETUp;
        bdt_METScaleUp = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, 0, 0);
	MVAVar_METScaleUp = getMVAVar(MVAVarType, passZZhinvSel, typePair, MET.Pt(), 0, dilepZll.M(), bdt_METScaleUp, xbins[nBinMVA]);
        MET = theFakeMETDown;
        bdt_METScaleDown = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, 0, 0);
	MVAVar_METScaleDown = getMVAVar(MVAVarType, passZZhinvSel, typePair, MET.Pt(), 0, dilepZll.M(), bdt_METScaleUp, xbins[nBinMVA]);
      }

      // begin event weighting
      vector<int>wBoson;
      vector<int>zBoson;
      vector<bool> isGenDupl;double bosonPtMin = 1000000000; bool isBosonFound = false;vector<bool> isNeuDupl;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) wBoson.push_back(ngen0);
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23) zBoson.push_back(ngen0);
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
        if(isNeuDupl[ngen] == 0) numberGoodGenLep[2]++;
	if(isGenDupl[ngen] == 1) continue;
	numberGoodGenLep[0]++;
	if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt() <= 10 ||
	   TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()) >= 2.5) continue;
	numberGoodGenLep[1]++;
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

      // luminosity
      double theLumi  = 1.0; if(infileCategory_[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight     = 1.0; if(infileCategory_[ifile] != 0) puWeight     = nPUScaleFactor(fhDPU    , (double)eventMonteCarlo.puTrueInt);
      double puWeightUp   = 1.0; if(infileCategory_[ifile] != 0) puWeightUp   = nPUScaleFactor(fhDPUUp  , (double)eventMonteCarlo.puTrueInt);
      double puWeightDown = 1.0; if(infileCategory_[ifile] != 0) puWeightDown = nPUScaleFactor(fhDPUDown, (double)eventMonteCarlo.puTrueInt);
      // lepton efficiency
      double effSF = 1.0;
      if(infileCategory_[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,true);
        }
      }

      // fake rate
      int theCategory = infileCategory_[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     ((infileCategory_[ifile] == 0 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add Z+jets from data
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    theCategory = 1;
          }
          if     (infileCategory_[ifile] != 0 && goodIsTight == idTight.size()-3) fakeSF = -1.0 * fakeSF; // triple fake, MC
          else if(infileCategory_[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF = +1.0 * fakeSF; // double fake, MC
          else if(infileCategory_[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
          else if(infileCategory_[ifile] == 0 && goodIsTight == idTight.size()-3) fakeSF = +1.0 * fakeSF; // triple fake, data
          else if(infileCategory_[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
          else if(infileCategory_[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF = +1.0 * fakeSF; // single fake, data
        }
        else if(infileCategory_[ifile] != 0 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infileCategory_[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infileCategory_[ifile] == 0){ // data or Z+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infileCategory_[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infileCategory_[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale;
      if(totalWeight == 0) continue;

      double the_rho = 0.0; if(the_rhoP4.P() > 0) the_rho = the_rhoP4.Pt()/the_rhoP4.P();
      double theZZCorr[2] {1,1};
      if(theCategory == 4 && infileName_[ifile].Contains("GluGlu") == kFALSE) {
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

      // end event weighting
      if((infileCategory_[ifile] != 0 || theCategory == 0) && passZZSel) sumEventsProcess[ifile] += totalWeight;

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[i+type4l*nSelTypes][theCategory] += totalWeight;
          weiDecay[i+type4l*nSelTypes][theCategory] += totalWeight*totalWeight;
        }
      }

      for(int thePlot=0; thePlot<allPlots; thePlot++){
	double theVar = 0.0;
	bool makePlot = false;
	if     (thePlot ==  0 && passZZLoose)             {makePlot = true;theVar = TMath::Min((double)mass4l,399.999);}
	else if(thePlot ==  1 && passZZSel)               {makePlot = true;theVar = type4l;}
	else if(thePlot ==  2 && passNMinusOne)           {makePlot = true;theVar = TMath::Min(minMassll,99.999);}
	else if(thePlot ==  3 && passZZLoose)             {makePlot = true;theVar = TMath::Min(minMassZ[0],199.999);}
	else if(thePlot ==  4 && passZZLoose)             {makePlot = true;theVar = TMath::Min(minMassZ[1],199.999);}
	else if(thePlot ==  5 && passZZSel)               {makePlot = true;theVar = TMath::Min(deltaRllMin,3.999);}
	else if(thePlot ==  6 && passZZSel)               {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot ==  7 && passZZSel)               {makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot ==  8 && passZZSel)               {makePlot = true;theVar = dPhiJetMET*180/TMath::Pi();}
	else if(thePlot ==  9 && passZZSel)               {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 10 && passZZSel)               {makePlot = true;theVar = TMath::Min((double)mass4l,399.999);}
	else if(thePlot == 11 && passZZhinvSel) 	  {makePlot = true;theVar = TMath::Min(TMath::Max(theZZllnnMET.Pt()     ,0.001),399.999);}
	else if(thePlot == 12 && passZZhinvSel) 	  {makePlot = true;theVar = TMath::Min(TMath::Max(theOtherZZllnnMET.Pt(),0.001),399.999);}
	else if(thePlot == 13 && passZZSel)               {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 14 && passZHWWSelNMinusOne[3]) {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 15 && passZHWWSelNMinusOne[0]) {makePlot = true;theVar = TMath::Min(minMassZ[0],199.999);}
	else if(thePlot == 16 && passZHWWSelNMinusOne[1]) {makePlot = true;theVar = TMath::Min(minMassZ[1],199.999);}
	else if(thePlot == 17 && passZHWWSelNMinusOne[2]) {makePlot = true;theVar = TMath::Min((double)mass4l,399.999);}
	else if(thePlot == 18 && passZHWWSel)             {makePlot = true;theVar = type4l;}
	else if(thePlot == 19 && passZHWWSelNMinusOne[4]) {makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot == 20 && passZHWWSelNMinusOne[5]) {makePlot = true;theVar = TMath::Min(fourLep.Pt(),199.999);}
	else if(thePlot == 21 && passZHWWSel)             {makePlot = true;theVar = TMath::Min(mtHWW,199.999);}
	else if(thePlot == 22 && passZHWWSel)             {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]))*180/TMath::Pi();}
	else if(thePlot == 23 && passZHWWSel)             {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]]))*180/TMath::Pi();}
	else if(thePlot == 24 && passZHWWSel)             {makePlot = true;theVar = TMath::Abs(dilepZ.DeltaPhi(dilepWW))*180/TMath::Pi();}
	else if(thePlot == 25 && passZHWWSelNMinusOne[4]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 26 && passZHWWSel)             {makePlot = true;theVar = TMath::Min(dilepWW.Pt()/dilepZ.Pt(),1.999);}
	else if(thePlot == 27 && passZZSel)               {makePlot = true;theVar = TMath::Min(theZZllnnMET.Pt(),xbins[nBinMVA]-0.001);}
	else if(thePlot == 28 && passZZSel)               {makePlot = true;theVar = TMath::Min(dilepZll.Pt(),xbins[nBinMVA]-0.001);}
	else if(thePlot == 29 && passZZSel)               {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
	else if(thePlot == 30 && passZZSel)               {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 31 && passZZSel)               {makePlot = true;theVar = dPhiDiLepMET;}
	else if(thePlot == 32 && passZZSel)               {makePlot = true;theVar = dPhiJetMET;}
	else if(thePlot == 33 && passZZSel)               {makePlot = true;theVar = TMath::Abs(dilepZll.DeltaPhi(dilepZnn));}

	if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
      }
      
      if(1) {
	double MVAVar, MVAVarMETSyst[2];
	if(MVAVarType==1) {
          MVAVar = TMath::Min(theZZllnnMET.Pt(),xbins[nBinMVA]-0.001);
          MVAVarMETSyst[0] = TMath::Min(theFakeMETUp.Pt(),xbins[nBinMVA]-0.001);
          MVAVarMETSyst[1] = TMath::Min(theFakeMETDown.Pt(),xbins[nBinMVA]-0.001);
        } else if(MVAVarType==3) {
          MVAVar = getMVAVar(MVAVarType, passZZhinvSel, typePair, theZZllnnMET.Pt(), 0, dilepZll.M(), bdt_value, xbins[nBinMVA]);
          MVAVarMETSyst[0] = MVAVar; 
          MVAVarMETSyst[1] = MVAVar; 
        }

	// Avoid QCD scale weights that are anomalous high
	double maxQCDscale = (TMath::Abs((double)eventMonteCarlo.r1f2)+TMath::Abs((double)eventMonteCarlo.r1f5)+TMath::Abs((double)eventMonteCarlo.r2f1)+
                              TMath::Abs((double)eventMonteCarlo.r2f2)+TMath::Abs((double)eventMonteCarlo.r5f1)+TMath::Abs((double)eventMonteCarlo.r5f5))/6.0;

        if     (theCategory == 0){
	  if(passZZhinvSel) histo_Data->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 1){
	  if(passZZhinvSel) histo_Fake->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 4){
	  if(passZZhinvSel) {
	     histo_ZZ->Fill(MVAVar,totalWeight);
	     if(infileName_[ifile].Contains("GluGlu") == kFALSE) {
	       if(the_rho <= 0.3) histo_ZZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*(1.0+TMath::Abs((theZZCorr[0]-1)*(15.99/9.89-1))));
	       else               histo_ZZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*(1.0+TMath::Abs((theZZCorr[0]-1)               )));
	     } else {
               histo_ZZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight);
	     }
	     if(infileName_[ifile].Contains("GluGlu") == kFALSE) histo_ZZ_CMS_ggCorrUp->Fill(MVAVar,totalWeight);
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
             if     (type4l == 0) histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(type4l == 1) histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else                {histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}

             if     (type4l == 0) histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02*1.02);
             else if(type4l == 1) histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02*1.02);
             else                {histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02*1.02);
	                          histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02*1.02);}

             if     (type4l == 0) histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98*0.98);
             else if(type4l == 1) histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98*0.98);
             else                {histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98*0.98);
	                          histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98*0.98);}

             histo_ZZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
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
	  if(passZZhinvSel) {
	     histo_VVV->Fill(MVAVar,totalWeight);
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
             if     (type4l == 0) histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(type4l == 1) histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else                {histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}

             if     (type4l == 0) histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02*1.02);
             else if(type4l == 1) histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02*1.02);
             else                {histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02*1.02);
	                          histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02*1.02);}

             if     (type4l == 0) histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98*0.98);
             else if(type4l == 1) histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98*0.98);
             else                {histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98*0.98);
	                          histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98*0.98);}

             histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	     if(useBDT) {
               histo_VVV_CMS_BDTMuonScaleBoundingUp      ->Fill(MVAVar_muonScaleUp      , totalWeight);
               histo_VVV_CMS_BDTMuonScaleBoundingDown    ->Fill(MVAVar_muonScaleDown    , totalWeight);
               histo_VVV_CMS_BDTElectronScaleBoundingUp  ->Fill(MVAVar_electronScaleUp  , totalWeight);
               histo_VVV_CMS_BDTElectronScaleBoundingDown->Fill(MVAVar_electronScaleDown, totalWeight);
               histo_VVV_CMS_BDTMETScaleBoundingUp       ->Fill(MVAVar_METScaleUp       , totalWeight);
               histo_VVV_CMS_BDTMETScaleBoundingDown     ->Fill(MVAVar_METScaleDown     , totalWeight);
             }
	  }
          if(passSystCuts[JESUP])  histo_VVV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VVV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);	  
        }
        else if(theCategory == 6){
	  if(passZZhinvSel) {
	     histo_Higgs->Fill(MVAVar,totalWeight);
	     histo_Higgs_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (type4l == 0) histo_Higgs_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(type4l == 1) histo_Higgs_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else                {histo_Higgs_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_Higgs_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}

             if     (type4l == 0) histo_Higgs_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02*1.02);
             else if(type4l == 1) histo_Higgs_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02*1.02);
             else                {histo_Higgs_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02*1.02);
	                          histo_Higgs_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02*1.02);}

             if     (type4l == 0) histo_Higgs_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98*0.98);
             else if(type4l == 1) histo_Higgs_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98*0.98);
             else                {histo_Higgs_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98*0.98);
	                          histo_Higgs_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98*0.98);}

             histo_Higgs_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Higgs_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
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
      }
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);
    the_input_file->Close();
  } // end of chain

  if(printMCEventList) eventList.close();
  printf("QCD Corr: ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) Higgs(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Higgs->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  for(int i=1; i<=histo_Fake->GetNbinsX(); i++) {
    histo_Fake->SetBinContent(i,TMath::Max(histo_Fake->GetBinContent(i),0.0));
  }

  histo[allPlots-1][0]->Add(histo_Data);
  histo[allPlots-1][1]->Add(histo_Fake);
  histo[allPlots-1][4]->Add(histo_ZZ);
  histo[allPlots-1][5]->Add(histo_VVV);
  histo[allPlots-1][6]->Add(histo_Higgs);  

  double theZZSF = 1.0;
  if(1) {
    double theData = bgdDecay[ZZSEL+nSelTypes*0][0]+bgdDecay[ZZSEL+nSelTypes*1][0]+bgdDecay[ZZSEL+nSelTypes*2][0]+bgdDecay[ZZSEL+nSelTypes*3][0];
    double theSig  = bgdDecay[ZZSEL+nSelTypes*0][4]+bgdDecay[ZZSEL+nSelTypes*1][4]+bgdDecay[ZZSEL+nSelTypes*2][4]+bgdDecay[ZZSEL+nSelTypes*3][4];
    double theBck  = bgdDecay[ZZSEL+nSelTypes*0][1]+bgdDecay[ZZSEL+nSelTypes*1][1]+bgdDecay[ZZSEL+nSelTypes*2][1]+bgdDecay[ZZSEL+nSelTypes*3][1]+
                     bgdDecay[ZZSEL+nSelTypes*0][5]+bgdDecay[ZZSEL+nSelTypes*1][5]+bgdDecay[ZZSEL+nSelTypes*2][5]+bgdDecay[ZZSEL+nSelTypes*3][5]+
                     bgdDecay[ZZSEL+nSelTypes*0][6]+bgdDecay[ZZSEL+nSelTypes*1][6]+bgdDecay[ZZSEL+nSelTypes*2][6]+bgdDecay[ZZSEL+nSelTypes*3][6];
    theZZSF = (theData-theBck)/theSig; 
    printf("SF: %f (%f/%f/%f)\n",theZZSF,theData,theBck,theSig); // histo[allPlots-1][4]->Scale(theZZSF);
    histo[ 0][4]->Scale(theZZSF);
    histo[ 1][4]->Scale(theZZSF);
    histo[ 2][4]->Scale(theZZSF);
    histo[ 3][4]->Scale(theZZSF);
    histo[ 4][4]->Scale(theZZSF);
    histo[ 5][4]->Scale(theZZSF);
    histo[ 6][4]->Scale(theZZSF);
    histo[ 7][4]->Scale(theZZSF);
    histo[ 8][4]->Scale(theZZSF);
    histo[ 9][4]->Scale(theZZSF);
    histo[10][4]->Scale(theZZSF);
    histo[13][4]->Scale(theZZSF);
    histo[27][4]->Scale(theZZSF);
    histo[28][4]->Scale(theZZSF);
    histo[29][4]->Scale(theZZSF);
    histo[30][4]->Scale(theZZSF);
    histo[31][4]->Scale(theZZSF);
    histo[32][4]->Scale(theZZSF);
    histo[33][4]->Scale(theZZSF);
  }

  for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {
    double factorUp = +1.0; double factorDown = -1.0;
    histo_Fake_CMS_MVAFakeStatBoundingUp       ->SetBinContent(i,TMath::Max(histo_Fake	 ->GetBinContent(i)+factorUp  *histo_Fake   ->GetBinError(i),0.000001));
    histo_Fake_CMS_MVAFakeStatBoundingDown     ->SetBinContent(i,TMath::Max(histo_Fake       ->GetBinContent(i)+factorDown*histo_Fake	    ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	       ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorUp  *histo_ZZ	->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorDown*histo_ZZ	->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorUp  *histo_VVV	->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorDown*histo_VVV	->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingUp     ->SetBinContent(i,TMath::Max(histo_Higgs	->GetBinContent(i)+factorUp  *histo_Higgs ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingDown   ->SetBinContent(i,TMath::Max(histo_Higgs       ->GetBinContent(i)+factorDown*histo_Higgs   ->GetBinError(i),0.000001));

    histo_Fake_CMS_MVAFakeStatBoundingBinUp[i-1]	 ->Add(histo_Fake     ); histo_Fake_CMS_MVAFakeStatBoundingBinUp[i-1]	      ->SetBinContent(i,TMath::Max(histo_Fake	    ->GetBinContent(i)+factorUp  *histo_Fake	   ->GetBinError(i),0.000001));
    histo_Fake_CMS_MVAFakeStatBoundingBinDown[i-1]       ->Add(histo_Fake     ); histo_Fake_CMS_MVAFakeStatBoundingBinDown[i-1]	  ->SetBinContent(i,TMath::Max(histo_Fake	->GetBinContent(i)+factorDown*histo_Fake       ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	         ->Add(histo_ZZ       ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]           ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorUp  *histo_ZZ       ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	         ->Add(histo_ZZ       ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]         ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorDown*histo_ZZ       ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	         ->Add(histo_VVV      ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]         ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorUp  *histo_VVV      ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]         ->Add(histo_VVV      ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]       ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorDown*histo_VVV      ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1]	 ->Add(histo_Higgs    ); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1] 	  ->SetBinContent(i,TMath::Max(histo_Higgs	 ->GetBinContent(i)+factorUp  *histo_Higgs	 ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]	 ->Add(histo_Higgs    ); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]	  ->SetBinContent(i,TMath::Max(histo_Higgs	 ->GetBinContent(i)+factorDown*histo_Higgs	 ->GetBinError(i),0.000001));
  }

  double mean,up,diff;
  for(int i=1; i<=histo_ZZ->GetNbinsX(); i++){
    mean = histo_ZZ              ->GetBinContent(i);
    up   = histo_ZZ_CMS_EWKCorrUp->GetBinContent(i);
    diff = mean-up;
    histo_ZZ_CMS_EWKCorrDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));

    mean = histo_ZZ             ->GetBinContent(i);
    up   = histo_ZZ_CMS_ggCorrUp->GetBinContent(i);
    diff = mean-up;
    histo_ZZ_CMS_ggCorrDown->SetBinContent(i,TMath::Max(mean+diff,0.000001));
  }

  printf("EWK Corr: ZZ(%f/%f/%f) ggZZ(%f/%f/%f)\n",
                     histo_ZZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_EWKCorrDown->GetSumOfWeights(),
		     histo_ZZ_CMS_ggCorrUp ->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_ggCorrDown ->GetSumOfWeights());

  if(1){
    char outputLimits[200];
    sprintf(outputLimits,"MitZHAnalysis/plots%s/zz4l%shinv%s_input_%s.root", subdirectory.c_str(), addChan.Data(), finalStateName, ECMsb.Data());
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
    
    if(verbose) {
      cout << histo_Data             ->GetSumOfWeights() << " ";
      cout << histo_Fake             ->GetSumOfWeights() << " ";
      cout << histo_ZZ               ->GetSumOfWeights() << " ";
      cout << histo_VVV              ->GetSumOfWeights() << " ";
      cout << histo_Higgs            ->GetSumOfWeights() << " ";
      cout << endl;
      printf("uncertainties Stat\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_Fake	->GetBinContent(i)>0)printf("%5.1f ",histo_Fake_CMS_MVAFakeStatBoundingUp		 ->GetBinContent(i)/histo_Fake       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_Fake	->GetBinContent(i)>0)printf("%5.1f ",histo_Fake_CMS_MVAFakeStatBoundingDown 	 ->GetBinContent(i)/histo_Fake	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp		 ->GetBinContent(i)/histo_ZZ       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown 	 ->GetBinContent(i)/histo_ZZ	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_VVV	->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp 	 ->GetBinContent(i)/histo_VVV	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_VVV	->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown	 ->GetBinContent(i)/histo_VVV	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_Higgs	->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingUp		 ->GetBinContent(i)/histo_Higgs       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ  ->GetNbinsX(); i++) {if(histo_Higgs	->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingDown 	 ->GetBinContent(i)/histo_Higgs	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffM\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffE\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties MET\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties JES\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties PU\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingUp   ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingDown ->GetBinContent(i)/histo_Higgs   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
    } 
    histo_Fake ->Write();
    histo_Data ->Write();
    histo_ZZ   ->Write();
    histo_VVV  ->Write();
    histo_Higgs->Write();

    histo_Fake_CMS_MVAFakeStatBoundingUp       ->Write();
    histo_Fake_CMS_MVAFakeStatBoundingDown     ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingUp	       ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingDown	       ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingUp	       ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingDown       ->Write();
    histo_Higgs_CMS_MVAHiggsStatBoundingUp     ->Write();
    histo_Higgs_CMS_MVAHiggsStatBoundingDown   ->Write();
    
    histo_ZZ_CMS_MVALepEffMBoundingUp	       ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingDown	       ->Write();
    histo_VVV_CMS_MVALepEffMBoundingUp	       ->Write();
    histo_VVV_CMS_MVALepEffMBoundingDown       ->Write();
    histo_Higgs_CMS_MVALepEffMBoundingUp       ->Write();
    histo_Higgs_CMS_MVALepEffMBoundingDown     ->Write();
    
    histo_ZZ_CMS_MVALepEffEBoundingUp	       ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingDown	       ->Write();
    histo_VVV_CMS_MVALepEffEBoundingUp	       ->Write();
    histo_VVV_CMS_MVALepEffEBoundingDown       ->Write();
    histo_Higgs_CMS_MVALepEffEBoundingUp       ->Write();
    histo_Higgs_CMS_MVALepEffEBoundingDown     ->Write();
    
    histo_ZZ_CMS_MVAMETBoundingUp	       ->Write();
    histo_ZZ_CMS_MVAMETBoundingDown	       ->Write();
    histo_VVV_CMS_MVAMETBoundingUp	       ->Write();
    histo_VVV_CMS_MVAMETBoundingDown	       ->Write();
    histo_Higgs_CMS_MVAMETBoundingUp	       ->Write();
    histo_Higgs_CMS_MVAMETBoundingDown	       ->Write();
    
    histo_ZZ_CMS_MVAJESBoundingUp 	       ->Write();
    histo_ZZ_CMS_MVAJESBoundingDown	       ->Write();
    histo_VVV_CMS_MVAJESBoundingUp	       ->Write();
    histo_VVV_CMS_MVAJESBoundingDown	       ->Write();
    histo_Higgs_CMS_MVAJESBoundingUp 	       ->Write();
    histo_Higgs_CMS_MVAJESBoundingDown	       ->Write();
    
    histo_ZZ_CMS_PUBoundingUp	               ->Write();
    histo_ZZ_CMS_PUBoundingDown	               ->Write();
    histo_VVV_CMS_PUBoundingUp	               ->Write();
    histo_VVV_CMS_PUBoundingDown	       ->Write();
    histo_Higgs_CMS_PUBoundingUp	       ->Write();
    histo_Higgs_CMS_PUBoundingDown	       ->Write();

    histo_ZZ_CMS_EWKCorrUp	               ->Write();
    histo_ZZ_CMS_EWKCorrDown	               ->Write();
    histo_ZZ_CMS_ggCorrUp	               ->Write();
    histo_ZZ_CMS_ggCorrDown	               ->Write();

    outFileLimits->Close();

    double lumiE = 1.030;
    double systLepResE[3] = {1.01,1.01,1.01};
    double systLepResM[3] = {1.01,1.01,1.01};

    for(int nb=1; nb<=nBinMVA; nb++){
      // QCD study
      double systQCDScale[3] = {TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0]   ->GetBinContent(nb)-histo_ZZ   ->GetBinContent(nb)),
                                TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]  ->GetBinContent(nb)-histo_VVV  ->GetBinContent(nb)),
                                TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)   -histo_ZZ     ->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)   -histo_ZZ	->GetBinContent(nb));
        if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)  -histo_VVV    ->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)  -histo_VVV	->GetBinContent(nb));
        if(TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs  ->GetBinContent(nb)) > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs	->GetBinContent(nb));
      }                 
      if(histo_ZZ->GetBinContent(nb) > 0)    systQCDScale[0] = 1 + systQCDScale[0]/histo_ZZ    ->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histo_VVV->GetBinContent(nb) > 0)   systQCDScale[1] = 1 + systQCDScale[1]/histo_VVV   ->GetBinContent(nb); else systQCDScale[1] = 1;
      if(histo_Higgs->GetBinContent(nb) > 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_Higgs ->GetBinContent(nb); else systQCDScale[2] = 1;

      for(int ntype=0; ntype<3; ntype++) if(systQCDScale[ntype] < 0) systQCDScale[ntype] = 1.0;
      if(verbose) printf("QCDScale(%d): %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2]);
  
      // PDF study
      double systPDF[3];
      histo_Diff->Reset();
      for(int npdf=0; npdf<102; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ZZ->GetBinContent(nb))/histo_ZZ->GetBinContent(nb));
      systPDF[0] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=0; npdf<102; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
      systPDF[1] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      for(int npdf=0; npdf<102; npdf++) histo_Diff->Fill((histo_Higgs_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))/histo_Higgs->GetBinContent(nb));
      systPDF[2] = 1.0+histo_Diff->GetRMS();
      histo_Diff->Reset();
      if(verbose) printf("PDF(%d): %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2]);
  
      double systLepEffM[3] = {1.0,1.0,1.0};
      if     (histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp          ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown        ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingUp         ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)      > 0 && histo_Higgs_CMS_MVALepEffMBoundingUp       ->GetBinContent(nb) > 0) systLepEffM[2] = histo_Higgs_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)      > 0 && histo_Higgs_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[2] = histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
  
      double systLepEffE[3] = {1.0,1.0,1.0};
      if     (histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffEBoundingUp          ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
      else if(histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffEBoundingDown        ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffEBoundingUp         ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
      else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)      > 0 && histo_Higgs_CMS_MVALepEffEBoundingUp       ->GetBinContent(nb) > 0) systLepEffE[2] = histo_Higgs_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
      else if(histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)      > 0 && histo_Higgs_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[2] = histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
  
      double systMetUp  [3] = {1.0,1.0,1.0};
      double systMetDown[3] = {1.0,1.0,1.0};
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_MVAMETBoundingUp          ->GetBinContent(nb) > 0) systMetUp  [0] = histo_ZZ_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_MVAMETBoundingDown        ->GetBinContent(nb) > 0) systMetDown[0] = histo_ZZ_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_MVAMETBoundingUp         ->GetBinContent(nb) > 0) systMetUp  [1] = histo_VVV_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMetDown[1] = histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb)     > 0 && histo_Higgs_CMS_MVAMETBoundingUp       ->GetBinContent(nb) > 0) systMetUp  [2] = histo_Higgs_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb)     > 0 && histo_Higgs_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMetDown[2] = histo_Higgs_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      for(int i=0; i<3; i++) if(systMetUp  [i] == 1) systMetUp  [i] = 0.998;
      for(int i=0; i<3; i++) if(systMetDown[i] == 1) systMetDown[i] = 1.002;
      for(int nmet=0; nmet<3; nmet++) if(systMetUp[nmet]   > 1.04) systMetUp[nmet]   = 1.04;
      for(int nmet=0; nmet<3; nmet++) if(systMetUp[nmet]   < 0.96) systMetUp[nmet]   = 0.96;
      for(int nmet=0; nmet<3; nmet++) if(systMetDown[nmet] > 1.04) systMetDown[nmet] = 1.04;
      for(int nmet=0; nmet<3; nmet++) if(systMetDown[nmet] < 0.96) systMetDown[nmet] = 0.96;
  
      double systJesUp  [3] = {1.0,1.0,1.0};
      double systJesDown[3] = {1.0,1.0,1.0};
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVAJESBoundingUp  	->GetBinContent(nb) > 0) systJesUp  [0] = histo_ZZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[0] = histo_ZZ_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVAJESBoundingUp 	->GetBinContent(nb) > 0) systJesUp  [1] = histo_VVV_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[1] = histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_MVAJESBoundingUp  	->GetBinContent(nb) > 0) systJesUp  [2] = histo_Higgs_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb)	  > 0 && histo_Higgs_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[2] = histo_Higgs_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      
      double systBDTMuonUp  [3] = {1.0,1.0,1.0};
      double systBDTMuonDown[3] = {1.0,1.0,1.0};
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMuonScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMuonUp  [0] = histo_ZZ_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMuonScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMuonDown[0] = histo_ZZ_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMuonScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMuonUp  [1] = histo_VVV_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMuonScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMuonDown[1] = histo_VVV_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_CMS_BDTMuonScaleBoundingUp	->GetBinContent(nb) > 0) systBDTMuonUp  [2] = histo_Higgs_CMS_BDTMuonScaleBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_CMS_BDTMuonScaleBoundingDown ->GetBinContent(nb) > 0) systBDTMuonDown[2] = histo_Higgs_CMS_BDTMuonScaleBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
   
      double systBDTElectronUp  [3] = {1.0,1.0,1.0};
      double systBDTElectronDown[3] = {1.0,1.0,1.0};
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTElectronScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTElectronUp  [0] = histo_ZZ_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTElectronScaleBoundingDown	->GetBinContent(nb) > 0) systBDTElectronDown[0] = histo_ZZ_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTElectronScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTElectronUp  [1] = histo_VVV_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTElectronScaleBoundingDown	->GetBinContent(nb) > 0) systBDTElectronDown[1] = histo_VVV_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_CMS_BDTElectronScaleBoundingUp	->GetBinContent(nb) > 0) systBDTElectronUp  [2] = histo_Higgs_CMS_BDTElectronScaleBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_CMS_BDTElectronScaleBoundingDown ->GetBinContent(nb) > 0) systBDTElectronDown[2] = histo_Higgs_CMS_BDTElectronScaleBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
   
      double systBDTMETUp  [3] = {1.0,1.0,1.0};
      double systBDTMETDown[3] = {1.0,1.0,1.0};
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMETScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMETUp  [0] = histo_ZZ_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_BDTMETScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMETDown[0] = histo_ZZ_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMETScaleBoundingUp  	->GetBinContent(nb) > 0) systBDTMETUp  [1] = histo_VVV_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_BDTMETScaleBoundingDown	->GetBinContent(nb) > 0) systBDTMETDown[1] = histo_VVV_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_CMS_BDTMETScaleBoundingUp	->GetBinContent(nb) > 0) systBDTMETUp  [2] = histo_Higgs_CMS_BDTMETScaleBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb) > 0 && histo_Higgs_CMS_BDTMETScaleBoundingDown ->GetBinContent(nb) > 0) systBDTMETDown[2] = histo_Higgs_CMS_BDTMETScaleBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
 
      double systPUUp  [3] = {1.0,1.0,1.0};
      double systPUDown[3] = {1.0,1.0,1.0};
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_PUBoundingUp          ->GetBinContent(nb) > 0) systPUUp  [0] = histo_ZZ_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)        > 0 && histo_ZZ_CMS_PUBoundingDown        ->GetBinContent(nb) > 0) systPUDown[0] = histo_ZZ_CMS_PUBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_PUBoundingUp         ->GetBinContent(nb) > 0) systPUUp  [1] = histo_VVV_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)       > 0 && histo_VVV_CMS_PUBoundingDown       ->GetBinContent(nb) > 0) systPUDown[1] = histo_VVV_CMS_PUBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb)     > 0 && histo_Higgs_CMS_PUBoundingUp       ->GetBinContent(nb) > 0) systPUUp  [2] = histo_Higgs_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      if(histo_Higgs->GetBinContent(nb)     > 0 && histo_Higgs_CMS_PUBoundingDown     ->GetBinContent(nb) > 0) systPUDown[2] = histo_Higgs_CMS_PUBoundingDown->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
      for(int npu=0; npu<3; npu++) if(systPUUp[npu]   > 1.02) systPUUp[npu]   = 1.02;
      for(int npu=0; npu<3; npu++) if(systPUUp[npu]   < 0.98) systPUUp[npu]   = 0.98;
      for(int npu=0; npu<3; npu++) if(systPUDown[npu] > 1.02) systPUDown[npu] = 1.02;
      for(int npu=0; npu<3; npu++) if(systPUDown[npu] < 0.98) systPUDown[npu] = 0.98;
   
      double syst_EWKCorrUp[2]   = {1.0,1.0};
      double syst_EWKCorrDown[2] = {1.0,1.0};
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [0] = histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[0] = histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb) > 0) syst_EWKCorrUp  [1] = histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb) > 0) syst_EWKCorrDown[1] = histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);

      char outputLimitsShape[200];                                            
      sprintf(outputLimitsShape,"MitZHAnalysis/datacards%s/histo_limits_zz4l%shinv%s_shape_%s_bin%d.txt",subdirectory.c_str(), addChan.Data(),finalStateName, ECMsb.Data(),nb-1);
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
      newcardShape << Form("bin hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process ZH_hinv ZZ VVV Higgs Fake\n");
      newcardShape << Form("process 0 1 2 3 4\n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f\n",0.0,histo_ZZ->GetBinContent(nb),TMath::Max(histo_VVV->GetBinContent(nb),0.0),TMath::Max(histo_Higgs->GetBinContent(nb),0.0),TMath::Max(histo_Fake->GetBinContent(nb),0.0));
      newcardShape << Form("lumi_%4s                               lnN  1.000 %7.5f %7.5f %7.5f   -\n",ECMsb.Data(),lumiE,lumiE,lumiE);		       
      newcardShape << Form("%s                                     lnN  1.000 %7.5f %7.5f %7.5f   -\n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2]);
      newcardShape << Form("%s                                     lnN  1.000 %7.5f %7.5f %7.5f   -\n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2]);
      newcardShape << Form("%s                                     lnN  1.000 %7.5f %7.5f %7.5f   -\n",momMName,systLepResM[0],systLepResM[1],systLepResM[2]);
      newcardShape << Form("%s                                     lnN  1.000 %7.5f %7.5f %7.5f   -\n",momEName,systLepResE[0],systLepResE[1],systLepResE[2]);
      newcardShape << Form("CMS_trigger2016                        lnN  1.000 %7.5f %7.5f %7.5f   -\n",1.01,1.01,1.01);
      newcardShape << Form("pdf_qqbar_ACCEPT                       lnN  1.000 %7.5f %7.5f %7.5f   -\n",TMath::Max(systPDF[0],1.01),TMath::Max(systPDF[1],1.01),TMath::Max(systPDF[2],1.01));
      newcardShape << Form("QCDscale_VVV                           lnN    -     -   %7.5f   -     -\n",systQCDScale[1]);
      newcardShape << Form("QCDscale_ZZ		                   lnN    -   %7.5f   -     -     -\n",systQCDScale[0]);
      newcardShape << Form("QCDscale_ggH		           lnN    -     -     -   %7.5f   -\n",systQCDScale[2]);

      newcardShape << Form("CMS_pu2016                             lnN  1.000 %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -\n",systPUUp[0] ,systPUDown[0] ,systPUUp[1] ,systPUDown[1] ,systPUUp[2] ,systPUDown[2] );
      newcardShape << Form("CMS_scale_met                          lnN  1.000 %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -\n",systMetUp[0],systMetDown[0],systMetUp[1],systMetDown[1],systMetUp[2],systMetDown[2]);
      newcardShape << Form("CMS_scale_j                            lnN  1.000 %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -\n",systJesUp[0],systJesDown[0],systJesUp[1],systJesDown[1],systJesUp[2],systJesDown[2]);   	     

      if(nb != 1){
      if(useZZWZEWKUnc){
      newcardShape << Form("CMS_zllhinv_ZZWZ_EWKCorr               lnN    -   %7.5f   -     -     -\n",1.+sqrt(0.02*0.02+(syst_EWKCorrUp[0]-1.0)*(syst_EWKCorrUp[0]-1.0)));		
      newcardShape << Form("CMS_hinv_vvnorm_bin%d rateParam  * ZZ 1 [0.1,10]\n",nb-1);	
      } else {
      newcardShape << Form("CMS_hinv_zznorm_bin%d rateParam  * ZZ 1 [0.1,10]\n",nb-1);	
      }
      }

      newcardShape << Form("CMS_zllhinv_ggZZCorr                   lnN    -     -     -     -   %7.5f/%7.5f   -     -  \n",syst_EWKCorrUp[1],syst_EWKCorrDown[1]);		
      if(MVAVarType == 3 || MVAVarType==4) {
      newcardShape << Form("CMS_BDT_scale_m                   lnN	- %7.5f/%7.5f  %7.5f/%7.5f   %7.5f/%7.5f   -  \n", systBDTMuonUp[0], systBDTMuonDown[0], systBDTMuonUp[1], systBDTMuonDown[1], systBDTMuonUp[2], systBDTMuonDown[2]);	    
      newcardShape << Form("CMS_BDT_scale_e                   lnN	- %7.5f/%7.5f  %7.5f/%7.5f   %7.5f/%7.5f   -  \n", systBDTElectronUp[0], systBDTElectronDown[0], systBDTElectronUp[1], systBDTElectronDown[1], systBDTElectronUp[2], systBDTElectronDown[2]);	    
      newcardShape << Form("CMS_BDT_scale_MET                 lnN	- %7.5f/%7.5f  %7.5f/%7.5f   %7.5f/%7.5f   -  \n", systBDTMETUp[0], systBDTMETDown[0], systBDTMETUp[1], systBDTMETDown[1], systBDTMETUp[2], systBDTMETDown[2]);	    
      }

      if(histo_ZZ   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding2016_%s_Bin%d      lnN     1.000 %7.5f   -	-    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ZZ   ->GetBinError(nb)/histo_ZZ   ->GetBinContent(nb));
      if(histo_VVV  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding2016_%s_Bin%d     lnN     1.000   -	%7.5f	-    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_VVV  ->GetBinError(nb)/histo_VVV  ->GetBinContent(nb));  
      if(histo_Higgs->GetBinContent(nb) > 0) newcardShape << Form("CMS_zllhinv%s_MVAHiggsStatBounding2016_%s_Bin%d   lnN     1.000   -	  -   %7.5f  -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_Higgs->GetBinError(nb)/histo_Higgs->GetBinContent(nb));
      if(histo_Fake ->GetBinContent(nb) > 0) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding2016_%s_Bin%d      lnN     1.000   -	  -	-  %7.5f\n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_Fake ->GetBinError(nb)/histo_Fake ->GetBinContent(nb));
    }
  }

  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsType[4] = {0,0,0,0}; double sumEventsTypeE[4] = {0,0,0,0};
    for(int np=0; np<histBins; np++) {       
      double theSF = 1.0;
      if((ns == ZZSEL) && np == 4) theSF = theZZSF;
       printf("(%6s): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f -> %8.2f +/- %6.2f\n",
       processName[np].Data(),bgdDecay[ns+nSelTypes*0][np]*theSF,sqrt(weiDecay[ns+nSelTypes*0][np])*theSF,bgdDecay[ns+nSelTypes*1][np]*theSF,sqrt(weiDecay[ns+nSelTypes*1][np])*theSF,
                              bgdDecay[ns+nSelTypes*2][np]*theSF,sqrt(weiDecay[ns+nSelTypes*2][np])*theSF,bgdDecay[ns+nSelTypes*3][np]*theSF,sqrt(weiDecay[ns+nSelTypes*3][np])*theSF,
                              (bgdDecay[ns+nSelTypes*0][np]+bgdDecay[ns+nSelTypes*1][np]+bgdDecay[ns+nSelTypes*2][np]+bgdDecay[ns+nSelTypes*3][np])*theSF,
                              sqrt(weiDecay[ns+nSelTypes*0][np]+weiDecay[ns+nSelTypes*1][np]+weiDecay[ns+nSelTypes*2][np]+weiDecay[ns+nSelTypes*3][np])*theSF);
       if(np != 0){
         sumEventsType[0] = sumEventsType[0] + bgdDecay[ns+nSelTypes*0][np]*theSF; sumEventsTypeE[0] = sumEventsTypeE[0] + weiDecay[ns+nSelTypes*0][np]*theSF;
         sumEventsType[1] = sumEventsType[1] + bgdDecay[ns+nSelTypes*1][np]*theSF; sumEventsTypeE[1] = sumEventsTypeE[1] + weiDecay[ns+nSelTypes*1][np]*theSF;
         sumEventsType[2] = sumEventsType[2] + bgdDecay[ns+nSelTypes*2][np]*theSF; sumEventsTypeE[2] = sumEventsTypeE[2] + weiDecay[ns+nSelTypes*2][np]*theSF;
         sumEventsType[3] = sumEventsType[3] + bgdDecay[ns+nSelTypes*3][np]*theSF; sumEventsTypeE[3] = sumEventsTypeE[3] + weiDecay[ns+nSelTypes*3][np]*theSF;
      }
    }
    printf("(...bkg): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f -> %8.2f +/- %6.2f\n",
           sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
	   sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]),
           sumEventsType[0]+sumEventsType[1]+sumEventsType[2]+sumEventsType[3],
           sqrt(sumEventsTypeE[0]+sumEventsTypeE[1]+sumEventsTypeE[2]+sumEventsTypeE[3]));
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"MitZHAnalysis/plots%s/histozz_nice_%d.root",subdirectory.c_str(), thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }
}
