#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"
#include "NeroProducer/Core/interface/BarePhotons.hpp"

#include "MitAnalysisRunII/macros/76x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

bool usePureMC = true; 
double mcPrescale = 1.0;
bool usingAllTriggers = false;
enum selType                     {SIGSEL=0,  ZSEL,   NOBTAG,   ZHPRESEL0,   ZHPRESEL1,   ZHSEL0,   ZHSEL1,   ZHSEL2,  ZHRECSEL0,   ZHRECSEL1,   ZHRECSEL2,  GJETSEL, nSelTypes};
TString selTypeName[nSelTypes]= {"SIGSEL",  "ZSEL", "NOBTAG", "ZHPRESEL0", "ZHPRESEL1", "ZHSEL0", "ZHSEL1", "ZHSEL2","ZHRECSEL0", "ZHRECSEL1", "ZHRECSEL2","GJETSEL"};

// nsel == 0 --> llg selection, pt>20/20/20, lepton triggers
//         1 --> llg selection, pt>20/20/60, photon trigger
//         2 --> g selection, pt>60, photon triggers
//         3 --> ll selection, pt> 20,20, lepton triggers, pt_ll > 60
//         4 --> lg selection, pt>20/100, photon trigger
//         5 --> g+jets(->fake lepton) selection, pt>60, photon triggers

void zgAnalysis(
 Int_t nsel = 0,
 Int_t typeSel = 3,
 TString typeLepSel = "medium",
 Int_t period = 1
 ){

  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/";
  Double_t lumi = 2.318;

  double denFRDAM[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double numFRDAM[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double denFRBGM[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double numFRBGM[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double denFRDAE[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double numFRDAE[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double denFRBGE[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double numFRBGE[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  if(nsel == 1 || nsel == 2 || nsel == 4 || nsel == 5) filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/";
  if(nsel == 1 || nsel == 2 || nsel == 4 || nsel == 5) filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/pho_";
  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath       = "";
  TString ptRatioPath  = "";
  TString etaRatioPath = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";
  if(nsel == 1 || nsel == 2 || nsel == 4 || nsel == 5) {
  infilenamev.push_back(Form("%sSinglePhoton+Run2015C_25ns-16Dec2015-v1+AOD.root",filesPathDA.Data())); 												 infilecatv.push_back(0);
  infilenamev.push_back(Form("%sSinglePhoton+Run2015D-16Dec2015-v1+AOD.root",filesPathDA.Data())); 												 infilecatv.push_back(0);
  } else {
  infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  }  

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));						   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(1);

  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                            infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(1);

  if(nsel == 2 || nsel == 5) {
  infilenamev.push_back(Form("%sGJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));          infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));          infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));          infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));          infilecatv.push_back(2);
  } else {
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                         infilecatv.push_back(2);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  }

  if(nsel == 1 || nsel == 2 || nsel == 4 || nsel == 5) {
  infilenamev.push_back(Form("%sDYJetsToNuNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));              infilecatv.push_back(3);
  }

  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(3);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(3);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(4);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);

  }
  else {assert(0);}
  
  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Medium_ele"));
  TH2D *fhDElTightSF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Tight_ele"));
  TH2D *fhDElMediumMVASF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_MediumMVA_ele"));
  TH2D *fhDElTightMVASF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_TightMVA_ele"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  assert(fhDElMediumMVASF);
  assert(fhDElTightMVASF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF ->SetDirectory(0);
  fhDElMediumMVASF->SetDirectory(0);
  fhDElTightMVASF ->SetDirectory(0);
  delete fElSF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Medium_mu"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Iso_mu"));
  assert(fhDMuMediumSF);
  assert(fhDMuIsoSF);
  fhDMuMediumSF->SetDirectory(0);
  fhDMuIsoSF->SetDirectory(0);
  delete fMuSF;

  double eventsTrg[4] = {0,0,0,0};
  double dataPrescale[4] = {21.479639,4.339259,2.163135,1};

  const int MVAVarType = 0; const int nBinMVA = 6; Float_t xbins[nBinMVA+1] = {200, 250, 300, 400, 600, 800, 1000};
  //const int MVAVarType = 1; const int nBinMVA = 6; Float_t xbins[nBinMVA+1] = {100, 125, 150, 175, 200, 225, 250};
  //const int MVAVarType = 2; const int nBinMVA = 17; Float_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 600, 800, 1000};
  //const int MVAVarType = 3; const int nBinMVA = 23; Float_t xbins[nBinMVA+1] = {45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 160, 170, 180, 190, 200};
  //const int MVAVarType = 4; const int nBinMVA = 20; Float_t xbins[nBinMVA+1] =  {0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();
  TH1D *histo_Zjets0 = (TH1D*) histoMVA->Clone("histo_Zjets0");
  TH1D *histo_Zjets1 = (TH1D*) histoMVA->Clone("histo_Zjets1");
  TH1D *histo_Zjets2 = (TH1D*) histoMVA->Clone("histo_Zjets2");
  TH1D *histo_ZjetsRec0 = (TH1D*) histoMVA->Clone("histo_ZjetsRec0");
  TH1D *histo_ZjetsRec1 = (TH1D*) histoMVA->Clone("histo_ZjetsRec1");
  TH1D *histo_ZjetsRec2 = (TH1D*) histoMVA->Clone("histo_ZjetsRec2");

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 40;
  const int histBins = 5;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "....VV", "....ZG", "...ZLL"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  2) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 100; xminPlot = -TMath::Pi();; xmaxPlot = TMath::Pi();}
    else if(thePlot >=  5 && thePlot <=  6) {nBinPlot = 1000; xminPlot =-100.0; xmaxPlot = 100.0;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >=  8 && thePlot <=  9) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot >= 11 && thePlot <= 11) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 12 && thePlot <= 12) {nBinPlot = 600; xminPlot = 0.0; xmaxPlot = 600.0;}
    else if(thePlot >= 13 && thePlot <= 13) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 14 && thePlot <= 15) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   2.5;}
    else if(thePlot >= 16 && thePlot <= 16) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 17 && thePlot <= 17) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   2.0;}
    else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =  0.05;}
    else if(thePlot >= 19 && thePlot <= 20) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 21 && thePlot <= 22) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot >= 23 && thePlot <= 24) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 25 && thePlot <= 26) {nBinPlot = 600; xminPlot = 0.0; xmaxPlot = 600.0;}
    else if(thePlot >= 27 && thePlot <= 28) {nBinPlot = 600; xminPlot = 0.0; xmaxPlot = 600.0;}
    else if(thePlot >= 29 && thePlot <= 30) {nBinPlot = 100; xminPlot = -TMath::Pi(); xmaxPlot = TMath::Pi();}
    else if(thePlot >= 31 && thePlot <= 32) {nBinPlot = 100; xminPlot = 0; xmaxPlot = TMath::Pi();}
    else if(thePlot >= 33 && thePlot <= 34) {nBinPlot = 600; xminPlot = 0.0; xmaxPlot = 600.0;}
    else if(thePlot >= 35 && thePlot <= 36) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >= 37 && thePlot <= 39) {nBinPlot =  10; xminPlot =-0.5; xmaxPlot =   9.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Clear();
  }

  double nSigCut[nSelTypes*4],nSigECut[nSelTypes*4];
  double bgdDecay[nSelTypes*4][histBins],weiDecay[nSelTypes*4][histBins];
  for(unsigned int i=0; i<nSelTypes*4; i++) {
    nSigCut[i] = 0.0; nSigECut[i] = 0.0;
    for(int j=0; j<histBins; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    eventJets.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BarePhotons eventPhotons;
    eventPhotons.SetExtend();
    eventPhotons.setBranchAddresses(the_input_tree);

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
    eventMonteCarlo.setBranchAddresses(the_input_tree);

    TNamed *triggerNames = (TNamed*)the_input_file.FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(infilecatv[ifile] == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }
    else {
      printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());
    }

    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_input_tree->GetEntry(i);

      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[4] = {kFALSE,kFALSE,kFALSE,kFALSE};
      vector<int> idLep; vector<int> idSoft;vector<int> idFake;  vector<int> idTight;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idLep.push_back(nlep);idFake.push_back(nlep);idTight.push_back(1);}
        else {
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP) {idSoft.push_back(nlep);}
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake)   {idFake.push_back(nlep);idTight.push_back(0);}
        }
      }

      if     ((nsel == 0 || nsel == 1 || nsel == 3) && idLep.size() == 2 && 
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 20 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20) passFilter[0] = kTRUE;
      else if(nsel == 2 && idLep.size() == 0)                        passFilter[0] = kTRUE;
      else if(nsel == 4 && idLep.size() == 1 && 
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 20) passFilter[0] = kTRUE;
      else if(nsel == 5)                                             passFilter[0] = kTRUE;
      if(passFilter[0] == kFALSE) continue;

      int signQ = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
      }
      if(signQ == 0 || nsel == 4 || nsel == 5) passFilter[1] = kTRUE;
      if(passFilter[1] == kFALSE) continue;

      vector<int> idPho;
      for(int npho=0; npho<eventPhotons.p4->GetEntriesFast(); npho++) {
        if(TMath::Abs(((TLorentzVector*)(*eventPhotons.p4)[npho])->Eta()) >= 1.479) continue;
	if(((TLorentzVector*)(*eventPhotons.p4)[npho])->Pt() <= 20) continue;
	if((double)(*eventPhotons.r9)[npho] <= 0.9) continue;
	//if((double)(*eventPhotons.sieie)[npho] >= 0.011) continue;
        bool isRecoLepton = false;
	for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventPhotons.p4)[npho])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3)
	    {isRecoLepton = true; break;}
        }
	if(isRecoLepton == true) continue;
        if(((int)(*eventPhotons.selBits)[npho] & BarePhotons::PhoTight)== BarePhotons::PhoTight){idPho.push_back(npho);}
      }
      if     (nsel == 0 && idPho.size() >= 1 &&
         ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt() > 20) passFilter[2] = kTRUE;
      else if((nsel == 1 || nsel == 2 || nsel == 5) && idPho.size() >= 1 &&
         ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt() > 60) passFilter[2] = kTRUE;
      else if(nsel == 4 && idPho.size() >= 1 &&
         ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt() > 130) passFilter[2] = kTRUE;
      else if(nsel == 3 && ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) ).Pt() > 60) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      double theDataPrescale = 1.0;
      for (int nt = 0; nt <(int)numtokens; nt++) {
        if((*eventTrigger.triggerFired)[nt] == 0) continue;
        if     (nsel == 0 || nsel == 3){
	  if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
             (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*") 	      == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu27_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu20_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoTkMu20_v*")				      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele23_WPLoose_Gsf_v*")			      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele22_eta2p1_WP75_Gsf_v*")			      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele27_WPLoose_Gsf_v*")			      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele27_WP85_Gsf_v*")  			      == 0)
             ) passFilter[3] = kTRUE;
          //if(infilecatv[ifile] != 0) passFilter[3] = kTRUE; // do not apply trigger filters to MC
        }
        else if(nsel == 2 || nsel == 5){
	  if(usingAllTriggers == false){
	  if     (strcmp(tokens[nt],"HLT_Photon50_R9Id90_HE10_IsoM_v*") == 0 && ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt() <=  85) {passFilter[3] = kTRUE; theDataPrescale = dataPrescale[0];}
	  else if(strcmp(tokens[nt],"HLT_Photon75_R9Id90_HE10_IsoM_v*") == 0 && ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt() <=  98) {passFilter[3] = kTRUE; theDataPrescale = dataPrescale[1];}
	  else if(strcmp(tokens[nt],"HLT_Photon90_R9Id90_HE10_IsoM_v*") == 0 && ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt() <= 130) {passFilter[3] = kTRUE; theDataPrescale = dataPrescale[2];}
	  else if(strcmp(tokens[nt],"HLT_Photon120_R9Id90_HE10_IsoM_v*") == 0 							  	    ) {passFilter[3] = kTRUE; theDataPrescale = dataPrescale[3];}
          } else {
	  if((strcmp(tokens[nt],"HLT_Photon50_R9Id90_HE10_IsoM_v*") == 0 ) ||
             (strcmp(tokens[nt],"HLT_Photon75_R9Id90_HE10_IsoM_v*") == 0 )||
             (strcmp(tokens[nt],"HLT_Photon90_R9Id90_HE10_IsoM_v*") == 0 )||
             (strcmp(tokens[nt],"HLT_Photon120_R9Id90_HE10_IsoM_v*") == 0)
             ) passFilter[3] = kTRUE;
          }
	}
        else if(nsel == 1 || nsel == 4){
	  if(strcmp(tokens[nt],"HLT_Photon120_R9Id90_HE10_IsoM_v*") == 0) passFilter[3] = kTRUE;
	}
      }
      if(passFilter[3] == kFALSE) continue;

      TLorentzVector dilep,dilepPho;
      if     (nsel == 0 || nsel == 1) {
        dilep = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) ); 
        dilepPho = dilep + ( *(TLorentzVector*)(eventPhotons.p4->At(idPho[0])) );
      }
      else if(nsel == 2 || nsel == 5) {
        dilep = ( *(TLorentzVector*)(eventPhotons.p4->At(idPho[0])) );
        dilepPho = dilep;
      }
      else if(nsel == 3) {
        dilep = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) ); 
        dilepPho = dilep;
      }
      else if(nsel == 4) {
        dilep = ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ); 
        dilepPho = dilep;
      }
      TLorentzVector recoil = dilep;

      // abs(eta) requirement
      if(TMath::Abs(dilep.Eta()) >= 1.479) continue;

      TLorentzVector theMET;
      theMET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px());
      theMET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py());
      //if(infilecatv[ifile] != 0) {
      //  theMET.SetPx(theMET.Px() - 2.614410);
      // theMET.SetPy(theMET.Py() + 1.727615);
      //}

      vector<int> idJet;
      vector<int> idGJet;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isPhoton = kFALSE;
        for(unsigned int np=0; np<idPho.size(); np++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventPhotons.p4)[idPho[np]])) < 0.3) isPhoton = kTRUE;
	}
	if(isPhoton == kTRUE) continue;

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 30) idGJet.push_back(nj);

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(theMET));

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 10 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;
        
	idJet.push_back(nj);
	recoil = recoil + ( *(TLorentzVector*)(eventJets.p4->At(nj)) );
      }

      int typePair = 0;
      if(nsel == 0 || nsel == 1 || nsel == 3){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typePair = 1;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typePair = 2;
      }
      else if(nsel == 2 || nsel == 5){
        typePair = 1;
      }
      else if(nsel == 4){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13) typePair = 1;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11) typePair = 2;
      }

      double dRLepPhoMin = 999.;
      double ptPho = 0;
      double etaPho = 0;
      double isoPho = 0;
      double sieiePho = 999.;
      double dRLepFakePhoMin = 999.;
      if(idPho.size() >= 1){
        ptPho = ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt();
	etaPho = TMath::Abs(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Eta());
	isoPho = (double)((*eventPhotons.iso)[idPho[0]]/((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->Pt());
	sieiePho = (double)(*eventPhotons.sieie)[idPho[0]];
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < dRLepPhoMin)
	    dRLepPhoMin = ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]]));
        }
        for(unsigned int nl=0; nl<idFake.size(); nl++){
          if(((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idFake[nl]])) < dRLepFakePhoMin)
	    dRLepFakePhoMin = ((TLorentzVector*)(*eventPhotons.p4)[idPho[0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idFake[nl]]));
        }
      }

      int numberGoodTaus = 0;
      for(int ntau=0; ntau<eventTaus.p4->GetEntriesFast(); ntau++) {
        if(((TLorentzVector*)(*eventTaus.p4)[ntau])->Pt() <= 20.0 ||
           TMath::Abs(((TLorentzVector*)(*eventTaus.p4)[ntau])->Eta()) >= 2.3) continue;
        bool isElMu = false;
        for(unsigned nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventTaus.p4)[ntau])) < 0.1) {
            isElMu = true;
            break;
          }
        }
        if(isElMu == false &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFinding	 ) == BareTaus::TauDecayModeFinding &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFindingNewDMs) == BareTaus::TauDecayModeFindingNewDMs &&
           (double)(*eventTaus.iso)[ntau] < 4.0){
          numberGoodTaus++;
        }
      }

      double dPhiDiLepMET = TMath::Abs(dilep.DeltaPhi(theMET));
      double dPhiRecoilMET = TMath::Abs(recoil.DeltaPhi(theMET));
      double ptFrac[3] = {TMath::Abs(dilep.Pt()   -theMET.Pt())/dilep.Pt(),
                          TMath::Abs(dilepPho.Pt()-theMET.Pt())/dilepPho.Pt(),
			  TMath::Abs(recoil.Pt()  -theMET.Pt())/recoil.Pt()};
      double deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(theMET));
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*theMET.Pt()*(1.0 - cos(deltaPhiDileptonMet)));

      bool passZGMass   = TMath::Abs(dilepPho.M()-91.1876) < 10.0;
      bool passZMass    = dilep.M() > 30.;
      bool passBtagVeto = true;
      if(nsel == 2 || nsel == 5) {
        passZGMass = true;
	passZMass = true;
      }
      else if(nsel == 4) {
        passZGMass = theMET.Pt() > 30.;
	passZMass = mtW > 30.;
        passBtagVeto = bDiscrMax < 0.800 && idSoft.size() == 0;
      }
      double metMIN = 100; double mtMIN = 200;
      if     (MVAVarType == 2 || MVAVarType == 3 || MVAVarType == 4) {metMIN = 45; mtMIN = 0;}
      bool passMET = theMET.Pt() > metMIN && mtW > mtMIN;

      bool passPTFrac    = ptFrac[0] < 0.4;
      if(MVAVarType == 4) passPTFrac = ptFrac[0] < 1.0;
      bool passDPhiZMET  = dPhiDiLepMET > 2.8;
      bool passPTLL      = dilep.Pt() > 60;

      bool passPTFracRecoil = ptFrac[2] < 0.4;
      if(MVAVarType == 4) passPTFracRecoil = ptFrac[2] < 1.0;
      bool passDPhiZMETRecoil = dPhiRecoilMET > 2.8;

      bool passDPhiJetMET = dPhiJetMET == -1 || dPhiJetMET >= 0.5;
      bool passTauVeto    = numberGoodTaus == 0;

      bool passAllCuts[nSelTypes] = {passZGMass && passZMass && passBtagVeto,
                                                   passZMass && passBtagVeto,
				     passZGMass && passZMass                ,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 0,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 1,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 0 && passMET &&  passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passTauVeto,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 1 && passMET &&  passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passTauVeto,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 2 && passMET &&  passPTFrac && passDPhiZMET && passPTLL && passDPhiJetMET && passTauVeto,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 0 && passMET &&  passPTFracRecoil && passDPhiZMETRecoil && passPTLL && passDPhiJetMET && passTauVeto,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 1 && passMET &&  passPTFracRecoil && passDPhiZMETRecoil && passPTLL && passDPhiJetMET && passTauVeto,
				     passZGMass && passZMass && passBtagVeto && idJet.size() == 2 && passMET &&  passPTFracRecoil && passDPhiZMETRecoil && passPTLL && passDPhiJetMET && passTauVeto,
				     passPTLL
				     };

      int theCategory = infilecatv[ifile];
      // begin event weighting
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = theDataPrescale;
      if(nsel == 2) {
	mcWeight = mcWeight * ratioFactor_gjets_zll(TMath::Min(dilep.Pt(),399.999));
	if(infilecatv[ifile] == 0 && theMET.Pt() > 200) mcWeight = mcWeight * 0.2;
      }
      // trigger efficiency
      double trigEff = 1.0;
      //if(infilecatv[ifile] != 0 && (nsel == 0 || nsel == 3)) {
      //  trigEff = trigLookup.GetExpectedTriggerEfficiency(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),
      //  						  ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),
      //  						 TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]));
      //}
      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight = 1.0; if(infilecatv[ifile] != 0) puWeight = nPUScaleFactor(fhDPU, (double)eventMonteCarlo.puTrueInt);
      // lepton efficiency
      double effSF = 1.0;
      if(infilecatv[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		typeLepSel.Data(),fhDMuMediumSF,fhDMuIsoSF,fhDElMediumSF,fhDElTightSF,fhDElMediumMVASF,fhDElTightMVASF);
        }
      }

      double totalWeight = mcWeight*theLumi*puWeight*effSF*theMCPrescale*trigEff;
      if(totalWeight == 0) continue;
      // end event weighting
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[SIGSEL]) sumEventsProcess[ifile] += totalWeight;

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        }
      }

      if(infilecatv[ifile] == 0 && passAllCuts[SIGSEL] && ptPho > 130){
        for (int nt = 0; nt <(int)numtokens; nt++) {
	  if(strcmp(tokens[nt],"HLT_Photon50_R9Id90_HE10_IsoM_v*") == 0  && (*eventTrigger.triggerFired)[nt] == 1) eventsTrg[0] = eventsTrg[0] + totalWeight;
	  if(strcmp(tokens[nt],"HLT_Photon75_R9Id90_HE10_IsoM_v*") == 0  && (*eventTrigger.triggerFired)[nt] == 1) eventsTrg[1] = eventsTrg[1] + totalWeight;
	  if(strcmp(tokens[nt],"HLT_Photon90_R9Id90_HE10_IsoM_v*") == 0  && (*eventTrigger.triggerFired)[nt] == 1) eventsTrg[2] = eventsTrg[2] + totalWeight;
	  if(strcmp(tokens[nt],"HLT_Photon120_R9Id90_HE10_IsoM_v*") == 0 && (*eventTrigger.triggerFired)[nt] == 1) eventsTrg[3] = eventsTrg[3] + totalWeight;
        }
      }

      if(nsel == 5 && passAllCuts[GJETSEL] && idFake.size() == 1 && dRLepFakePhoMin > 0.3){
	int iPt = -1;
	if     (((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Pt() < 15){
          iPt = 0;
	}
	else if(((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Pt() < 20){
          iPt = 1;
	}
	else if(((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Pt() < 25){
          iPt = 2;
	}
	else if(((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Pt() < 30){
          iPt = 3;
	}
	else {
          iPt = 4;
	}
	int iEta = -1;
	if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Eta()) < 0.5){
           iEta = 0;
	}
	else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Eta()) < 1.0){
           iEta = 1;
	}
	else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Eta()) < 1.5){
           iEta = 2;
	}
	else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idFake[0]])->Eta()) < 2.0){
           iEta = 3;
	}
	else {
           iEta = 4;
	}
	if(infilecatv[ifile] == 0) {
	  if     (TMath::Abs((int)(*eventLeptons.pdgId)[idFake[0]]) == 13){
                                denFRDAM[iPt][iEta] = denFRDAM[iPt][iEta] + totalWeight;
            if(idTight[0] == 1) numFRDAM[iPt][iEta] = numFRDAM[iPt][iEta] + totalWeight;
	  }
	  else if(TMath::Abs((int)(*eventLeptons.pdgId)[idFake[0]]) == 11){
                                denFRDAE[iPt][iEta] = denFRDAE[iPt][iEta] + totalWeight;
            if(idTight[0] == 1) numFRDAE[iPt][iEta] = numFRDAE[iPt][iEta] + totalWeight;
	  }
	}
	else {
	  if     (TMath::Abs((int)(*eventLeptons.pdgId)[idFake[0]]) == 13){
                                denFRBGM[iPt][iEta] = denFRBGM[iPt][iEta] + totalWeight;
            if(idTight[0] == 1) numFRBGM[iPt][iEta] = numFRBGM[iPt][iEta] + totalWeight;
	  }
	  else if(TMath::Abs((int)(*eventLeptons.pdgId)[idFake[0]]) == 11){
                                denFRBGE[iPt][iEta] = denFRBGE[iPt][iEta] + totalWeight;
            if(idTight[0] == 1) numFRBGE[iPt][iEta] = numFRBGE[iPt][iEta] + totalWeight;
	  }
	}
      }

      if((typeSel == typePair) || (typeSel == 3 && (typePair == 1 || typePair == 2))) {
	double MVAVar = 0.0;
	if     (MVAVarType == 0) MVAVar = TMath::Max(TMath::Min(mtW,999.999),200.001);
	else if(MVAVarType == 1) MVAVar = TMath::Min((double)theMET.Pt(),249.999);
	else if(MVAVarType == 2) MVAVar = TMath::Min(mtW,999.999);
	else if(MVAVarType == 3) MVAVar = TMath::Min((double)theMET.Pt(),199.999);
	else if(MVAVarType == 4) MVAVar = TMath::Min(ptFrac[0],0.999);
	else {assert(0); return;}
	
	double theTemplateWeight = totalWeight;
	if     (theCategory == 0) theTemplateWeight = theTemplateWeight;
	else if(theCategory == 3) theTemplateWeight = -1.0 * theTemplateWeight;
	else                      theTemplateWeight = 0.0;
	if(theTemplateWeight != 0){
	  if(passAllCuts[ZHSEL0]   ) {histo_Zjets0   ->Fill(MVAVar,theTemplateWeight);}
	  if(passAllCuts[ZHSEL1]   ) {histo_Zjets1   ->Fill(MVAVar,theTemplateWeight);}
	  if(passAllCuts[ZHSEL2]   ) {histo_Zjets2   ->Fill(MVAVar,theTemplateWeight);}
	  if(passAllCuts[ZHRECSEL0]) {histo_ZjetsRec0->Fill(MVAVar,theTemplateWeight);}
	  if(passAllCuts[ZHRECSEL1]) {histo_ZjetsRec1->Fill(MVAVar,theTemplateWeight);}
	  if(passAllCuts[ZHRECSEL2]) {histo_ZjetsRec2->Fill(MVAVar,theTemplateWeight);}
        }

	for(int thePlot=0; thePlot<allPlots; thePlot++){
	  double theVar = 0.0;
	  bool makePlot = false;
	  if     (thePlot ==  0 && passAllCuts[ZSEL])      {makePlot = true;theVar = TMath::Min(dilep.M(),199.999);}
	  else if(thePlot ==  1 && passAllCuts[ZSEL])      {makePlot = true;theVar = TMath::Min(dilepPho.M(),199.999);}
	  else if(thePlot ==  2 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min((double)theMET.Pt(),199.999);}
	  else if(thePlot ==  3 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(dilep.Pt(),499.999);}
	  else if(thePlot ==  4 && passAllCuts[SIGSEL])    {makePlot = true;theVar = (double)theMET.Phi();}
	  else if(thePlot ==  5 && passAllCuts[SIGSEL])    {makePlot = true;theVar = (double)theMET.Px();}
	  else if(thePlot ==  6 && passAllCuts[SIGSEL])    {makePlot = true;theVar = (double)theMET.Py();}
	  else if(thePlot ==  7 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	  else if(thePlot ==  8 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(ptFrac[0],0.999);}
	  else if(thePlot ==  9 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(ptFrac[1],0.999);}
	  else if(thePlot == 10 && passAllCuts[SIGSEL])    {makePlot = true;theVar = dPhiDiLepMET;}
	  else if(thePlot == 11 && passAllCuts[NOBTAG])    {makePlot = true;theVar = TMath::Max(TMath::Min(bDiscrMax,0.999),0.001);}
	  else if(thePlot == 12 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(mtW,599.999);}
	  else if(thePlot == 13 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(ptPho,199.999);}
	  else if(thePlot == 14 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(etaPho,2.499);}
	  else if(thePlot == 15 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.Eta()),2.499);}
	  else if(thePlot == 16 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(isoPho,0.999);}
	  else if(thePlot == 17 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(dRLepPhoMin,1.999);}
	  else if(thePlot == 18 && passAllCuts[SIGSEL])    {makePlot = true;theVar = TMath::Min(sieiePho,0.04999);}
	  else if(thePlot == 19 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = TMath::Min((double)theMET.Pt(),199.999);}
	  else if(thePlot == 20 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = TMath::Min((double)theMET.Pt(),199.999);}
	  else if(thePlot == 21 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = dPhiDiLepMET;}
	  else if(thePlot == 22 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = dPhiDiLepMET;}
	  else if(thePlot == 23 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = TMath::Min(ptFrac[0],0.999);}
	  else if(thePlot == 24 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = TMath::Min(ptFrac[0],0.999);}
	  else if(thePlot == 25 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = TMath::Min(dilep.Pt(),499.999);}
	  else if(thePlot == 26 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = TMath::Min(dilep.Pt(),499.999);}
	  else if(thePlot == 27 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = TMath::Min(mtW,599.999);}
	  else if(thePlot == 28 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = TMath::Min(mtW,599.999);}
	  else if(thePlot == 29 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = (double)theMET.Phi();}
	  else if(thePlot == 30 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = (double)theMET.Phi();}
	  else if(thePlot == 31 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = TMath::Abs((double)dilep.Phi());}
	  else if(thePlot == 32 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = TMath::Abs((double)dilep.Phi());}
	  else if(thePlot == 33 && passAllCuts[ZHSEL0])    {makePlot = true;theVar = TMath::Min(mtW,599.999);}
	  else if(thePlot == 34 && passAllCuts[ZHSEL1])    {makePlot = true;theVar = TMath::Min(mtW,599.999);}
	  else if(thePlot == 35 && passAllCuts[ZHPRESEL0]) {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	  else if(thePlot == 36 && passAllCuts[ZHPRESEL1]) {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	  else if(thePlot == 37 && passAllCuts[GJETSEL])   {makePlot = true;theVar = TMath::Min((double)idGJet.size(),4.499);}
	  else if(thePlot == 38 && passAllCuts[GJETSEL])   {makePlot = true;theVar = TMath::Min((double)idFake.size(),4.499);}
	  else if(thePlot == 39 && passAllCuts[GJETSEL])   {makePlot = true;theVar = TMath::Min((double)idLep.size(),4.499);}

	  if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
	}
      }
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);

  } // end of chain

  if(eventsTrg[0] > 0 && eventsTrg[1] > 0 && eventsTrg[2] > 0 && eventsTrg[3] > 0){
    printf("gamma130 yields: %f %f %f %f --> %f %f %f\n",eventsTrg[0],eventsTrg[1],eventsTrg[2],eventsTrg[3],eventsTrg[3]/eventsTrg[0],eventsTrg[3]/eventsTrg[1],eventsTrg[3]/eventsTrg[2]);
  }

  if(nsel == 5){
    double sumTot[2] = {0.,0.};
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
	sumTot[0] = sumTot[0] + denFRDAM[iPt][iEta];
	sumTot[1] = sumTot[1] + denFRBGM[iPt][iEta];
	printf("M(%d,%d): (%6.1f-%6.1f)/(%6.1f-%6.1f)=%4.3f | ",iPt,iEta,numFRDAM[iPt][iEta],numFRBGM[iPt][iEta] , denFRDAM[iPt][iEta],denFRBGM[iPt][iEta],
                                                                	(numFRDAM[iPt][iEta]-numFRBGM[iPt][iEta])/(denFRDAM[iPt][iEta]-denFRBGM[iPt][iEta]));
	if(iPt==4) printf("\n");
      }
    }
    printf("sumTotM(da/bg) = %f / %f = %f\n",sumTot[0],sumTot[1],sumTot[0]/sumTot[1]);

    sumTot[0] = 0.; sumTot[1] = 0.;
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
	sumTot[0] = sumTot[0] + denFRDAE[iPt][iEta];
	sumTot[1] = sumTot[1] + denFRBGE[iPt][iEta];
	printf("E(%d,%d): (%6.1f-%6.1f)/(%6.1f-%6.1f)=%4.3f | ",iPt,iEta,numFRDAE[iPt][iEta],numFRBGE[iPt][iEta] , denFRDAE[iPt][iEta],denFRBGE[iPt][iEta],
                                                                	(numFRDAE[iPt][iEta]-numFRBGE[iPt][iEta])/(denFRDAE[iPt][iEta]-denFRBGE[iPt][iEta]));
	if(iPt==4) printf("\n");
      }
    }
    printf("sumTotE(da/bg) = %f / %f = %f\n",sumTot[0],sumTot[1],sumTot[0]/sumTot[1]);


    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
	printf("(%d,%d): (%6.1f)/(%6.1f)=%4.3f | ",iPt,iEta,numFRDAM[iPt][iEta] , denFRDAM[iPt][iEta],
  							   (numFRDAM[iPt][iEta])/(denFRDAM[iPt][iEta]));
	if(iPt==4) printf("\n");
      }
    }
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
	printf("(%d,%d): (%6.1f)/(%6.1f)=%4.3f | ",iPt,iEta,numFRDAE[iPt][iEta] , denFRDAE[iPt][iEta],
  							   (numFRDAE[iPt][iEta])/(denFRDAE[iPt][iEta]));
	if(iPt==4) printf("\n");
      }
    }

    printf("double fake_rate_m[%d][%d] = {\n",5,5);
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
	printf("%4.3f",TMath::Abs((numFRDAM[iPt][iEta]-numFRBGM[iPt][iEta])/(denFRDAM[iPt][iEta]-denFRBGM[iPt][iEta])));
	if(iPt!=4||iEta!=4) printf(",");
	if(iPt==4) printf("\n");
      }
    }
    printf("};\n");

    printf("double fake_rate_e[%d][%d] = {\n",5,5);
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
	printf("%4.3f",TMath::Abs((numFRDAE[iPt][iEta]-numFRBGE[iPt][iEta])/(denFRDAE[iPt][iEta]-denFRBGE[iPt][iEta])));
	if(iPt!=4||iEta!=4) printf(",");
	if(iPt==4) printf("\n");
      }
    }
    printf("};\n");

  }

  for(int i=1; i<=histo_Zjets0->GetNbinsX(); i++){
    histo_Zjets0   ->SetBinContent(i,TMath::Max(histo_Zjets0   ->GetBinContent(i),0.000001));
    histo_Zjets1   ->SetBinContent(i,TMath::Max(histo_Zjets1   ->GetBinContent(i),0.000001));
    histo_Zjets2   ->SetBinContent(i,TMath::Max(histo_Zjets2   ->GetBinContent(i),0.000001));
    histo_ZjetsRec0->SetBinContent(i,TMath::Max(histo_ZjetsRec0->GetBinContent(i),0.000001));
    histo_ZjetsRec1->SetBinContent(i,TMath::Max(histo_ZjetsRec1->GetBinContent(i),0.000001));
    histo_ZjetsRec2->SetBinContent(i,TMath::Max(histo_ZjetsRec2->GetBinContent(i),0.000001));
  }

  char output[200];
  sprintf(output,"zjets_13TeV_25ns.root");	  
  TFile* outFilePlotsZjets = new TFile(output,"recreate");
  outFilePlotsZjets->cd();
  histo_Zjets0->Write();
  histo_Zjets1->Write();
  histo_Zjets2->Write();
  histo_ZjetsRec0->Write();
  histo_ZjetsRec1->Write();
  histo_ZjetsRec2->Write();
  outFilePlotsZjets->Close();

  double sumEvents = 0;
  for(int np=1; np<histBins; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    sprintf(output,"histozg_nice_%d_%d.root",nsel,thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }

  printf("                    em                      mm                      ee                      ll\n");
  printf("-----------------------------------------------------------------------------------------------------------\n");
  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsType[4] = {0,0,0,0}; double sumEventsTypeE[4] = {0,0,0,0};
    for(int np=0; np<histBins; np++) {       
       bgdDecay[ns+nSelTypes*3][np] = bgdDecay[ns+nSelTypes*1][np] + bgdDecay[ns+nSelTypes*2][np];
       weiDecay[ns+nSelTypes*3][np] = weiDecay[ns+nSelTypes*1][np] + weiDecay[ns+nSelTypes*2][np];
       printf("(%6s): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
       processName[np].Data(),bgdDecay[ns+nSelTypes*0][np],sqrt(weiDecay[ns+nSelTypes*0][np]),bgdDecay[ns+nSelTypes*1][np],sqrt(weiDecay[ns+nSelTypes*1][np]),
                              bgdDecay[ns+nSelTypes*2][np],sqrt(weiDecay[ns+nSelTypes*2][np]),bgdDecay[ns+nSelTypes*3][np],sqrt(weiDecay[ns+nSelTypes*3][np]));
       if(np!=0 && np!=6){
         sumEventsType[0] = sumEventsType[0] + bgdDecay[ns+nSelTypes*0][np]; sumEventsTypeE[0] = sumEventsTypeE[0] + weiDecay[ns+nSelTypes*0][np];
         sumEventsType[1] = sumEventsType[1] + bgdDecay[ns+nSelTypes*1][np]; sumEventsTypeE[1] = sumEventsTypeE[1] + weiDecay[ns+nSelTypes*1][np];
         sumEventsType[2] = sumEventsType[2] + bgdDecay[ns+nSelTypes*2][np]; sumEventsTypeE[2] = sumEventsTypeE[2] + weiDecay[ns+nSelTypes*2][np];
         sumEventsType[3] = sumEventsType[3] + bgdDecay[ns+nSelTypes*3][np]; sumEventsTypeE[3] = sumEventsTypeE[3] + weiDecay[ns+nSelTypes*3][np];
       }
    }
    printf("(...bkg): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
	   sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]));
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
}
