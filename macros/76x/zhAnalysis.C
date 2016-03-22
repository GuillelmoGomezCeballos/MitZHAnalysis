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

#include "MitAnalysisRunII/macros/76x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

bool useZjetsTemplate = true;
bool usePureMC = true; 
bool useEMFromData = true;
double mcPrescale = 1.0;
enum selType                     {ZSEL=0,  SIGSEL,   WWSEL,   WWLOOSESEL,   BTAGSEL,   WZSEL,   PRESEL, nSelTypes};
TString selTypeName[nSelTypes]= {"ZSEL",  "SIGSEL", "WWSEL", "WWLOOSESEL", "BTAGSEL", "WZSEL", "PRESEL"};
enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};
const TString typeLepSel = "medium";

void zhAnalysis(
 int mH = 125,
 unsigned int nJetsType = 0,
 bool useGGZH = false,
 bool isBlinded = false,
 Int_t typeSel = 3
 ){

  Int_t period = 1;
  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/met_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/met_";
  Double_t lumi = 2.318;
  TString processTag = "";

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;
  vector<Int_t> infilecatv;  

  TString puPath = "";
  TString zjetsTemplatesPath = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";

  infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));						   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(1);

  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(1);
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

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(2);

  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);

  //infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                         infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);

  if(mH>100){
  processTag = Form("mh%d",mH);
  infilenamev.push_back(Form("%sZH_ZToMM_HToInvisible_M%d_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data(),mH)); infilecatv.push_back(6);
  infilenamev.push_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data(),mH)); infilecatv.push_back(6);
  }
  else {assert(0);}
  
  if(mH==125 && useGGZH == true){
  infilenamev.push_back(Form("%sggZH_HToInv_ZToLL_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); infilecatv.push_back(7);
  }
  }

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}
  
  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Medium_ele"));
  TH2D *fhDElTightSF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Tight_ele"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF ->SetDirectory(0);
  delete fElSF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Medium_mu"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Iso_mu"));
  assert(fhDMuMediumSF);
  assert(fhDMuIsoSF);
  fhDMuMediumSF->SetDirectory(0);
  fhDMuIsoSF->SetDirectory(0);
  delete fMuSF;

  TString ECMsb  = "13TeV2015";
  const int MVAVarType = 0; const int nBinMVA = 6; Float_t xbins[nBinMVA+1] = {200, 250, 300, 400, 600, 800, 1000}; TString addChan = "";
  //const int MVAVarType = 1; const int nBinMVA = 6; Float_t xbins[nBinMVA+1] = {100, 125, 150, 175, 200, 225, 250}; TString addChan = "1";
  //const int MVAVarType = 2; const int nBinMVA = 17; Float_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 600, 800, 1000}; TString addChan = "2";
  //const int MVAVarType = 3; const int nBinMVA = 23; Float_t xbins[nBinMVA+1] = {45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150, 160, 170, 180, 190, 200}; TString addChan = "3";
  //const int MVAVarType = 4; const int nBinMVA = 20; Float_t xbins[nBinMVA+1] =  {0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00}; TString addChan = "4";
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  if     (MVAVarType == 0) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt100_mt.root";
  else if(MVAVarType == 1) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt100_met.root";
  else if(MVAVarType == 2) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt45_mt.root";
  else if(MVAVarType == 3) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt45_met.root";
  else if(MVAVarType == 4) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt45_ptfraclt1_ptfrac.root";
  else {printf("PROBLEM with MVAVarType\n");}

  TFile *fZjetsTemplatesFile = TFile::Open(Form("%s",zjetsTemplatesPath.Data()));

  TH1D *fhDZjets;
  if     (nJetsType == 0) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
  else if(nJetsType == 1) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
  else if(nJetsType == 2) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets2"));
  assert(fhDZjets);
  fhDZjets->SetDirectory(0);

  TH1D *fhDZjetsSyst;
  if     (nJetsType == 0) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets0"));
  else if(nJetsType == 1) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets0"));
  else if(nJetsType == 2) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
  assert(fhDZjetsSyst);
  fhDZjetsSyst->SetDirectory(0);

  delete fZjetsTemplatesFile;

  const int numberCuts = 12;
  TH1D* histoZHSEL[4];
  histoZHSEL[0] = new TH1D("histoZHSEL_0", "histoZHSEL_0", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[1] = new TH1D("histoZHSEL_1", "histoZHSEL_1", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[2] = new TH1D("histoZHSEL_2", "histoZHSEL_2", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[3] = new TH1D("histoZHSEL_3", "histoZHSEL_3", numberCuts+1, -0.5, numberCuts+0.5);
  TString cutName[numberCuts+1] = {"ptl>20/20","ptll>60","MET>100","Njets=0","Z mass","3rd lepton veto","btag-veto","dPhi(Z-MET)>2.8","dPhi(l-l)<pt/2","|ptll-MET|/ptll<0.4","MT>200","dPhiJetMet>0.5","tauVeto"};

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 25;
  const int histBins = 8;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "...DY", "...WZ", "....ZZ", "...VVV", "....ZH", "..ggZH"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 400; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >=  1 && thePlot <=  1) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >=  8 && thePlot <=  9) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot >= 11 && thePlot <= 11) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >= 12 && thePlot <= 14) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot >= 15 && thePlot <= 15) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 16 && thePlot <= 16) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 17 && thePlot <= 17) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 500; xminPlot =-0.5; xmaxPlot = 499.5;}
    else if(thePlot >= 19 && thePlot <= 19) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 20 && thePlot <= 20) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   2.5;}
    else if(thePlot >= 21 && thePlot <= 21) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot =1000.0;}
    else if(thePlot >= 22 && thePlot <= 22) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot == allPlots-2)          {nBinPlot =  numberCuts+1; xminPlot =-0.5; xmaxPlot =  numberCuts+0.5;}
    TH1D* histos;
    if(thePlot != allPlots-1) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    else                      histos = new TH1D("histos", "histos", nBinMVA, xbins);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  TH1D *histo_Data     = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_ZH_hinv  = (TH1D*) histoMVA->Clone("histo_ZH_hinv"); 
  TH1D *histo_Zjets    = (TH1D*) histoMVA->Clone("histo_Zjets");	 
  TH1D *histo_VVV      = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZ       = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZ       = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_EM       = (TH1D*) histoMVA->Clone("histo_EM");	 
  TH1D *histo_ggZH_hinv= (TH1D*) histoMVA->Clone("histo_ggZH_hinv"); 
  TH1D *histo_ZH_hinvNoW  = (TH1D*) histoMVA->Clone("histo_ZH_hinv"); 
  TH1D *histo_ZjetsNoW    = (TH1D*) histoMVA->Clone("histo_Zjets");	 
  TH1D *histo_VVVNoW      = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZNoW       = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZNoW       = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_EMNoW       = (TH1D*) histoMVA->Clone("histo_EM");	 
  TH1D *histo_ggZH_hinvNoW= (TH1D*) histoMVA->Clone("histo_ggZH_hinv"); 

  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff_m");sprintf(momMName,"CMS_scale_m");
  sprintf(effEName,"CMS_eff_e");sprintf(momEName,"CMS_scale_e");
  if     (typeSel == 3 && nJetsType == 0) {sprintf(finalStateName,"ll0j");}
  else if(typeSel == 1 && nJetsType == 0) {sprintf(finalStateName,"mm0j");}
  else if(typeSel == 2 && nJetsType == 0) {sprintf(finalStateName,"ee0j");}
  else if(typeSel == 3 && nJetsType == 1) {sprintf(finalStateName,"ll1j");}
  else if(typeSel == 1 && nJetsType == 1) {sprintf(finalStateName,"mm1j");}
  else if(typeSel == 2 && nJetsType == 1) {sprintf(finalStateName,"ee1j");}
  else if(typeSel == 3 && nJetsType == 2) {sprintf(finalStateName,"ll2j");}
  else if(typeSel == 1 && nJetsType == 2) {sprintf(finalStateName,"mm2j");}
  else if(typeSel == 2 && nJetsType == 2) {sprintf(finalStateName,"ee2j");}
  else {printf("Wrong lSel/nJetsType: %d/%d\n",typeSel,nJetsType); assert(0); return;}

  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingUp      = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingDown    = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingDown->Sumw2();
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

  TH1D* histo_ZH_hinv_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ggZH_hinv_CMS_QCDScaleBounding[6];
  for(int nb=0; nb<6; nb++){
    histo_ZH_hinv_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_ZH_hinv_QCDScale_f%d",nb), Form("histo_ZH_hinv_QCDScale_f%d",nb),nBinMVA, xbins); histo_ZH_hinv_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]       = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding[nb]	     = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),      Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_WZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding[nb]	     = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ggZH_hinv_CMS_QCDScaleBounding[nb] = new TH1D(Form("histo_ggZH_hinv_QCDScale_f%d",nb), Form("histo_ggZH_hinv_QCDScale_f%d",nb),nBinMVA, xbins); histo_ggZH_hinv_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_ZH_hinv_CMS_PDFBounding[102];
  TH1D* histo_VVV_CMS_PDFBounding[102];
  TH1D* histo_WZ_CMS_PDFBounding[102];
  TH1D* histo_ZZ_CMS_PDFBounding[102];
  TH1D* histo_ggZH_hinv_CMS_PDFBounding[102];
  for(int nb=0; nb<102; nb++){
    histo_ZH_hinv_CMS_PDFBounding[nb]   = new TH1D(Form("histo_ZH_hinv_PDF_f%d",nb), Form("histo_ZH_hinv_PDF_f%d",nb),nBinMVA, xbins); histo_ZH_hinv_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]       = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WZ_CMS_PDFBounding[nb]        = new TH1D(Form("histo_WZ_PDF_f%d",nb),      Form("histo_WZ_PDF_f%d",nb),     nBinMVA, xbins); histo_WZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ZZ_CMS_PDFBounding[nb]        = new TH1D(Form("histo_ZZ_PDF_f%d",nb),      Form("histo_ZZ_PDF_f%d",nb),     nBinMVA, xbins); histo_ZZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ggZH_hinv_CMS_PDFBounding[nb] = new TH1D(Form("histo_ggZH_hinv_PDF_f%d",nb), Form("histo_ggZH_hinv_PDF_f%d",nb),nBinMVA, xbins); histo_ggZH_hinv_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nBinMVA];
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
  for(int nb=0; nb<nBinMVA; nb++){
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]        = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dUp"       ,finalStateName,ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb]      = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dDown"     ,finalStateName,ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nb] ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]       = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]  ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]     = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	        = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"          ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	        = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"        ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	        = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"            ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	        = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"          ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	        = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"            ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	        = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"          ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	        = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"            ,finalStateName,ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingBinUp[nb]	      ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]	        = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"          ,finalStateName,ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingBinDown[nb]	  ->Sumw2();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[nb]    = new TH1D(Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[nb]  = new TH1D(Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[nb] ->Sumw2();
  }

  TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingUp     = new TH1D( Form("histo_ZH_hinv_%sUp",effMName)  , Form("histo_ZH_hinv_%sUp",effMName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingDown   = new TH1D( Form("histo_ZH_hinv_%sDown",effMName), Form("histo_ZH_hinv_%sDown",effMName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	   = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	   = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingUp    	   = new TH1D( Form("histo_WZ_%sUp",effMName)  , Form("histo_WZ_%sUp",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingDown  	   = new TH1D( Form("histo_WZ_%sDown",effMName), Form("histo_WZ_%sDown",effMName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp    	   = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown  	   = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingUp   = new TH1D( Form("histo_ggZH_hinv_%sUp",effMName)  , Form("histo_ggZH_hinv_%sUp",effMName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingDown = new TH1D( Form("histo_ggZH_hinv_%sDown",effMName), Form("histo_ggZH_hinv_%sDown",effMName), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_ZH_hinv_%sAvg",effMName)  , Form("histo_ZH_hinv_%sAvg",effMName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	   = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg    	   = new TH1D( Form("histo_WZ_%sAvg",effMName)  , Form("histo_WZ_%sAvg",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    	   = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg  = new TH1D( Form("histo_ggZH_hinv_%sAvg",effMName)  , Form("histo_ggZH_hinv_%sAvg",effMName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingUp     = new TH1D( Form("histo_ZH_hinv_%sUp",effEName)  , Form("histo_ZH_hinv_%sUp",effEName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingDown   = new TH1D( Form("histo_ZH_hinv_%sDown",effEName), Form("histo_ZH_hinv_%sDown",effEName), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	   = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	   = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingUp    	   = new TH1D( Form("histo_WZ_%sUp",effEName)  , Form("histo_WZ_%sUp",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingDown  	   = new TH1D( Form("histo_WZ_%sDown",effEName), Form("histo_WZ_%sDown",effEName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp    	   = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown  	   = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingUp   = new TH1D( Form("histo_ggZH_hinv_%sUp",effEName)  , Form("histo_ggZH_hinv_%sUp",effEName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingDown = new TH1D( Form("histo_ggZH_hinv_%sDown",effEName), Form("histo_ggZH_hinv_%sDown",effEName), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_ZH_hinv_%sAvg",effEName)  , Form("histo_ZH_hinv_%sAvg",effEName)  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	   = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg    	   = new TH1D( Form("histo_WZ_%sAvg",effEName)  , Form("histo_WZ_%sAvg",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    	   = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg  = new TH1D( Form("histo_ggZH_hinv_%sAvg",effEName)  , Form("histo_ggZH_hinv_%sAvg",effEName)  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVAMETBoundingUp      = new TH1D( Form("histo_ZH_hinv_CMS_scale_metUp")  , Form("histo_ZH_hinv_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAMETBoundingDown    = new TH1D( Form("histo_ZH_hinv_CMS_scale_metDown"), Form("histo_ZH_hinv_CMS_scale_metDown"), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown"), nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAMETBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_scale_metUp")  , Form("histo_ggZH_hinv_CMS_scale_metUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAMETBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_scale_metDown"), Form("histo_ggZH_hinv_CMS_scale_metDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingUp      = new TH1D( Form("histo_ZH_hinv_CMS_scale_jUp")  , Form("histo_ZH_hinv_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_MVAJESBoundingDown    = new TH1D( Form("histo_ZH_hinv_CMS_scale_jDown"), Form("histo_ZH_hinv_CMS_scale_jDown"), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp      	= new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown    	= new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_WZ_CMS_scale_jUp")  , Form("histo_WZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_WZ_CMS_scale_jDown"), Form("histo_WZ_CMS_scale_jDown"), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingUp       	= new TH1D( Form("histo_ZZ_CMS_scale_jUp")  , Form("histo_ZZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown     	= new TH1D( Form("histo_ZZ_CMS_scale_jDown"), Form("histo_ZZ_CMS_scale_jDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jUp")  , Form("histo_ggZH_hinv_CMS_scale_jUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jDown"), Form("histo_ggZH_hinv_CMS_scale_jDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_ZH_hinv_CMS_PUBoundingUp          = new TH1D( Form("histo_ZH_hinv_CMS_puUp")  , Form("histo_ZH_hinv_CMS_puUp")  , nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ZH_hinv_CMS_PUBoundingDown        = new TH1D( Form("histo_ZH_hinv_CMS_puDown"), Form("histo_ZH_hinv_CMS_puDown"), nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingUp           	= new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown         	= new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingUp            	= new TH1D( Form("histo_WZ_CMS_puUp")  , Form("histo_WZ_CMS_puUp")  , nBinMVA, xbins); histo_WZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingDown          	= new TH1D( Form("histo_WZ_CMS_puDown"), Form("histo_WZ_CMS_puDown"), nBinMVA, xbins); histo_WZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingUp            	= new TH1D( Form("histo_ZZ_CMS_puUp")  , Form("histo_ZZ_CMS_puUp")  , nBinMVA, xbins); histo_ZZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingDown  	        = new TH1D( Form("histo_ZZ_CMS_puDown"), Form("histo_ZZ_CMS_puDown"), nBinMVA, xbins); histo_ZZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_PUBoundingUp        = new TH1D( Form("histo_ggZH_hinv_CMS_puUp")  , Form("histo_ggZH_hinv_CMS_puUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ggZH_hinv_CMS_PUBoundingDown      = new TH1D( Form("histo_ggZH_hinv_CMS_puDown"), Form("histo_ggZH_hinv_CMS_puDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_PUBoundingDown->Sumw2();

  TH1D* histo_WZ_CMS_EWKCorrUp    	        = new TH1D( Form("histo_WZ_EWKCorrUp")  , Form("histo_WZ_EWKCorrUp")  , nBinMVA, xbins); histo_WZ_CMS_EWKCorrUp  ->Sumw2();
  TH1D* histo_WZ_CMS_EWKCorrDown                = new TH1D( Form("histo_WZ_EWKCorrDown"), Form("histo_WZ_EWKCorrDown"), nBinMVA, xbins); histo_WZ_CMS_EWKCorrDown->Sumw2();
  TH1D* histo_ZZ_CMS_EWKCorrUp                  = new TH1D( Form("histo_ZZ_EWKCorrUp")  , Form("histo_ZZ_EWKCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_EWKCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_EWKCorrDown                = new TH1D( Form("histo_ZZ_EWKCorrDown"), Form("histo_ZZ_EWKCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_EWKCorrDown->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrUp                   = new TH1D( Form("histo_ZZ_ggCorrUp")  , Form("histo_ZZ_ggCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_ggCorrUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_ggCorrDown                 = new TH1D( Form("histo_ZZ_ggCorrDown"), Form("histo_ZZ_ggCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_ggCorrDown->Sumw2();
  TH1D* histo_Zjets_CMS_ZjetsSystUp    	        = new TH1D( Form("histo_Zjets_ZjetsSystUp")  , Form("histo_Zjets_ZjetsSystUp")  , nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_ZjetsSystDown           = new TH1D( Form("histo_Zjets_ZjetsSystDown"), Form("histo_Zjets_ZjetsSystDown"), nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystDown->Sumw2();

  double bgdDecay[nSelTypes*4][histBins],weiDecay[nSelTypes*4][histBins];
  for(unsigned int i=0; i<nSelTypes*4; i++) {
    for(int j=0; j<histBins; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }

  unsigned int numberOfLeptons = 2;

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");
    TTree *the_PDF_tree   = (TTree*)the_input_file.FindObjectAny("pdfReweight");

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
	  printf("PDF(%d): %s\n",i,weightDef);
        }
      }
      else {
        printf("PDFTree not available\n");
      }
      //return;
    }

    histoZHSEL[0]->Scale(0.0);
    histoZHSEL[1]->Scale(0.0);
    histoZHSEL[2]->Scale(0.0);
    histoZHSEL[3]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_input_tree->GetEntry(i);

      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[4] = {kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      for (int nt = 0; nt <(int)numtokens; nt++) {
        if((*eventTrigger.triggerFired)[nt] == 0) continue;
        if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*") 	    == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")	    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu27_v*") 				    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu20_v*") 				    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoTkMu20_v*")				    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele23_WPLoose_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele22_eta2p1_WP75_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WPLoose_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WP85_Gsf_v*")  			    == 0)
           ) passFilter[1] = kTRUE;
      }

      //if(infilecatv[ifile] != 0) passFilter[1] = kTRUE; // do not apply trigger filters to MC
      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()>=numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      if(idLep.size()==goodIsTight) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
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

      vector<int> idJet,idJetUp,idJetDown;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 15) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1 && ((TLorentzVector*)(*eventJets.p4)[nj])->Pt()> 30) dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 15 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()      > 30) {idJet.push_back(nj);}
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*1.05 > 30) idJetUp.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*0.95 > 30) idJetDown.push_back(nj);
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

      int typePair = 0;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typePair = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typePair = 2;

      double dPhiDiLepMET = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double ptFrac = TMath::Abs(dilep.Pt()-((TLorentzVector*)(*eventMet.p4)[0])->Pt())/dilep.Pt();
      double deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));

      bool passZMass     = TMath::Abs(dilep.M()-91.1876) < 15.0;
      bool passNjets     = idJet.size() == nJetsType;

      double metMIN = 100; double mtMIN = 200;
      if     (MVAVarType == 2 || MVAVarType == 3 || MVAVarType == 4) {metMIN = 45; mtMIN = 0;}
      bool passMET = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > metMIN;
      bool passMT = mtW > mtMIN;
      if(infilecatv[ifile] == 0 && isBlinded) passMET = passMET && ((TLorentzVector*)(*eventMet.p4)[0])->Pt() < 120;

      bool passPTFrac    = ptFrac < 0.4;
      if(MVAVarType == 4) passPTFrac = ptFrac < 1.0;
      bool passDPhiZMET  = dPhiDiLepMET > 2.8;
      bool passBtagVeto  = bDiscrMax < 0.800 && idSoft.size() == 0;
      bool passPTLL      = dilep.Pt() > 60;
      bool passMETMin    = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 45.;
      bool pass3rdLVeto  = idLep.size() == numberOfLeptons && TMath::Abs(signQ) == 0;
      bool passDelphiLL  = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]])) < TMath::Pi()/2.;

      bool passZMassLarge = TMath::Abs(dilep.M()-91.1876) < 30.0;
      bool passZMassSB    = (dilep.M() > 110.0 && dilep.M() < 200.0);

      bool passDPhiJetMET = dPhiJetMET == -1 || dPhiJetMET >= 0.5;
      bool passTauVeto    = numberGoodTaus == 0;

      bool passNMinusOne[11] = {           passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                 passMT &&              passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                 passMT && passZMass &&              passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                 passMT && passZMass && passNjets &&            passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                 passMT && passZMass && passNjets && passMET &&               passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                 passMT && passZMass && passNjets && passMET && passPTFrac                 && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                 passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET                 && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                 passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto             &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
			         passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto                 && passDPhiJetMET && passTauVeto,
				 passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL &&                   passTauVeto,
				 passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET};

      bool passAllCuts[nSelTypes] = {                   passZMass && passNjets                                                                                 &&  pass3rdLVeto                ,
                                                        passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                     passZMassLarge && !passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                     passZMassSB    && !passZMass && passNjets && passMETMin                                      && !passBtagVeto             &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                                        passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET && !passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
                                                        passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL && !pass3rdLVeto,
							passZMass && passNjets && passMETMin                      && passDPhiZMET                  && passPTLL &&  pass3rdLVeto && ((TLorentzVector*)(*eventMet.p4)[0])->Pt() < 100.};

     bool passEvolFilter[numberCuts] = {passPTLL,passMET,passNjets,passZMass,pass3rdLVeto,passBtagVeto,passDPhiZMET,passDelphiLL,passPTFrac,passMT,passDPhiJetMET,passTauVeto};

     int sumEvol = 0;
     bool totalSel = kTRUE;
     for(int isel=0; isel<numberCuts; isel++) {
       totalSel = totalSel && passEvolFilter[isel];
       if(totalSel == kTRUE) sumEvol++;
     }

     double mtWSyst[2] = {TMath::Sqrt(2.0*dilep.Pt()*(double)(*eventMet.ptJESUP)[0]  *(1.0 - cos(deltaPhiDileptonMet))),
                          TMath::Sqrt(2.0*dilep.Pt()*(double)(*eventMet.ptJESDOWN)[0]*(1.0 - cos(deltaPhiDileptonMet)))};
     bool passSystCuts[nSystTypes] = {
          passZMass && idJetUp.size() == nJetsType  && passMET && passMT                            && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && idJetDown.size()== nJetsType && passMET && passMT                            && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && passNjets && (double)(*eventMet.ptJESUP)[0]   > metMIN && mtWSyst[0] > mtMIN && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && passNjets && (double)(*eventMet.ptJESDOWN)[0] > metMIN && mtWSyst[1] > mtMIN && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto
     };

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

      // trigger efficiency
      double trigEff = 1.0;
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
		period,typeLepSel.Data(),fhDMuMediumSF,fhDMuIsoSF,fhDElMediumSF,fhDElTightSF);
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

      if(totalWeight == 0) continue;

      double the_rho = 0.0; if(the_rhoP4.P() > 0) the_rho = the_rhoP4.Pt()/the_rhoP4.P();
      double theZZCorr[2] {1,1};
      if(theCategory == 4 && infilenamev[ifile].Contains("GluGlu") == kFALSE) {
        float GENmZZ = 0.0;
	if(zzBoson.size() >= 2) GENmZZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[1])) ) ).M();
	theZZCorr[0] = weightEWKCorr(bosonPtMin,1);
	theZZCorr[1] = kfactor_qqZZ_qcd_M(GENmZZ);
        totalWeight = totalWeight * (theZZCorr[0]*theZZCorr[1]);
      }
      // end event weighting

      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[SIGSEL]) sumEventsProcess[ifile] += totalWeight;

      for(int nl=0; nl <=sumEvol; nl++) histo[allPlots-2][theCategory]->Fill((double)nl,totalWeight);
      for(int nl=0; nl <=sumEvol; nl++) histoZHSEL[typePair ]         ->Fill((double)nl,totalWeight);
      if(typePair == 1 || typePair == 2)
      for(int nl=0; nl <=sumEvol; nl++) histoZHSEL[3]                 ->Fill((double)nl,totalWeight);

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        }
      }

      if((typeSel == typePair) || (typeSel == 3 && (typePair == 1 || typePair == 2))) {
	for(int thePlot=0; thePlot<allPlots-2; thePlot++){
	  double theVar = 0.0;
	  bool makePlot = false;
	  if     (thePlot ==  0 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(mtW,999.999);}
	  else if(thePlot ==  1 && passNMinusOne[1])   {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.M()-91.1876),99.999);}
	  else if(thePlot ==  2 && passNMinusOne[2])   {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	  else if(thePlot ==  3 && passNMinusOne[3])   {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	  else if(thePlot ==  4 && passNMinusOne[4])   {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
	  else if(thePlot ==  5 && passNMinusOne[5])   {makePlot = true;theVar = dPhiDiLepMET;}
	  else if(thePlot ==  6 && passNMinusOne[6])   {makePlot = true;theVar = TMath::Max(TMath::Min(bDiscrMax,0.999),0.001);}
	  else if(thePlot ==  7 && passNMinusOne[7])   {makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
	  else if(thePlot ==  8 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
	  else if(thePlot ==  9 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
	  else if(thePlot == 10 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min((double)eventEvent.rho,39.999);}
	  else if(thePlot == 11 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	  else if(thePlot == 12 && passNMinusOne[9])   {makePlot = true;theVar = dPhiJetMET;}
	  else if(thePlot == 13 && passAllCuts[SIGSEL]){makePlot = true;theVar = dPhiLepMETMin;}
	  else if(thePlot == 14 && passNMinusOne[8])   {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]));}
	  else if(thePlot == 15 && passAllCuts[PRESEL]){makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	  else if(thePlot == 16 && passAllCuts[PRESEL]){makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
	  else if(thePlot == 17 && passAllCuts[PRESEL]){makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
	  else if(thePlot == 18 && passAllCuts[SIGSEL]){makePlot = true;theVar = (double)(numberGoodGenLep[0]+10*numberGoodGenLep[1]+100*numberGoodGenLep[2]);}
	  else if(thePlot == 19 && passNMinusOne[10])  {makePlot = true;theVar = TMath::Min((double)numberGoodTaus,3.499);}
	  else if(thePlot == 20 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.Eta()),2.499);}
	  else if(thePlot == 21 && passNMinusOne[0])   {makePlot = true;theVar = TMath::Min(mtW,999.999);}
	  else if(thePlot == 22 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(the_rho,0.999);}

	  if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
	}

	double MVAVar = 0.0;
	if     (MVAVarType == 0) MVAVar = TMath::Max(TMath::Min(mtW,999.999),200.001);
	else if(MVAVarType == 1) MVAVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),249.999);
	else if(MVAVarType == 2) MVAVar = TMath::Min(mtW,999.999);
	else if(MVAVarType == 3) MVAVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);
	else if(MVAVarType == 4) MVAVar = TMath::Min(ptFrac,0.999);
	else {assert(0); return;}

        if     (theCategory == 0){
	  if(passAllCuts[SIGSEL]) histo_Data->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 1){
	  if(passAllCuts[SIGSEL]) histo_EM   ->Fill(MVAVar,totalWeight);
	  if(passAllCuts[SIGSEL]) histo_EMNoW->Fill(MVAVar,1.);
        }
        else if(theCategory == 2){
	  if(passAllCuts[SIGSEL]) histo_Zjets   ->Fill(MVAVar,totalWeight);
	  if(passAllCuts[SIGSEL]) histo_ZjetsNoW->Fill(MVAVar,1.);
        }
        else if(theCategory == 3){
	  if(passAllCuts[SIGSEL]) {
	     histo_WZ              ->Fill(MVAVar,totalWeight);
	     histo_WZNoW           ->Fill(MVAVar,1.);
	     histo_WZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*1.03);
	     histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
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
          if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
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
	     histo_ZZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_ZZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_ZZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_ZZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_ZZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_ZZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
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
          if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 5){
	  if(passAllCuts[SIGSEL]) {
	     histo_VVV   ->Fill(MVAVar,totalWeight);
	     histo_VVVNoW->Fill(MVAVar,1.);
	     histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
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
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 6){
	  if(passAllCuts[SIGSEL]) {
	     histo_ZH_hinv   ->Fill(MVAVar,totalWeight);
	     histo_ZH_hinvNoW->Fill(MVAVar,1.);
	     histo_ZH_hinv_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_ZH_hinv_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_ZH_hinv_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_ZH_hinv_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_ZH_hinv_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_ZH_hinv_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_ZH_hinv_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_ZH_hinv_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_ZH_hinv_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZH_hinv_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_ZH_hinv_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZH_hinv_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZH_hinv_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_ZH_hinv_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 7){
	  if(passAllCuts[SIGSEL]) {
	     histo_ggZH_hinv   ->Fill(MVAVar,totalWeight);
	     histo_ggZH_hinvNoW->Fill(MVAVar,1.0);
	     histo_ggZH_hinv_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_ggZH_hinv_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_ggZH_hinv_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_ggZH_hinv_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_ggZH_hinv_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_ggZH_hinv_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
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
          if(passSystCuts[METUP])  histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_ggZH_hinv_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
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

  } // end of chain

  // "-1" to remove the Higgs contribution
  double sumEvents = 0;
  for(int np=1; np<histBins-1; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);

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
       if(np!=0 && np!=6 && np!=7){
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

  // This is the alternative method
  double kEoverM  = sqrt(bgdDecay[ZSEL+nSelTypes*2][0]/bgdDecay[ZSEL+nSelTypes*1][0]);
  double NemFact[3] = {(1.0/kEoverM)*0.5, (kEoverM)*0.5, (kEoverM+1.0/kEoverM)*0.5};

  // This is the default method
  double NemFact_FromMLLSB[3] = {NemFact[0], NemFact[1], NemFact[2]}; double NemFact_FromMLLSBE[3] = {1.0, 1.0, 1.0};
  if(bgdDecay[WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[WWLOOSESEL+nSelTypes*1][0]) > 0) {
    NemFact_FromMLLSB[0]  = (bgdDecay[WWLOOSESEL+nSelTypes*1][0])/bgdDecay[WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[0] = sqrt(1./bgdDecay[WWLOOSESEL+nSelTypes*1][0]+1./bgdDecay[WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[0];
  }
  if(bgdDecay[WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[WWLOOSESEL+nSelTypes*2][0]) > 0) {
    NemFact_FromMLLSB[1]  = (bgdDecay[WWLOOSESEL+nSelTypes*2][0])/bgdDecay[WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[1] = sqrt(1./bgdDecay[WWLOOSESEL+nSelTypes*2][0]+1./bgdDecay[WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[1];
  }
  if(bgdDecay[WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[WWLOOSESEL+nSelTypes*3][0]) > 0) {
    NemFact_FromMLLSB[2]  = (bgdDecay[WWLOOSESEL+nSelTypes*3][0])/bgdDecay[WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[2] = sqrt(1./bgdDecay[WWLOOSESEL+nSelTypes*3][0]+1./bgdDecay[WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[2];
  }

  printf("(mm) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[0],NemFact_FromMLLSB[0],NemFact_FromMLLSBE[0]);
  printf("(ee) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[1],NemFact_FromMLLSB[1],NemFact_FromMLLSBE[1]);
  printf("(ll) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[2],NemFact_FromMLLSB[2],NemFact_FromMLLSBE[2]);

  // There uncertainties: closure test, different between default and alternative method, and data statistics
  // Only the second item is used as systematic uncertainty
  double EMSystTotal[3] = {1.0,1.0,1.0}; double EMSyst[2][3] = {0.0,0.0,0.0,0.0,0.0,0.0};

  EMSyst[0][0] = bgdDecay[SIGSEL+nSelTypes*1][1]/bgdDecay[SIGSEL+nSelTypes*0][1]/NemFact_FromMLLSB[0];
  if(EMSyst[0][0] < 1.0) EMSyst[0][0] = 1/EMSyst[0][0]; EMSyst[0][0] = EMSyst[0][0] - 1.0;

  EMSyst[0][1] = bgdDecay[SIGSEL+nSelTypes*2][1]/bgdDecay[SIGSEL+nSelTypes*0][1]/NemFact_FromMLLSB[1];
  if(EMSyst[0][1] < 1.0) EMSyst[0][1] = 1/EMSyst[0][1]; EMSyst[0][1] = EMSyst[0][1] - 1.0;

  EMSyst[0][2] = bgdDecay[SIGSEL+nSelTypes*3][1]/bgdDecay[SIGSEL+nSelTypes*0][1]/NemFact_FromMLLSB[2];
  if(EMSyst[0][2] < 1.0) EMSyst[0][2] = 1/EMSyst[0][2]; EMSyst[0][2] = EMSyst[0][2] - 1.0;

  EMSyst[1][0] = NemFact_FromMLLSB[0]/NemFact[0];
  if(EMSyst[1][0] < 1.0) EMSyst[1][0] = 1/EMSyst[1][0]; EMSyst[1][0] = EMSyst[1][0] - 1.0;

  EMSyst[1][1] = NemFact_FromMLLSB[1]/NemFact[1];
  if(EMSyst[1][1] < 1.0) EMSyst[1][1] = 1/EMSyst[1][1]; EMSyst[1][1] = EMSyst[1][1] - 1.0;

  EMSyst[1][2] = NemFact_FromMLLSB[2]/NemFact[2];
  if(EMSyst[1][2] < 1.0) EMSyst[1][2] = 1/EMSyst[1][2]; EMSyst[1][2] = EMSyst[1][2] - 1.0;

  if(bgdDecay[SIGSEL][0] > 0) EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1/bgdDecay[SIGSEL][0]);
  else                        EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1.0);

  if(bgdDecay[SIGSEL][0] > 0) EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1/bgdDecay[SIGSEL][0]);
  else                        EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1.0);

  if(bgdDecay[SIGSEL][0] > 0) EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1/bgdDecay[SIGSEL][0]);
  else                        EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1.0);

  double EMbkg = bgdDecay[SIGSEL][2]+bgdDecay[SIGSEL][3]+bgdDecay[SIGSEL][4]+bgdDecay[SIGSEL][5];

  printf("(mm) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[SIGSEL+nSelTypes*1][1],sqrt(weiDecay[SIGSEL+nSelTypes*1][4]),
	 bgdDecay[SIGSEL][1]*NemFact_FromMLLSB[0]  ,sqrt(weiDecay[SIGSEL][4])*NemFact_FromMLLSB[0],bgdDecay[SIGSEL][0],EMbkg,EMSystTotal[0],EMSyst[0][0],EMSyst[1][0],sqrt(EMSystTotal[0]*EMSystTotal[0] - EMSyst[0][0]*EMSyst[0][0] - EMSyst[1][0]*EMSyst[1][0]));

  printf("(ee) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[SIGSEL+nSelTypes*2][1],sqrt(weiDecay[SIGSEL+nSelTypes*2][4]),
	 bgdDecay[SIGSEL][1]*NemFact_FromMLLSB[1]  ,sqrt(weiDecay[SIGSEL][4])*NemFact_FromMLLSB[1],bgdDecay[SIGSEL][0],EMbkg,EMSystTotal[1],EMSyst[0][1],EMSyst[1][1],sqrt(EMSystTotal[1]*EMSystTotal[1] - EMSyst[0][1]*EMSyst[0][1] - EMSyst[1][1]*EMSyst[1][1]));

  printf("(ll) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[SIGSEL+nSelTypes*3][1],sqrt(weiDecay[SIGSEL+nSelTypes*3][4]),
	 bgdDecay[SIGSEL][1]*NemFact_FromMLLSB[2]  ,sqrt(weiDecay[SIGSEL][4])*NemFact_FromMLLSB[2],bgdDecay[SIGSEL][0],EMbkg,EMSystTotal[2],EMSyst[0][2],EMSyst[1][2],sqrt(EMSystTotal[2]*EMSystTotal[2] - EMSyst[0][2]*EMSyst[0][2] - EMSyst[1][2]*EMSyst[1][2]));

  // DY background estimation
  histo_Zjets_CMS_ZjetsSystUp->Add(fhDZjets);
  if(useZjetsTemplate){
    histo_Zjets->Scale(0.0);
    histo_Zjets->Add(fhDZjets);
    histo_ZjetsNoW->Scale(0.0);
    histo_ZjetsNoW->Add(fhDZjets);
    double ZJetsNorm[3] = {0.56, 0.52, 1.};
    if     (MVAVarType == 2 || MVAVarType == 3 || MVAVarType == 4) {
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

  double mean,up,diff;
  for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++){
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
  printf("EWK Corr: WZ(%f/%f/%f) ZZ(%f/%f/%f) ggZZ(%f/%f/%f)\n",
                     histo_WZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_WZ->GetSumOfWeights(),histo_WZ_CMS_EWKCorrDown->GetSumOfWeights(),
                     histo_ZZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_EWKCorrDown->GetSumOfWeights(),
		     histo_ZZ_CMS_ggCorrUp ->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_ggCorrDown ->GetSumOfWeights());

  printf("QCD Corr: WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) ZH(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZH_hinv->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  double systEM[2] = {1.0, 1.0};
  if(histo_EM->GetSumOfWeights() > 1) systEM[0] = 1. + 1./sqrt(histo_EM->GetSumOfWeights());
  else  			      systEM[0] = 2.;
  if(useEMFromData == true){
    double EMNormFact[3] = {((bgdDecay[SIGSEL][0]-EMbkg)*NemFact[0])/bgdDecay[SIGSEL+nSelTypes*(1)][1],
                            ((bgdDecay[SIGSEL][0]-EMbkg)*NemFact[1])/bgdDecay[SIGSEL+nSelTypes*(2)][1],
                            ((bgdDecay[SIGSEL][0]-EMbkg)*NemFact[2])/bgdDecay[SIGSEL+nSelTypes*(3)][1]};
                            
    printf("EM(1): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[SIGSEL+nSelTypes*1][1],bgdDecay[SIGSEL][0],EMbkg,NemFact[0],bgdDecay[SIGSEL+nSelTypes*(1)][1],bgdDecay[SIGSEL+nSelTypes*1][1]*EMNormFact[0],bgdDecay[SIGSEL+nSelTypes*1][1]*EMNormFact[0]*EMSystTotal[0]);
    printf("EM(2): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[SIGSEL+nSelTypes*2][1],bgdDecay[SIGSEL][0],EMbkg,NemFact[1],bgdDecay[SIGSEL+nSelTypes*(2)][1],bgdDecay[SIGSEL+nSelTypes*2][1]*EMNormFact[1],bgdDecay[SIGSEL+nSelTypes*2][1]*EMNormFact[1]*EMSystTotal[1]);
    printf("EM(3): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[SIGSEL+nSelTypes*3][1],bgdDecay[SIGSEL][0],EMbkg,NemFact[2],bgdDecay[SIGSEL+nSelTypes*(3)][1],bgdDecay[SIGSEL+nSelTypes*3][1]*EMNormFact[2],bgdDecay[SIGSEL+nSelTypes*3][1]*EMNormFact[2]*EMSystTotal[2]);

    //systEM[0] = 1. + EMSystTotal[typeSel-1];
    systEM[0] = 1. + EMSyst[1][typeSel-1];
    systEM[1] = TMath::Max(bgdDecay[SIGSEL][0],1.0);
    if(bgdDecay[SIGSEL][0] > 0) histo_EM->Scale(EMNormFact[typeSel-1]);
  }

  histo[allPlots-1][0]->Add(histo_Data);
  histo[allPlots-1][1]->Add(histo_EM);
  histo[allPlots-1][2]->Add(histo_Zjets);
  histo[allPlots-1][3]->Add(histo_WZ);
  histo[allPlots-1][4]->Add(histo_ZZ);
  histo[allPlots-1][5]->Add(histo_VVV);
  histo[allPlots-1][6]->Add(histo_ZH_hinv);
  histo[allPlots-1][7]->Add(histo_ggZH_hinv);
  
  double qcdScaleTotal[2] = {0.035, 0.231};
  double pdfTotal[2] = {0.016, 0.051};
  
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histo%szh%s_nice_%s_%d.root",addChan.Data(),finalStateName,processTag.Data(),thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }

  for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_ZH_hinv  ->GetBinContent(i)+factorUp  *histo_ZH_hinv  ->GetBinError(i),0.000001));
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_ZH_hinv  ->GetBinContent(i)+factorDown*histo_ZH_hinv  ->GetBinError(i),0.000001));
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

    histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[i-1]        ->Add(histo_ZH_hinv  ); histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[i-1]      ->SetBinContent(i,TMath::Max(histo_ZH_hinv  ->GetBinContent(i)+factorUp  *histo_ZH_hinv  ->GetBinError(i),0.000001));
    histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[i-1]      ->Add(histo_ZH_hinv  ); histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[i-1]    ->SetBinContent(i,TMath::Max(histo_ZH_hinv  ->GetBinContent(i)+factorDown*histo_ZH_hinv  ->GetBinError(i),0.000001));
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
  }
  char outputLimits[200];
  sprintf(outputLimits,"zll%shinv%s_%s_input_%s.root",addChan.Data(),finalStateName,processTag.Data(),ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data     ->Write();
  histo_ZH_hinv  ->Write();
  histo_Zjets    ->Write();
  histo_VVV      ->Write();
  histo_WZ       ->Write();
  histo_ZZ       ->Write();
  histo_EM       ->Write();
  histo_ggZH_hinv->Write();
  cout << histo_Data     ->GetSumOfWeights() << " ";
  cout << histo_ZH_hinv  ->GetSumOfWeights() << " ";
  cout << histo_Zjets    ->GetSumOfWeights() << " ";
  cout << histo_VVV      ->GetSumOfWeights() << " ";
  cout << histo_WZ       ->GetSumOfWeights() << " ";
  cout << histo_ZZ       ->GetSumOfWeights() << " ";
  cout << histo_EM       ->GetSumOfWeights() << " ";
  cout << histo_ggZH_hinv->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_ZH_hinv_CMS_MVAZHStatBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_ZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingUp	 ->GetBinContent(i)/histo_ZH_hinv  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAZHStatBoundingDown     ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_ZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingDown	 ->GetBinContent(i)/histo_ZH_hinv  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_VVV	->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp 	 ->GetBinContent(i)/histo_VVV	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown        ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_VVV	->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown	 ->GetBinContent(i)/histo_VVV	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp		 ->GetBinContent(i)/histo_WZ       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingDown	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown 	 ->GetBinContent(i)/histo_WZ	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp		 ->GetBinContent(i)/histo_ZZ       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingDown	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown 	 ->GetBinContent(i)/histo_ZZ	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_EM_CMS_MVAEMStatBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_EM	->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingUp		 ->GetBinContent(i)/histo_EM       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_EM_CMS_MVAEMStatBoundingDown	      ->Write(); for(int i=1; i<=histo_ZH_hinv  ->GetNbinsX(); i++) {if(histo_EM	->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingDown 	 ->GetBinContent(i)/histo_EM	   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp   ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp   ->GetBinContent(i)/histo_ggZH_hinv->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv	->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown ->GetBinContent(i)/histo_ggZH_hinv->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffM\n");
  histo_ZH_hinv_CMS_MVALepEffMBoundingUp  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVALepEffMBoundingDown->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffE\n");
  histo_ZH_hinv_CMS_MVALepEffEBoundingUp  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVALepEffEBoundingDown->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties MET\n");
  histo_ZH_hinv_CMS_MVAMETBoundingUp      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_ZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAMETBoundingDown    ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETBoundingDown   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVAMETBoundingUp      ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVAMETBoundingDown    ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAMETBoundingDown   ->GetBinContent(i)/histo_ggZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_ZH_hinv_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_ZH_hinv     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_MVAJESBoundingDown    ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingDown	 ->GetBinContent(i)/histo_ZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingUp 	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingUp 	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_ggZH_hinv     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_MVAJESBoundingDown    ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAJESBoundingDown	 ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties PU\n");
  histo_ZH_hinv_CMS_PUBoundingUp      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_PUBoundingUp      ->GetBinContent(i)/histo_ZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZH_hinv_CMS_PUBoundingDown    ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_PUBoundingDown   ->GetBinContent(i)/histo_ZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingDown	      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_PUBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_PUBoundingDown	      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_PUBoundingUp	      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_PUBoundingDown	      ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_PUBoundingUp    ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_PUBoundingUp	  ->GetBinContent(i)/histo_ggZH_hinv	  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggZH_hinv_CMS_PUBoundingDown  ->Write(); for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_PUBoundingDown   ->GetBinContent(i)/histo_ggZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties EWK Corr\n");
  histo_WZ_CMS_EWKCorrUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_EWKCorrUp  ->GetBinContent(i)/histo_WZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_EWKCorrDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_WZ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_EWKCorrDown->GetBinContent(i)/histo_WZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_EWKCorrUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_EWKCorrDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_EWKCorrDown->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_ggCorrUp	          ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ggCorrUp  ->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_ggCorrDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ggCorrDown->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties Zjets\n");
  histo_Zjets_CMS_ZjetsSystUp	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystUp  ->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zjets_CMS_ZjetsSystDown	  ->Write(); for(int i=1; i<=histo_ZH_hinv->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystDown->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  outFileLimits->Close();

  double lumiE = 1.027;
  double systLepResE[5] = {1.01,1.01,1.01,1.01,1.01};
  double systLepResM[5] = {1.01,1.01,1.01,1.01,1.01};
  double syst_btag = 1.02;
  
  for(int nb=1; nb<=nBinMVA; nb++){
     // QCD study
    double systQCDScale[5] = {TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[0]  ->GetBinContent(nb)-histo_ZH_hinv  ->GetBinContent(nb)),
                              TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]      ->GetBinContent(nb)-histo_VVV      ->GetBinContent(nb)),
                              TMath::Abs(histo_WZ_CMS_QCDScaleBounding[0]       ->GetBinContent(nb)-histo_WZ       ->GetBinContent(nb)),
                              TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0]       ->GetBinContent(nb)-histo_ZZ       ->GetBinContent(nb)),
			      TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_ggZH_hinv->GetBinContent(nb))};
    for(int nqcd=1; nqcd<6; nqcd++) {
      if(TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)  -histo_ZH_hinv->GetBinContent(nb))   > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)  -histo_ZH_hinv  ->GetBinContent(nb));
      if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)      -histo_VVV    ->GetBinContent(nb))   > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)      -histo_VVV      ->GetBinContent(nb));
      if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_WZ     ->GetBinContent(nb))   > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_WZ       ->GetBinContent(nb));
      if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_ZZ     ->GetBinContent(nb))   > systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_ZZ       ->GetBinContent(nb));
      if(TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggZH_hinv->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggZH_hinv->GetBinContent(nb));
    }                 
    systQCDScale[0] = 1 + sqrt(TMath::Max(systQCDScale[0]*systQCDScale[0]/histo_ZH_hinv->GetBinContent(nb)/histo_ZH_hinv->GetBinContent(nb) - qcdScaleTotal[0]*qcdScaleTotal[0],0.0));
    systQCDScale[1] = 1 + systQCDScale[1]/histo_VVV      ->GetBinContent(nb);
    systQCDScale[2] = 1 + systQCDScale[2]/histo_WZ       ->GetBinContent(nb);
    systQCDScale[3] = 1 + systQCDScale[3]/histo_ZZ       ->GetBinContent(nb);
    if(histo_ggZH_hinv->GetBinContent(nb) > 0)
    systQCDScale[4] = 1 + sqrt(TMath::Max(systQCDScale[4]*systQCDScale[4]/histo_ggZH_hinv->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb) - qcdScaleTotal[1]*qcdScaleTotal[1],0.0));
    for(int ntype=0; ntype<5; ntype++) if(systQCDScale[ntype] < 0) systQCDScale[ntype] = 1.0;
    printf("QCDScale(%d): %f %f %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[4]);

    // PDF study
    double systPDF[5];
    histo_Diff->Reset();
    for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_ZH_hinv_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ZH_hinv->GetBinContent(nb))/histo_ZH_hinv->GetBinContent(nb));
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

    printf("PDF(%d): %f %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4]);

    double systLepEffM[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffMBoundingUp     ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZH_hinv_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffMBoundingDown   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingUp         ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffMBoundingUp          ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffMBoundingDown        ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp          ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown        ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[4] = histo_ggZH_hinv_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[4] = histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(nb);

    double systLepEffE[5] = {1.0,1.0,1.0,1.0,1.0};
    if     (histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffEBoundingUp     ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZH_hinv_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffEBoundingDown   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
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
    if(histo_ZH_hinv->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAMETBoundingUp     ->GetBinContent(nb) > 0) systMetUp  [0] = histo_ZH_hinv_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_ZH_hinv->GetBinContent(nb);
    if(histo_ZH_hinv->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAMETBoundingDown   ->GetBinContent(nb) > 0) systMetDown[0] = histo_ZH_hinv_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_ZH_hinv->GetBinContent(nb);
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

    double systJesUp  [5] = {1.0,1.0,1.0,1.0,1.0};
    double systJesDown[5] = {1.0,1.0,1.0,1.0,1.0};
    if(histo_ZH_hinv->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAJESBoundingUp	->GetBinContent(nb) > 0) systJesUp  [0] = histo_ZH_hinv_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ZH_hinv->GetBinContent(nb);
    if(histo_ZH_hinv->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[0] = histo_ZH_hinv_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ZH_hinv->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVAJESBoundingUp 	->GetBinContent(nb) > 0) systJesUp  [1] = histo_VVV_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_VVV->GetBinContent(nb)	  > 0 && histo_VVV_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[1] = histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_MVAJESBoundingUp  	->GetBinContent(nb) > 0) systJesUp  [2] = histo_WZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb)	  > 0 && histo_WZ_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[2] = histo_WZ_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVAJESBoundingUp  	->GetBinContent(nb) > 0) systJesUp  [3] = histo_ZZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb)	  > 0 && histo_ZZ_CMS_MVAJESBoundingDown	->GetBinContent(nb) > 0) systJesDown[3] = histo_ZZ_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAJESBoundingUp	->GetBinContent(nb) > 0) systJesUp  [4] = histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
    if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJesDown[4] = histo_ggZH_hinv_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);

    double systPUUp  [5] = {1.0,1.0,1.0,1.0,1.0};
    double systPUDown[5] = {1.0,1.0,1.0,1.0,1.0};
    if(histo_ZH_hinv->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_PUBoundingUp     ->GetBinContent(nb) > 0) systPUUp  [0] = histo_ZH_hinv_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_ZH_hinv->GetBinContent(nb);
    if(histo_ZH_hinv->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_PUBoundingDown   ->GetBinContent(nb) > 0) systPUDown[0] = histo_ZH_hinv_CMS_PUBoundingDown->GetBinContent(nb)/histo_ZH_hinv->GetBinContent(nb);
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

    double syst_EWKCorrUp[3]   = {1.0,1.0,1.0};
    double syst_EWKCorrDown[3] = {1.0,1.0,1.0};
    if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [0] = histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[0] = histo_WZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [1] = histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[1] = histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb) > 0) syst_EWKCorrUp  [2] = histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb) > 0) syst_EWKCorrDown[2] = histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);

    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_zll%shinv%s_%s_shape_%s_bin%d.txt",addChan.Data(),finalStateName,processTag.Data(),ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
    newcardShape << Form("bin hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    newcardShape << Form("process ZH_hinv Zjets VVV WZ ZZ EM ggZH_hinv\n");
    newcardShape << Form("process 0 1 2 3 4 5 -1\n");
    newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",histo_ZH_hinv->GetBinContent(nb),histo_Zjets->GetBinContent(nb),TMath::Max(histo_VVV->GetBinContent(nb),0.0),histo_WZ->GetBinContent(nb),histo_ZZ->GetBinContent(nb),histo_EM->GetBinContent(nb),histo_ggZH_hinv->GetBinContent(nb));
    newcardShape << Form("lumi_%4s                               lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE);		     
    newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3],systLepEffM[4]);
    newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3],systLepEffE[4]);
    newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3],systLepResM[4]);
    newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3],systLepResE[4]);
    newcardShape << Form("CMS_pu                                 lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systPUUp[0],systPUDown[0],systPUUp[1],systPUDown[1],systPUUp[2],systPUDown[2],systPUUp[3],systPUDown[3],systPUUp[0],systPUDown[0]); // 0 --> 4
    newcardShape << Form("CMS_scale_met                          lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systMetUp[0],systMetDown[0],systMetUp[1],systMetDown[1],systMetUp[2],systMetDown[2],systMetUp[3],systMetDown[3],systMetUp[0],systMetDown[0]); // 0 --> 4
    newcardShape << Form("CMS_scale_j                            lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systJesUp[0],systJesDown[0],systJesUp[1],systJesDown[1],systJesUp[2],systJesDown[2],systJesUp[3],systJesDown[3],systJesUp[0],systJesDown[0]); // 0 --> 4		 
    newcardShape << Form("UEPS			                 lnN  1.030   -     -     -     -     -   1.030\n");
    newcardShape << Form("CMS_eff_b                              lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",syst_btag,syst_btag,syst_btag,syst_btag,syst_btag);
    newcardShape << Form("pdf_qqbar_ACCEPT                       lnN  %7.5f   -   %7.5f %7.5f %7.5f   -     -  \n",systPDF[0],systPDF[1],systPDF[2],systPDF[3]);
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
    newcardShape << Form("QCDscale_VV		                 lnN    -     -     -   %7.5f %7.5f   -     -  \n",systQCDScale[2],systQCDScale[3]);		
    newcardShape << Form("CMS_zllhinv_WZ_EWKCorr                 lnN    -     -     -   %7.5f/%7.5f   -      -    -  \n",syst_EWKCorrUp[0],syst_EWKCorrDown[0]);		
    newcardShape << Form("CMS_zllhinv_ZZ_EWKCorr                 lnN    -     -     -     -   %7.5f/%7.5f   -     -  \n",syst_EWKCorrUp[1],syst_EWKCorrDown[1]);		
    newcardShape << Form("CMS_zllhinv_ggZZCorr                   lnN    -     -     -     -   %7.5f/%7.5f   -     -  \n",syst_EWKCorrUp[2],syst_EWKCorrDown[2]);		
    newcardShape << Form("CMS_zllhinv_ZLLNorm_%s_%s              lnN	-   %7.5f   -	  -     -     -     -  \n",finalStateName,ECMsb.Data(),2.0);	    
    newcardShape << Form("CMS_zllhinv_ZLLShape_%s_%s             lnN	-   %7.5f/%7.5f   -	-     -     -     -  \n",finalStateName,ECMsb.Data(),systZjetsUp[0],systZjetsDown[0]);	    
    newcardShape << Form("CMS_zllhinv_EMSyst_%s_%s               lnN	-     -     -	  -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),systEM[0]);	       
    if(useEMFromData == true)
    newcardShape << Form("CMS_zllhinv_EMNorm_%s_%s      gmN %d  	-     -     -	  -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),(int)systEM[1],histo_EM->GetBinContent(nb)/systEM[1]);	       

    if     (histo_ZH_hinv->GetBinContent(nb)   > 0 && histo_ZH_hinvNoW->GetBinContent(nb) < 20  ) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%d	  gmN %d  %7.5f -      -    -	 -    -      -  \n",finalStateName,ECMsb.Data(),nb-1,(int)histo_ZH_hinvNoW  ->GetBinContent(nb)  ,histo_ZH_hinv  ->GetBinContent(nb)/histo_ZH_hinvNoW  ->GetBinContent(nb));
    else if(histo_ZH_hinv->GetBinContent(nb)   > 0                                              ) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding_%s_Bin%d	    lnN    %7.5f -      -    -    -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ZH_hinv->GetBinError(nb)  /histo_ZH_hinv->GetBinContent(nb)  );

    if     (histo_Zjets->GetBinContent(nb)     > 0 && histo_ZjetsNoW->GetBinContent(nb) < 20    ) newcardShape << Form("CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%d  gmN %d    -  %7.5f   -    -	 -    -      -  \n",finalStateName,ECMsb.Data(),nb-1,(int)histo_ZjetsNoW    ->GetBinContent(nb)  ,histo_Zjets	 ->GetBinContent(nb)/histo_ZjetsNoW    ->GetBinContent(nb));
    else if(histo_Zjets->GetBinContent(nb)     > 0                                              ) newcardShape << Form("CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%d  lnN      -  %7.5f   -    -    -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_Zjets->GetBinError(nb)    /histo_Zjets->GetBinContent(nb)    );

    if     (histo_VVV->GetBinContent(nb)       > 0 && histo_VVVNoW->GetBinContent(nb) < 20      ) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%d    gmN %d    -	  -  %7.5f   -    -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,(int)histo_VVVNoW      ->GetBinContent(nb)  ,histo_VVV	 ->GetBinContent(nb)/histo_VVVNoW      ->GetBinContent(nb));
    else if(histo_VVV->GetBinContent(nb)       > 0                                              ) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%d    lnN      -	-  %7.5f   -    -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_VVV->GetBinError(nb)      /histo_VVV->GetBinContent(nb)      );

    if     (histo_WZ->GetBinContent(nb)        > 0 && histo_WZNoW->GetBinContent(nb) < 20       ) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%d	  gmN %d    -	  -	-  %7.5f  -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,(int)histo_WZNoW	    ->GetBinContent(nb)  ,histo_WZ	 ->GetBinContent(nb)/histo_WZNoW       ->GetBinContent(nb));
    else if(histo_WZ->GetBinContent(nb)        > 0                                              ) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%d	    lnN      -	-     -  %7.5f  -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_WZ->GetBinError(nb)       /histo_WZ->GetBinContent(nb)       );

    if      (histo_ZZ->GetBinContent(nb)       > 0 && histo_ZZNoW->GetBinContent(nb) < 20       ) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%d	  gmN %d    -	  -	-    -  %7.5f  -     -  \n",finalStateName,ECMsb.Data(),nb-1,(int)histo_ZZNoW	    ->GetBinContent(nb)  ,histo_ZZ	 ->GetBinContent(nb)/histo_ZZNoW       ->GetBinContent(nb));
    else if(histo_ZZ->GetBinContent(nb)	       > 0                                              ) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%d	    lnN      -	-     -    -  %7.5f  -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ZZ->GetBinError(nb)       /histo_ZZ->GetBinContent(nb)       );

    if      (histo_EM->GetBinContent(nb)       > 0                                              ) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%d	    lnN      -	-     -    -    -  %7.5f   -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_EM->GetBinError(nb)       /histo_EM->GetBinContent(nb)	  );

    if     (histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinvNoW->GetBinContent(nb) < 20) newcardShape << Form("CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%d	  gmN %d    -	  -	-    -    -    -   %7.5f\n",finalStateName,ECMsb.Data(),nb-1,(int)histo_ggZH_hinvNoW->GetBinContent(nb)  ,histo_ggZH_hinv->GetBinContent(nb)/histo_ggZH_hinvNoW->GetBinContent(nb));
    else if(histo_ggZH_hinv->GetBinContent(nb) > 0                                              ) newcardShape << Form("CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%d   lnN      -    -     -    -    -    -   %7.5f\n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ggZH_hinv->GetBinError(nb)/histo_ggZH_hinv->GetBinContent(nb));

  }
}
