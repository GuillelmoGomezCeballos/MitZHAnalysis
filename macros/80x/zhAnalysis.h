#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <iostream>
#include <fstream>
#include "TMVA/Reader.h"

#include "MitZHAnalysis/macros/80x/FlatFile.h"
#include "PandaAnalysis/Flat/interface/GeneralLeptonicTree.h"
#include "MitZHAnalysis/macros/80x/zhMVA.h"

// MVA variable types:
// 1: MET only
// 2: MET x mll
// 3: classifier only
// 4: MET x classifier

//const int MVAVarType = 0; const int nBinMVA = 8; Double_t xbins[nBinMVA+1] = {0, 50, 200, 250, 300, 400, 600, 800, 1000}; TString addChan = "";
const int MVAVarType = 1; const int nBinMVA = 12; Double_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600}; TString addChan = "1";
//const int MVAVarType = 2; const int nBinMVA = 20; Double_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 350,
//                                                                                         1125,1150,1175,1200,1250,1350,
//                                                                                             2125,2150,2175,2200,2250,2350}; TString addChan = "2";
//const int MVAVarType = 3; const int nBinMVA = 12; Double_t xbins[nBinMVA+1] =  {-2, -1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}; TString addChan = "3";
//const int MVAVarType = 4; const int nBinMVA = 26; Double_t xbins[nBinMVA+1] = {0, 50, 100, 125, 150, 175, 200, 250, 350,
//                                                                                         1125,1150,1175,1200,1250,1350,
//                                                                                         2125,2150,2175,2200,2250,2350,
//                                                                                         3125,3150,3175,3200,3250,3350}; TString addChan = "4";
const int      allPlots     = 37; //number of variable plots
const int      processTypes = 8; // number of process categories
const int      numberCuts   = 11; //number of cuts for the cutflow numbers
const unsigned num_bdt_toys = 1000;

class zhAnalysis {
  public:
    // Settings
    bool       isMIT                   = true;
    bool       isBlinded               = false;
    int        plotModel               = 0;
    double     lumi                    = 3.8;
    bool       isMINIAOD               = true;
    bool       useZjetsTemplate        = false;
    bool       usePureMC               = true; 
    bool       useEMFromData           = true;
    bool       useVVFromData           = true;
    bool       useZZWZEWKUnc           = true;
    double     mcPrescale              = 1.;
    bool       useBDT                  = false;
    bool       useCachedBDTSystematics = false;
    string     the_BDT_weights         = "";

    // kluge fix because aclic does not want to compile CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h
    // which is included in PandaLeptonicAnalyzer.h
    enum SelectionBit {
      kLoose   =(1<<0),
      kFake    =(1<<1),
      kMedium  =(1<<2),
      kTight   =(1<<3)
    };

    enum TriggerBits {
      kMETTrig       =(1<<0),
      kSinglePhoTrig =(1<<1),
      kMuEGTrig      =(1<<2),
      kMuMuTrig      =(1<<3),
      kMuTrig        =(1<<4),
      kEGEGTrig      =(1<<5),
      kEGTrig        =(1<<6)
    };
    enum selType {
      ZSEL=0,
      SIGSEL,
      ZHGSEL,
      WWLOOSESEL,
      BTAGSEL,
      WZSEL,
      PRESEL,
      CR1SEL,
      CR2SEL,
      CR12SEL,
      TIGHTSEL,
      DYSANESEL1,
      DYSANESEL2,
      nSelTypes
    };
    enum systType {
      JESUP=0,
      JESDOWN,
      METUP,
      METDOWN,
      nSystTypes
    };
    zhAnalysis(string subdirectory_="");
    ~zhAnalysis();
    bool LoadFlatFiles(bool doDM=false);
    void Run(bool verbose_=false, int nJetsType_=1, int typeSel_=3, int plotModel_=0);
  private:
    // Hardcoded settings
    
    TString selTypeName[nSelTypes]={"ZSEL",  "SIGSEL", "ZHGSEL", "WWLOOSESEL", "BTAGSEL", "WZSEL", "PRESEL", "CR1SEL", "CR2SEL", "CR12SEL", "TIGHTSEL", "DYSANESEL1", "DYSANESEL2"};
    TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};
    const TString  typeLepSel   = "medium";
    const bool     useDYPT      = true;
    const double   bTagCut      = 0.8484;
    const int period            = 1;
    const TString ECMsb         = "13TeV2017";
    char finalStateName[4];
    string         subdirectory = "";
    TString effMName="CMS_eff2017_m",   effEName="CMS_eff2017_e",
            momMName="CMS_scale2017_m", momEName="CMS_scale2017_e";

    // Systematics numbers
    double systEM[2] = {1.0, 1.0};
    double qcdScaleTotal[2] = {0.035, 0.231};
    double pdfTotal[2] = {0.016, 0.051};
    double quantileProbs[3]={0.159,0.5,0.841};

    // Physics quantities
    GeneralLeptonicTree gltEvent;
    float normalizedWeight, sf_btag0, sf_btag0BUp, sf_btag0BDown, sf_btag0MUp, sf_btag0MDown;
    TTree *inputTree;
    
    // Switches, counters, other storage
    unsigned       randomToySeed;
    int nSigModels;
    unsigned nInputFiles=0;
    bool verbose=false;
    bool madeHistos=false;
    float bdt_toy_scale[num_bdt_toys];
    unsigned   nJetsType;
    int        typeSel;
    

    // File prefixes/paths
    TString filesPathDA   = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathMC	  = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathMC2  = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathDMMC = "/data/t3home000/ceballos/panda/v_003_0/";
    TString zjetsTemplatesPath = "";
    char filenameBDTSysts[200]; 

    // Files, signals, category/process names
    vector<FlatFile> inputFlatFiles; vector<TString> signalName_;
    TString plotName_[allPlots];
    TString categoryName_[processTypes]={"Data","Nonresonant","Z+jets","WZ","ZZ","VVV","qqZH","ggZH"};
    TString processName[processTypes] = {"..Data", "....EM", "...DY", "...WZ", "....ZZ", "...VVV", "....ZH", "..ggZH"};
    TFile* cachedSystFile;

    // Allocate memory for histograms. Support only 500 signal models (hard-coded).
    TH1D* histoZHSEL[4];
    TH1D* histoMVA;
    TH1D *histo_Data;
    TH1D *histo_Zjets;         
    TH1D *histo_VVV;         
    TH1D *histo_WZ;         
    TH1D *histo_ZZ;
    TH1D *histo_EM;         
    TH1D *histo_ggZH_hinv; 
    TH1D *histo_ZjetsNoW;         
    TH1D *histo_VVVNoW;         
    TH1D *histo_WZNoW;         
    TH1D *histo_ZZNoW;
    TH1D *histo_EMNoW;         
    TH1D *histo_ggZH_hinvNoW; 
    TH1D *histo_ZH_hinv[500];
    TH1D *histo_ZH_hinvNoW[500];
    TH1D* histo[allPlots][processTypes];
    TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingDown[500];
    TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingUp;
    TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingDown;
    TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp;
    TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown;
    TH1D* histo_WZ_CMS_MVAWZStatBoundingUp;
    TH1D* histo_WZ_CMS_MVAWZStatBoundingDown;
    TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp;
    TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown;
    TH1D* histo_EM_CMS_MVAEMStatBoundingUp;
    TH1D* histo_EM_CMS_MVAEMStatBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown;
    TH1D* histo_ZH_hinv_CMS_QCDScaleBounding[500][6];
    TH1D* histo_VVV_CMS_QCDScaleBounding[6];
    TH1D* histo_WZ_CMS_QCDScaleBounding[6];
    TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
    TH1D* histo_ggZH_hinv_CMS_QCDScaleBounding[6];
    TH1D* histo_ZH_hinv_CMS_PDFUp[500];
    TH1D* histo_ZH_hinv_CMS_PDFDown[500];
    TH1D* histo_VVV_CMS_PDFUp;
    TH1D* histo_VVV_CMS_PDFDown;
    TH1D* histo_WZ_CMS_PDFUp;
    TH1D* histo_WZ_CMS_PDFDown;
    TH1D* histo_ZZ_CMS_PDFUp;
    TH1D* histo_ZZ_CMS_PDFDown;
    TH1D* histo_ggZH_hinv_CMS_PDFUp;
    TH1D* histo_ggZH_hinv_CMS_PDFDown;
 
    TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[500][nBinMVA];
    TH1D* histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[500][nBinMVA];
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
    TH1D* histo_VVV_CMS_MVALepEffMBoundingUp;
    TH1D* histo_VVV_CMS_MVALepEffMBoundingDown;
    TH1D* histo_WZ_CMS_MVALepEffMBoundingUp;
    TH1D* histo_WZ_CMS_MVALepEffMBoundingDown;
    TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp;
    TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingDown;
    TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingUp[500]; 
    TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingDown[500];

    TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg;
    TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg;
    TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg;
    TH1D* histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg;
    TH1D* histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[500];

    TH1D* histo_VVV_CMS_MVALepEffEBoundingUp;
    TH1D* histo_VVV_CMS_MVALepEffEBoundingDown;
    TH1D* histo_WZ_CMS_MVALepEffEBoundingUp;
    TH1D* histo_WZ_CMS_MVALepEffEBoundingDown;
    TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp;
    TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingDown;
    TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingDown[500];

    TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg;
    TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg;
    TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg;
    TH1D* histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg;
    TH1D* histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[500];

    TH1D* histo_VVV_CMS_MVAMETBoundingUp;
    TH1D* histo_VVV_CMS_MVAMETBoundingDown;
    TH1D* histo_WZ_CMS_MVAMETBoundingUp;
    TH1D* histo_WZ_CMS_MVAMETBoundingDown;
    TH1D* histo_ZZ_CMS_MVAMETBoundingUp;
    TH1D* histo_ZZ_CMS_MVAMETBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_MVAMETBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_MVAMETBoundingDown;
    TH1D* histo_ZH_hinv_CMS_MVAMETBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_MVAMETBoundingDown[500];

    TH1D* histo_VVV_CMS_MVAJESBoundingUp;
    TH1D* histo_VVV_CMS_MVAJESBoundingDown;
    TH1D* histo_WZ_CMS_MVAJESBoundingUp;
    TH1D* histo_WZ_CMS_MVAJESBoundingDown;
    TH1D* histo_ZZ_CMS_MVAJESBoundingUp;
    TH1D* histo_ZZ_CMS_MVAJESBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_MVAJESBoundingDown;
    TH1D* histo_ZH_hinv_CMS_MVAJESBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_MVAJESBoundingDown[500];

    TH1D* histo_VVV_CMS_MVABTAGBoundingUp;
    TH1D* histo_VVV_CMS_MVABTAGBoundingDown;
    TH1D* histo_WZ_CMS_MVABTAGBoundingUp;
    TH1D* histo_WZ_CMS_MVABTAGBoundingDown;
    TH1D* histo_ZZ_CMS_MVABTAGBoundingUp;
    TH1D* histo_ZZ_CMS_MVABTAGBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_MVABTAGBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_MVABTAGBoundingDown;
    TH1D* histo_ZH_hinv_CMS_MVABTAGBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_MVABTAGBoundingDown[500];

    TH1D* histo_VVV_CMS_BDTMuonScaleBoundingUp;
    TH1D* histo_VVV_CMS_BDTMuonScaleBoundingDown;
    TH1D* histo_WZ_CMS_BDTMuonScaleBoundingUp;
    TH1D* histo_WZ_CMS_BDTMuonScaleBoundingDown;
    TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingUp;
    TH1D* histo_ZZ_CMS_BDTMuonScaleBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_BDTMuonScaleBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_BDTMuonScaleBoundingDown;
    TH1D* histo_ZH_hinv_CMS_BDTMuonScaleBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_BDTMuonScaleBoundingDown[500];

    TH1D* histo_VVV_CMS_BDTElectronScaleBoundingUp;
    TH1D* histo_VVV_CMS_BDTElectronScaleBoundingDown;
    TH1D* histo_WZ_CMS_BDTElectronScaleBoundingUp;
    TH1D* histo_WZ_CMS_BDTElectronScaleBoundingDown;
    TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingUp;
    TH1D* histo_ZZ_CMS_BDTElectronScaleBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_BDTElectronScaleBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_BDTElectronScaleBoundingDown;
    TH1D* histo_ZH_hinv_CMS_BDTElectronScaleBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_BDTElectronScaleBoundingDown[500];

    TH1D* histo_VVV_CMS_BDTMETScaleBoundingUp;
    TH1D* histo_VVV_CMS_BDTMETScaleBoundingDown;
    TH1D* histo_WZ_CMS_BDTMETScaleBoundingUp;
    TH1D* histo_WZ_CMS_BDTMETScaleBoundingDown;
    TH1D* histo_ZZ_CMS_BDTMETScaleBoundingUp;
    TH1D* histo_ZZ_CMS_BDTMETScaleBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_BDTMETScaleBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_BDTMETScaleBoundingDown;
    TH1D* histo_ZH_hinv_CMS_BDTMETScaleBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_BDTMETScaleBoundingDown[500];

    TH1D* histo_VVV_CMS_BDTJetScaleBoundingUp;
    TH1D* histo_VVV_CMS_BDTJetScaleBoundingDown;
    TH1D* histo_WZ_CMS_BDTJetScaleBoundingUp;
    TH1D* histo_WZ_CMS_BDTJetScaleBoundingDown;
    TH1D* histo_ZZ_CMS_BDTJetScaleBoundingUp;
    TH1D* histo_ZZ_CMS_BDTJetScaleBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_BDTJetScaleBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_BDTJetScaleBoundingDown;
    TH1D* histo_ZH_hinv_CMS_BDTJetScaleBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_BDTJetScaleBoundingDown[500];

    TH1D* histo_VVV_CMS_PUBoundingUp;
    TH1D* histo_VVV_CMS_PUBoundingDown;
    TH1D* histo_WZ_CMS_PUBoundingUp;
    TH1D* histo_WZ_CMS_PUBoundingDown;
    TH1D* histo_ZZ_CMS_PUBoundingUp;
    TH1D* histo_ZZ_CMS_PUBoundingDown;
    TH1D* histo_ggZH_hinv_CMS_PUBoundingUp;
    TH1D* histo_ggZH_hinv_CMS_PUBoundingDown;
    TH1D* histo_ZH_hinv_CMS_PUBoundingUp[500];  
    TH1D* histo_ZH_hinv_CMS_PUBoundingDown[500];

    TH1D* histo_ZH_hinv_CMS_EWKCorrUp[500];
    TH1D* histo_ZH_hinv_CMS_EWKCorrDown[500];
    TH1D* histo_WZ_CMS_EWKCorrUp;
    TH1D* histo_WZ_CMS_EWKCorrDown;
    TH1D* histo_ZZ_CMS_EWKCorrUp;
    TH1D* histo_ZZ_CMS_EWKCorrDown;
    TH1D* histo_ZZ_CMS_ggCorrUp;
    TH1D* histo_ZZ_CMS_ggCorrDown;
    TH1D* histo_Zjets_CMS_ZjetsSystUp;
    TH1D* histo_Zjets_CMS_ZjetsSystDown;

    // Pointers for TH2 objects to store the toy BDT shapes
    TH2F* histo_bdt_toys_electronScale_VVV, *histo_bdt_toys_electronScale_WZ, *histo_bdt_toys_electronScale_ZZ, *histo_bdt_toys_electronScale_ggZH_hinv, *histo_bdt_toys_electronScale_ZH_hinv[500];
    TH2F* histo_bdt_toys_muonScale_VVV, *histo_bdt_toys_muonScale_WZ, *histo_bdt_toys_muonScale_ZZ, *histo_bdt_toys_muonScale_ggZH_hinv, *histo_bdt_toys_muonScale_ZH_hinv[500];
    TH2F* histo_bdt_toys_METScale_VVV, *histo_bdt_toys_METScale_WZ, *histo_bdt_toys_METScale_ZZ, *histo_bdt_toys_METScale_ggZH_hinv, *histo_bdt_toys_METScale_ZH_hinv[500];
    // Pointers for TH1 arrays to store the bin yields
    TH1F *bdt_toy_binyields_electronScale_VVV[nBinMVA], *bdt_toy_binyields_electronScale_WZ[nBinMVA], *bdt_toy_binyields_electronScale_ZZ[nBinMVA], *bdt_toy_binyields_electronScale_ZH_hinv[500][nBinMVA], *bdt_toy_binyields_electronScale_ggZH_hinv[nBinMVA];
    TH1F *bdt_toy_binyields_muonScale_VVV[nBinMVA], *bdt_toy_binyields_muonScale_WZ[nBinMVA], *bdt_toy_binyields_muonScale_ZZ[nBinMVA], *bdt_toy_binyields_muonScale_ZH_hinv[500][nBinMVA], *bdt_toy_binyields_muonScale_ggZH_hinv[nBinMVA];
    TH1F *bdt_toy_binyields_METScale_VVV[nBinMVA], *bdt_toy_binyields_METScale_WZ[nBinMVA], *bdt_toy_binyields_METScale_ZZ[nBinMVA], *bdt_toy_binyields_METScale_ZH_hinv[500][nBinMVA], *bdt_toy_binyields_METScale_ggZH_hinv[nBinMVA];
    // Pointers for TH1's to store the relative up/down systematics (this is what will be cached)
    TH1F *bdt_syst_electronScaleUp_VVV, *bdt_syst_electronScaleUp_WZ, *bdt_syst_electronScaleUp_ZZ, *bdt_syst_electronScaleUp_ZH_hinv[500], *bdt_syst_electronScaleUp_ggZH_hinv;
    TH1F *bdt_syst_muonScaleUp_VVV, *bdt_syst_muonScaleUp_WZ, *bdt_syst_muonScaleUp_ZZ, *bdt_syst_muonScaleUp_ZH_hinv[500], *bdt_syst_muonScaleUp_ggZH_hinv;
    TH1F *bdt_syst_METScaleUp_VVV, *bdt_syst_METScaleUp_WZ, *bdt_syst_METScaleUp_ZZ, *bdt_syst_METScaleUp_ZH_hinv[500], *bdt_syst_METScaleUp_ggZH_hinv;
    TH1F *bdt_syst_electronScaleDown_VVV, *bdt_syst_electronScaleDown_WZ, *bdt_syst_electronScaleDown_ZZ, *bdt_syst_electronScaleDown_ZH_hinv[500], *bdt_syst_electronScaleDown_ggZH_hinv;
    TH1F *bdt_syst_muonScaleDown_VVV, *bdt_syst_muonScaleDown_WZ, *bdt_syst_muonScaleDown_ZZ, *bdt_syst_muonScaleDown_ZH_hinv[500], *bdt_syst_muonScaleDown_ggZH_hinv;
    TH1F *bdt_syst_METScaleDown_VVV, *bdt_syst_METScaleDown_WZ, *bdt_syst_METScaleDown_ZZ, *bdt_syst_METScaleDown_ZH_hinv[500], *bdt_syst_METScaleDown_ggZH_hinv;
      // Pointers for TH2's which hold the 2D BDT toy envelope maps
    TH2F *bdt_toy_envelope_electronScale_VVV, *bdt_toy_envelope_electronScale_WZ, *bdt_toy_envelope_electronScale_ZZ, *bdt_toy_envelope_electronScale_ZH_hinv[500], *bdt_toy_envelope_electronScale_ggZH_hinv;
    TH2F *bdt_toy_envelope_muonScale_VVV, *bdt_toy_envelope_muonScale_WZ, *bdt_toy_envelope_muonScale_ZZ, *bdt_toy_envelope_muonScale_ZH_hinv[500], *bdt_toy_envelope_muonScale_ggZH_hinv;
    TH2F *bdt_toy_envelope_METScale_VVV, *bdt_toy_envelope_METScale_WZ, *bdt_toy_envelope_METScale_ZZ, *bdt_toy_envelope_METScale_ZH_hinv[500], *bdt_toy_envelope_METScale_ggZH_hinv;

    // Private member functions
    bool  ComputeSystematics();
    bool  LoadCachedBDTSystematics();
    TH1D* MakeHisto(unsigned int thePlot, TString &plotName);
    bool  MakeHistos();
    bool  SaveBDTSystematics(int nModel);
    bool  SaveHistos();
    void  SetFinalStateName();
    void  SetBranchAddresses(TTree *theInputTree, bool isData=false);
    bool  SetupBDTSystematics();
};
