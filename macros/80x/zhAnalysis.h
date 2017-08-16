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

#include "PandaAnalysis/Flat/interface/GeneralLeptonicTree.h"
//#include "PandaAnalysis/Flat/interface/PandaLeptonicAnalyzer.h"
#include "MitZHAnalysis/macros/80x/zhMVA.h"
#include "MitZHAnalysis/macros/80x/FlatFile.h"

class zhAnalysis {
  public:
    // Settings
    int        nJetsType;
    int        typeSel;
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
    void LoadFlatFiles(bool doDM=false);
    void Run(bool verbose=true, int nJetsType_=1, int typeSel_=3);
  private:
    // Hardcoded settings
    
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
    TString selTypeName[nSelTypes]={"ZSEL",  "SIGSEL", "ZHGSEL", "WWLOOSESEL", "BTAGSEL", "WZSEL", "PRESEL", "CR1SEL", "CR2SEL", "CR12SEL", "TIGHTSEL", "DYSANESEL1", "DYSANESEL2"};
    TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};
    const TString typeLepSel = "medium";
    const bool useDYPT = true;
    const unsigned int  num_bdt_toys = 1000;
    string subdirectory="";
    unsigned int randomToySeed=0;
    TString effMName="CMS_eff2017_m",   effEName="CMS_eff2017_e",
            momMName="CMS_scale2017_m", momEName="CMS_scale2017_e";
    char finalStateName[4];

    // Switches and counters
    int nSigModels;
    unsigned nInputFiles=0;

    // File prefixes
    TString filesPathDA   = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathMC	  = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathMC2  = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathDMMC = "/data/t3home000/ceballos/panda/v_003_0/";

    vector<FlatFile> inputFlatFiles; vector<TString> signalName_;


    void SetFinalStateName();
};
