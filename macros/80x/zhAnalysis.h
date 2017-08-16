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
    bool       isMIT                   = true;
    unsigned   nJetsType               = 1;
    bool       isBlinded               = false;
    int        typeSel                 = 3;
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
    void Run(bool verbose=true);
  private:
    TString selTypeName[nSelTypes]={"ZSEL",  "SIGSEL", "ZHGSEL", "WWLOOSESEL", "BTAGSEL", "WZSEL", "PRESEL", "CR1SEL", "CR2SEL", "CR12SEL", "TIGHTSEL", "DYSANESEL1", "DYSANESEL2"};
    TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};
    const TString typeLepSel = "medium";
    const double bTagCuts[1] = {0.8484}; // 0.5426/0.8484/0.9535 (check BTagCalibration2Reader!)
    const bool useDYPT = true;
    const unsigned int  num_bdt_toys = 1000;
    string subdirectory="";
    unsigned int randomToySeed=0;
    int nSigModels;
    unsigned nInputFiles;

    // File prefixes
    TString filesPathDA   = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathMC	  = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathMC2  = "/data/t3home000/ceballos/panda/v_003_0/";
    TString filesPathDMMC = "/data/t3home000/ceballos/panda/v_003_0/";

    vector<FlatFile> inputFlatFiles; vector<TString> signalName_;

    void loadFlatFiles(bool doDM=false);
};
