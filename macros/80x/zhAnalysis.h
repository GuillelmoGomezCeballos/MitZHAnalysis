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
#include "PandaAnalysis/Flat/interface/PandaLeptonicAnalyzer.h"

#include "MitZHAnalysis/macros/80x/zhMVA.h"

Double_t lumi = 3.8;
bool isMIT = true;
// File instances on EOS
// These paths are OUT OF DATE! ~DGH Aug 1 2017
TString filesPathDA   = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/met_";
TString filesPathMC   = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/met_";
TString filesPathMC2  = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/met_";
TString filesPathDMMC = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/";
// File instances on T3 hadoop
if(isMIT){
  filesPathDA   = "/data/t3home000/ceballos/panda/v_003_0/";
  filesPathMC	  = "/data/t3home000/ceballos/panda/v_003_0/";
  filesPathMC2  = "/data/t3home000/ceballos/panda/v_003_0/";
  filesPathDMMC = "/data/t3home000/ceballos/panda/v_003_0/";
}

