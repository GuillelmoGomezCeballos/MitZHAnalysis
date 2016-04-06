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

const TString typeLepSel = "medium";

void zhAnalysis_DGH(
  unsigned int nJets=0

) {
  Double_t lumi = 2.318;
  int mH=125;
  TString processTag = "";
  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
  
  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  

  double events_0jet_ee, events_0jet_mm, events_1jet_ee, events_1jet_mm;

  // Input files
  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/met_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/met_";
  vector<TString> infilenamev;
  vector<Int_t> infilecatv;
  TString puPath = "";
  TString zjetsTemplatesPath = "";
  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  //TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  //TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
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

  infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));                                                infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));                                                infilecatv.push_back(0);

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));                    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                        infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(1);

  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));        infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));             infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                            infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));       infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));               infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3+AODSIM.root",filesPathMC.Data()));                        infilecatv.push_back(1);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));      infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(2);

  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));               infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));            infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                    infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));               infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                       infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));         infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));              infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));             infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));            infilecatv.push_back(4);

  //infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));             infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                         infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));          infilecatv.push_back(5);

  processTag = Form("mh%d",mH);
  infilenamev.push_back(Form("%sZH_ZToMM_HToInvisible_M%d_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data(),mH)); infilecatv.push_back(6);
  infilenamev.push_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data(),mH)); infilecatv.push_back(6);

  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    TFile input_file(infilenamev[ifile]);
    assert(input_file.IsOpen() && !input_file.IsZombie());
    TTree *input_tree = (TTree*)input_file.FindObjectAny("events");
    TTree *input_all  = (TTree*)input_file.FindObjectAny("all");
    TTree *PDF_tree   = (TTree*)input_file.FindObjectAny("pdfReweight");

    BareEvent eventEvent;
    BareJets eventJets;
    BareLeptons eventLeptons;
    BareMet eventMet;
    BareMonteCarlo eventMonteCarlo;
    BareTaus eventTaus;
    BareTrigger eventTrigger;
    BareVertex eventVertex;
    eventEvent.setBranchAddresses(input_tree);
    eventJets.setBranchAddresses(input_tree);
    eventLeptons.setBranchAddresses(input_tree);
    eventMet.SetExtend();
    eventMet.setBranchAddresses(input_tree);
    eventMonteCarlo.setBranchAddresses(input_tree);
    eventTaus.SetExtend();
    eventTaus.setBranchAddresses(input_tree);
    eventTrigger.setBranchAddresses(input_tree);
    eventVertex.setBranchAddresses(input_tree);

    TNamed *triggerNamesCommaDelimited = (TNamed*)input_file.FindObjectAny("triggerNames");
    char **triggerNames;
    size_t numTriggers;
    triggerNames = strsplit(triggerNamesCommaDelimited->GetTitle(), ",", &numTriggers);

    //if(infilecatv[ifile] != 2) continue; // only look at data for now
    
    //if(infilecatv[ifile] == 0){
    //  for (int i = 0; i < (int)numTriggers; i++) {
    //    printf("triggerNames(%2d): \"%s\"\n",(int)i, triggerNames[i]);
    //  }
    //}
    //else {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());
    //}
    
    double numFoundLeptonCandidates, numPassTrigger;
    double numTwoGoodLeptons                [5]={0,0,0,0,0},
           numPassDileptonPtCut             [5]={0,0,0,0,0},
           numPassMetPreselection           [5]={0,0,0,0,0},
           numPassPreselection              [5]={0,0,0,0,0},
           numPassMetSelection              [5]={0,0,0,0,0},
           numPassNJetsCut                  [5]={0,0,0,0,0},
           numPassZMassWindowCut            [5]={0,0,0,0,0},
           numPass3rdLeptonVeto             [5]={0,0,0,0,0},
           numPassBJetVeto                  [5]={0,0,0,0,0},
           numPassSoftMuonVeto              [5]={0,0,0,0,0},
           numPassDeltaPhiDileptonMetCut    [5]={0,0,0,0,0},
           numPassDeltaPhiLeptonsCut        [5]={0,0,0,0,0},
           numPassMetBalanceCut             [5]={0,0,0,0,0},
           numPassTransverseMassCut         [5]={0,0,0,0,0},
           numPassDeltaPhiJetMetCut         [5]={0,0,0,0,0},
           numPassTauVeto                   [5]={0,0,0,0,0},
           numPassEventSelection            [5]={0,0,0,0,0},
           badBoostedNonZLeptons            [5]={0,0,0,0,0},
           goodBoostedNonZLeptons           [5]={0,0,0,0,0};
    vector<ULong64_t> eventNums;
    for(int i=0; i<(int)input_tree->GetEntries(); i++) {
      input_tree->GetEntry(i); 
      
      // Look for lepton candidates
      bool foundLeptonCandidates = (eventLeptons.p4->GetEntriesFast() >= 2 &&
         ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 &&
         ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10);
      
      // Check if the event passes the trigger
      bool passTrigger=false;
      for(int ntrigger=0; ntrigger < (int)numTriggers && !passTrigger; ntrigger++) {
        if( (*eventTrigger.triggerFired)[ntrigger] != 0 && (
          (strcmp( triggerNames[ntrigger], "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*")         == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")       == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_IsoMu27_v*")                     == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_IsoMu20_v*")                     == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_IsoTkMu20_v*")                   == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")       == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Ele23_WPLoose_Gsf_v*")               == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Ele22_eta2p1_WP75_Gsf_v*")               == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Ele27_WPLoose_Gsf_v*")               == 0) ||
          (strcmp( triggerNames[ntrigger], "HLT_Ele27_WP85_Gsf_v*")                  == 0)
        )) passTrigger=true;
      }
      
      if(foundLeptonCandidates) ++numFoundLeptonCandidates; else continue;
      if(passTrigger) ++numPassTrigger; else continue;
      double met = ((TLorentzVector*)(*eventMet.p4)[0])->Pt();
      // Apply the pre-selection cuts
      bool passMetPreselection = met > 45.;
      bool passDileptonPtCut=false;
      bool passPreselection=false;
      
      // Record how many well identified leptons we have
      vector<int> goodLeptons, softLeptons, fakeLeptons, goodOrFakeLeptons;
      int totalLeptonCharge=0;
      for(int il1=0; il1 < eventLeptons.p4->GetEntriesFast(); il1++) {
        if(selectIdIsoCut(
            typeLepSel.Data(),
            TMath::Abs((int)(*eventLeptons.pdgId)[il1]),
            ((TLorentzVector*)(*eventLeptons.p4)[il1])->Pt(),
            TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[il1])->Eta()),
            (double)(*eventLeptons.iso)[il1],
            (int)(*eventLeptons.selBits)[il1],
            (double)(*eventLeptons.mva)[il1]
          ) 
        ) { 
          goodLeptons.push_back(il1); 
          goodOrFakeLeptons.push_back(il1); 
          totalLeptonCharge += (int)(*eventLeptons.pdgId)[il1] / TMath::Abs((int)(*eventLeptons.pdgId)[il1]);
        } else if(((int)(*eventLeptons.selBits)[il1] & BareLeptons::LepFake)== BareLeptons::LepFake ) {
          fakeLeptons.push_back(il1);
          goodOrFakeLeptons.push_back(il1);
          totalLeptonCharge += (int)(*eventLeptons.pdgId)[il1] / TMath::Abs((int)(*eventLeptons.pdgId)[il1]);
        } else if( ((int)(*eventLeptons.selBits)[il1] & BareLeptons::LepSoftIP) == BareLeptons::LepSoftIP ) {
          softLeptons.push_back(il1);
          totalLeptonCharge += (int)(*eventLeptons.pdgId)[il1] / TMath::Abs((int)(*eventLeptons.pdgId)[il1]);
        }

      }
      if(goodOrFakeLeptons.size() < 2 || goodOrFakeLeptons.size() != goodLeptons.size() ) continue;
      // Now try to associate pairs of leptons
      int flavor=-1; // 0 for e-mu, 1 for ee, 2 for mm.
      int lepton1=goodOrFakeLeptons[0], lepton2=goodOrFakeLeptons[1];
      // record event flavor
      if(TMath::Abs((int)(*eventLeptons.pdgId)[lepton1]) == 13 &&
         TMath::Abs((int)(*eventLeptons.pdgId)[lepton2]) == 13) flavor=1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[lepton1]) == 11 &&
              TMath::Abs((int)(*eventLeptons.pdgId)[lepton2]) == 11) flavor=2;
      else flavor=0;
      TLorentzVector dileptonSystem;
//      for(int igl1=0; igl1 < (int)goodLeptons.size() && !passPreselection && foundLeptonCandidates && !passDileptonPtCut; igl1++) {
//        for(int igl2=0; igl2 < (int)goodLeptons.size() && !passPreselection && !passDileptonPtCut; igl2++) {
//          if(igl1==igl2) continue;
//          int il1=goodLeptons[igl1], il2=goodLeptons[igl2];
      if( ((TLorentzVector*)(*eventLeptons.p4)[lepton1])->Pt() <= 20 || ((TLorentzVector*)(*eventLeptons.p4)[lepton2])->Pt() <= 20) continue;
      
      // Now that we have the two leptons and the flavor, handle event weight here
      double totalWeight=1;
      {
        double mcWeight = (infilecatv[ifile] == 0) ? 1 : eventMonteCarlo.mcWeight;
        
        //lepton efficiency
        double leptonScaleFactor = 1.;
        if(infilecatv[ifile] != 0) { for(unsigned int igl1=0; igl1 < (int)goodOrFakeLeptons.size(); igl1++) {
          int il1=goodOrFakeLeptons[igl1];
          leptonScaleFactor *= effhDScaleFactor(
            ((TLorentzVector*)(*eventLeptons.p4)[il1])->Pt(),
            ((TLorentzVector*)(*eventLeptons.p4)[il1])->Eta(),
            TMath::Abs((int)(*eventLeptons.pdgId)[il1]),
            typeLepSel.Data(),
            fhDMuMediumSF,
            fhDMuIsoSF,
            fhDElMediumSF,
            fhDElTightSF,
            fhDElMediumMVASF,
            fhDElTightMVASF          
          );
        }}

        //fake rate
        double fakeScaleFactor = 1.; //not implemented for now

        //trigger efficiency
        double triggerScaleFactor = 1.; // assume 1 for now

        //luminosity
        double luminosityRescale = (infilecatv[ifile] != 0) ? lumi : 1.; 
        
        // pileup
        double puWeight     = (infilecatv[ifile] != 0) ? nPUScaleFactor(fhDPU    , (double)eventMonteCarlo.puTrueInt) : 1.0;

        // electroweak corrections
        vector<int>zzBoson;
        vector<bool> isGenDupl; double bosonPtMin = 1000000000; bool isBosonFound = false; vector<bool> isNeuDupl;
        if(infilecatv[ifile] == 4 && infilenamev[ifile].Contains("GluGlu") == kFALSE) { for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
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
        }}
        
        //double the_rho = 0.0; if(the_rhoP4.P() > 0) the_rho = the_rhoP4.Pt()/the_rhoP4.P();
        double zzCorrection[2] {1,1};
        if(infilecatv[ifile] == 4 && infilenamev[ifile].Contains("GluGlu") == kFALSE) {
          float GENmZZ = 0.0; 
          if(zzBoson.size() >= 2) GENmZZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[1])) ) ).M();
          zzCorrection[0] = weightEWKCorr(bosonPtMin,1);
          zzCorrection[1] = kfactor_qqZZ_qcd_M(GENmZZ);
        }
        totalWeight *= mcWeight * leptonScaleFactor * fakeScaleFactor * triggerScaleFactor * luminosityRescale * zzCorrection[0] * zzCorrection[1];
    

      }
      numTwoGoodLeptons[flavor]+=totalWeight;
      // Ptll cut
      dileptonSystem = ( (* (TLorentzVector*)(eventLeptons.p4->At(lepton1))) + (* (TLorentzVector*)(eventLeptons.p4->At(lepton2))) );
      if(dileptonSystem.Pt() > 60) passDileptonPtCut=true;
      if(passDileptonPtCut) numPassDileptonPtCut[flavor]+=totalWeight; else continue;
      if(passMetPreselection) numPassMetPreselection[flavor]+=totalWeight; else continue;
      if(passMetPreselection && passDileptonPtCut) passPreselection=true;
      if(passPreselection) numPassPreselection[flavor]+=totalWeight; else continue;

      // Start applying event selection
      bool passMetSelection = met > 100.;
      bool passEventSelection=false;
      if(passMetSelection) numPassMetSelection[flavor] += totalWeight; else continue;
      
      // Look for jets in this event and also find the max b discriminator value
      vector<int> goodJets;
      double maxBDiscriminatorValue=0;
      vector<double> deltaPhiJetMet;
      for(int njet=0; njet<(int)eventJets.p4->GetEntriesFast(); njet++) {
        //if(eventEvent.eventNum==77208166) printf("event 77208166 jet #%d has pT %f and B discriminator value %f\n", njet, (((TLorentzVector*)(*eventJets.p4)[njet])->Pt()), (float)(*eventJets.bDiscr)[njet]);
        //if(
        //  !passJetId(fMVACut, (float)(*eventJets.puId)[njet], ((TLorentzVector*)(*eventJets.p4)[njet])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[njet])->Eta()))
        //) continue;
        bool jetMatchedToLepton=false;
        for( int igl1=0; igl1 < (int)goodOrFakeLeptons.size() && !jetMatchedToLepton; igl1++) {
          int il1=goodOrFakeLeptons[igl1];
          // See if jet can be delta-R matched to a good lepton
          if(((TLorentzVector*)(*eventJets.p4)[njet])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[il1])) < 0.3) { jetMatchedToLepton=true; break; }
          //if(eventEvent.eventNum==77208166) printf("\tFailed to delta R match this jet with (pT, eta)=(%f, %f) to lepton with (pT, eta)=(%f,%f)\n", ((TLorentzVector*)(*eventJets.p4)[njet])->Pt(), ((TLorentzVector*)(*eventJets.p4)[njet])->Eta(), ((TLorentzVector*)(*eventLeptons.p4)[il1])->Pt(), ((TLorentzVector*)(*eventLeptons.p4)[il1])->Eta());
        }
        if(jetMatchedToLepton) continue;
        if((((TLorentzVector*)(*eventJets.p4)[njet])->Pt() > 15) && (float)(*eventJets.bDiscr)[njet] > maxBDiscriminatorValue)
          maxBDiscriminatorValue = (float)(*eventJets.bDiscr)[njet];
        if(!(((TLorentzVector*)(*eventJets.p4)[njet])->Pt() > 30)) continue;
        goodJets.push_back(njet);
        deltaPhiJetMet.push_back(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[njet])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))));
      }

      // Cut on number of jets (0 for now)
      if(goodJets.size()==nJets) numPassNJetsCut[flavor]+=totalWeight; else continue;

      // Z mass window
      double dileptonMass = dileptonSystem.M();
      if(TMath::Abs(dileptonMass - 91.1876) < 15.) numPassZMassWindowCut[flavor]+=totalWeight; else continue;
      
      // 3rd lepton veto
      if((int)goodOrFakeLeptons.size() == 2 && totalLeptonCharge == 0) numPass3rdLeptonVeto[flavor]+=totalWeight; else { if(eventEvent.eventNum==77208166) printf("event 77208166 killed by 3rd lepton veto\n"); continue; }
      
      // b-jet and soft muon veto
      if(maxBDiscriminatorValue < .8) numPassBJetVeto[flavor]+=totalWeight; //else { if(eventEvent.eventNum==77208166) printf("event 77208166 killed by bjet veto (max discriminator value %f)\n\n", maxBDiscriminatorValue); continue; }
      else continue;
      if(softLeptons.size() == 0) numPassSoftMuonVeto[flavor]+=totalWeight; //else { if(eventEvent.eventNum==77208166) printf("event 77208166 killed by soft muon veto\n"); continue; }
      else continue;

      // delta phi cut between dilepton system and the MET
      double deltaPhiDileptonMet = TMath::Abs(dileptonSystem.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      if(deltaPhiDileptonMet > 2.8) numPassDeltaPhiDileptonMetCut[flavor]+=totalWeight; else continue;
      
      // delta phi cut between the two leptons
      double deltaPhiLeptons = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[lepton1])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[lepton2]));
      if(deltaPhiLeptons < TMath::Pi()/2.) numPassDeltaPhiLeptonsCut[flavor]+=totalWeight; else continue;

      // MET balance cut
      double metBalance = TMath::Abs(dileptonSystem.Pt()-((TLorentzVector*)(*eventMet.p4)[0])->Pt() ) / dileptonSystem.Pt();
      if(metBalance < 0.4) numPassMetBalanceCut[flavor]+=totalWeight; else continue;

      // transverse mass cut
      double transverseMass = TMath::Sqrt(2.0*dileptonSystem.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - TMath::Cos(deltaPhiDileptonMet)));
      if(transverseMass > 200.) numPassTransverseMassCut[flavor]+=totalWeight; else continue;

      // delta phi cut between any jets and the MET (not implemented yet)
      if(true) numPassDeltaPhiJetMetCut[flavor]+=totalWeight; else continue;

      // veto on taus
      int numberOfGoodTaus=0;
      for(int ntau=0; ntau < (int)eventTaus.p4->GetEntriesFast(); ntau++) {
        if(((TLorentzVector*)(*eventTaus.p4)[ntau])->Pt() <= 20.0 ||
           TMath::Abs(((TLorentzVector*)(*eventTaus.p4)[ntau])->Eta()) >= 2.3) continue;
        bool isLeptonicTau = false;
        for(int igl1=0; igl1<(int)goodOrFakeLeptons.size(); igl1++) {
          int il1 = goodOrFakeLeptons[igl1];
          if(((TLorentzVector*)(*eventLeptons.p4)[il1])->DeltaR(*((TLorentzVector*)(*eventTaus.p4)[ntau])) < 0.1) {
            isLeptonicTau = true;
            break;
          }
        }
        if(isLeptonicTau == false &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFinding  ) == BareTaus::TauDecayModeFinding &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFindingNewDMs) == BareTaus::TauDecayModeFindingNewDMs &&
           (double)(*eventTaus.iso)[ntau] < 4.0) numberOfGoodTaus++;
      }
      if(numberOfGoodTaus==0) numPassTauVeto[flavor]+=totalWeight; else continue;

      passEventSelection=true;
      if((flavor==1 || flavor==2) && infilecatv[ifile]==2) { 

        /*printf("event # %lld \n\tlepton1: pT=%f, eta=%f, phi=%f, selBits=%d\n\tlepton2: pT=%f, eta=%f, phi=%f, selBits=%d\n\tptll = %f, MET = %f\n\t%d jets\n\tmll = %f\n\tmax B discriminator value = %f\n\tdPhi(ll,MET) = %f, dPhi(l1,l2) = %f\n\t|ptll - MET|/ptll = %f, transverse mass = %f\n\n",
          eventEvent.eventNum,
          ((TLorentzVector*)(*eventLeptons.p4)[lepton1])->Pt(),
          ((TLorentzVector*)(*eventLeptons.p4)[lepton1])->Eta(),
          ((TLorentzVector*)(*eventLeptons.p4)[lepton1])->Phi(),
          (int)(*eventLeptons.selBits)[lepton1],
          ((TLorentzVector*)(*eventLeptons.p4)[lepton2])->Pt(),
          ((TLorentzVector*)(*eventLeptons.p4)[lepton2])->Eta(),
          ((TLorentzVector*)(*eventLeptons.p4)[lepton2])->Phi(),
          (int)(*eventLeptons.selBits)[lepton2],
          dileptonSystem.Pt(),
          met,
          (int)goodJets.size(),
          dileptonMass,
          maxBDiscriminatorValue,
          deltaPhiDileptonMet,
          deltaPhiLeptons,
          metBalance,
          transverseMass);*/
        eventNums.push_back(eventEvent.eventNum);
      }

      if(passEventSelection) { for(int il1=0; il1 < TMath::Max(lepton1, lepton2); il1++) {
          if(il1==lepton1 || il1==lepton2) continue;
          if(
            !selectIdIsoCut(
              typeLepSel.Data(),
              TMath::Abs((int)(*eventLeptons.pdgId)[il1]),
              ((TLorentzVector*)(*eventLeptons.p4)[il1])->Pt(),
              TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[il1])->Eta()),
              (double)(*eventLeptons.iso)[il1],
              (int)(*eventLeptons.selBits)[il1],
              (double)(*eventLeptons.mva)[il1]
            )) badBoostedNonZLeptons[flavor]+=totalWeight;
        //else goodBoostedNonZLeptons[flavor]+=totalWeight;
      }}

      if(passEventSelection) numPassEventSelection[flavor]+=totalWeight;
    }
    numTwoGoodLeptons[3] = numTwoGoodLeptons[1] + numTwoGoodLeptons[2]; numTwoGoodLeptons[4] = numTwoGoodLeptons[0] + numTwoGoodLeptons[3];
    numPassDileptonPtCut[3] = numPassDileptonPtCut[1] + numPassDileptonPtCut[2]; numPassDileptonPtCut[4] = numPassDileptonPtCut[0] + numPassDileptonPtCut[3];
    numPassMetPreselection[3] = numPassMetPreselection[1] + numPassMetPreselection[2]; numPassMetPreselection[4] = numPassMetPreselection[0] + numPassMetPreselection[3];
    numPassPreselection[3] = numPassPreselection[1] + numPassPreselection[2]; numPassPreselection[4] = numPassPreselection[0] + numPassPreselection[3];
    numPassMetSelection[3] = numPassMetSelection[1] + numPassMetSelection[2]; numPassMetSelection[4] = numPassMetSelection[0] + numPassMetSelection[3];
    numPassNJetsCut[3] = numPassNJetsCut[1] + numPassNJetsCut[2]; numPassNJetsCut[4] = numPassNJetsCut[0] + numPassNJetsCut[3];
    numPassZMassWindowCut[3] = numPassZMassWindowCut[1] + numPassZMassWindowCut[2]; numPassZMassWindowCut[4] = numPassZMassWindowCut[0] + numPassZMassWindowCut[3];
    numPass3rdLeptonVeto[3] = numPass3rdLeptonVeto[1] + numPass3rdLeptonVeto[2]; numPass3rdLeptonVeto[4] = numPass3rdLeptonVeto[0] + numPass3rdLeptonVeto[3];
    numPassEventSelection[3] = numPassEventSelection[1] + numPassEventSelection[2]; numPassEventSelection[4] = numPassEventSelection[0] + numPassEventSelection[3];
    numPassBJetVeto[3] = numPassBJetVeto[1] + numPassBJetVeto[2]; numPassBJetVeto[4] = numPassBJetVeto[0] + numPassBJetVeto[3];
    numPassSoftMuonVeto[3] = numPassSoftMuonVeto[1] + numPassSoftMuonVeto[2]; numPassSoftMuonVeto[4] = numPassSoftMuonVeto[0] + numPassSoftMuonVeto[3];
    numPassDeltaPhiDileptonMetCut[3] = numPassDeltaPhiDileptonMetCut[1] + numPassDeltaPhiDileptonMetCut[2]; numPassDeltaPhiDileptonMetCut[4] = numPassDeltaPhiDileptonMetCut[0] + numPassDeltaPhiDileptonMetCut[3];
    numPassDeltaPhiLeptonsCut[3] = numPassDeltaPhiLeptonsCut[1] + numPassDeltaPhiLeptonsCut[2]; numPassDeltaPhiLeptonsCut[4] = numPassDeltaPhiLeptonsCut[0] + numPassDeltaPhiLeptonsCut[3];
    numPassMetBalanceCut[3] = numPassMetBalanceCut[1] + numPassMetBalanceCut[2]; numPassMetBalanceCut[4] = numPassMetBalanceCut[0] + numPassMetBalanceCut[3];
    numPassTransverseMassCut[3] = numPassTransverseMassCut[1] + numPassTransverseMassCut[2]; numPassTransverseMassCut[4] = numPassTransverseMassCut[0] + numPassTransverseMassCut[3];
    numPassDeltaPhiJetMetCut[3] = numPassDeltaPhiJetMetCut[1] + numPassDeltaPhiJetMetCut[2]; numPassDeltaPhiJetMetCut[4] = numPassDeltaPhiJetMetCut[0] + numPassDeltaPhiJetMetCut[3];
    numPassTauVeto[3] = numPassTauVeto[1] + numPassTauVeto[2]; numPassTauVeto[4] = numPassTauVeto[0] + numPassTauVeto[3];
    badBoostedNonZLeptons[3] = badBoostedNonZLeptons[1] + badBoostedNonZLeptons[2]; badBoostedNonZLeptons[4] = badBoostedNonZLeptons[0] + badBoostedNonZLeptons[3];
    goodBoostedNonZLeptons[3] = goodBoostedNonZLeptons[1] + goodBoostedNonZLeptons[2]; goodBoostedNonZLeptons[4] = goodBoostedNonZLeptons[0] + goodBoostedNonZLeptons[3];
    printf("\nTable of Event Yields\n"); printf("----------------------------------------------------------------------------------------------------------------------\n");
    printf("\t\t\t%15.15s %15.15s %15.15s %15.15s %15.15s \n", "em", "mm", "ee", "ee+mm", "all");
    printf("%23.23s %15.15s %15.15s %15.15s %15.15s %15.2f \n", "2+ lep.cands.", "-", "-", "-", "-", numFoundLeptonCandidates);
    printf("%23.23s %15.15s %15.15s %15.15s %15.15s %15.2f \n", "trigger", "-", "-", "-", "-", numPassTrigger);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "l1,l2 pt>20 &no fakes", numTwoGoodLeptons[0], numTwoGoodLeptons[1], numTwoGoodLeptons[2], numTwoGoodLeptons[3], numTwoGoodLeptons[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "pre: ptll>60", numPassDileptonPtCut[0], numPassDileptonPtCut[1], numPassDileptonPtCut[2], numPassDileptonPtCut[3], numPassDileptonPtCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "pre: MET>45", numPassMetPreselection[0], numPassMetPreselection[1], numPassMetPreselection[2], numPassMetPreselection[3], numPassMetPreselection[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "MET>100", numPassMetSelection[0], numPassMetSelection[1], numPassMetSelection[2], numPassMetSelection[3], numPassMetSelection[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "# of jets", numPassNJetsCut[0], numPassNJetsCut[1], numPassNJetsCut[2], numPassNJetsCut[3], numPassNJetsCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "|mll - mZ| < 15", numPassZMassWindowCut[0], numPassZMassWindowCut[1], numPassZMassWindowCut[2], numPassZMassWindowCut[3], numPassZMassWindowCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "Q=0 and 2 leptons", numPass3rdLeptonVeto[0], numPass3rdLeptonVeto[1], numPass3rdLeptonVeto[2], numPass3rdLeptonVeto[3], numPass3rdLeptonVeto[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "b jet veto", numPassBJetVeto[0], numPassBJetVeto[1], numPassBJetVeto[2], numPassBJetVeto[3], numPassBJetVeto[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "soft muon veto", numPassSoftMuonVeto[0], numPassSoftMuonVeto[1], numPassSoftMuonVeto[2], numPassSoftMuonVeto[3], numPassSoftMuonVeto[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "dPhi(ll,MET) > 2.8", numPassDeltaPhiDileptonMetCut[0], numPassDeltaPhiDileptonMetCut[1], numPassDeltaPhiDileptonMetCut[2], numPassDeltaPhiDileptonMetCut[3], numPassDeltaPhiDileptonMetCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "dPhi(l1,l2) < pi/2", numPassDeltaPhiLeptonsCut[0], numPassDeltaPhiLeptonsCut[1], numPassDeltaPhiLeptonsCut[2], numPassDeltaPhiLeptonsCut[3], numPassDeltaPhiLeptonsCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "|ptll - MET|/ptll < 0.4", numPassMetBalanceCut[0], numPassMetBalanceCut[1], numPassMetBalanceCut[2], numPassMetBalanceCut[3], numPassMetBalanceCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "transverse mass > 200", numPassTransverseMassCut[0], numPassTransverseMassCut[1], numPassTransverseMassCut[2], numPassTransverseMassCut[3], numPassTransverseMassCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "dPhi(jet_i, MET)>0.5", numPassDeltaPhiJetMetCut[0], numPassDeltaPhiJetMetCut[1], numPassDeltaPhiJetMetCut[2], numPassDeltaPhiJetMetCut[3], numPassDeltaPhiJetMetCut[4]);
    printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "tau veto", numPassTauVeto[0], numPassTauVeto[1], numPassTauVeto[2], numPassTauVeto[3], numPassTauVeto[4]);
    printf("\nDebugging Yields\n"); printf("----------------------------------------------------------------------------------------------------------------------\n");
    //printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "bad boosted nonZ lep", badBoostedNonZLeptons[0], badBoostedNonZLeptons[1], badBoostedNonZLeptons[2], badBoostedNonZLeptons[3], badBoostedNonZLeptons[4]);
    //printf("%23.23s %15.2f %15.2f %15.2f %15.2f %15.2f \n", "nonZ lep[0] good", goodBoostedNonZLeptons[0], goodBoostedNonZLeptons[1], goodBoostedNonZLeptons[2], goodBoostedNonZLeptons[3], goodBoostedNonZLeptons[4]);
    //printf("----------------------------------------------------------------------------------------------------------------------\n"); 
    //printf("Selected event numbers:\n");
    //for(int iEvent=0; iEvent < (int) eventNums.size(); iEvent++) {
    //  printf("\t%lld\n", eventNums[iEvent]);
    //}
    //printf("----------------------------------------------------------------------------------------------------------------------\n"); 
  }
}
