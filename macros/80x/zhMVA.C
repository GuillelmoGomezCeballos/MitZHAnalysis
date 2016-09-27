#include <TMVA/Factory.h>
#include <TMVA/Types.h>
#include <TFile.h>
#include <TCut.h>
#include <TTree.h>

void zhMVA(string signal_model) {
  TFile *mva_input_trees = TFile::Open("MitZHAnalysis/mva/mva_input_trees.root", "READ");
  TTree *Zjets_mva_tree    = (TTree*) mva_input_trees->Get("bkg_mva_tree_Zjets");
  TTree *EM_mva_tree       = (TTree*) mva_input_trees->Get("bkg_mva_tree_EM"); 
  TTree *WZ_mva_tree       = (TTree*) mva_input_trees->Get("bkg_mva_tree_WZ");
  TTree *ZZ_mva_tree       = (TTree*) mva_input_trees->Get("bkg_mva_tree_ZZ");
  TTree *VVV_mva_tree      = (TTree*) mva_input_trees->Get("bkg_mva_tree_VVV");
  char signal_tree_name[512];
  sprintf(signal_tree_name, "signal_mva_tree_%s", signal_model.c_str());
  TTree *signal_mva_tree  = (TTree*) mva_input_trees->Get(signal_tree_name);
  
  TFile *output_file=TFile::Open(("MitZHAnalysis/mva/training_result_"+signal_model+".root").c_str(), "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("bdt", output_file, "AnalysisType=Classification");
  factory->AddSignalTree(signal_mva_tree ,    1);
  factory->AddBackgroundTree(Zjets_mva_tree , 1);
  factory->AddBackgroundTree(EM_mva_tree    , 1);
  factory->AddBackgroundTree(WZ_mva_tree    , 1);
  factory->AddBackgroundTree(ZZ_mva_tree    , 1);
  factory->AddBackgroundTree(VVV_mva_tree   , 1);
  
  factory->AddVariable( "mva_njets"             , "n_{jets}"                          , ""    , 'I');
  factory->AddVariable( "mva_ntaus"             , "n_{taus}"                          , ""    , 'I');
  factory->AddVariable( "mva_btag_veto"         , "b-tag veto (0|1)"                  , ""    , 'I');
  factory->AddVariable( "mva_balance"           , "Balance"                           , ""    , 'F');
  factory->AddVariable( "mva_cos_theta_star_l1" , "cos #theta^{*}_{l1}"               , ""    , 'F');
  factory->AddVariable( "mva_cos_theta_star_l2" , "cos #theta^{*}_{l2}"               , ""    , 'F');
  factory->AddVariable( "mva_cos_theta_CS_l1"   , "cos #theta^{CS}_{l1}"              , ""    , 'F');
  factory->AddVariable( "mva_cos_theta_CS_l2"   , "cos #theta^{CS}_{l2}"              , ""    , 'F');
  factory->AddVariable( "mva_delphi_ptll_MET"   , "#Delta#Phi(p^{ll}, p^{miss})"      , "rad" , 'F');
  factory->AddVariable( "mva_delphi_ll"         , "#Delta#Phi^{ll}"                   , "rad" , 'F');
  factory->AddVariable( "mva_delphi_jet_MET"    , "#Delta#Phi(jet, p^{miss})"         , "rad" , 'F');
  factory->AddVariable( "mva_deltaR_ll"         , "#DeltaR(l1, l2)"                   , ""    , 'F');
  factory->AddVariable( "mva_etall"             , "#eta_{ll}"                         , ""    , 'F');
  factory->AddVariable( "mva_etal1"             , "#eta_{l1}"                         , ""    , 'F');
  factory->AddVariable( "mva_etal2"             , "#eta_{l2}"                         , ""    , 'F');
  factory->AddVariable( "mva_MET"               , "E_{T}^{miss}"                      , "GeV" , 'F');
  factory->AddVariable( "mva_mll_minus_mZ"      , "|m_{ll} - m_{Z}|"                  , "GeV" , 'F');
  factory->AddVariable( "mva_mTjetMET"          , "m_{T}(jet, p_{T}^{miss})"          , "GeV" , 'F');
  factory->AddVariable( "mva_mTll"              , "m_{T}(p_{T}^{ll}, p^{miss})"       , "GeV" , 'F');
  factory->AddVariable( "mva_mTl1MET"           , "m_{T}(p_{T}^{l1}, p_{T}^{miss})"   , "GeV" , 'F');
  factory->AddVariable( "mva_mTl2MET"           , "m_{T}(p_{T}^{l2}, p_{T}^{miss})"   , "GeV" , 'F');
  factory->AddVariable( "mva_ptll"              , "p_{T}^{ll}"                        , "GeV" , 'F');
  factory->AddVariable( "mva_ptl1"              , "p_{T}^{l1}"                        , "GeV" , 'F');
  factory->AddVariable( "mva_ptl2"              , "p_{T}^{l2}"                        , "GeV" , 'F');
  factory->AddVariable( "ptl1mptl2_over_ptll"   , "Lepton Balance"                    , ""    , 'F');
  //factory->AddVariable( "mva_3lveto"          , "3rd lep. veto (0|1)", "", 'I');

  factory->SetWeightExpression( "mva_weight");
  TCut preselectionCut= "mva_MET>80 && mva_balance<5 && mva_3lveto==1 && mva_mll_minus_mZ <= 30";
  factory->PrepareTrainingAndTestTree(preselectionCut, "");
  factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=7:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();



}

