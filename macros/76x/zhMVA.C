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
  
  factory->AddVariable( "mva_balance"         , "#zeta_{ll}", "", 'F');
  factory->AddVariable( "mva_delphi_ptll_MET" , "#Delta#Phi(p^{ll}, p^{miss})", "rad", 'F');
  factory->AddVariable( "mva_delphi_ll"       , "#Delta#Phi^{ll}", "rad", 'F');
  factory->AddVariable( "mva_delphi_jet_MET"  , "#Delta#Phi(jet, p^{miss})", "rad", 'F');
  factory->AddVariable( "mva_MET"             , "E_T^{miss}", "GeV", 'F');
  factory->AddVariable( "mva_mll"             , "m_{ll}", "GeV", 'F');
  factory->AddVariable( "mva_mTll"            , "m_{T} (p^{ll}, p^{miss})", "GeV", 'F');
  factory->AddVariable( "mva_mTl1MET"         , "m_{T} (p^{l_{1}}, p^{miss})", "GeV", 'F');
  factory->AddVariable( "mva_mTl2MET"         , "m_{T} (p^{l_{2}}, p^{miss})", "GeV", 'F');
  factory->AddVariable( "mva_ptll"            , "p_{T}^{ll}", "GeV", 'F');
  factory->AddVariable( "mva_ptl1"            , "p_{T}^{l_{1}}", "GeV", 'F');
  factory->AddVariable( "mva_ptl2"            , "p_{T}^{l_{2}}", "GeV", 'F');
  factory->AddVariable( "mva_response"        , "Response", "", 'F');
  factory->AddVariable( "mva_njets"           , "n_{jets}", "", 'I');
  factory->AddVariable( "mva_ntaus"           , "n_{taus}", "", 'I');
  factory->AddVariable( "mva_btag"            , "b-tag (0|1)", "", 'I');

  factory->SetWeightExpression( "mva_weight");
  TCut preselectionCut= "mva_MET>80 && mva_balance<5";
  factory->PrepareTrainingAndTestTree(preselectionCut, "");
  factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();



}

