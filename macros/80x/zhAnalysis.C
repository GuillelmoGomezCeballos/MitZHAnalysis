#include "MitAnalysisRunII/macros/LeptonScaleLookup.cc"
#include "MitZHAnalysis/macros/80x/zhAnalysis.h"

// zhAnalysis Constructor
// sets up the subdirectory output and the random toy seed
zhAnalysis::zhAnalysis(string subdirectory_) {
  subdirectory=subdirectory_;
  if(subdirectory!="" && subdirectory.c_str()[0]!='/') subdirectory = "/"+subdirectory;
  system(("mkdir -p MitZHAnalysis/datacards"+subdirectory).c_str());
  system(("mkdir -p MitZHAnalysis/plots"+subdirectory).c_str());
  // Initialize the random seed based on Dylan's age in seconds
  std::time_t t = std::time(0);
  unsigned long int time_now = static_cast<unsigned long int>(time(NULL));
  randomToySeed=(time_now-731178000);
  if     (MVAVarType == 0) zjetsTemplatesPath = "MitZHAnalysis/data/76x/zjets_13TeV_25ns_metgt50_mt.root";
  else if(MVAVarType == 1) zjetsTemplatesPath = "MitZHAnalysis/data/80x/zjets_13TeV_25ns_metgt50_met.root";
  else useZjetsTemplate = false;
  return;                                                                                                             
}

zhAnalysis::~zhAnalysis() {}

void zhAnalysis::Run(
  bool verbose_/*=false*/, int nJetsType_/*=1*/, int typeSel_/*=3*/, int plotModel_/*=0*/
) {
  printf("zhAnalysis::Run : Begin running the whole analysis\n");
  unsigned long int t1 = static_cast<unsigned long int>(time(NULL));
  
  verbose=verbose_;
  nJetsType=nJetsType_;
  typeSel=typeSel_;
  plotModel=plotModel_;

  if(nInputFiles==0) assert(LoadFlatFiles());
  assert(MakeHistos());
  if(nInputFiles==0) { printf("Error in zhAnalysis::Run: no input files loaded (problem with zhAnalysis::LoadFlatFiles)\n"); return; }
  SetFinalStateName();

  TString processTag = "";

  TH1D *fhDZjets;
  TH1D *fhDZjetsSyst;
  if(useZjetsTemplate){
    TFile *fZjetsTemplatesFile = TFile::Open(Form("%s",zjetsTemplatesPath.Data()));

    if     (nJetsType == 0) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets0"));
    else if(nJetsType == 1) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
    else if(nJetsType == 2) fhDZjets = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets2"));
    assert(fhDZjets);
    fhDZjets->SetDirectory(0);

    if     (nJetsType == 0) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
    else if(nJetsType == 1) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets0"));
    else if(nJetsType == 2) fhDZjetsSyst = (TH1D*)(fZjetsTemplatesFile->Get("histo_Zjets1"));
    assert(fhDZjetsSyst);
    fhDZjetsSyst->SetDirectory(0);

    delete fZjetsTemplatesFile;
  }

  histoZHSEL[0] = new TH1D("histoZHSEL_0", "histoZHSEL_0", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[1] = new TH1D("histoZHSEL_1", "histoZHSEL_1", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[2] = new TH1D("histoZHSEL_2", "histoZHSEL_2", numberCuts+1, -0.5, numberCuts+0.5);
  histoZHSEL[3] = new TH1D("histoZHSEL_3", "histoZHSEL_3", numberCuts+1, -0.5, numberCuts+0.5);
  TString cutName[numberCuts+1] = {"ptl>20/20","3rd lepton veto","btag-veto","tauVeto","Njets","Z mass","ptll>60","MET>100","dPhi(Z-MET)>2.8","|ptll-MET|/ptll<0.4","dPhiJetMet>0.5","all"};


  if (MVAVarType==3 || MVAVarType==4) useBDT=true;

  double bgdDecay[nSigModels][nSelTypes*4][processTypes],weiDecay[nSigModels][nSelTypes*4][processTypes];
  for(int nModel=0; nModel<nSigModels; nModel++) { for(unsigned int i=0; i<nSelTypes*4; i++) { for(int j=0; j<processTypes; j++) {       
    bgdDecay[nModel][i][j] = 0.0; weiDecay[nModel][i][j] = 0.0; 
  }}}
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
           mva_weight;//,
           //aux_PDFscale[102],
           //aux_PUscale,
           //aux_QCDscale_r1f2,
           //aux_QCDscale_r1f5,
           //aux_QCDscale_r2f1,
           //aux_QCDscale_r2f2,
           //aux_QCDscale_r5f1,
           //aux_QCDscale_r5f5,
           //aux_MET_JESup,
           //aux_MET_JESdown;
  UChar_t  mva_njets,
           mva_ntaus,
           aux_njets_JESup,
           aux_njets_JESdown;
  Bool_t   mva_btag_veto,
           mva_3lveto;
  if(useBDT) {
    reader=new TMVA::Reader();
    if(MVAVarType==3) {
      //reader->AddVariable( "mva_balance"                       , &mva_balance            );
      //reader->AddVariable( "mva_cos_theta_star_l1"             , &mva_cos_theta_star_l1  );
      reader->AddVariable( "TMath::Abs(mva_cos_theta_CS_l1)"   , &mva_cos_theta_CS_l1    );
      reader->AddVariable( "mva_delphi_ptll_MET"               , &mva_delphi_ptll_MET    );
      //reader->AddVariable( "mva_delphi_ll"                     , &mva_delphi_ll          );
      reader->AddVariable( "mva_deltaR_ll"                     , &mva_deltaR_ll          );
      //reader->AddVariable( "TMath::Abs(mva_etall)"             , &mva_etall              );
      reader->AddVariable( "TMath::Abs(mva_etal1)"             , &mva_etal1              );
      reader->AddVariable( "TMath::Abs(mva_etal2)"             , &mva_etal2              );
      reader->AddVariable( "mva_MET"                           , &mva_MET                );
      reader->AddVariable( "mva_mll_minus_mZ"                  , &mva_mll_minus_mZ       );
      //reader->AddVariable( "mva_mTll"                          , &mva_mTll               );
      reader->AddVariable( "mva_mTl1MET"                       , &mva_mTl1MET            );
      reader->AddVariable( "mva_mTl2MET"                       , &mva_mTl2MET            );
      reader->AddVariable( "mva_ptll"                          , &mva_ptll               );
      reader->AddVariable( "mva_ptl1"                          , &mva_ptl1               );
      reader->AddVariable( "mva_ptl2"                          , &mva_ptl2               );
      //reader->AddVariable( "mva_ptl1mptl2_over_ptll"           , &mva_ptl1mptl2_over_ptll);
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
  
  // Either set up the random BDT toys or load the cached systematics here
  if (useBDT) { 
    // Filename for the cached bdt systematics, we will either read them from here or write fresh ones to here 
    sprintf(filenameBDTSysts,"MitZHAnalysis/plots%s/zll%shinv%s_%s_BDTsyst_%s.root", subdirectory.c_str(), addChan.Data(), finalStateName, signalName_[plotModel].Data(), ECMsb.Data());
    if(useCachedBDTSystematics) LoadCachedBDTSystematics();
    else                        SetupBDTSystematics();
  }

  unsigned int numberOfLeptons = 2;
  TString signalName="";
  double totalEventsProcess[50];

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(unsigned ifile=0; ifile<nInputFiles; ifile++) {
    FlatFile inputFlatFile = inputFlatFiles[ifile];
    printf("zhAnalysis::Run : Reading sample #%d (%s)\n",ifile, inputFlatFile.name.Data());
    TFile *inputTFile = inputFlatFile.Open();
    int nModel = (inputFlatFile.category==6 || inputFlatFile.category==7) ? inputFlatFile.signalIndex:-1;
    if(nModel>=0) signalName=signalName_[nModel];
    if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
    TTree *inputTree = (TTree*)inputTFile->FindObjectAny("events");
    bool isData=(inputFlatFile.category==0);
    SetBranchAddresses(inputTree,isData);

    int initPDFTag = -1;
    bool errorMsgQCDscale = false;
    histoZHSEL[0]->Scale(0.0);
    histoZHSEL[1]->Scale(0.0);
    histoZHSEL[2]->Scale(0.0);
    histoZHSEL[3]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(inputFlatFile.category == 0) theMCPrescale = 1.0;
    // Catch all of the possible ways that someone could name a glu-glu ntuple
    bool isGluonInduced=
      inputFlatFile.name.Contains("gg") ||
      inputFlatFile.name.Contains("GG") ||
      inputFlatFile.name.Contains("GluGlu") ||
      inputFlatFile.name.Contains("Cont");

    for (int i=0; i<int(inputTree->GetEntries()/theMCPrescale); ++i) {
      if(i%1000000==0) printf("zhAnalysis::Run : Processed events %d/%d for sample #%d\n",i,(int)inputTree->GetEntries(),ifile);
      inputTree->GetEntry(i);

      Bool_t passFilter[4] = {kFALSE,kFALSE,kFALSE,kFALSE};
      if(gltEvent.looseLep1Pt>20 && gltEvent.looseLep2Pt>10) passFilter[0] = kTRUE;
      if(passFilter[0] == kFALSE) continue;
      
      // To do: get the trigger bits and selection bits defined outside of zhAnalysis and PandaLeptonicAnalyzer namespaces? ~DGH
      passFilter[1] = (!isData) || ( 
        (gltEvent.trigger & zhAnalysis::TriggerBits::kEGTrig  ) != 0 ||
        (gltEvent.trigger & zhAnalysis::TriggerBits::kMuTrig  ) != 0 ||
        (gltEvent.trigger & zhAnalysis::TriggerBits::kEGEGTrig) != 0 ||
        (gltEvent.trigger & zhAnalysis::TriggerBits::kMuMuTrig) != 0 ||
        (gltEvent.trigger & zhAnalysis::TriggerBits::kMuEGTrig) != 0
      ); // pass filter if it's a MC file or if it passes the trigger soup
      //if(passFilter[1] == kFALSE) continue;
      // Begin the offline leptonic selection 
      vector<int> idLep;
      vector<bool> idTight, idSoft; unsigned int goodIsTight = 0;
      vector<float*> idLepPts, idLepEtas, idLepPhis;            // Vectors of addresses to the kinematics of chosen leptons
      vector<float*> idLepIdSfs, idLepTrkSfs;                    // Vectors of addresses to the scale factors of chosen leptons
      vector<int*> idLepSelBits, idLepGenPdgIds, idLepPdgIds; // Vectors of addresses to the integer properties of chosen leptons


      // Implement the (flavor -> selection) correspondence as a map for now
      // Currently the ID working point is hardcoded, maybe there is a better solution ~DGH
      std::map<unsigned, unsigned> fsMap;
      fsMap[-11] = zhAnalysis::SelectionBit::kMedium;
      fsMap[11] = zhAnalysis::SelectionBit::kMedium;
      fsMap[-13] = zhAnalysis::SelectionBit::kTight;
      fsMap[13] = zhAnalysis::SelectionBit::kTight;

      // Create the lepton container. Store "true" in idTight for the leptons that pass the selection.
      // Note: idLep is currently redundant because the minimum fakeable object definition is equivalent to that of PandaLeptonicAnalyzer ~DGH
      // It might be faster to allocate the size of these vectors once based on gltEvent.looseLep*Pt, then just set using [] ~DGH
      { bool isTight;
        // Not storing soft muons for now, need to support it in PandaLeptonicAnalyzer if we want to do it here! ~DGH
        if(gltEvent.looseLep1Pt>=10 && fsMap.find(gltEvent.looseLep1PdgId)!=fsMap.end()) {
          isTight=((gltEvent.looseLep1SelBit & fsMap[gltEvent.looseLep1PdgId])!=0); idTight.push_back(isTight); goodIsTight+=isTight;
          idLepPts.push_back(&gltEvent.looseLep1Pt); idLepEtas.push_back(&gltEvent.looseLep1Eta); idLepPhis.push_back(&gltEvent.looseLep1Phi); idLepPdgIds.push_back(&gltEvent.looseLep1PdgId); idLep.push_back((int)idLep.size());
          idLepIdSfs.push_back( (gltEvent.looseLep1PdgId==13 || gltEvent.looseLep1PdgId==-13)? &gltEvent.sf_tight1 : &gltEvent.sf_medium1 );
          idLepTrkSfs.push_back( &gltEvent.sf_trk1 );
        } if(gltEvent.looseLep2Pt>=10 && fsMap.find(gltEvent.looseLep2PdgId)!=fsMap.end()) {
          isTight=((gltEvent.looseLep2SelBit & fsMap[gltEvent.looseLep2PdgId])!=0); idTight.push_back(isTight); goodIsTight+=isTight;
          idLepPts.push_back(&gltEvent.looseLep2Pt); idLepEtas.push_back(&gltEvent.looseLep2Eta); idLepPhis.push_back(&gltEvent.looseLep2Phi); idLepPdgIds.push_back(&gltEvent.looseLep2PdgId); idLep.push_back((int)idLep.size());
          idLepIdSfs.push_back( (gltEvent.looseLep2PdgId==13 || gltEvent.looseLep2PdgId==-13)? &gltEvent.sf_tight2 : &gltEvent.sf_medium2 );
          idLepTrkSfs.push_back( &gltEvent.sf_trk2 );
        } if(gltEvent.looseLep3Pt>=10 && fsMap.find(gltEvent.looseLep3PdgId)!=fsMap.end()) {
          isTight=((gltEvent.looseLep3SelBit & fsMap[gltEvent.looseLep3PdgId])!=0); idTight.push_back(isTight); goodIsTight+=isTight;
          idLepPts.push_back(&gltEvent.looseLep3Pt); idLepEtas.push_back(&gltEvent.looseLep3Eta); idLepPhis.push_back(&gltEvent.looseLep3Phi); idLepPdgIds.push_back(&gltEvent.looseLep3PdgId); idLep.push_back((int)idLep.size());
          idLepIdSfs.push_back( (gltEvent.looseLep3PdgId==13 || gltEvent.looseLep3PdgId==-13)? &gltEvent.sf_tight3 : &gltEvent.sf_medium3 );
          idLepTrkSfs.push_back( &gltEvent.sf_trk3 );
        } if(gltEvent.looseLep4Pt>=10 && fsMap.find(gltEvent.looseLep4PdgId)!=fsMap.end()) {
          isTight=((gltEvent.looseLep4SelBit & fsMap[gltEvent.looseLep4PdgId])!=0); idTight.push_back(isTight); goodIsTight+=isTight;
          idLepPts.push_back(&gltEvent.looseLep4Pt); idLepEtas.push_back(&gltEvent.looseLep4Eta); idLepPhis.push_back(&gltEvent.looseLep4Phi); idLepPdgIds.push_back(&gltEvent.looseLep4PdgId); idLep.push_back((int)idLep.size());
          idLepIdSfs.push_back( (gltEvent.looseLep4PdgId==13 || gltEvent.looseLep4PdgId==-13)? &gltEvent.sf_tight4 : &gltEvent.sf_medium4 );
          idLepTrkSfs.push_back( &gltEvent.sf_trk4 );
        }          
      }
      if(idTight.size()>=numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue; 
      if(idTight.size()==goodIsTight) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(gltEvent.looseLep1Pt<=25 || gltEvent.looseLep2Pt<=20) continue;

      double dPhiLepMETMin = 999.;
      int signQ = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        signQ = signQ + (int)*idLepPdgIds[idLep[nl]]/TMath::Abs((int)*idLepPdgIds[idLep[nl]]);
        double dPhiLepMETOption = TMath::Abs(TVector2::Phi_mpi_pi( gltEvent.pfmetphi - *idLepPhis[idLep[nl]]));
        if(dPhiLepMETOption < dPhiLepMETMin) dPhiLepMETMin = dPhiLepMETOption;
      }
      double minMET  = TMath::Min(gltEvent.pfmet, gltEvent.trkmetphi);
      double minPMET = TMath::Min(gltEvent.pfmet, gltEvent.trkmetphi);
      if(dPhiLepMETMin < TMath::Pi()/2.) minPMET = minPMET * sin(dPhiLepMETMin);

      // Gen jet flavor tagging already handled by PandaLeptonicAnalyzer ~DGH
      // Photon-lepton cleaning handled already by PandaLeptonicAnalyzer ~DGH
      bool isGoodPhoton=(gltEvent.loosePho1Pt>30 && TMath::Abs(gltEvent.loosePho1Eta)<2.5);
      // The equivalent of Nero BarePhotons::PhoElectronVeto is the panda::Photon::csafeVeto
      // This is currently not implemented in PandaLeptonicAnalyzer, need to add it ~DGH
      // Old code:
      //vector<int> idPho;
      //  if(((int)(*eventPhotons.selBits)[npho] & BarePhotons::PhoTight)== BarePhotons::PhoTight &&
      //   ((int)(*eventPhotons.selBits)[npho] & BarePhotons::PhoElectronVeto)== BarePhotons::PhoElectronVeto){idPho.push_back(npho);}

      TLorentzVector idLep1P4, idLep2P4, metP4;
      // Use the precise values of the lepton masses for now, even though that is not consistent with CMSSW. ~DGH
      idLep1P4.SetPtEtaPhiM(*idLepPts[idLep[0]],*idLepEtas[idLep[0]],*idLepPhis[idLep[0]],(*idLepPdgIds[idLep[0]]==13||*idLepPdgIds[idLep[0]]==-13)? 0.105658 : 0.000511);
      idLep2P4.SetPtEtaPhiM(*idLepPts[idLep[1]],*idLepEtas[idLep[1]],*idLepPhis[idLep[1]],(*idLepPdgIds[idLep[1]]==13||*idLepPdgIds[idLep[1]]==-13)? 0.105658 : 0.000511);
      metP4.SetPtEtaPhiM(gltEvent.pfmet,0,gltEvent.pfmetphi,0);
      // Form dilepton and dilepton+MET system
      TLorentzVector dilep(idLep1P4+idLep2P4); TLorentzVector dilepMET(dilep+metP4);
      TLorentzVector dilepg(dilep); //dilepton+y system
      if(isGoodPhoton){
        TLorentzVector idPho1P4;
        idPho1P4.SetPtEtaPhiM(gltEvent.loosePho1Pt,gltEvent.loosePho1Eta,gltEvent.loosePho1Phi,0);
        dilepg = dilepg + idPho1P4;
      }

      // Begin the offline jet selection 
      vector<int> idJet, idJetUp, idJetDown, idBJet, idJetNoPh; 
      vector<float*> jetPts, jetEtas, jetPhis, jetPtsUp, jetPtsDown, jetEtasUp, jetEtasDown; // Vectors of addresses to the kinematics of jets
      vector<float*> jetBTags, jetGenPts;
      vector<int*> jetSelBits, jetFlavs; // Vectors of addresses to the integer properties of jets

      // Create the jet container.
      if(gltEvent.jet1Pt>=30) {
        jetPts.push_back(&gltEvent.jet1Pt); jetEtas.push_back(&gltEvent.jet1Eta); jetPhis.push_back(&gltEvent.jet1Phi); jetBTags.push_back(&gltEvent.jet1BTag); jetSelBits.push_back(&gltEvent.jet1SelBit);
        jetPtsUp.push_back(&gltEvent.jet1PtUp); jetPtsDown.push_back(&gltEvent.jet1PtDown); jetEtasUp.push_back(&gltEvent.jet1EtaUp); jetEtasDown.push_back(&gltEvent.jet1EtaDown);
        idJet.push_back((int)idJet.size());
        if(inputFlatFile.category != 0) { jetGenPts.push_back(&gltEvent.jet1GenPt); jetFlavs.push_back(&gltEvent.jet1Flav); }
      } if(gltEvent.jet2Pt>=30) {
        jetPts.push_back(&gltEvent.jet2Pt); jetEtas.push_back(&gltEvent.jet2Eta); jetPhis.push_back(&gltEvent.jet2Phi); jetBTags.push_back(&gltEvent.jet2BTag); jetSelBits.push_back(&gltEvent.jet2SelBit);
        jetPtsUp.push_back(&gltEvent.jet2PtUp); jetPtsDown.push_back(&gltEvent.jet2PtDown); jetEtasUp.push_back(&gltEvent.jet2EtaUp); jetEtasDown.push_back(&gltEvent.jet2EtaDown);
        idJet.push_back((int)idJet.size());
        if(inputFlatFile.category != 0) { jetGenPts.push_back(&gltEvent.jet2GenPt); jetFlavs.push_back(&gltEvent.jet2Flav); }
      } if(gltEvent.jet3Pt>=30) {
        jetPts.push_back(&gltEvent.jet3Pt); jetEtas.push_back(&gltEvent.jet3Eta); jetPhis.push_back(&gltEvent.jet3Phi); jetBTags.push_back(&gltEvent.jet3BTag); jetSelBits.push_back(&gltEvent.jet3SelBit);
        jetPtsUp.push_back(&gltEvent.jet3PtUp); jetPtsDown.push_back(&gltEvent.jet3PtDown); jetEtasUp.push_back(&gltEvent.jet3EtaUp); jetEtasDown.push_back(&gltEvent.jet3EtaDown);
        idJet.push_back((int)idJet.size());
        if(inputFlatFile.category != 0) { jetGenPts.push_back(&gltEvent.jet3GenPt); jetFlavs.push_back(&gltEvent.jet3Flav); }
      } if(gltEvent.jet4Pt>=30) {
        jetPts.push_back(&gltEvent.jet4Pt); jetEtas.push_back(&gltEvent.jet4Eta); jetPhis.push_back(&gltEvent.jet4Phi); jetBTags.push_back(&gltEvent.jet4BTag); jetSelBits.push_back(&gltEvent.jet4SelBit);
        jetPtsUp.push_back(&gltEvent.jet4PtUp); jetPtsDown.push_back(&gltEvent.jet4PtDown); jetEtasUp.push_back(&gltEvent.jet4EtaUp); jetEtasDown.push_back(&gltEvent.jet4EtaDown);
        idJet.push_back((int)idJet.size());
        if(inputFlatFile.category != 0) { jetGenPts.push_back(&gltEvent.jet4GenPt); jetFlavs.push_back(&gltEvent.jet4Flav); }
      }  
      // End create jet container

      bool isBtag = kFALSE;
      double sumPtJets = 0.0;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0; double dPhiJetNoPhMET = -1.0;
      double mTJetMET = -1;
      double dPhiJetDiLep = -1.0;
      TLorentzVector dilepJet = dilep;
      for(unsigned nj=0; nj<(unsigned)jetPts.size(); nj++){
        // B-tagging stuff removed, we are covered already by PandaLeptonicAnalyzer ~DGH
        // Removed photon-jet cleaning here, already covered by PandaLeptonicAnalyzer ~DGH

        // WARNING: Previously computed the b-veto using 20 GeV jets, but PandaLeptonicAnalyzer only considers 30 GeV+ jets! ~DGH
        if(*jetPts[nj]      > 20) {
          if( *jetBTags[nj] > bDiscrMax ) bDiscrMax = *jetBTags[nj];
          if( *jetBTags[nj] > 0.8       ) idBJet.push_back(nj);
        } if(*jetPts[nj]      > 30) {
          sumPtJets += *jetPts[nj];
          idJet.push_back(nj);
          // WARNING: jet energy/mass not available in PandaLeptonicAnalyzer, need to add it. ~DGH
          TLorentzVector jetP4; jetP4.SetPtEtaPhiM( *jetPts[nj], *jetEtas[nj], *jetPhis[nj], 0.0); dilepJet += jetP4;
        }
        if(*jetPtsUp[nj]   > 30.) idJetUp.push_back(nj);
        if(*jetPtsDown[nj] > 30.) idJetDown.push_back(nj);
      }

      // Removed finding "numberGoodTaus" section, as this is already done in PandaLeptonicAnalyzer and stored in gltEvent.nTau ~DGH
      
      TLorentzVector theCaloMET; theCaloMET.SetPtEtaPhiM(gltEvent.calomet, 0, gltEvent.calometphi, 0);
      for(unsigned int nl=0; nl<idLep.size(); nl++) { // Calo MET No Mu calculation
        if(TMath::Abs((int)*idLepPdgIds[idLep[nl]]) == 13) {
          // It's redundant to make this 4vector again after having made idLep1P4, idLep2P4, but the code is more portable to other selections this way. ~DGH
          TLorentzVector muP4; muP4.SetPtEtaPhiM(*idLepPts[idLep[nl]], *idLepEtas[idLep[nl]], *idLepPhis[idLep[nl]], 0.105658);
          theCaloMET.SetPx(theCaloMET.Px()-muP4.Px());
          theCaloMET.SetPy(theCaloMET.Py()-muP4.Py());
        }
      }

      // Determine flavor of the pair (0 means e-mu pair, 1 means mu-mu, 2 means e-e)
      int typePair = 0;
      if     (TMath::Abs((int)*idLepPdgIds[idLep[0]])==13 && TMath::Abs((int)*idLepPdgIds[idLep[1]])==13) typePair = 1;
      else if(TMath::Abs((int)*idLepPdgIds[idLep[0]])==11 && TMath::Abs((int)*idLepPdgIds[idLep[1]])==11) typePair = 2;

      // Calculate a lot of physics quantities used for the rectangular selection

      double dPhiDiLepMET = TMath::Abs(dilep.DeltaPhi(metP4));
      double ptFrac = TMath::Abs(dilep.Pt()-metP4.Pt())/dilep.Pt(); // TMath::Abs(dilepJet.Pt()-metP4.Pt())/dilepJet.Pt();
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*metP4.Pt()*(1.0 - TMath::Cos(dPhiDiLepMET)));

      double dPhiDiLepGMET = TMath::Abs(dilepg.DeltaPhi(metP4));
      double ptFracG = TMath::Abs(dilepg.Pt()-metP4.Pt())/dilepg.Pt();

      double caloMinusPFMETRel = TMath::Abs( gltEvent.calomet - gltEvent.pfmet ) / gltEvent.pfmet;
      
      TVector2 metv(metP4.Px(), metP4.Py());
      TVector2 dilv(dilep.Px(), dilep.Py());
      TVector2 utv = -1.*(metv+dilv);
      double phiv = utv.DeltaPhi(dilv);
      double the_upara = TMath::Abs(utv.Mod()*TMath::Cos(phiv))/dilep.Pt();
      
      bool passZMass = dilep.M() > 76.1876 && dilep.M() < 106.1876;
      bool passNjets = idJet.size() <= nJetsType;
      // idJetNoPh is meaningless/empty right now, because PandaLeptonicAnalyzer only stores a single photon. Need to check jet dR < 0.3 with all the photons if we want to use it. ~DGH
      bool passNjetsG = idJetNoPh.size() <= 1;

      double metMIN = 100; double mtMIN = 200; double metTIGHT = 100;
      double bdtMIN = xbins[2]; //boundary after the drell yan bin [-1, ?]
      
      bool passMETTight  = metP4.Pt() > metTIGHT;

      if(MVAVarType == 0)                                            {metMIN = 50; mtMIN = 200;}
      else if(MVAVarType == 1 || MVAVarType == 2 || MVAVarType == 4) {metMIN = 50; mtMIN = 0;}
      else if(MVAVarType ==3) { metMIN = 80; mtMIN = 0;}
      bool passMET = metP4.Pt() > metMIN;
      bool passMT = mtW > mtMIN || !passMETTight;
      if(inputFlatFile.category == 0 && isBlinded) passMET = passMET && metP4.Pt() < 100;

      bool passPTFracG    = ptFracG < 0.5;
      bool passDPhiZGMET  = dPhiDiLepGMET > 2.6;

      bool passPTFrac    = ptFrac < 0.4;
      bool passDPhiZMET  = dPhiDiLepMET > 2.6;
      bool passBtagVeto  = bDiscrMax < bTagCut;
      bool passPTLL      = dilep.Pt() > 60;
      bool pass3rdLVeto  = idLep.size() == numberOfLeptons && TMath::Abs(signQ) == 0;
      double dphill = TMath::Abs(idLep1P4.DeltaPhi(idLep2P4));
      double detall = TMath::Abs(idLep1P4.Eta()-idLep2P4.Eta());
      double drll = sqrt(dphill*dphill+detall*detall);
      bool passDelphiLL  = drll < 1.8;//dphill < TMath::Pi()/2.;

      bool passZMassLarge = TMath::Abs(dilep.M()-91.1876) < 30.0;
      bool passZMassSB    = (dilep.M() > 110.0 && dilep.M() < 200.0);

      bool passDPhiJetMET     = dPhiJetMET     == -1 || dPhiJetMET     >= 0.5;
      bool passDPhiJetNoPhMET = dPhiJetNoPhMET == -1 || dPhiJetNoPhMET >= 0.5;
      bool passTauVeto    = gltEvent.nTau == 0;

      bool passNMinusOne[11] = {
                  passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT &&              passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass &&              passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
                  passZMass && passNjets &&            passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
        passMT && passZMass && passNjets && passMET &&               passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac                 && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET                 && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto             &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto                 && passDPhiJetMET && passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL &&                   passTauVeto && passMETTight,
        passMT && passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET                && passMETTight
      };

      bool passAllCuts[nSelTypes] = {                   
                             passZMass && passNjets                                                                                && pass3rdLVeto                                                   ,         // ZSEL
                             passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,         // SIGSEL
        isGoodPhoton      && passZMass && passNjetsG&& passMETTight      && passPTFracG&& passDPhiZGMET&& passBtagVeto && passPTLL && pass3rdLVeto &&             passDPhiJetNoPhMET && passTauVeto,         // ZHGSEL
        passZMassSB       &&!passZMass && passNjets && passMET                                         &&!passBtagVeto             && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,         // WWLOOSESEL
                             passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET &&!passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,         // BTAGSEL
                             passZMass && passNjets && passMT && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL &&!pass3rdLVeto,                                                         // WZSEL
                             passZMass && passNjets && passMET                         && passDPhiZMET                 && passPTLL && pass3rdLVeto && metP4.Pt() < 100., // PRESEL
                             passZMass && passNjets && passMT && passMET &&!passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,         // CR1SEL
                             passZMass && passNjets && passMT && passMET && passPTFrac &&!passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,         // CR2SEL
                             passZMass && passNjets && passMT && passMET &&!passPTFrac &&!passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,         // CR12SEL
                             passZMass && passNjets && passMETTight      && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto, // TIGHTSEL
                             passZMass && passNjets && passMETTight      &&!passPTFrac                                 && passPTLL && pass3rdLVeto,                                                  // DYSANESEL1
                             passZMass && passNjets && passMETTight      &&!passDPhiZMET                               && passPTLL && pass3rdLVeto                                                   // DYSANESEL2
      };
     // Evaluate nominal BDT value
     double bdt_value=-1;
     TLorentzVector jet1P4(0,0,0,0);
     if(useBDT) {
       if(idJet.size()>0) jet1P4.SetPtEtaPhiM(*jetPts[0], *jetEtas[0], *jetPhis[0], 0); // Wrong jet energy, need to fix ~DGH
       TLorentzVector lepton1 = idLep1P4,
                      lepton2 = idLep2P4,
                      MET     = metP4,
                      jet1    = jet1P4;
        bdt_value = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0,0,0,0);
     }
     bool passBDT = bdt_value>bdtMIN;
     if(MVAVarType==3) {
       passAllCuts[WZSEL]  = passZMassLarge && passNjets && passMET && passBtagVeto  && !pass3rdLVeto;
       passAllCuts[PRESEL] = passZMassLarge && passNjets && passMET && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET;
       passAllCuts[SIGSEL] = (passAllCuts[PRESEL] && passBDT && passMETTight) || (passZMass && passNjets && passMT && passMET && !passMETTight && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL &&  pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto);
       passAllCuts[TIGHTSEL] = passAllCuts[PRESEL] && passBDT && passMETTight;
     }
     bool passEvolFilter[numberCuts] = {pass3rdLVeto,passBtagVeto,passTauVeto,passNjets,passZMass,passPTLL,passMETTight,passDPhiZMET,passPTFrac,passDPhiJetMET,passDelphiLL&&passMT};

     int sumEvol = 0;
     bool totalSel = kTRUE;
     for(int isel=0; isel<numberCuts; isel++) {
       totalSel = totalSel && passEvolFilter[isel];
       if(totalSel == kTRUE) sumEvol++;
     }
     double mtWSyst[2] = {TMath::Sqrt(2.0*dilep.Pt()*gltEvent.pfmetUp  *(1.0 - TMath::Cos(dPhiDiLepMET))),
                          TMath::Sqrt(2.0*dilep.Pt()*gltEvent.pfmetDown*(1.0 - TMath::Cos(dPhiDiLepMET)))};
     // Syst cuts for MVA
     bool passSystCuts[nSystTypes] = {
          passZMass && idJetUp.size() <= nJetsType  && passMET && passMT                                                                                                                                         && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && idJetDown.size()<= nJetsType && passMET && passMT                                                                                                                                         && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && passNjets && gltEvent.pfmetUp > metMIN && (mtWSyst[0] > mtMIN || gltEvent.pfmetUp < metTIGHT) && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto,
          passZMass && passNjets && gltEvent.pfmetDown > metMIN && (mtWSyst[1] > mtMIN || gltEvent.pfmetDown < metTIGHT) && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL && pass3rdLVeto && passDelphiLL && passDPhiJetMET && passTauVeto
     };
     if(MVAVarType==3) {
         // (Cuts in common between DY bin and signal region) && ((Exclusive BDT signal region cuts || Exclusive DY bin cuts))  
         passSystCuts[0] = (idJetUp.size()   <= nJetsType && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && passMETTight && passBDT) || (passZMass && passMET && !passMETTight && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
         passSystCuts[1] = (idJetDown.size() <= nJetsType && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && passMETTight && passBDT) || (passZMass && passMET && !passMETTight && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
         passSystCuts[2] = (passNjets && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && gltEvent.pfmetUp > metTIGHT && passBDT) || (passZMass && gltEvent.pfmetUp > metMIN && gltEvent.pfmetUp <= metTIGHT && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
         passSystCuts[3] = (passNjets && passBtagVeto && pass3rdLVeto && passTauVeto && passDPhiJetMET) && ((passZMassLarge && gltEvent.pfmetDown > metTIGHT && passBDT) || (passZMass && gltEvent.pfmetDown > metMIN && gltEvent.pfmetDown <= metTIGHT && passPTFrac && passDPhiZMET && passPTLL && passDelphiLL));
     }
      
      // begin event weighting
      

      // Deleted code for the generator leptons. There is some support in PLA using the generator dressed leptons, if we want to use those for the fake rate study. ~DGH

      // trigger efficiency
      double trigEff = 1.0;
      //if(inputFlatFile.category != 0) {
      //  trigEff = trigLookup.GetExpectedTriggerEfficiency(idLep1P4.Eta(),idLep1P4.Pt(),
      //                                                    idLep2P4.Eta(),idLep2P4.Pt(),
      //                                                   TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]));
      //}
      // luminosity
      double theLumi  = 1.0; if(inputFlatFile.category != 0) theLumi  = lumi;
      // pile-up
      double puWeight     = 1.0; if(inputFlatFile.category != 0) puWeight     = gltEvent.sf_pu; 
      double puWeightUp   = 1.0; if(inputFlatFile.category != 0) puWeightUp   = gltEvent.sf_puUp;
      double puWeightDown = 1.0; if(inputFlatFile.category != 0) puWeightDown = gltEvent.sf_puDown;
      // lepton efficiency
      double effSF = 1.0;
      if(inputFlatFile.category != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF *= (*idLepTrkSfs[idLep[nl]]) * (*idLepIdSfs[idLep[nl]]); // too many asterisks :(
        }
      }

      // fake rate
      int theCategory = inputFlatFile.category;
      double fakeSF = 1.0;
      if(usePureMC == false){
        printf("NEED TO WORK ON IT IF WE WANT TO USE IT\n");return;
        // What can I say... we need to work on it, if we want to use it. ~DGH

        //if     (theCategory == 5){ // remove W+jets from MC
        //  fakeSF = 0.0;
        //}
        //else if(theCategory == 2 && goodIsTight != idTight.size()){ // remove Z+jets from MC as fakeable objects
        //  fakeSF = 0.0;
        //}
        //else if((inputFlatFile.category == 0 || inputFlatFile.category == 6 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add W+jets from data
        //  for(unsigned int nl=0; nl<idLep.size(); nl++){
        //    if(idTight[nl] == 1) continue;
        //    effSF = effSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
        //    theCategory = 5;
        //  }
        //  if     (inputFlatFile.category != 0 && goodIsTight == idTight.size()-2) effSF =  1.0 * effSF; // double fake, MC
        //  else if(inputFlatFile.category != 0 && goodIsTight == idTight.size()-1) effSF = -1.0 * effSF; // single fake, MC
        //  else if(inputFlatFile.category == 0 && goodIsTight == idTight.size()-2) effSF = -1.0 * effSF; // double fake, data
        //  else if(inputFlatFile.category == 0 && goodIsTight == idTight.size()-1) effSF =  1.0 * effSF; // single fake, data
        //}
        //else if(inputFlatFile.category != 0 && inputFlatFile.category != 6 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
        //  fakeSF = 0.0;
        //}
        //else if(inputFlatFile.category != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
        //  fakeSF = 1.0;
        //}
        //else if(inputFlatFile.category == 0 || inputFlatFile.category == 6){ // data or W+gamma with all good leptons
        //  fakeSF = 1.0;
        //}
        //else {
        //  printf("PROBLEM: %d %d %d %d %d\n",inputFlatFile.category,goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
        //  assert(0);
        //}
      }
                          
      // Begin event weighting                     
      double totalWeight = isData? 1.0 : normalizedWeight; if(!isData) {
        totalWeight *= theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff; // lumi is in inverse femtobarns
        // Btag scale factor
        totalWeight *= sf_btag0; if(totalWeight == 0) continue;
        // ZH EWK correction (only for SM case)
        if(theCategory==6 && signalName_[nModel]=="sm") totalWeight *= gltEvent.sf_zh;
        // ZZ corrections for only the quark induced datasets
        if( theCategory==4 && !isGluonInduced) totalWeight *= gltEvent.sf_zz;
        // WZ
        if(theCategory == 3) totalWeight *= gltEvent.sf_wz;
      } // End event weighting

      for(int nl=0; nl <=sumEvol; nl++) histo[allPlots-2][theCategory]->Fill((double)nl,totalWeight);
      for(int nl=0; nl <=sumEvol; nl++) histoZHSEL[typePair ]         ->Fill((double)nl,totalWeight);
      if(typePair == 1 || typePair == 2)
      for(int nl=0; nl <=sumEvol; nl++) histoZHSEL[3]                 ->Fill((double)nl,totalWeight);
      
      // Need to remind myself what these histos do ~DGH
      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i] && (theCategory == 6 || theCategory == 7)) {
          bgdDecay[nModel][i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[nModel][i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        } else if(passAllCuts[i]) { for(int mModel=0; mModel<nSigModels; mModel++) { 
          bgdDecay[mModel][i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[mModel][i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        }}
      }

      if((typeSel == typePair) || (typeSel == 3 && (typePair == 1 || typePair == 2))) {
        if(nModel<0 || nModel==plotModel){
          for(int thePlot=0; thePlot<allPlots-2; thePlot++){
            double theVar = 0.0;
            bool makePlot = false;
                 if(thePlot ==  0 && passZMass)              {makePlot = true;theVar = idJet.size();}
            else if(thePlot ==  1 && passZMass)              {makePlot = true;theVar = idBJet.size();}
            else if(thePlot ==  2 && passZMass)              {makePlot = true;theVar = idLep.size();;}
            else if(thePlot ==  3 && passZMass)              {makePlot = true;theVar = TMath::Min(TMath::Max(dPhiDiLepMET,-0.05),3.099);}
            else if(thePlot ==  4 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min((double)gltEvent.pfmet,199.999);}
            else if(thePlot ==  5 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
            else if(thePlot ==  6 && passAllCuts[PRESEL])    {makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
            else if(thePlot ==  7 && passNMinusOne[10])      {makePlot = true;theVar = TMath::Min((double)gltEvent.nTau,3.499);}
            else if(thePlot ==  8 && passNMinusOne[2])       {makePlot = true;theVar = TMath::Min(sumPtJets/(sumPtJets+gltEvent.pfmet+dilep.Pt()),0.999);}
            else if(thePlot ==  9 && passNMinusOne[0])       {makePlot = true;theVar = TMath::Min(mtW,999.999);}
            else if(thePlot == 10 && passNMinusOne[1])       {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.M()-91.1876),99.999);}
            else if(thePlot == 11 && passNMinusOne[2])       {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
            else if(thePlot == 12 && passNMinusOne[3])       {makePlot = true;theVar = TMath::Min((double)gltEvent.pfmet,399.999);}
            else if(thePlot == 13 && passNMinusOne[4])       {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
            else if(thePlot == 14 && passNMinusOne[5])       {makePlot = true;theVar = dPhiDiLepMET;}
            else if(thePlot == 15 && passNMinusOne[6])       {makePlot = true;theVar = TMath::Max(TMath::Min(bDiscrMax,0.999),0.001);}
            else if(thePlot == 16 && passNMinusOne[7])       {makePlot = true;theVar = TMath::Min(dilep.Pt(),249.999);}
            else if(thePlot == 17 && passNMinusOne[9])       {makePlot = true;theVar = dPhiJetMET;}
            else if(thePlot == 18 && passNMinusOne[8])       {makePlot = true;theVar = dphill;}
            else if(thePlot == 19 && passNMinusOne[3])       {makePlot = true;theVar = gltEvent.pfmet;}
            else if(thePlot == 20 && passNMinusOne[8])       {makePlot = true;theVar = TMath::Min(drll,2.999);}
            else if(thePlot == 21 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(mtW,999.999);}
            else if(thePlot == 22 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min((double)*idLepPts[idLep[0]],199.999);}
            else if(thePlot == 23 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min((double)*idLepPts[idLep[1]],199.999);}
            else if(thePlot == 24 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min((double)gltEvent.npv,39.499);}
            else if(thePlot == 25 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(dilep.Pt()/mtW,0.999);}
            else if(thePlot == 26 && passAllCuts[TIGHTSEL] && gltEvent.pfmet>150) {makePlot=true;theVar = TMath::Min(1., TMath::Max(-1.,bdt_value));}
            else if(thePlot == 27 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = dPhiLepMETMin;}
            else if(thePlot == 28 && passAllCuts[TIGHTSEL])  {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.Eta()),2.499);}
            else if(thePlot == 29 && passAllCuts[DYSANESEL1]){makePlot = true;theVar = TMath::Min(TMath::Max(dPhiJetMET,-0.05),3.099);}
            else if(thePlot == 30 && passAllCuts[DYSANESEL1]){makePlot = true;theVar = TMath::Min(caloMinusPFMETRel,1.999);}
            else if(thePlot == 31 && passAllCuts[DYSANESEL1]){makePlot = true;theVar = TMath::Min((double)gltEvent.pfmet,499.999);}
            else if(thePlot == 32 && passAllCuts[DYSANESEL2]){makePlot = true;theVar = TMath::Min(TMath::Max(dPhiJetMET,-0.05),3.099);}
            else if(thePlot == 33 && passAllCuts[DYSANESEL2]){makePlot = true;theVar = TMath::Min(caloMinusPFMETRel,1.999);}
            else if(thePlot == 34 && passAllCuts[DYSANESEL2]){makePlot = true;theVar = TMath::Min((double)gltEvent.pfmet,499.999);}
            else if(thePlot == 35 && passAllCuts[ZHGSEL])    {makePlot = true;theVar = TMath::Min((double)gltEvent.loosePho1Eta,2.499);}
            if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
          }
        }
      }
      if(typeSel == typePair || typeSel == 3) {
        double MVAVar = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, metP4.Pt(), mtW, dilep.M(), bdt_value, xbins[nBinMVA]);
        double MVAVarMETSyst[2] = {MVAVar, MVAVar};        
        if(MVAVarType==1) {
          MVAVarMETSyst[0] = TMath::Min((double)gltEvent.pfmetUp,xbins[nBinMVA]-0.001);
          MVAVarMETSyst[1] = TMath::Min((double)gltEvent.pfmetDown,xbins[nBinMVA]-0.001);
        }

        // Throw O(1000) toys for the BDT nuisance evaluations
        if(useBDT && passAllCuts[SIGSEL] && !useCachedBDTSystematics && theCategory >= 3) {
          double bdt_toy_value_muonScale, bdt_toy_value_electronScale, bdt_toy_value_METScale, MVAVar_toy_muonScale, MVAVar_toy_electronScale, MVAVar_toy_METScale;
          TLorentzVector lepton1 = idLep1P4,
                         lepton2 = idLep2P4,
                         MET     = metP4,
                         jet1    = jet1P4;
          double lepton1_scale_variation, lepton2_scale_variation;
          bool lep1IsMuon=(*idLepPdgIds[idLep[0]])==13 || (*idLepPdgIds[idLep[0]])==-13;
          bool lep2IsMuon=(*idLepPdgIds[idLep[1]])==13 || (*idLepPdgIds[idLep[1]])==-13;
          for(unsigned int i_toy=0; i_toy<num_bdt_toys; i_toy++) {
            // BDT variation with the muon scale variation (flat 1%)
            bdt_toy_value_muonScale = bdt_value; MVAVar_toy_muonScale = MVAVar;
            if(lep1IsMuon||lep2IsMuon) { // don't perform expensive BDT evaluation unless there is a muon
              lepton1_scale_variation = lep1IsMuon? 0.01 * bdt_toy_scale[i_toy] : 0;
              lepton2_scale_variation = lep2IsMuon? 0.01 * bdt_toy_scale[i_toy] : 0;
              bdt_toy_value_muonScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation, 0, 0);
                MVAVar_toy_muonScale = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, metP4.Pt(), mtW, dilep.M(), bdt_toy_value_muonScale, xbins[nBinMVA]);
            }
            // BDT variation with the electron scale variation (flat 1%)
            bdt_toy_value_electronScale = bdt_value; MVAVar_toy_electronScale = MVAVar;
            if(!lep1IsMuon || !lep2IsMuon) { // don't perform expensive BDT evaluation unless there is an electron
              lepton1_scale_variation = !lep1IsMuon? 0.01 * bdt_toy_scale[i_toy] : 0;
              lepton2_scale_variation = !lep1IsMuon? 0.01 * bdt_toy_scale[i_toy] : 0;
              bdt_toy_value_electronScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, lepton1_scale_variation, lepton2_scale_variation, 0, 0);
                  MVAVar_toy_electronScale = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, metP4.Pt(), mtW, dilep.M(), bdt_toy_value_electronScale, xbins[nBinMVA]);
            }
            // BDT variation with the MET scale variation (from Nero)
            if(bdt_toy_scale[i_toy] >= 0) 
              bdt_toy_value_METScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, bdt_toy_scale[i_toy] * (gltEvent.pfmetUp / MET.Pt() - 1.), 0);
            else
              bdt_toy_value_METScale = mvaNuisances(reader, lepton1, lepton2, MET, jet1, mva_balance, mva_cos_theta_star_l1, mva_cos_theta_CS_l1, mva_delphi_ptll_MET, mva_delphi_ll, mva_delphi_jet_MET, mva_deltaR_ll, mva_etall, mva_etal1, mva_etal2, mva_MET, mva_mll_minus_mZ, mva_mTjetMET, mva_mTll, mva_mTl1MET, mva_mTl2MET, mva_ptll, mva_ptl1, mva_ptl2, mva_ptl1mptl2_over_ptll, 0, 0, -bdt_toy_scale[i_toy] * (gltEvent.pfmetDown / MET.Pt() - 1.), 0);
                  MVAVar_toy_METScale = getMVAVar(MVAVarType, passAllCuts[TIGHTSEL], typePair, metP4.Pt(), mtW, dilep.M(), bdt_toy_value_METScale, xbins[nBinMVA]);

            if(theCategory == 3)      { histo_bdt_toys_muonScale_WZ       ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_WZ       ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_WZ       ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 4) { histo_bdt_toys_muonScale_ZZ       ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_ZZ       ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_ZZ       ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 5) { histo_bdt_toys_muonScale_VVV      ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_VVV      ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_VVV      ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 6) { histo_bdt_toys_muonScale_ZH_hinv[nModel]  ->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_ZH_hinv[nModel]  ->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_ZH_hinv[nModel]  ->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
            else if(theCategory == 7) { histo_bdt_toys_muonScale_ggZH_hinv->Fill(MVAVar_toy_muonScale, i_toy, totalWeight); histo_bdt_toys_electronScale_ggZH_hinv->Fill(MVAVar_toy_electronScale, i_toy, totalWeight); histo_bdt_toys_METScale_ggZH_hinv->Fill(MVAVar_toy_METScale, i_toy, totalWeight); }
          }
        }

        if     (theCategory == 0){
          if(passAllCuts[SIGSEL]) histo_Data->Fill(MVAVar,totalWeight);
          if(passAllCuts[TIGHTSEL] && verbose) {
            if      (typePair==1) printf("mm data event: ");
            else if (typePair==2) printf("ee data event: " );
            else if (typePair==3) printf("em data event: " );
            printf("runnumber %d lumisection %d eventnumber %lld ptll %f MET %f njets %d l1_pt %f l1_eta %f l2_pt %f l2_eta %f balance %f\n",
              gltEvent.runNumber,
              gltEvent.lumiNumber,
              gltEvent.eventNumber,
              dilep.Pt(),
              metP4.Pt(),
              (int)idJet.size(),
              idLep1P4.Pt(), idLep1P4.Eta(),
              idLep2P4.Pt(), idLep2P4.Eta(),
              ptFrac
            );
          }
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
             histo_WZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*1.02);
             histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[0]));
             histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[1]));
             histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[2]));
             histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[3]));
             histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[4]));
             histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[5]));
             histo_WZ_CMS_PDFUp  ->Fill(MVAVar, totalWeight*gltEvent.pdfUp  );
             histo_WZ_CMS_PDFDown->Fill(MVAVar, totalWeight*gltEvent.pdfDown);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_WZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_WZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*sf_btag0MUp  /sf_btag0);
             histo_WZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*sf_btag0MDown/sf_btag0);
          }
          if(passSystCuts[JESUP])  histo_WZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 4){
          if(passAllCuts[SIGSEL]) {
             histo_ZZ              ->Fill(MVAVar,totalWeight);
             histo_ZZNoW           ->Fill(MVAVar,1.);
             histo_ZZ_CMS_EWKCorrUp->Fill(MVAVar,totalWeight*gltEvent.sf_zzUnc);
             if(true)  histo_ZZ_CMS_ggCorrUp->Fill(MVAVar,totalWeight*1.30); // WARNING: Check that this is only for the gluon induced ZZ samples! ~DGH
             histo_ZZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[0]));
             histo_ZZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[1]));
             histo_ZZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[2]));
             histo_ZZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[3]));
             histo_ZZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[4]));
             histo_ZZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[5]));
             histo_ZZ_CMS_PDFUp  ->Fill(MVAVar, totalWeight*gltEvent.pdfUp  );
             histo_ZZ_CMS_PDFDown->Fill(MVAVar, totalWeight*gltEvent.pdfDown);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_ZZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ZZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*sf_btag0MUp/sf_btag0);
             histo_ZZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*sf_btag0MDown/sf_btag0);
          }
          if(passSystCuts[JESUP])  histo_ZZ_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZZ_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 5){
          if(passAllCuts[SIGSEL]) {
             histo_VVV   ->Fill(MVAVar,totalWeight);
             histo_VVVNoW->Fill(MVAVar,1.);
             histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[0]));
             histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[1]));
             histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[2]));
             histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[3]));
             histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[4]));
             histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[5]));
             histo_VVV_CMS_PDFUp  ->Fill(MVAVar, totalWeight*gltEvent.pdfUp  );
             histo_VVV_CMS_PDFDown->Fill(MVAVar, totalWeight*gltEvent.pdfDown);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_VVV_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*sf_btag0MUp/sf_btag0);
             histo_VVV_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*sf_btag0MDown/sf_btag0);
          }
          if(passSystCuts[JESUP])  histo_VVV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VVV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 6){
          if(passAllCuts[SIGSEL]) {
             histo_ZH_hinv[nModel]   ->Fill(MVAVar,totalWeight);
             histo_ZH_hinvNoW[nModel]->Fill(MVAVar,1.);
             if(nModel == 0) {
               histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->Fill(MVAVar,totalWeight*gltEvent.sf_zhUp);
               histo_ZH_hinv_CMS_EWKCorrDown[nModel]->Fill(MVAVar,totalWeight*gltEvent.sf_zhDown);
             } else {
               histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->Fill(MVAVar,totalWeight);
               histo_ZH_hinv_CMS_EWKCorrDown[nModel]->Fill(MVAVar,totalWeight);
             }
             histo_ZH_hinv_CMS_QCDScaleBounding[nModel][0]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[0]));
             histo_ZH_hinv_CMS_QCDScaleBounding[nModel][1]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[1]));
             histo_ZH_hinv_CMS_QCDScaleBounding[nModel][2]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[2]));
             histo_ZH_hinv_CMS_QCDScaleBounding[nModel][3]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[3]));
             histo_ZH_hinv_CMS_QCDScaleBounding[nModel][4]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[4]));
             histo_ZH_hinv_CMS_QCDScaleBounding[nModel][5]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[5]));
             histo_ZH_hinv_CMS_PDFUp  [nModel]->Fill(MVAVar, totalWeight*gltEvent.pdfUp  );
             histo_ZH_hinv_CMS_PDFDown[nModel]->Fill(MVAVar, totalWeight*gltEvent.pdfDown);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingAvg [nModel]->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingAvg [nModel]->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingUp  [nModel]->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingUp  [nModel]->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Fill(MVAVar,totalWeight*0.98);
             histo_ZH_hinv_CMS_PUBoundingUp  [nModel]->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZH_hinv_CMS_PUBoundingDown[nModel]->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]  ->Fill(MVAVar,totalWeight*sf_btag0MUp/sf_btag0);
             histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]->Fill(MVAVar,totalWeight*sf_btag0MDown/sf_btag0);
          }
          if(passSystCuts[JESUP])  histo_ZH_hinv_CMS_MVAJESBoundingUp  [nModel]->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ZH_hinv_CMS_MVAMETBoundingUp  [nModel]->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else if(theCategory == 7){
          if(passAllCuts[SIGSEL]) {
             histo_ggZH_hinv   ->Fill(MVAVar,totalWeight);
             histo_ggZH_hinvNoW->Fill(MVAVar,1.0);
             histo_ggZH_hinv_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[0]));
             histo_ggZH_hinv_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[1]));
             histo_ggZH_hinv_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[2]));
             histo_ggZH_hinv_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[3]));
             histo_ggZH_hinv_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[4]));
             histo_ggZH_hinv_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*(1.+gltEvent.scale[5]));
             histo_ggZH_hinv_CMS_PDFUp  ->Fill(MVAVar, totalWeight*gltEvent.pdfUp  );
             histo_ggZH_hinv_CMS_PDFDown->Fill(MVAVar, totalWeight*gltEvent.pdfDown);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
             if     (typePair == 1) histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
             else if(typePair == 2) histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
             histo_ggZH_hinv_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ggZH_hinv_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ggZH_hinv_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*sf_btag0MUp/sf_btag0);
             histo_ggZH_hinv_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*sf_btag0MDown/sf_btag0);
          }
          if(passSystCuts[JESUP])  histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ggZH_hinv_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->Fill(MVAVarMETSyst[0],totalWeight);
          if(passSystCuts[METDOWN])histo_ggZH_hinv_CMS_MVAMETBoundingDown->Fill(MVAVarMETSyst[1],totalWeight);
        }
        else {
          printf("CATEGORY PROBLEM!\n"); return;
        }
      }
    }
    inputFlatFile.Close();
  } // end of chain

  // "-1" to remove the Higgs contribution
  double sumEvents = 0;
  for(int np=1; np<processTypes-1; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  //printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<processTypes; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  //printf(" = %.3f\n",sumEvents);

  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Printing yields and statistical uncertainties for all signal models\n\n");
  for(int nModel=0; nModel<nSigModels; nModel++) {
    if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
    printf("Model: %s\n", signalName_[nModel].Data()); 
    printf("                    em                      mm                      ee                      ll\n");
    printf("-----------------------------------------------------------------------------------------------------------\n");
    for(int ns=0; ns<nSelTypes; ns++) {
      printf("Selection: %s\n",selTypeName[ns].Data());
      double sumEventsType[4] = {0,0,0,0}; double sumEventsTypeE[4] = {0,0,0,0};
      for(int np=0; np<processTypes; np++) {       
         bgdDecay[nModel][ns+nSelTypes*3][np] = bgdDecay[nModel][ns+nSelTypes*1][np] + bgdDecay[nModel][ns+nSelTypes*2][np];
         weiDecay[nModel][ns+nSelTypes*3][np] = weiDecay[nModel][ns+nSelTypes*1][np] + weiDecay[nModel][ns+nSelTypes*2][np];
         printf("(%6s): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
         processName[np].Data(),bgdDecay[nModel][ns+nSelTypes*0][np],sqrt(weiDecay[nModel][ns+nSelTypes*0][np]),bgdDecay[nModel][ns+nSelTypes*1][np],sqrt(weiDecay[nModel][ns+nSelTypes*1][np]),
                                bgdDecay[nModel][ns+nSelTypes*2][np],sqrt(weiDecay[nModel][ns+nSelTypes*2][np]),bgdDecay[nModel][ns+nSelTypes*3][np],sqrt(weiDecay[nModel][ns+nSelTypes*3][np]));
         if(np!=0 && np!=6 && np!=7){
           sumEventsType[0] = sumEventsType[0] + bgdDecay[nModel][ns+nSelTypes*0][np]; sumEventsTypeE[0] = sumEventsTypeE[0] + weiDecay[nModel][ns+nSelTypes*0][np];
           sumEventsType[1] = sumEventsType[1] + bgdDecay[nModel][ns+nSelTypes*1][np]; sumEventsTypeE[1] = sumEventsTypeE[1] + weiDecay[nModel][ns+nSelTypes*1][np];
           sumEventsType[2] = sumEventsType[2] + bgdDecay[nModel][ns+nSelTypes*2][np]; sumEventsTypeE[2] = sumEventsTypeE[2] + weiDecay[nModel][ns+nSelTypes*2][np];
           sumEventsType[3] = sumEventsType[3] + bgdDecay[nModel][ns+nSelTypes*3][np]; sumEventsTypeE[3] = sumEventsTypeE[3] + weiDecay[nModel][ns+nSelTypes*3][np];
         }
      }
      printf("(...bkg): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
             sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
               sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]));
    }
  }

  // This is the alternative method
  double kEoverM  = sqrt(bgdDecay[0][ZSEL+nSelTypes*2][0]/bgdDecay[0][ZSEL+nSelTypes*1][0]);
  double NemFact[3] = {(1.0/kEoverM)*0.5, (kEoverM)*0.5, (kEoverM+1.0/kEoverM)*0.5};

  // This is the default method
  double NemFact_FromMLLSB[3] = {NemFact[0], NemFact[1], NemFact[2]}; double NemFact_FromMLLSBE[3] = {1.0, 1.0, 1.0};
  if(bgdDecay[0][WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWLOOSESEL+nSelTypes*1][0]) > 0) {
    NemFact_FromMLLSB[0]  = (bgdDecay[0][WWLOOSESEL+nSelTypes*1][0])/bgdDecay[0][WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[0] = sqrt(1./bgdDecay[0][WWLOOSESEL+nSelTypes*1][0]+1./bgdDecay[0][WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[0];
  }
  if(bgdDecay[0][WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWLOOSESEL+nSelTypes*2][0]) > 0) {
    NemFact_FromMLLSB[1]  = (bgdDecay[0][WWLOOSESEL+nSelTypes*2][0])/bgdDecay[0][WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[1] = sqrt(1./bgdDecay[0][WWLOOSESEL+nSelTypes*2][0]+1./bgdDecay[0][WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[1];
  }
  if(bgdDecay[0][WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[0][WWLOOSESEL+nSelTypes*3][0]) > 0) {
    NemFact_FromMLLSB[2]  = (bgdDecay[0][WWLOOSESEL+nSelTypes*3][0])/bgdDecay[0][WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[2] = sqrt(1./bgdDecay[0][WWLOOSESEL+nSelTypes*3][0]+1./bgdDecay[0][WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[2];
  }
  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Computing the flavor k-factors used for the EM background\n\n");

  printf("(mm) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[0],NemFact_FromMLLSB[0],NemFact_FromMLLSBE[0]);
  printf("(ee) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[1],NemFact_FromMLLSB[1],NemFact_FromMLLSBE[1]);
  printf("(ll) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[2],NemFact_FromMLLSB[2],NemFact_FromMLLSBE[2]);

  // There uncertainties: closure test, different between default and alternative method, and data statistics
  // Only the second item is used as systematic uncertainty
  double EMSystTotal[3] = {1.0,1.0,1.0}; double EMSyst[2][3] = {0.0,0.0,0.0,0.0,0.0,0.0};

  EMSyst[0][0] = bgdDecay[0][SIGSEL+nSelTypes*1][1]/bgdDecay[0][SIGSEL+nSelTypes*0][1]/NemFact[0];
  if(EMSyst[0][0] < 1.0) EMSyst[0][0] = 1/EMSyst[0][0]; EMSyst[0][0] = EMSyst[0][0] - 1.0;

  EMSyst[0][1] = bgdDecay[0][SIGSEL+nSelTypes*2][1]/bgdDecay[0][SIGSEL+nSelTypes*0][1]/NemFact[1];
  if(EMSyst[0][1] < 1.0) EMSyst[0][1] = 1/EMSyst[0][1]; EMSyst[0][1] = EMSyst[0][1] - 1.0;

  EMSyst[0][2] = bgdDecay[0][SIGSEL+nSelTypes*3][1]/bgdDecay[0][SIGSEL+nSelTypes*0][1]/NemFact[2];
  if(EMSyst[0][2] < 1.0) EMSyst[0][2] = 1/EMSyst[0][2]; EMSyst[0][2] = EMSyst[0][2] - 1.0;

  EMSyst[1][0] = NemFact_FromMLLSB[0]/NemFact[0];
  if(EMSyst[1][0] < 1.0) EMSyst[1][0] = 1/EMSyst[1][0]; EMSyst[1][0] = EMSyst[1][0] - 1.0;

  EMSyst[1][1] = NemFact_FromMLLSB[1]/NemFact[1];
  if(EMSyst[1][1] < 1.0) EMSyst[1][1] = 1/EMSyst[1][1]; EMSyst[1][1] = EMSyst[1][1] - 1.0;

  EMSyst[1][2] = NemFact_FromMLLSB[2]/NemFact[2];
  if(EMSyst[1][2] < 1.0) EMSyst[1][2] = 1/EMSyst[1][2]; EMSyst[1][2] = EMSyst[1][2] - 1.0;

  if(bgdDecay[0][SIGSEL][0] > 0) EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1/bgdDecay[0][SIGSEL][0]);
  else                        EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1.0);

  if(bgdDecay[0][SIGSEL][0] > 0) EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1/bgdDecay[0][SIGSEL][0]);
  else                        EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1.0);

  if(bgdDecay[0][SIGSEL][0] > 0) EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1/bgdDecay[0][SIGSEL][0]);
  else                        EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1.0);

  double EMbkg = bgdDecay[0][SIGSEL][2]+bgdDecay[0][SIGSEL][3]+bgdDecay[0][SIGSEL][4]+bgdDecay[0][SIGSEL][5];

  printf("(mm) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][SIGSEL+nSelTypes*1][1],sqrt(weiDecay[0][SIGSEL+nSelTypes*1][4]),
         bgdDecay[0][SIGSEL][1]*NemFact[0]  ,sqrt(weiDecay[0][SIGSEL][4])*NemFact[0],bgdDecay[0][SIGSEL][0],EMbkg,EMSystTotal[0],EMSyst[0][0],EMSyst[1][0],sqrt(EMSystTotal[0]*EMSystTotal[0] - EMSyst[0][0]*EMSyst[0][0] - EMSyst[1][0]*EMSyst[1][0]));

  printf("(ee) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][SIGSEL+nSelTypes*2][1],sqrt(weiDecay[0][SIGSEL+nSelTypes*2][4]),
         bgdDecay[0][SIGSEL][1]*NemFact[1]  ,sqrt(weiDecay[0][SIGSEL][4])*NemFact[1],bgdDecay[0][SIGSEL][0],EMbkg,EMSystTotal[1],EMSyst[0][1],EMSyst[1][1],sqrt(EMSystTotal[1]*EMSystTotal[1] - EMSyst[0][1]*EMSyst[0][1] - EMSyst[1][1]*EMSyst[1][1]));

  printf("(ll) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[0][SIGSEL+nSelTypes*3][1],sqrt(weiDecay[0][SIGSEL+nSelTypes*3][4]),
         bgdDecay[0][SIGSEL][1]*NemFact[2]  ,sqrt(weiDecay[0][SIGSEL][4])*NemFact[2],bgdDecay[0][SIGSEL][0],EMbkg,EMSystTotal[2],EMSyst[0][2],EMSyst[1][2],sqrt(EMSystTotal[2]*EMSystTotal[2] - EMSyst[0][2]*EMSyst[0][2] - EMSyst[1][2]*EMSyst[1][2]));

  if(histo_EM->GetSumOfWeights() > 1) systEM[0] = 1. + 1./sqrt(histo_EM->GetSumOfWeights());
  else                                systEM[0] = 2.;
  if(useEMFromData == true && histo_Data->GetSumOfWeights() > 0.0){
    EMbkg = bgdDecay[0][TIGHTSEL][2]+bgdDecay[0][TIGHTSEL][3]+bgdDecay[0][TIGHTSEL][4]+bgdDecay[0][TIGHTSEL][5];
    double EMNormFact[4] = {((bgdDecay[0][TIGHTSEL][0]-EMbkg)*NemFact[0])/bgdDecay[0][TIGHTSEL+nSelTypes*(1)][1],
                            ((bgdDecay[0][TIGHTSEL][0]-EMbkg)*NemFact[1])/bgdDecay[0][TIGHTSEL+nSelTypes*(2)][1],
                            ((bgdDecay[0][TIGHTSEL][0]-EMbkg)*NemFact[2])/bgdDecay[0][TIGHTSEL+nSelTypes*(3)][1],
                            ((bgdDecay[0][TIGHTSEL][0]-EMbkg)*1.00000000)/bgdDecay[0][TIGHTSEL+nSelTypes*(0)][1]};

    printf("EM(1): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*1][1],bgdDecay[0][TIGHTSEL][0],EMbkg,NemFact[0],bgdDecay[0][TIGHTSEL+nSelTypes*(1)][1],bgdDecay[0][TIGHTSEL+nSelTypes*1][1]*EMNormFact[0],bgdDecay[0][TIGHTSEL+nSelTypes*1][1]*EMNormFact[0]*EMSystTotal[0]);
    printf("EM(2): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*2][1],bgdDecay[0][TIGHTSEL][0],EMbkg,NemFact[1],bgdDecay[0][TIGHTSEL+nSelTypes*(2)][1],bgdDecay[0][TIGHTSEL+nSelTypes*2][1]*EMNormFact[1],bgdDecay[0][TIGHTSEL+nSelTypes*2][1]*EMNormFact[1]*EMSystTotal[1]);
    printf("EM(3): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*3][1],bgdDecay[0][TIGHTSEL][0],EMbkg,NemFact[2],bgdDecay[0][TIGHTSEL+nSelTypes*(3)][1],bgdDecay[0][TIGHTSEL+nSelTypes*3][1]*EMNormFact[2],bgdDecay[0][TIGHTSEL+nSelTypes*3][1]*EMNormFact[2]*EMSystTotal[2]);
    printf("EM(4): %f * (%f-%f)*%f/%f = %f +/- %f\n",bgdDecay[0][TIGHTSEL+nSelTypes*0][1],bgdDecay[0][TIGHTSEL][0],EMbkg,1.00000000,bgdDecay[0][TIGHTSEL+nSelTypes*(0)][1],bgdDecay[0][TIGHTSEL+nSelTypes*0][1]*EMNormFact[3],bgdDecay[0][TIGHTSEL+nSelTypes*0][1]*EMNormFact[3]*0.000000000000);

    //systEM[0] = 1. + EMSystTotal[typeSel-1];
    //systEM[0] = 1. + EMSyst[1][typeSel-1];
    systEM[0] = 1. + 0.20;
    systEM[1] = TMath::Max(bgdDecay[0][TIGHTSEL][0],1.0);
    if(bgdDecay[0][TIGHTSEL][0] > 0) {
      histo_EM->Scale(EMNormFact[typeSel-1]);
      histo_EM->SetBinContent(1,histo_EM->GetBinContent(1)*EMNormFact[3]/EMNormFact[typeSel-1]);
    }
  }

  // DY background estimation
  if(useZjetsTemplate){
    histo_Zjets_CMS_ZjetsSystUp->Add(fhDZjets);
    histo_Zjets->Scale(0.0);
    histo_Zjets->Add(fhDZjets);
    histo_ZjetsNoW->Scale(0.0);
    histo_ZjetsNoW->Add(fhDZjets);
    double ZJetsNorm[3] = {0.56, 0.52, 1.};
    if     (MVAVarType == 0 || MVAVarType == 1 ||MVAVarType == 2 || MVAVarType == 3) {
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

  // computing DY scale factor using the first bin
  if((MVAVarType == 0 || MVAVarType == 1 || MVAVarType == 2 || MVAVarType==3 || MVAVarType==4) && histo_Zjets->GetSumOfWeights() > 0.0 && histo_Data->GetSumOfWeights() > 0.0){
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("Computing the Drell-Yan data/MC scale factor using the first shape bin\n\n");
    int theBin = 2;
    double the_bck = histo_EM->GetBinContent(theBin) + histo_WZ->GetBinContent(theBin) + histo_ZZ->GetBinContent(theBin) + histo_VVV->GetBinContent(theBin);
    double the_data = histo_Data->GetBinContent(theBin);
    double the_sf;
    //if(MVAVarType!=3) { 
    the_sf = (the_data-the_bck)/histo_Zjets->GetBinContent(theBin);
    printf("DY SF: data/bck/DY = %f %f %f ==> %f\n",the_data,the_bck,histo_Zjets->GetBinContent(theBin),the_sf);
    //} else {
    //  the_sf=1.068889;
    //}
    histo_Zjets->Scale(the_sf);
    histo_Zjets_CMS_ZjetsSystUp->Scale(the_sf);
  }

  double nZjetsNorm = histo_Zjets->GetSumOfWeights();
  if(nZjetsNorm > 0){
    for(int i=1; i<=histo_Zjets->GetNbinsX(); i++){
      histo_Zjets->SetBinContent(i,TMath::Max(histo_Zjets->GetBinContent(i),0.000001));
    }
    histo_Zjets->Scale(nZjetsNorm/histo_Zjets->GetSumOfWeights());
  }

  double mean,up,diff;
  for(int i=1; i<=histo_ZH_hinv[0]->GetNbinsX(); i++){
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
  printf("-----------------------------------------------------------------------------------------------------------\n");
  printf("Computing diboson corrections\n\n");
  
  printf("EWK Corr: WZ(%f/%f/%f) ZZ(%f/%f/%f) ggZZ(%f/%f/%f) ZH(%f/%f/%f)\n",
                     histo_WZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_WZ->GetSumOfWeights(),histo_WZ_CMS_EWKCorrDown->GetSumOfWeights(),
                     histo_ZZ_CMS_EWKCorrUp->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_EWKCorrDown->GetSumOfWeights(),
                     histo_ZZ_CMS_ggCorrUp ->GetSumOfWeights(),histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_ggCorrDown ->GetSumOfWeights(),
                     histo_ZH_hinv_CMS_EWKCorrUp[0]->GetSumOfWeights(),histo_ZH_hinv[0]->GetSumOfWeights(),histo_ZH_hinv_CMS_EWKCorrDown[0]->GetSumOfWeights());

  printf("QCD Corr: WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) ZH(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZH_hinv[0]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][0]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][1]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][2]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][3]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][4]->GetSumOfWeights(),histo_ZH_hinv_CMS_QCDScaleBounding[0][5]->GetSumOfWeights());

  histo[allPlots-1][0]->Add(histo_Data);
  histo[allPlots-1][1]->Add(histo_EM);
  histo[allPlots-1][2]->Add(histo_Zjets);
  histo[allPlots-1][3]->Add(histo_WZ);
  histo[allPlots-1][4]->Add(histo_ZZ);
  histo[allPlots-1][5]->Add(histo_VVV);
  histo[allPlots-1][6]->Add(histo_ZH_hinv[plotModel]);
  histo[allPlots-1][7]->Add(histo_ggZH_hinv);
  
  

  for(int i=1; i<=histo_ZH_hinv[0]->GetNbinsX(); i++) {
    double factorUp = +1.0; double factorDown = -1.0;
    histo_VVV_CMS_MVAVVVStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorUp  *histo_VVV        ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorDown*histo_VVV        ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp               ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorUp  *histo_WZ        ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorDown*histo_WZ        ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp               ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorUp  *histo_ZZ        ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorDown*histo_ZZ        ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingUp               ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorUp  *histo_EM        ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorDown*histo_EM       ->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorUp  *histo_ggZH_hinv->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorDown*histo_ggZH_hinv->GetBinError(i),0.000001));

    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]                 ->Add(histo_VVV      ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]         ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorUp  *histo_VVV      ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]         ->Add(histo_VVV      ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]       ->SetBinContent(i,TMath::Max(histo_VVV      ->GetBinContent(i)+factorDown*histo_VVV      ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]                 ->Add(histo_WZ       ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]           ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorUp  *histo_WZ       ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]                 ->Add(histo_WZ       ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]         ->SetBinContent(i,TMath::Max(histo_WZ       ->GetBinContent(i)+factorDown*histo_WZ       ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]                 ->Add(histo_ZZ       ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]           ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorUp  *histo_ZZ       ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]                 ->Add(histo_ZZ       ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]         ->SetBinContent(i,TMath::Max(histo_ZZ       ->GetBinContent(i)+factorDown*histo_ZZ       ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]                 ->Add(histo_EM       ); histo_EM_CMS_MVAEMStatBoundingBinUp[i-1]           ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorUp  *histo_EM       ->GetBinError(i),0.000001));
    histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]           ->Add(histo_EM       ); histo_EM_CMS_MVAEMStatBoundingBinDown[i-1]         ->SetBinContent(i,TMath::Max(histo_EM       ->GetBinContent(i)+factorDown*histo_EM       ->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[i-1]    ->Add(histo_ggZH_hinv); histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorUp  *histo_ggZH_hinv->GetBinError(i),0.000001));
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[i-1]  ->Add(histo_ggZH_hinv); histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_ggZH_hinv->GetBinContent(i)+factorDown*histo_ggZH_hinv->GetBinError(i),0.000001));

    for(int nModel=0; nModel<nSigModels; nModel++) { 
      if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
      histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]      ->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorUp  *histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
      histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]    ->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorDown*histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][i-1]  ->Add(histo_ZH_hinv[nModel]);
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][i-1]->Add(histo_ZH_hinv[nModel]);
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][i-1]  ->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorUp  *histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][i-1]->SetBinContent(i,TMath::Max(histo_ZH_hinv[nModel]->GetBinContent(i)+factorDown*histo_ZH_hinv[nModel]->GetBinError(i),0.000001));
    }

  }
  assert(SaveHistos());  
  unsigned long int t2 = static_cast<unsigned long int>(time(NULL));
  printf("zhAnalysis::Run : Finished running the whole analysis (%lu seconds)\n", t2-t1);
}
void zhAnalysis::SetFinalStateName() {
  if     (typeSel == 3 && nJetsType == 0) {sprintf(finalStateName,"ll0j");}
  else if(typeSel == 1 && nJetsType == 0) {sprintf(finalStateName,"mm0j");}
  else if(typeSel == 2 && nJetsType == 0) {sprintf(finalStateName,"ee0j");}
  else if(typeSel == 3 && nJetsType == 1) {sprintf(finalStateName,"ll1j");}
  else if(typeSel == 1 && nJetsType == 1) {sprintf(finalStateName,"mm1j");}
  else if(typeSel == 2 && nJetsType == 1) {sprintf(finalStateName,"ee1j");}
  else if(typeSel == 3 && nJetsType == 2) {sprintf(finalStateName,"ll2j");}
  else if(typeSel == 1 && nJetsType == 2) {sprintf(finalStateName,"mm2j");}
  else if(typeSel == 2 && nJetsType == 2) {sprintf(finalStateName,"ee2j");}
  else if(typeSel == 0 && nJetsType == 0) {sprintf(finalStateName,"em0j");}
  else if(typeSel == 0 && nJetsType == 1) {sprintf(finalStateName,"em1j");}
  else if(typeSel == 0 && nJetsType == 2) {sprintf(finalStateName,"em2j");}
  else {printf("Wrong lSel/nJetsType: %d/%d\n",typeSel,nJetsType); assert(0); return;}
}

bool zhAnalysis::LoadFlatFiles(bool doDM) {
  printf("zhAnalysis::LoadFlatFiles : Begin loading flat input files\n");
  unsigned long int t1 = static_cast<unsigned long int>(time(NULL));
  //*******************************************************
  //Input Files
  //*******************************************************
  // Data files - Macro Category 0
  inputFlatFiles.emplace_back(Form("%sdata.root",filesPathDA.Data()),0,-1);
  /*
  // Begin MC backgrounds 
      
    // Combo / flavor-symmetric / non-resonant Backgrounds - Macro Category 1
    inputFlatFiles.emplace_back(Form("%sggWW.root",filesPathMC.Data()),1,-1);
    inputFlatFiles.emplace_back(Form("%sTT2L.root",filesPathMC.Data()),1,-1);
    // Previously counted TTZToLLNuNu as VVV background? ~DGH
    inputFlatFiles.emplace_back(Form("%sTTV.root",filesPathMC.Data()),1,-1);
    inputFlatFiles.emplace_back(Form("%sTW.root",filesPathMC.Data()),1,-1);
    inputFlatFiles.emplace_back(Form("%sqqWW.root",filesPathMC.Data()),1,-1);
    inputFlatFiles.emplace_back(Form("%sggWW.root",filesPathMC.Data()),1,-1);
    inputFlatFiles.emplace_back(Form("%sWGstar.root",filesPathMC.Data()),1,-1);
    inputFlatFiles.emplace_back(Form("%sWWdps.root",filesPathMC.Data()),1,-1);
    // Include this one? Does it include (WGToLNuG) ? ~DGH
    inputFlatFiles.emplace_back(Form("%sVG.root",filesPathMC.Data()),1,-1);
    inputFlatFiles.emplace_back(Form("%sH125.root",filesPathMC.Data()),1,-1);
    
    // Drell-Yan / Z backgrounds - Macro Category 2
    if(useDYPT==false){
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_M-10to50.root",filesPathMC.Data()),2,-1);
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_M-50_LO.root",filesPathMC.Data()),2,-1);
    }
    else {
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_Pt0To50.root",filesPathMC.Data()),2,-1);
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_Pt50To100.root",filesPathMC.Data()),2,-1);
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_Pt100To250.root",filesPathMC.Data()),2,-1);
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_Pt250To400.root",filesPathMC.Data()),2,-1);
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_Pt400To650.root",filesPathMC.Data()),2,-1);
      inputFlatFiles.emplace_back(Form("%sDYJetsToLL_Pt650ToInf.root",filesPathMC.Data()),2,-1);
      //inputFlatFiles.emplace_back(Form("%s",filesPathMC.Data()),2,-1);
    }
    inputFlatFiles.emplace_back(Form("%sDYJetsToTauTau.root",filesPathMC.Data()),2,-1);
    // WZ backgrounds - Macro Category 3
    inputFlatFiles.emplace_back(Form("%sWZ.root",filesPathMC.Data()),3,-1);
    // ZZ Backgrounds - Macro Category 4
    inputFlatFiles.emplace_back(Form("%sqqZZ.root",filesPathMC.Data()),4,-1);
    inputFlatFiles.emplace_back(Form("%sggZZ.root",filesPathMC.Data()),4,-1);
      
    // Triboson / VVV Backgrounds - Macro Category 5
    // Does this include tZq? ~DGH
    inputFlatFiles.emplace_back(Form("%sVVV.root",filesPathMC.Data()),5,-1);
      
  // End MC backgrounds
  */
  // Monte Carlo signals
  if(false){ // Model 0: standard model Higgs (125) with glu-glu
    int mH=125;
    signalName_.push_back("sm");
    inputFlatFiles.emplace_back(Form("%sZH_ZToMM_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH),6,0); 
    inputFlatFiles.emplace_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH),6,0); 
    inputFlatFiles.emplace_back(Form("%sggZH_HToInv_ZToLL_M125_13TeV_powheg_pythia8.root",filesPathDMMC.Data()),7,0);       
  }  // Models 1 thru 8: standard-model-like Higgs mass points without glu-glu (8 models)
  else signalName_.push_back("sm"); // catch for now, since we have no models
  if(false){ int mH_[10]={110, 125, 150, 200, 300, 400, 500, 600, 800, 1000}; int iH=0; for(int i=1; i<=10; i++) { int mH = mH_[iH]; 
    signalName_.push_back(Form("mh%d", mH));
    inputFlatFiles.emplace_back(Form("%sZH_ZToMM_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH),6,iH+1); 
    inputFlatFiles.emplace_back(Form("%sZH_ZToEE_HToInvisible_M%d_13TeV_powheg_pythia8.root",filesPathDMMC.Data(),mH),6,iH+1);
    iH++;
  }}

  if(doDM){ 
    int i=signalName_.size();
    signalName_.push_back("ADDMonoZ_ZToLL_MD-3_d-2_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZToLL_MD-3_d-2_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZToLL_MD-3_d-3_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZToLL_MD-3_d-3_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZToLL_MD-3_d-4_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZToLL_MD-3_d-4_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZToLL_MD-3_d-5_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZToLL_MD-3_d-5_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZToLL_MD-3_d-6_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZToLL_MD-3_d-6_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZToLL_MD-3_d-7_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZToLL_MD-3_d-7_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-1_d-2_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-1_d-2_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-1_d-3_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-1_d-3_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-1_d-4_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-1_d-4_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-1_d-5_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-1_d-5_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-1_d-6_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-1_d-6_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-1_d-7_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-1_d-7_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-2_d-2_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-2_d-2_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-2_d-3_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-2_d-3_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-2_d-4_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-2_d-4_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-2_d-5_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-2_d-5_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-2_d-6_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-2_d-6_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("ADDMonoZ_ZtoLL_MD-2_d-7_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sADDMonoZ_ZtoLL_MD-2_d-7_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-1995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-1995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1000_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1000_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-10_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-10_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-150_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-150_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-1_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-1_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-500_Mv-995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-500_Mv-995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_A_Mx-50_Mv-95_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_A_Mx-50_Mv-95_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-0_Mv-20_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-0_Mv-20_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1000_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1000_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1000_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1000_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-100_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-100_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-100_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-100_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-10_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-140_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-140_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-150_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-1_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-200_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-200_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-300_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-300_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-40_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-40_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-490_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-490_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-500_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-500_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-500_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-500_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-500_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-500_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-50_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-75_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-75_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-75_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-75_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Axial_Mx-990_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Axial_Mx-990_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-100_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-10_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-10_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-10_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-10_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-10_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-10_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-20_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-20_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-1_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-200_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-200_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-350_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-350_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-350_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-350_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-40_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-40_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-200_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-200_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Pseudo_Mx-50_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-100_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-1_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-200_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-200_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-200_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-200_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-200_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-200_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-300_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-300_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-300_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-300_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_SMM_Mx-400_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_SMM_Mx-400_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-0_Mv-20_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-0_Mv-20_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-100_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-10_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-10_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-10_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-10_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-10_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-10_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-200_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-200_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-20_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-20_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-1_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-200_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-200_Mv-500_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-350_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-350_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-350_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-350_Mv-750_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-40_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-40_Mv-100_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-450_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-450_Mv-1000_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-10_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-200_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-200_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-300_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-350_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-400_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Scalar_Mx-50_Mv-50_gDM1_gQ1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-0_Mv-20_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-0_Mv-20_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1000_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1000_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-100_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-100_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-10_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-150_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-150_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-150_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-150_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-150_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-150_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-10_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-50_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-1_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-300_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-300_Mv-750_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-40_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-40_Mv-100_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-490_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-490_Mv-1000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-500_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-500_Mv-10000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-500_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-500_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-50_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-50_Mv-200_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-50_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-50_Mv-300_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-50_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-50_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-75_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-75_Mv-500_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_NLO_Vector_Mx-990_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_NLO_Vector_Mx-990_Mv-2000_gDM1_gQ0p25_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-1995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-1995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1000_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1000_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-10_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-10_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-295_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-150_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-150_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-1000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-100_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-10_1-gDMgQ_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-10_1-gDMgQ_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-20_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-1_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-1_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-2000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-500_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-500_Mv-995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-500_Mv-995_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-10_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-200_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-300_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-5000_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-50_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("DarkMatter_MonoZToLL_V_Mx-50_Mv-95_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph"); inputFlatFiles.emplace_back(Form("%sDarkMatter_MonoZToLL_V_Mx-50_Mv-95_gDMgQ-1_TuneCUETP8M1_13TeV-madgraph.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p01_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p01_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p02_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p02_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p04_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p04_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p06_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p06_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p09_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p09_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p10_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p10_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p20_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p20_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p30_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p30_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p40_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p40_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p50_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p50_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p60_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p60_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p70_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p70_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p80_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p80_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-1p90_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-1p90_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-2p00_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-2p00_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
    signalName_.push_back("Unpart_ZToLL_SU-0_dU-2p20_LU-15_TuneCUETP8M1_13TeV-pythia8"); inputFlatFiles.emplace_back(Form("%sUnpart_ZToLL_SU-0_dU-2p20_LU-15_TuneCUETP8M1_13TeV-pythia8.root", filesPathDMMC.Data()), 6, i); i++;
  }

  nSigModels=(unsigned)signalName_.size();
  nInputFiles=(unsigned)inputFlatFiles.size();
  assert(nSigModels<=500);
  unsigned long int t2 = static_cast<unsigned long int>(time(NULL));
  printf("zhAnalysis::LoadFlatFiles : Finished getting list of flat input files (%lu seconds)\n", t2-t1);
  return true;
}

bool zhAnalysis::MakeHistos() {
  printf("zhAnalysis::MakeHistos : Start making histograms...\n");
  unsigned long int t1 = static_cast<unsigned long int>(time(NULL));
  // Make plotting histos
  TString plotName;
  for(int thePlot=0; thePlot<allPlots; thePlot++) {
    TH1D* histos = MakeHisto(thePlot, plotName);
    plotName_[thePlot]=plotName;
    for(int iType=0; iType<processTypes; iType++) histo[thePlot][iType] = (TH1D*) histos->Clone(Form("histo %s %s",categoryName_[iType].Data(), plotName.Data()));
    histos->Reset();histos->Clear();
  }
  // Make MVA shape histos
  histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins); histoMVA->Sumw2();
  histo_Data     = (TH1D*) histoMVA->Clone("histo_Data");
  histo_Zjets    = (TH1D*) histoMVA->Clone("histo_Zjets");         
  histo_VVV      = (TH1D*) histoMVA->Clone("histo_VVV");         
  histo_WZ       = (TH1D*) histoMVA->Clone("histo_WZ");         
  histo_ZZ       = (TH1D*) histoMVA->Clone("histo_ZZ");
  histo_EM       = (TH1D*) histoMVA->Clone("histo_EM");         
  histo_ggZH_hinv= (TH1D*) histoMVA->Clone("histo_ggZH_hinv"); 
  histo_ZjetsNoW    = (TH1D*) histoMVA->Clone("histo_Zjets");         
  histo_VVVNoW      = (TH1D*) histoMVA->Clone("histo_VVV");         
  histo_WZNoW       = (TH1D*) histoMVA->Clone("histo_WZ");         
  histo_ZZNoW       = (TH1D*) histoMVA->Clone("histo_ZZ");
  histo_EMNoW       = (TH1D*) histoMVA->Clone("histo_EM");         
  histo_ggZH_hinvNoW= (TH1D*) histoMVA->Clone("histo_ggZH_hinv"); 
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histo_ZH_hinv[nModel]     = (TH1D*) histoMVA->Clone(Form("histo_ZH_hinv_%s",   signalName_[nModel].Data())); 
    histo_ZH_hinvNoW[nModel]  = (TH1D*) histoMVA->Clone(Form("histo_ZH_hinvNoW_%s",signalName_[nModel].Data())); 
  }

  // Make stat. bounding histos
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]      = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sUp"  ,finalStateName, signalName_[nModel].Data(), ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sUp"  ,finalStateName, signalName_[nModel].Data(), ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]    = new TH1D( Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sDown",finalStateName, signalName_[nModel].Data(), ECMsb.Data()), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%sDown",finalStateName, signalName_[nModel].Data(), ECMsb.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]->Sumw2();
  }
  histo_Zjets_CMS_MVAZjetsStatBoundingUp     = new TH1D( Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingUp  ->Sumw2();
  histo_Zjets_CMS_MVAZjetsStatBoundingDown   = new TH1D( Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingDown->Sumw2();
  histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  histo_WZ_CMS_MVAWZStatBoundingUp           = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  histo_WZ_CMS_MVAWZStatBoundingDown         = new TH1D( Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  histo_EM_CMS_MVAEMStatBoundingUp           = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingUp  ->Sumw2();
  histo_EM_CMS_MVAEMStatBoundingDown         = new TH1D( Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EM_CMS_MVAEMStatBoundingDown->Sumw2();
  histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp  = new TH1D( Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown= new TH1D( Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown->Sumw2();

  // Make QCD scale bounding histos
  for(int nb=0; nb<6; nb++){
    for(int nModel=0; nModel<nSigModels; nModel++) {
      histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nb]   = new TH1D(Form("histo_ZH_hinv_%s_QCDScale_f%d", signalName_[nModel].Data(), nb), Form("histo_ZH_hinv_%s_QCDScale_f%d", signalName_[nModel].Data(), nb),nBinMVA, xbins); histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nb]->Sumw2();
    }
    histo_VVV_CMS_QCDScaleBounding[nb]       = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding[nb]             = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),      Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_WZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding[nb]             = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ggZH_hinv_CMS_QCDScaleBounding[nb] = new TH1D(Form("histo_ggZH_hinv_QCDScale_f%d",nb), Form("histo_ggZH_hinv_QCDScale_f%d",nb),nBinMVA, xbins); histo_ggZH_hinv_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  
  // Make PDF bounding histos
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histo_ZH_hinv_CMS_PDFUp[nModel]   = new TH1D(Form("histo_ZH_hinv_%s_PDFUp", signalName_[nModel].Data()), Form("histo_ZH_hinv_%s_PDFUp", signalName_[nModel].Data()),nBinMVA, xbins); histo_ZH_hinv_CMS_PDFUp[nModel]->Sumw2();
    histo_ZH_hinv_CMS_PDFDown[nModel] = new TH1D(Form("histo_ZH_hinv_%s_PDFDown", signalName_[nModel].Data()), Form("histo_ZH_hinv_%s_PDFDown", signalName_[nModel].Data()),nBinMVA, xbins); histo_ZH_hinv_CMS_PDFDown[nModel]->Sumw2();
  }
  histo_VVV_CMS_PDFUp         = new TH1D( "histo_VVV_PDFUp", "histo_VVV_PDFUp", nBinMVA, xbins); histo_VVV_CMS_PDFUp  ->Sumw2();
  histo_VVV_CMS_PDFDown       = new TH1D( "histo_VVV_PDFDown", "histo_VVV_PDFDown", nBinMVA, xbins); histo_VVV_CMS_PDFDown->Sumw2();
  histo_WZ_CMS_PDFUp          = new TH1D( "histo_WZ_PDFUp", "histo_WZ_PDFUp", nBinMVA, xbins); histo_WZ_CMS_PDFUp  ->Sumw2();
  histo_WZ_CMS_PDFDown        = new TH1D( "histo_WZ_PDFDown", "histo_WZ_PDFDown", nBinMVA, xbins); histo_WZ_CMS_PDFDown->Sumw2();
  histo_ZZ_CMS_PDFUp          = new TH1D( "histo_ZZ_PDFUp", "histo_ZZ_PDFUp", nBinMVA, xbins); histo_ZZ_CMS_PDFUp  ->Sumw2();
  histo_ZZ_CMS_PDFDown        = new TH1D( "histo_ZZ_PDFDown", "histo_ZZ_PDFDown", nBinMVA, xbins); histo_ZZ_CMS_PDFDown->Sumw2();
  histo_ggZH_hinv_CMS_PDFUp   = new TH1D( "histo_ggZH_hinv_PDFUp", "histo_ggZH_hinv_PDFUp", nBinMVA, xbins); histo_ggZH_hinv_CMS_PDFUp  ->Sumw2();
  histo_ggZH_hinv_CMS_PDFDown = new TH1D( "histo_ggZH_hinv_PDFDown", "histo_ggZH_hinv_PDFDown", nBinMVA, xbins); histo_ggZH_hinv_CMS_PDFDown->Sumw2();
  
  // Make stat bounding bin histos 
  for(int nb=0; nb<nBinMVA; nb++) {
    for(int nModel=0; nModel<nSigModels; nModel++) { 
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][nb]        = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dUp"       ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dUp"       ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb),nBinMVA, xbins); 
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][nb]      = new TH1D(Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dDown"     ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb), Form("histo_ZH_hinv_CMS_zllhinv%s_%s_MVAZHStatBounding_%s_Bin%dDown"     ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb),nBinMVA, xbins); 
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinUp[nModel][nb]   ->Sumw2();
      histo_ZH_hinv_CMS_MVAZHStatBoundingBinDown[nModel][nb] ->Sumw2();
    }
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]       = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"      ,finalStateName,  ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dUp"      ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]     = new TH1D(Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown"    ,finalStateName,  ECMsb.Data(),nb), Form("histo_Zjets_CMS_zllhinv%s_MVAZjetsStatBounding_%s_Bin%dDown"    ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]            = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dUp"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]            = new TH1D(Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"        ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_zllhinv%s_MVAVVVStatBounding_%s_Bin%dDown"        ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]                = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]            = new TH1D(Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_WZ_CMS_zllhinv%s_MVAWZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]                = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]            = new TH1D(Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_zllhinv%s_MVAZZStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]                = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dUp"            ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]            = new TH1D(Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb), Form("histo_EM_CMS_zllhinv%s_MVAEMStatBounding_%s_Bin%dDown"          ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[nb]    = new TH1D(Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dUp"   ,finalStateName,  ECMsb.Data(),nb), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dUp"   ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[nb]  = new TH1D(Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dDown" ,finalStateName,  ECMsb.Data(),nb), Form("histo_ggZH_hinv_CMS_zllhinv%s_MVAggZHStatBounding_%s_Bin%dDown" ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_Zjets_CMS_MVAZjetsStatBoundingBinUp[nb]  ->Sumw2();
    histo_Zjets_CMS_MVAZjetsStatBoundingBinDown[nb]->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]              ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]              ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinUp[nb]              ->Sumw2();
    histo_EM_CMS_MVAEMStatBoundingBinDown[nb]          ->Sumw2();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingBinDown[nb] ->Sumw2();
  }

  // Make lepton efficiency bounding histos
  histo_VVV_CMS_MVALepEffMBoundingUp              = new TH1D( Form("histo_VVV_%sUp",effMName.Data())  , Form("histo_VVV_%sUp",effMName.Data())  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  histo_VVV_CMS_MVALepEffMBoundingDown            = new TH1D( Form("histo_VVV_%sDown",effMName.Data()), Form("histo_VVV_%sDown",effMName.Data()), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  histo_WZ_CMS_MVALepEffMBoundingUp               = new TH1D( Form("histo_WZ_%sUp",effMName.Data())  , Form("histo_WZ_%sUp",effMName.Data())  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  histo_WZ_CMS_MVALepEffMBoundingDown             = new TH1D( Form("histo_WZ_%sDown",effMName.Data()), Form("histo_WZ_%sDown",effMName.Data()), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();
  histo_ZZ_CMS_MVALepEffMBoundingUp               = new TH1D( Form("histo_ZZ_%sUp",effMName.Data())  , Form("histo_ZZ_%sUp",effMName.Data())  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  histo_ZZ_CMS_MVALepEffMBoundingDown             = new TH1D( Form("histo_ZZ_%sDown",effMName.Data()), Form("histo_ZZ_%sDown",effMName.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
  histo_ggZH_hinv_CMS_MVALepEffMBoundingUp   = new TH1D( Form("histo_ggZH_hinv_%sUp",effMName.Data())  , Form("histo_ggZH_hinv_%sUp",effMName.Data())  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_MVALepEffMBoundingDown = new TH1D( Form("histo_ggZH_hinv_%sDown",effMName.Data()), Form("histo_ggZH_hinv_%sDown",effMName.Data()), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->Sumw2();

  histo_VVV_CMS_MVALepEffMBoundingAvg              = new TH1D( Form("histo_VVV_%sAvg",effMName.Data())  , Form("histo_VVV_%sAvg",effMName.Data())  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  histo_WZ_CMS_MVALepEffMBoundingAvg               = new TH1D( Form("histo_WZ_%sAvg",effMName.Data())  , Form("histo_WZ_%sAvg",effMName.Data())  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  histo_ZZ_CMS_MVALepEffMBoundingAvg               = new TH1D( Form("histo_ZZ_%sAvg",effMName.Data())  , Form("histo_ZZ_%sAvg",effMName.Data())  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg  = new TH1D( Form("histo_ggZH_hinv_%sAvg",effMName.Data())  , Form("histo_ggZH_hinv_%sAvg",effMName.Data())  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  histo_VVV_CMS_MVALepEffEBoundingUp              = new TH1D( Form("histo_VVV_%sUp",effEName.Data())  , Form("histo_VVV_%sUp",effEName.Data())  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  histo_VVV_CMS_MVALepEffEBoundingDown            = new TH1D( Form("histo_VVV_%sDown",effEName.Data()), Form("histo_VVV_%sDown",effEName.Data()), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  histo_WZ_CMS_MVALepEffEBoundingUp               = new TH1D( Form("histo_WZ_%sUp",effEName.Data())  , Form("histo_WZ_%sUp",effEName.Data())  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  histo_WZ_CMS_MVALepEffEBoundingDown             = new TH1D( Form("histo_WZ_%sDown",effEName.Data()), Form("histo_WZ_%sDown",effEName.Data()), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  histo_ZZ_CMS_MVALepEffEBoundingUp               = new TH1D( Form("histo_ZZ_%sUp",effEName.Data())  , Form("histo_ZZ_%sUp",effEName.Data())  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  histo_ZZ_CMS_MVALepEffEBoundingDown             = new TH1D( Form("histo_ZZ_%sDown",effEName.Data()), Form("histo_ZZ_%sDown",effEName.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  histo_ggZH_hinv_CMS_MVALepEffEBoundingUp   = new TH1D( Form("histo_ggZH_hinv_%sUp",effEName.Data())  , Form("histo_ggZH_hinv_%sUp",effEName.Data())  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_MVALepEffEBoundingDown = new TH1D( Form("histo_ggZH_hinv_%sDown",effEName.Data()), Form("histo_ggZH_hinv_%sDown",effEName.Data()), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->Sumw2();

  histo_VVV_CMS_MVALepEffEBoundingAvg              = new TH1D( Form("histo_VVV_%sAvg",effEName.Data())  , Form("histo_VVV_%sAvg",effEName.Data())  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  histo_WZ_CMS_MVALepEffEBoundingAvg               = new TH1D( Form("histo_WZ_%sAvg",effEName.Data())  , Form("histo_WZ_%sAvg",effEName.Data())  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  histo_ZZ_CMS_MVALepEffEBoundingAvg               = new TH1D( Form("histo_ZZ_%sAvg",effEName.Data())  , Form("histo_ZZ_%sAvg",effEName.Data())  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg  = new TH1D( Form("histo_ggZH_hinv_%sAvg",effEName.Data())  , Form("histo_ggZH_hinv_%sAvg",effEName.Data())  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  // Make MET scale bounding histos
  histo_VVV_CMS_MVAMETBoundingUp           = new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  histo_VVV_CMS_MVAMETBoundingDown         = new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  histo_WZ_CMS_MVAMETBoundingUp            = new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  histo_WZ_CMS_MVAMETBoundingDown          = new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown"), nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();
  histo_ZZ_CMS_MVAMETBoundingUp            = new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  histo_ZZ_CMS_MVAMETBoundingDown          = new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  histo_ggZH_hinv_CMS_MVAMETBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_scale_metUp")  , Form("histo_ggZH_hinv_CMS_scale_metUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAMETBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_MVAMETBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_scale_metDown"), Form("histo_ggZH_hinv_CMS_scale_metDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAMETBoundingDown->Sumw2();

  // Make JES bounding histos
  histo_VVV_CMS_MVAJESBoundingUp              = new TH1D( Form("histo_VVV_CMS_eff_b_2016Up")  , Form("histo_VVV_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  histo_VVV_CMS_MVAJESBoundingDown            = new TH1D( Form("histo_VVV_CMS_eff_b_2016Down"), Form("histo_VVV_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  histo_WZ_CMS_MVAJESBoundingUp               = new TH1D( Form("histo_WZ_CMS_eff_b_2016Up")  , Form("histo_WZ_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  histo_WZ_CMS_MVAJESBoundingDown             = new TH1D( Form("histo_WZ_CMS_eff_b_2016Down"), Form("histo_WZ_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  histo_ZZ_CMS_MVAJESBoundingUp               = new TH1D( Form("histo_ZZ_CMS_eff_b_2016Up")  , Form("histo_ZZ_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  histo_ZZ_CMS_MVAJESBoundingDown             = new TH1D( Form("histo_ZZ_CMS_eff_b_2016Down"), Form("histo_ZZ_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();
  histo_ggZH_hinv_CMS_MVAJESBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_eff_b_2016Up")  , Form("histo_ggZH_hinv_CMS_eff_b_2016Up")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_MVAJESBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_eff_b_2016Down"), Form("histo_ggZH_hinv_CMS_eff_b_2016Down"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVAJESBoundingDown->Sumw2();

  // Make Btag bounding histos
  histo_VVV_CMS_MVABTAGBoundingUp              = new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingUp  ->Sumw2();
  histo_VVV_CMS_MVABTAGBoundingDown            = new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingDown->Sumw2();
  histo_WZ_CMS_MVABTAGBoundingUp               = new TH1D( Form("histo_WZ_CMS_scale_jUp")  , Form("histo_WZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  histo_WZ_CMS_MVABTAGBoundingDown             = new TH1D( Form("histo_WZ_CMS_scale_jDown"), Form("histo_WZ_CMS_scale_jDown"), nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingDown->Sumw2();
  histo_ZZ_CMS_MVABTAGBoundingUp               = new TH1D( Form("histo_ZZ_CMS_scale_jUp")  , Form("histo_ZZ_CMS_scale_jUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  histo_ZZ_CMS_MVABTAGBoundingDown             = new TH1D( Form("histo_ZZ_CMS_scale_jDown"), Form("histo_ZZ_CMS_scale_jDown"), nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingDown->Sumw2();
  histo_ggZH_hinv_CMS_MVABTAGBoundingUp   = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jUp")  , Form("histo_ggZH_hinv_CMS_scale_jUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_MVABTAGBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_MVABTAGBoundingDown = new TH1D( Form("histo_ggZH_hinv_CMS_scale_jDown"), Form("histo_ggZH_hinv_CMS_scale_jDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_MVABTAGBoundingDown->Sumw2();

  // Make BDT muons scale bounding histos
  histo_VVV_CMS_BDTMuonScaleBoundingUp           = new TH1D( Form("histo_VVV_CMS_bdt_muonUp")  , Form("histo_VVV_CMS_bdt_muonUp")  , nBinMVA, xbins);histo_VVV_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  histo_VVV_CMS_BDTMuonScaleBoundingDown         = new TH1D( Form("histo_VVV_CMS_bdt_muonDown"), Form("histo_VVV_CMS_bdt_muonDown"), nBinMVA, xbins);histo_VVV_CMS_BDTMuonScaleBoundingDown->Sumw2();
  histo_WZ_CMS_BDTMuonScaleBoundingUp            = new TH1D( Form("histo_WZ_CMS_bdt_muonUp")  , Form("histo_WZ_CMS_bdt_muonUp")  , nBinMVA, xbins);  histo_WZ_CMS_BDTMuonScaleBoundingUp   ->Sumw2();
  histo_WZ_CMS_BDTMuonScaleBoundingDown          = new TH1D( Form("histo_WZ_CMS_bdt_muonDown"), Form("histo_WZ_CMS_bdt_muonDown"), nBinMVA, xbins);  histo_WZ_CMS_BDTMuonScaleBoundingDown ->Sumw2();
  histo_ZZ_CMS_BDTMuonScaleBoundingUp            = new TH1D( Form("histo_ZZ_CMS_bdt_muonUp")  , Form("histo_ZZ_CMS_bdt_muonUp")  , nBinMVA, xbins);  histo_ZZ_CMS_BDTMuonScaleBoundingUp   ->Sumw2();
  histo_ZZ_CMS_BDTMuonScaleBoundingDown          = new TH1D( Form("histo_ZZ_CMS_bdt_muonDown"), Form("histo_ZZ_CMS_bdt_muonDown"), nBinMVA, xbins);  histo_ZZ_CMS_BDTMuonScaleBoundingDown ->Sumw2();
  histo_ggZH_hinv_CMS_BDTMuonScaleBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_muonUp")  , Form("histo_ggZH_hinv_CMS_bdt_muonUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMuonScaleBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_BDTMuonScaleBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_muonDown"), Form("histo_ggZH_hinv_CMS_bdt_muonDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMuonScaleBoundingDown->Sumw2();

  // Make BDT electron scale bounding histos
  histo_VVV_CMS_BDTElectronScaleBoundingUp           = new TH1D( Form("histo_VVV_CMS_bdt_electronUp")  , Form("histo_VVV_CMS_bdt_electronUp")  , nBinMVA, xbins);             histo_VVV_CMS_BDTElectronScaleBoundingUp         ->Sumw2();
  histo_VVV_CMS_BDTElectronScaleBoundingDown         = new TH1D( Form("histo_VVV_CMS_bdt_electronDown"), Form("histo_VVV_CMS_bdt_electronDown"), nBinMVA, xbins);             histo_VVV_CMS_BDTElectronScaleBoundingDown         ->Sumw2();
  histo_WZ_CMS_BDTElectronScaleBoundingUp            = new TH1D( Form("histo_WZ_CMS_bdt_electronUp")  , Form("histo_WZ_CMS_bdt_electronUp")  , nBinMVA, xbins);               histo_WZ_CMS_BDTElectronScaleBoundingUp         ->Sumw2();
  histo_WZ_CMS_BDTElectronScaleBoundingDown          = new TH1D( Form("histo_WZ_CMS_bdt_electronDown"), Form("histo_WZ_CMS_bdt_electronDown"), nBinMVA, xbins);               histo_WZ_CMS_BDTElectronScaleBoundingDown         ->Sumw2();
  histo_ZZ_CMS_BDTElectronScaleBoundingUp            = new TH1D( Form("histo_ZZ_CMS_bdt_electronUp")  , Form("histo_ZZ_CMS_bdt_electronUp")  , nBinMVA, xbins);               histo_ZZ_CMS_BDTElectronScaleBoundingUp         ->Sumw2();
  histo_ZZ_CMS_BDTElectronScaleBoundingDown          = new TH1D( Form("histo_ZZ_CMS_bdt_electronDown"), Form("histo_ZZ_CMS_bdt_electronDown"), nBinMVA, xbins);               histo_ZZ_CMS_BDTElectronScaleBoundingDown         ->Sumw2();
  histo_ggZH_hinv_CMS_BDTElectronScaleBoundingUp  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_electronUp")  , Form("histo_ggZH_hinv_CMS_bdt_electronUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTElectronScaleBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_BDTElectronScaleBoundingDown= new TH1D( Form("histo_ggZH_hinv_CMS_bdt_electronDown"), Form("histo_ggZH_hinv_CMS_bdt_electronDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTElectronScaleBoundingDown->Sumw2();

  // Make BDT met scale bounding histos
  histo_VVV_CMS_BDTMETScaleBoundingUp           = new TH1D( Form("histo_VVV_CMS_bdt_METUp")  , Form("histo_VVV_CMS_bdt_METUp")  , nBinMVA, xbins); histo_VVV_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  histo_VVV_CMS_BDTMETScaleBoundingDown         = new TH1D( Form("histo_VVV_CMS_bdt_METDown"), Form("histo_VVV_CMS_bdt_METDown"), nBinMVA, xbins); histo_VVV_CMS_BDTMETScaleBoundingDown->Sumw2();
  histo_WZ_CMS_BDTMETScaleBoundingUp            = new TH1D( Form("histo_WZ_CMS_bdt_METUp")  , Form("histo_WZ_CMS_bdt_METUp")  , nBinMVA, xbins);   histo_WZ_CMS_BDTMETScaleBoundingUp        ->Sumw2();
  histo_WZ_CMS_BDTMETScaleBoundingDown          = new TH1D( Form("histo_WZ_CMS_bdt_METDown"), Form("histo_WZ_CMS_bdt_METDown"), nBinMVA, xbins);   histo_WZ_CMS_BDTMETScaleBoundingDown ->Sumw2();
  histo_ZZ_CMS_BDTMETScaleBoundingUp            = new TH1D( Form("histo_ZZ_CMS_bdt_METUp")  , Form("histo_ZZ_CMS_bdt_METUp")  , nBinMVA, xbins);   histo_ZZ_CMS_BDTMETScaleBoundingUp        ->Sumw2();
  histo_ZZ_CMS_BDTMETScaleBoundingDown          = new TH1D( Form("histo_ZZ_CMS_bdt_METDown"), Form("histo_ZZ_CMS_bdt_METDown"), nBinMVA, xbins);   histo_ZZ_CMS_BDTMETScaleBoundingDown ->Sumw2();
  histo_ggZH_hinv_CMS_BDTMETScaleBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_METUp")  , Form("histo_ggZH_hinv_CMS_bdt_METUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMETScaleBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_BDTMETScaleBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_METDown"), Form("histo_ggZH_hinv_CMS_bdt_METDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTMETScaleBoundingDown->Sumw2();

  // Make BDT jet scale bounding histos
  histo_VVV_CMS_BDTJetScaleBoundingUp           = new TH1D( Form("histo_VVV_CMS_bdt_JESUp")  , Form("histo_VVV_CMS_bdt_JESUp")  , nBinMVA, xbins); histo_VVV_CMS_BDTJetScaleBoundingUp  ->Sumw2();
  histo_VVV_CMS_BDTJetScaleBoundingDown         = new TH1D( Form("histo_VVV_CMS_bdt_JESDown"), Form("histo_VVV_CMS_bdt_JESDown"), nBinMVA, xbins); histo_VVV_CMS_BDTJetScaleBoundingDown->Sumw2();
  histo_WZ_CMS_BDTJetScaleBoundingUp            = new TH1D( Form("histo_WZ_CMS_bdt_JESUp")  , Form("histo_WZ_CMS_bdt_JESUp")  , nBinMVA, xbins);   histo_WZ_CMS_BDTJetScaleBoundingUp        ->Sumw2();
  histo_WZ_CMS_BDTJetScaleBoundingDown          = new TH1D( Form("histo_WZ_CMS_bdt_JESDown"), Form("histo_WZ_CMS_bdt_JESDown"), nBinMVA, xbins);   histo_WZ_CMS_BDTJetScaleBoundingDown ->Sumw2();
  histo_ZZ_CMS_BDTJetScaleBoundingUp            = new TH1D( Form("histo_ZZ_CMS_bdt_JESUp")  , Form("histo_ZZ_CMS_bdt_JESUp")  , nBinMVA, xbins);   histo_ZZ_CMS_BDTJetScaleBoundingUp        ->Sumw2();
  histo_ZZ_CMS_BDTJetScaleBoundingDown          = new TH1D( Form("histo_ZZ_CMS_bdt_JESDown"), Form("histo_ZZ_CMS_bdt_JESDown"), nBinMVA, xbins);   histo_ZZ_CMS_BDTJetScaleBoundingDown ->Sumw2();
  histo_ggZH_hinv_CMS_BDTJetScaleBoundingUp    = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_JESUp")  , Form("histo_ggZH_hinv_CMS_bdt_JESUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTJetScaleBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_BDTJetScaleBoundingDown  = new TH1D( Form("histo_ggZH_hinv_CMS_bdt_JESDown"), Form("histo_ggZH_hinv_CMS_bdt_JESDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_BDTJetScaleBoundingDown->Sumw2();

  // Make pileup bounding histos
  histo_VVV_CMS_PUBoundingUp                   = new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  histo_VVV_CMS_PUBoundingDown                 = new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  histo_WZ_CMS_PUBoundingUp                    = new TH1D( Form("histo_WZ_CMS_puUp")  , Form("histo_WZ_CMS_puUp")  , nBinMVA, xbins); histo_WZ_CMS_PUBoundingUp  ->Sumw2();
  histo_WZ_CMS_PUBoundingDown                  = new TH1D( Form("histo_WZ_CMS_puDown"), Form("histo_WZ_CMS_puDown"), nBinMVA, xbins); histo_WZ_CMS_PUBoundingDown->Sumw2();
  histo_ZZ_CMS_PUBoundingUp                    = new TH1D( Form("histo_ZZ_CMS_puUp")  , Form("histo_ZZ_CMS_puUp")  , nBinMVA, xbins); histo_ZZ_CMS_PUBoundingUp  ->Sumw2();
  histo_ZZ_CMS_PUBoundingDown                  = new TH1D( Form("histo_ZZ_CMS_puDown"), Form("histo_ZZ_CMS_puDown"), nBinMVA, xbins); histo_ZZ_CMS_PUBoundingDown->Sumw2();
  histo_ggZH_hinv_CMS_PUBoundingUp        = new TH1D( Form("histo_ggZH_hinv_CMS_puUp")  , Form("histo_ggZH_hinv_CMS_puUp")  , nBinMVA, xbins); histo_ggZH_hinv_CMS_PUBoundingUp  ->Sumw2();
  histo_ggZH_hinv_CMS_PUBoundingDown      = new TH1D( Form("histo_ggZH_hinv_CMS_puDown"), Form("histo_ggZH_hinv_CMS_puDown"), nBinMVA, xbins); histo_ggZH_hinv_CMS_PUBoundingDown->Sumw2();
  
  // Make electroweak correction bounding histos
  histo_WZ_CMS_EWKCorrUp                    = new TH1D( Form("histo_WZ_EWKCorrUp")  , Form("histo_WZ_EWKCorrUp")  , nBinMVA, xbins); histo_WZ_CMS_EWKCorrUp  ->Sumw2();
  histo_WZ_CMS_EWKCorrDown                = new TH1D( Form("histo_WZ_EWKCorrDown"), Form("histo_WZ_EWKCorrDown"), nBinMVA, xbins); histo_WZ_CMS_EWKCorrDown->Sumw2();
  histo_ZZ_CMS_EWKCorrUp                  = new TH1D( Form("histo_ZZ_EWKCorrUp")  , Form("histo_ZZ_EWKCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_EWKCorrUp  ->Sumw2();
  histo_ZZ_CMS_EWKCorrDown                = new TH1D( Form("histo_ZZ_EWKCorrDown"), Form("histo_ZZ_EWKCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_EWKCorrDown->Sumw2();
  histo_ZZ_CMS_ggCorrUp                   = new TH1D( Form("histo_ZZ_ggCorrUp")  , Form("histo_ZZ_ggCorrUp")  , nBinMVA, xbins); histo_ZZ_CMS_ggCorrUp  ->Sumw2();
  histo_ZZ_CMS_ggCorrDown                 = new TH1D( Form("histo_ZZ_ggCorrDown"), Form("histo_ZZ_ggCorrDown"), nBinMVA, xbins); histo_ZZ_CMS_ggCorrDown->Sumw2();
  histo_Zjets_CMS_ZjetsSystUp                    = new TH1D( Form("histo_Zjets_ZjetsSystUp")  , Form("histo_Zjets_ZjetsSystUp")  , nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystUp  ->Sumw2();
  histo_Zjets_CMS_ZjetsSystDown           = new TH1D( Form("histo_Zjets_ZjetsSystDown"), Form("histo_Zjets_ZjetsSystDown"), nBinMVA, xbins); histo_Zjets_CMS_ZjetsSystDown->Sumw2();
  
  // Make the above bounding histos for each uncertainty for each signal model
  for(int nModel=0; nModel<nSigModels; nModel++) { 
    histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]          = new TH1D( Form("histo_ZH_hinv_%s_%sUp",   signalName_[nModel].Data(), effMName.Data()), Form("histo_ZH_hinv_%s_%sUp",  signalName_[nModel].Data(), effMName.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]        = new TH1D( Form("histo_ZH_hinv_%s_%sDown", signalName_[nModel].Data(), effMName.Data()), Form("histo_ZH_hinv_%s_%sDown",signalName_[nModel].Data(), effMName.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffMBoundingAvg [nModel]        = new TH1D( Form("histo_ZH_hinv_%s_%sAvg",             signalName_[nModel].Data(), effMName.Data()), Form("histo_ZH_hinv_%s_%sAvg" ,           signalName_[nModel].Data(), effMName.Data())  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingUp [nModel]         = new TH1D( Form("histo_ZH_hinv_%s_%sUp",                   signalName_[nModel].Data(), effEName.Data()), Form("histo_ZH_hinv_%s_%sUp"  ,           signalName_[nModel].Data(), effEName.Data())  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingDown [nModel]       = new TH1D( Form("histo_ZH_hinv_%s_%sDown",            signalName_[nModel].Data(), effEName.Data()), Form("histo_ZH_hinv_%s_%sDown",           signalName_[nModel].Data(), effEName.Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVALepEffEBoundingAvg [nModel]        = new TH1D( Form("histo_ZH_hinv_%s_%sAvg",             signalName_[nModel].Data(), effEName.Data()), Form("histo_ZH_hinv_%s_%sAvg" ,           signalName_[nModel].Data(), effEName.Data())  , nBinMVA, xbins); histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAMETBoundingUp [nModel]             = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_metUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_metUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAMETBoundingDown [nModel]           = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_metDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_metDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVAJESBoundingUp [nModel]             = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_jUp"         , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_jUp"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVAJESBoundingDown [nModel]           = new TH1D( Form("histo_ZH_hinv_%s_CMS_scale_jDown"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_scale_jDown"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_MVABTAGBoundingUp [nModel]            = new TH1D( Form("histo_ZH_hinv_%s_CMS_eff_b_2016Up"        , signalName_[nModel].Data()),  Form("histo_ZH_hinv_%s_CMS_eff_b_2016Up"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_MVABTAGBoundingDown [nModel]          = new TH1D( Form("histo_ZH_hinv_%s_CMS_eff_b_2016Down"  , signalName_[nModel].Data()),          Form("histo_ZH_hinv_%s_CMS_eff_b_2016Down"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTMuonScaleBoundingUp [nModel]       = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_muonUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_muonUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMuonScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTMuonScaleBoundingDown [nModel]     = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_muonDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_muonDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMuonScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTElectronScaleBoundingUp [nModel]   = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_electronUp"  , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_electronUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTElectronScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTElectronScaleBoundingDown [nModel] = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_electronDown", signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_bdt_electronDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTElectronScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTMETScaleBoundingUp [nModel]        = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_METUp"  , signalName_[nModel].Data()),                Form("histo_ZH_hinv_%s_CMS_bdt_METUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMETScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTMETScaleBoundingDown [nModel]      = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_METDown", signalName_[nModel].Data()),                Form("histo_ZH_hinv_%s_CMS_bdt_METDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTMETScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_BDTJetScaleBoundingUp [nModel]        = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_JESUp"  , signalName_[nModel].Data()),                Form("histo_ZH_hinv_%s_CMS_bdt_JESUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTJetScaleBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_BDTJetScaleBoundingDown [nModel]      = new TH1D( Form("histo_ZH_hinv_%s_CMS_bdt_JESDown", signalName_[nModel].Data()),                Form("histo_ZH_hinv_%s_CMS_bdt_JESDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_BDTJetScaleBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_PUBoundingUp [nModel]                 = new TH1D( Form("histo_ZH_hinv_%s_CMS_puUp"         , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_puUp"           , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_PUBoundingDown [nModel]               = new TH1D( Form("histo_ZH_hinv_%s_CMS_puDown"         , signalName_[nModel].Data()),           Form("histo_ZH_hinv_%s_CMS_puDown"           , signalName_[nModel].Data()), nBinMVA, xbins); histo_ZH_hinv_CMS_PUBoundingDown[nModel]->Sumw2();
    histo_ZH_hinv_CMS_EWKCorrUp[nModel]                     = new TH1D( Form("histo_ZH_hinv_%s_%sUp",   signalName_[nModel].Data(), "CMS_EWKCorr"), Form("histo_ZH_hinv_%s_%sUp",  signalName_[nModel].Data(), "CMS_EWKCorr"), nBinMVA, xbins); histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->Sumw2();
    histo_ZH_hinv_CMS_EWKCorrDown[nModel]                   = new TH1D( Form("histo_ZH_hinv_%s_%sDown", signalName_[nModel].Data(), "CMS_EWKCorr"), Form("histo_ZH_hinv_%s_%sDown",signalName_[nModel].Data(), "CMS_EWKCorr"), nBinMVA, xbins); histo_ZH_hinv_CMS_EWKCorrDown[nModel]->Sumw2();

  }
  unsigned long int t2 = static_cast<unsigned long int>(time(NULL));
  printf("zhAnalysis::MakeHistos : Complete! (%lu seconds)\n", t2-t1);
  madeHistos=true;
  return true;
}

TH1D* zhAnalysis::MakeHisto(unsigned int thePlot, TString &plotName) {
  TH1D *theHisto; unsigned nBinPlot; float xminPlot, xmaxPlot, pi=TMath::Pi();
       if(thePlot== 0) {nBinPlot=   7; xminPlot=-0.5; xmaxPlot=   6.5; plotName="passZMass nJet";} 
  else if(thePlot== 1) {nBinPlot=   7; xminPlot=-0.5; xmaxPlot=   6.5; plotName="passZMass nBJet";} 
  else if(thePlot== 2) {nBinPlot=   7; xminPlot=-0.5; xmaxPlot=   6.5; plotName="passZMass nLep";} 
  else if(thePlot== 3) {nBinPlot=  32; xminPlot=-0.1; xmaxPlot=   3.1; plotName="passZMass dPhi dilep MET";} 
  else if(thePlot== 4) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot= 200.0; plotName="PRESEL MET";}
  else if(thePlot== 5) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot=   1.0; plotName="PRESEL MET balance";}
  else if(thePlot== 6) {nBinPlot= 100; xminPlot=50.0; xmaxPlot= 250.0; plotName="PRESEL pTll";}
  else if(thePlot== 7) {nBinPlot=   4; xminPlot=-0.5; xmaxPlot=   3.5; plotName="N-1 nTau";}
  else if(thePlot== 8) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot=   1.0; plotName="N-1 jet ET balance";}
  else if(thePlot== 9) {nBinPlot= 500; xminPlot= 0.0; xmaxPlot=1000.0; plotName="N-1 mT";}
  else if(thePlot==10) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot= 100.0; plotName="N-1 dilep mass";}
  else if(thePlot==11) {nBinPlot=   7; xminPlot=-0.5; xmaxPlot=   6.5; plotName="N-1 nJet";}
  else if(thePlot==12) {nBinPlot= 200; xminPlot= 0.0; xmaxPlot= 400.0; plotName="N-1 MET";}
  else if(thePlot==13) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot=   1.0; plotName="N-1 MET balance";}
  else if(thePlot==14) {nBinPlot= 200; xminPlot= 0.0; xmaxPlot=    pi; plotName="N-1 dPhi dilep MET";}
  else if(thePlot==15) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot=   1.0; plotName="N-1 max CSV2";}
  else if(thePlot==16) {nBinPlot= 100; xminPlot=50.0; xmaxPlot= 250.0; plotName="N-1 dilep pT";}
  else if(thePlot==17) {nBinPlot= 200; xminPlot= 0.0; xmaxPlot=    pi; plotName="N-1 dPhi jet1 MET";}
  else if(thePlot==18) {nBinPlot= 200; xminPlot= 0.0; xmaxPlot=    pi; plotName="N-1 dPhi lep1 lep2";}
  else if(thePlot==19) {nBinPlot=  60; xminPlot=40.0; xmaxPlot= 100.0; plotName="N-1 MET";}
  else if(thePlot==20) {nBinPlot=  60; xminPlot= 0.0; xmaxPlot=   3.0; plotName="N-1 dR lep1 lep2";}
  else if(thePlot==21) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot=1000.0; plotName="TIGHTSEL mTll";}
  else if(thePlot==22) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot= 200.0; plotName="TIGHTSEL lep1 pT";}
  else if(thePlot==23) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot= 200.0; plotName="TIGHTSEL lep2 pT";}
  else if(thePlot==24) {nBinPlot=  40; xminPlot=-0.5; xmaxPlot=  39.5; plotName="TIGHTSEL NPV";}
  else if(thePlot==25) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot=   1.0; plotName="TIGHTSEL pTll/mTll";}
  else if(thePlot==26) {nBinPlot= 100; xminPlot=-1.0; xmaxPlot=   1.0; plotName="TIGHTSEL MET150 BDT";}
  else if(thePlot==27) {nBinPlot= 200; xminPlot= 0.0; xmaxPlot=    pi; plotName="TIGHTSEL min dPhi lep MET";}
  else if(thePlot==28) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot=   2.5; plotName="TIGHTSEL dilep eta";}
  else if(thePlot==29) {nBinPlot=  32; xminPlot=-0.1; xmaxPlot=   3.1; plotName="DYSANESEL1 dPhi jet1 MET";} 
  else if(thePlot==30) {nBinPlot=  40; xminPlot= 0.0; xmaxPlot=   2.0; plotName="DYSANESEL1 calo PF balance";} 
  else if(thePlot==31) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot= 500.0; plotName="DYSANESEL1 MET";} 
  else if(thePlot==32) {nBinPlot=  32; xminPlot=-0.1; xmaxPlot=   3.1; plotName="DYSANESEL2 dPhi jet1 MET";} 
  else if(thePlot==33) {nBinPlot=  40; xminPlot= 0.0; xmaxPlot=   2.0; plotName="DYSANESEL2 calo PF balance";} 
  else if(thePlot==34) {nBinPlot= 100; xminPlot= 0.0; xmaxPlot= 500.0; plotName="DYSANESEL2 MET";} 
  else if(thePlot==35) {nBinPlot=  50; xminPlot= 0.0; xmaxPlot=   2.5; plotName="ZHGSEL pho abs eta";}
  else if(thePlot==allPlots-2) {nBinPlot =  numberCuts+1; xminPlot =-0.5; xmaxPlot =  numberCuts+0.5; plotName="Cut flow";}
  else if(thePlot==allPlots-1) plotName="Shape analysis";
  else { printf("error with zhAnalysis::MakeHisto: thePlot out of bounds (%d). check allPlots in zhAnalysis.h\n", thePlot); assert(0); theHisto=new TH1D; return theHisto;}
  if(thePlot != allPlots-1) theHisto = new TH1D("theHisto", "theHisto", nBinPlot, xminPlot, xmaxPlot);
  else                      theHisto = new TH1D("theHisto", "theHisto", nBinMVA, xbins);
  theHisto->Sumw2();
  return theHisto;
}
bool zhAnalysis::SaveHistos() {
  printf("zhAnalysis::SaveHistos : Begin writing the plot and shape histograms to files\n");
  unsigned long int t1 = static_cast<unsigned long int>(time(NULL));
  char output[200];
  
  // Save plotting histos
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    sprintf(output,"MitZHAnalysis/plots%s/histo%szh%s_nice_%s_%d.root",subdirectory.c_str(),addChan.Data(),finalStateName, signalName_[plotModel].Data(),thePlot);          
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    if(!outFilePlotsNote || !outFilePlotsNote->IsOpen()) {
      printf("Error in zhAnalysis::SaveHistos : could not open \"%s\" for writing\n", output);
      return false;
    }
    outFilePlotsNote->cd();
    for(int np=0; np<processTypes; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }
  
  // Save the shapes and the datacards

  double process_syst[nSigModels][7], yield_processTypes[nSigModels][processTypes+1],
         stat_processTypes[nSigModels][processTypes+1], syst_processTypes[nSigModels][processTypes+1],
         syst_types_allBackground[28];
  char outputLimits[200];
  for(int nModel=0; nModel<nSigModels; nModel++) { // Output the limits for all the models
    if(nModel > 0 && nModel != plotModel && MVAVarType==3) continue;
    sprintf(outputLimits,"MitZHAnalysis/plots%s/zll%shinv%s_%s_input_%s.root", subdirectory.c_str(), addChan.Data(), finalStateName, signalName_[nModel].Data(), ECMsb.Data());
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
    
    if(verbose) {
      cout << histo_Data             ->GetSumOfWeights() << " ";
      cout << histo_ZH_hinv[nModel]  ->GetSumOfWeights() << " ";
      cout << histo_Zjets            ->GetSumOfWeights() << " ";
      cout << histo_VVV              ->GetSumOfWeights() << " ";
      cout << histo_WZ               ->GetSumOfWeights() << " ";
      cout << histo_ZZ               ->GetSumOfWeights() << " ";
      cout << histo_EM               ->GetSumOfWeights() << " ";
      cout << histo_ggZH_hinv        ->GetSumOfWeights() << " ";
      cout << endl;
      printf("uncertainties Stat\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]         ->GetBinContent(i)/histo_ZH_hinv[nModel]  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel]         ->GetBinContent(i)/histo_ZH_hinv[nModel]  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_VVV        ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp          ->GetBinContent(i)/histo_VVV           ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_VVV        ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown         ->GetBinContent(i)/histo_VVV           ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_WZ        ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp                 ->GetBinContent(i)/histo_WZ       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_WZ        ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown          ->GetBinContent(i)/histo_WZ           ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZZ        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp                 ->GetBinContent(i)/histo_ZZ       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_ZZ        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown          ->GetBinContent(i)/histo_ZZ           ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_EM        ->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingUp                 ->GetBinContent(i)/histo_EM       ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]  ->GetNbinsX(); i++) {if(histo_EM        ->GetBinContent(i)>0)printf("%5.1f ",histo_EM_CMS_MVAEMStatBoundingDown          ->GetBinContent(i)/histo_EM           ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp   ->GetBinContent(i)/histo_ggZH_hinv->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown ->GetBinContent(i)/histo_ggZH_hinv->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffM\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]  ->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffE\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]  ->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->GetBinContent(i)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_ggZH_hinv_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties MET\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]      ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]   ->GetBinContent(i)/histo_ZH_hinv[nModel]   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAMETBoundingDown   ->GetBinContent(i)/histo_ggZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties JES\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]    ->GetBinContent(i)/histo_ZH_hinv[nModel]     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]         ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAJESBoundingUp    ->GetBinContent(i)/histo_ggZH_hinv     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVAJESBoundingDown         ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties BTAG\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]    ->GetBinContent(i)/histo_ZH_hinv[nModel]     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]         ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVABTAGBoundingUp    ->GetBinContent(i)/histo_ggZH_hinv     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_MVABTAGBoundingDown         ->GetBinContent(i)/histo_ggZH_hinv      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties PU\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_PUBoundingUp[nModel]      ->GetBinContent(i)/histo_ZH_hinv[nModel]      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]        ->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_PUBoundingDown[nModel]   ->GetBinContent(i)/histo_ZH_hinv[nModel]   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_PUBoundingUp          ->GetBinContent(i)/histo_ggZH_hinv          ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ggZH_hinv->GetNbinsX(); i++) {if(histo_ggZH_hinv        ->GetBinContent(i)>0)printf("%5.1f ",histo_ggZH_hinv_CMS_PUBoundingDown   ->GetBinContent(i)/histo_ggZH_hinv   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties EWK Corr\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_EWKCorrUp  ->GetBinContent(i)/histo_WZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_WZ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_EWKCorrDown->GetBinContent(i)/histo_WZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_EWKCorrDown->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ggCorrUp  ->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZZ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_ggCorrDown->GetBinContent(i)/histo_ZZ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->GetBinContent(i)/histo_ZH_hinv[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_ZH_hinv[nModel]->GetBinContent(i)>0)printf("%5.1f ",histo_ZH_hinv_CMS_EWKCorrDown[nModel]->GetBinContent(i)/histo_ZH_hinv[nModel]->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties Zjets\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystUp  ->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZH_hinv[nModel]->GetNbinsX(); i++) {if(histo_Zjets->GetBinContent(i)>0)printf("%5.1f ",histo_Zjets_CMS_ZjetsSystDown->GetBinContent(i)/histo_Zjets->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
    } 

    if(useBDT) {
      for(int nb=1; nb<=nBinMVA; nb++) {
        for(unsigned int i_toy=0; i_toy<num_bdt_toys; i_toy++) {
          int nbin2d = histo_bdt_toys_electronScale_VVV->GetBin(nb, i_toy);
          if(histo_VVV            ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_VVV      ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_VVV            ->GetBinContent(nbin2d) / histo_VVV            ->GetBinContent(nb));
          if(histo_WZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_WZ       ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_WZ             ->GetBinContent(nbin2d) / histo_WZ             ->GetBinContent(nb));
          if(histo_ZZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_ZZ       ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_ZZ             ->GetBinContent(nbin2d) / histo_ZZ             ->GetBinContent(nb));
          if(histo_ZH_hinv[nModel]->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_ZH_hinv[nModel]  ->Fill( xbins[nb-1], histo_bdt_toys_electronScale_ZH_hinv[nModel]->GetBinContent(nbin2d) / histo_ZH_hinv[nModel]->GetBinContent(nb));
          if(histo_ggZH_hinv      ->GetBinContent(nb) != 0 ) bdt_toy_envelope_electronScale_ggZH_hinv->Fill( xbins[nb-1], histo_bdt_toys_electronScale_ggZH_hinv      ->GetBinContent(nbin2d) / histo_ggZH_hinv      ->GetBinContent(nb));
          if(histo_VVV            ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_VVV      ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_VVV            ->GetBinContent(nbin2d) / histo_VVV            ->GetBinContent(nb));
          if(histo_WZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_WZ       ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_WZ             ->GetBinContent(nbin2d) / histo_WZ             ->GetBinContent(nb));
          if(histo_ZZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_ZZ       ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_ZZ             ->GetBinContent(nbin2d) / histo_ZZ             ->GetBinContent(nb));
          if(histo_ZH_hinv[nModel]->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_ZH_hinv[nModel]  ->Fill( xbins[nb-1], histo_bdt_toys_muonScale_ZH_hinv[nModel]->GetBinContent(nbin2d) / histo_ZH_hinv[nModel]->GetBinContent(nb));
          if(histo_ggZH_hinv      ->GetBinContent(nb) != 0 ) bdt_toy_envelope_muonScale_ggZH_hinv->Fill( xbins[nb-1], histo_bdt_toys_muonScale_ggZH_hinv      ->GetBinContent(nbin2d) / histo_ggZH_hinv      ->GetBinContent(nb));
          if(histo_VVV            ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_VVV      ->Fill( xbins[nb-1], histo_bdt_toys_METScale_VVV            ->GetBinContent(nbin2d) / histo_VVV            ->GetBinContent(nb));
          if(histo_WZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_WZ       ->Fill( xbins[nb-1], histo_bdt_toys_METScale_WZ             ->GetBinContent(nbin2d) / histo_WZ             ->GetBinContent(nb));
          if(histo_ZZ             ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_ZZ       ->Fill( xbins[nb-1], histo_bdt_toys_METScale_ZZ             ->GetBinContent(nbin2d) / histo_ZZ             ->GetBinContent(nb));
          if(histo_ZH_hinv[nModel]->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_ZH_hinv[nModel]  ->Fill( xbins[nb-1], histo_bdt_toys_METScale_ZH_hinv[nModel]->GetBinContent(nbin2d) / histo_ZH_hinv[nModel]->GetBinContent(nb));
          if(histo_ggZH_hinv      ->GetBinContent(nb) != 0 ) bdt_toy_envelope_METScale_ggZH_hinv->Fill( xbins[nb-1], histo_bdt_toys_METScale_ggZH_hinv      ->GetBinContent(nbin2d) / histo_ggZH_hinv      ->GetBinContent(nb));
        }
        bdt_toy_binyields_electronScale_VVV[nb-1]             = (TH1F*) bdt_toy_envelope_electronScale_VVV            ->ProjectionY(Form("bdt_toy_binyields_electronScale_VVV_bin%d", nb)      , nb, nb);
        bdt_toy_binyields_electronScale_WZ[nb-1]              = (TH1F*) bdt_toy_envelope_electronScale_WZ             ->ProjectionY(Form("bdt_toy_binyields_electronScale_WZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_electronScale_ZZ[nb-1]              = (TH1F*) bdt_toy_envelope_electronScale_ZZ             ->ProjectionY(Form("bdt_toy_binyields_electronScale_ZZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_electronScale_ZH_hinv[nModel][nb-1] = (TH1F*) bdt_toy_envelope_electronScale_ZH_hinv[nModel]->ProjectionY(Form("bdt_toy_binyields_electronScale_ZH_hinv_%s_bin%d",signalName_[nModel].Data(),nb), nb, nb);
        bdt_toy_binyields_electronScale_ggZH_hinv[nb-1]       = (TH1F*) bdt_toy_envelope_electronScale_ggZH_hinv      ->ProjectionY(Form("bdt_toy_binyields_electronScale_ggZH_hinv_bin%d", nb), nb, nb);
        bdt_toy_binyields_muonScale_VVV[nb-1]                 = (TH1F*) bdt_toy_envelope_muonScale_VVV                ->ProjectionY(Form("bdt_toy_binyields_muonScale_VVV_bin%d", nb)      , nb, nb);
        bdt_toy_binyields_muonScale_WZ[nb-1]                  = (TH1F*) bdt_toy_envelope_muonScale_WZ                 ->ProjectionY(Form("bdt_toy_binyields_muonScale_WZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_muonScale_ZZ[nb-1]                  = (TH1F*) bdt_toy_envelope_muonScale_ZZ                 ->ProjectionY(Form("bdt_toy_binyields_muonScale_ZZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_muonScale_ZH_hinv[nModel][nb-1]     = (TH1F*) bdt_toy_envelope_muonScale_ZH_hinv[nModel]    ->ProjectionY(Form("bdt_toy_binyields_muonScale_ZH_hinv_%s_bin%d", signalName_[nModel].Data(),nb), nb, nb);
        bdt_toy_binyields_muonScale_ggZH_hinv[nb-1]           = (TH1F*) bdt_toy_envelope_muonScale_ggZH_hinv          ->ProjectionY(Form("bdt_toy_binyields_muonScale_ggZH_hinv_bin%d", nb), nb, nb);
        bdt_toy_binyields_METScale_VVV[nb-1]                  = (TH1F*) bdt_toy_envelope_METScale_VVV                 ->ProjectionY(Form("bdt_toy_binyields_METScale_VVV_bin%d", nb)      , nb, nb);
        bdt_toy_binyields_METScale_WZ[nb-1]                   = (TH1F*) bdt_toy_envelope_METScale_WZ                  ->ProjectionY(Form("bdt_toy_binyields_METScale_WZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_METScale_ZZ[nb-1]                   = (TH1F*) bdt_toy_envelope_METScale_ZZ                  ->ProjectionY(Form("bdt_toy_binyields_METScale_ZZ_bin%d", nb)       , nb, nb);
        bdt_toy_binyields_METScale_ZH_hinv[nModel][nb-1]      = (TH1F*) bdt_toy_envelope_METScale_ZH_hinv[nModel]     ->ProjectionY(Form("bdt_toy_binyields_METScale_ZH_hinv_%s_bin%d", signalName_[nModel].Data(),nb), nb, nb);
        bdt_toy_binyields_METScale_ggZH_hinv[nb-1]            = (TH1F*) bdt_toy_envelope_METScale_ggZH_hinv           ->ProjectionY(Form("bdt_toy_binyields_METScale_ggZH_hinv_bin%d", nb), nb, nb);
        bdt_toy_binyields_electronScale_VVV[nb-1]             ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_electronScale_WZ[nb-1]              ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_electronScale_ZZ[nb-1]              ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_electronScale_ZH_hinv[nModel][nb-1] ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_electronScale_ggZH_hinv[nb-1]       ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_muonScale_VVV[nb-1]                 ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_muonScale_WZ[nb-1]                  ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_muonScale_ZZ[nb-1]                  ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_muonScale_ZH_hinv[nModel][nb-1]     ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_muonScale_ggZH_hinv[nb-1]           ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_METScale_VVV[nb-1]                  ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_METScale_WZ[nb-1]                   ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_METScale_ZZ[nb-1]                   ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_METScale_ZH_hinv[nModel][nb-1]      ->GetXaxis()->SetRangeUser(0.8,1.199);
        bdt_toy_binyields_METScale_ggZH_hinv[nb-1]            ->GetXaxis()->SetRangeUser(0.8,1.199);
 
        // Since not using cached results, finally calculate the bin-by-bin % systematics (up and down) for the BDT 
        double theQuantilesElectron[3], theQuantilesMuon[3], theQuantilesMET[3];
        if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 ) { 
          bdt_toy_binyields_electronScale_ZH_hinv[nModel][nb-1]->GetQuantiles(3, theQuantilesElectron, quantileProbs);
          bdt_toy_binyields_muonScale_ZH_hinv[nModel][nb-1]    ->GetQuantiles(3, theQuantilesMuon    , quantileProbs);
          bdt_toy_binyields_METScale_ZH_hinv[nModel][nb-1]     ->GetQuantiles(3, theQuantilesMET     , quantileProbs);
          bdt_syst_electronScaleUp_ZH_hinv[nModel]   ->SetBinContent(nb, theQuantilesElectron[2]);
          bdt_syst_electronScaleDown_ZH_hinv[nModel] ->SetBinContent(nb, theQuantilesElectron[0]);
          bdt_syst_muonScaleUp_ZH_hinv[nModel]       ->SetBinContent(nb, theQuantilesMuon[2]);
          bdt_syst_muonScaleDown_ZH_hinv[nModel]     ->SetBinContent(nb, theQuantilesMuon[0]);
          bdt_syst_METScaleUp_ZH_hinv[nModel]        ->SetBinContent(nb, theQuantilesMET[2]);
          bdt_syst_METScaleDown_ZH_hinv[nModel]      ->SetBinContent(nb, theQuantilesMET[0]);
        }
        if(histo_VVV->GetBinContent(nb) > 0 ) {
          bdt_toy_binyields_electronScale_VVV[nb-1]->GetQuantiles(3, theQuantilesElectron, quantileProbs);
          bdt_toy_binyields_muonScale_VVV[nb-1]    ->GetQuantiles(3, theQuantilesMuon    , quantileProbs);
          bdt_toy_binyields_METScale_VVV[nb-1]     ->GetQuantiles(3, theQuantilesMET     , quantileProbs);
          bdt_syst_electronScaleUp_VVV   ->SetBinContent(nb, theQuantilesElectron[2]);
          bdt_syst_electronScaleDown_VVV ->SetBinContent(nb, theQuantilesElectron[0]);
          bdt_syst_muonScaleUp_VVV       ->SetBinContent(nb, theQuantilesMuon[2]);
          bdt_syst_muonScaleDown_VVV     ->SetBinContent(nb, theQuantilesMuon[0]);
          bdt_syst_METScaleUp_VVV        ->SetBinContent(nb, theQuantilesMET[2]);
          bdt_syst_METScaleDown_VVV      ->SetBinContent(nb, theQuantilesMET[0]);
        }
        if(histo_WZ->GetBinContent(nb) > 0 ) {
          bdt_toy_binyields_electronScale_WZ[nb-1]->GetQuantiles(3, theQuantilesElectron, quantileProbs);
          bdt_toy_binyields_muonScale_WZ[nb-1]    ->GetQuantiles(3, theQuantilesMuon    , quantileProbs);
          bdt_toy_binyields_METScale_WZ[nb-1]     ->GetQuantiles(3, theQuantilesMET     , quantileProbs);
          bdt_syst_electronScaleUp_WZ   ->SetBinContent(nb, theQuantilesElectron[2]);
          bdt_syst_electronScaleDown_WZ ->SetBinContent(nb, theQuantilesElectron[0]);
          bdt_syst_muonScaleUp_WZ       ->SetBinContent(nb, theQuantilesMuon[2]);
          bdt_syst_muonScaleDown_WZ     ->SetBinContent(nb, theQuantilesMuon[0]);
          bdt_syst_METScaleUp_WZ        ->SetBinContent(nb, theQuantilesMET[2]);
          bdt_syst_METScaleDown_WZ      ->SetBinContent(nb, theQuantilesMET[0]);
        }
        if(histo_ZZ->GetBinContent(nb) > 0 ) {
          bdt_toy_binyields_electronScale_ZZ[nb-1]->GetQuantiles(3, theQuantilesElectron, quantileProbs);
          bdt_toy_binyields_muonScale_ZZ[nb-1]    ->GetQuantiles(3, theQuantilesMuon    , quantileProbs);
          bdt_toy_binyields_METScale_ZZ[nb-1]     ->GetQuantiles(3, theQuantilesMET     , quantileProbs);
          bdt_syst_electronScaleUp_ZZ   ->SetBinContent(nb, theQuantilesElectron[2]);
          bdt_syst_electronScaleDown_ZZ ->SetBinContent(nb, theQuantilesElectron[0]);
          bdt_syst_muonScaleUp_ZZ       ->SetBinContent(nb, theQuantilesMuon[2]);
          bdt_syst_muonScaleDown_ZZ     ->SetBinContent(nb, theQuantilesMuon[0]);
          bdt_syst_METScaleUp_ZZ        ->SetBinContent(nb, theQuantilesMET[2]);
          bdt_syst_METScaleDown_ZZ      ->SetBinContent(nb, theQuantilesMET[0]);
        }
        if(histo_ggZH_hinv->GetBinContent(nb) > 0 ) {
          bdt_toy_binyields_electronScale_ggZH_hinv[nb-1]->GetQuantiles(3, theQuantilesElectron, quantileProbs);
          bdt_toy_binyields_muonScale_ggZH_hinv[nb-1]    ->GetQuantiles(3, theQuantilesMuon    , quantileProbs);
          bdt_toy_binyields_METScale_ggZH_hinv[nb-1]     ->GetQuantiles(3, theQuantilesMET     , quantileProbs);
          bdt_syst_electronScaleUp_ggZH_hinv   ->SetBinContent(nb, theQuantilesElectron[2]);
          bdt_syst_electronScaleDown_ggZH_hinv ->SetBinContent(nb, theQuantilesElectron[0]);
          bdt_syst_muonScaleUp_ggZH_hinv       ->SetBinContent(nb, theQuantilesMuon[2]);
          bdt_syst_muonScaleDown_ggZH_hinv     ->SetBinContent(nb, theQuantilesMuon[0]);
          bdt_syst_METScaleUp_ggZH_hinv        ->SetBinContent(nb, theQuantilesMET[2]);
          bdt_syst_METScaleDown_ggZH_hinv      ->SetBinContent(nb, theQuantilesMET[0]);
        }
      }
    }

    histo_Data             ->Write();
    histo_ZH_hinv[nModel]  ->Write();
    histo_Zjets            ->Write();
    histo_VVV              ->Write();
    histo_WZ               ->Write();
    histo_ZZ               ->Write();
    histo_EM               ->Write();
    histo_ggZH_hinv        ->Write();
    
    histo_ZH_hinv_CMS_MVAZHStatBoundingUp[nModel]   ->Write();
    histo_ZH_hinv_CMS_MVAZHStatBoundingDown[nModel] ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingUp                    ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingDown            ->Write();
    histo_WZ_CMS_MVAWZStatBoundingUp                    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingDown                    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingUp                    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingDown                    ->Write();
    histo_EM_CMS_MVAEMStatBoundingUp                    ->Write();
    histo_EM_CMS_MVAEMStatBoundingDown                    ->Write();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingUp       ->Write();
    histo_ggZH_hinv_CMS_MVAggZHStatBoundingDown     ->Write();
    
    histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]  ->Write();
    histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->Write();
    histo_VVV_CMS_MVALepEffMBoundingUp                    ->Write();
    histo_VVV_CMS_MVALepEffMBoundingDown            ->Write();
    histo_WZ_CMS_MVALepEffMBoundingUp                    ->Write();
    histo_WZ_CMS_MVALepEffMBoundingDown                    ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingUp                    ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingDown                    ->Write();
    histo_ggZH_hinv_CMS_MVALepEffMBoundingUp        ->Write();
    histo_ggZH_hinv_CMS_MVALepEffMBoundingDown      ->Write();
    
    histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]  ->Write();
    histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->Write();
    histo_VVV_CMS_MVALepEffEBoundingUp                    ->Write();
    histo_VVV_CMS_MVALepEffEBoundingDown            ->Write();
    histo_WZ_CMS_MVALepEffEBoundingUp                    ->Write();
    histo_WZ_CMS_MVALepEffEBoundingDown                    ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingUp                    ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingDown                    ->Write();
    histo_ggZH_hinv_CMS_MVALepEffEBoundingUp        ->Write();
    histo_ggZH_hinv_CMS_MVALepEffEBoundingDown      ->Write();
    
    histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]      ->Write();
    histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]    ->Write();
    histo_VVV_CMS_MVAMETBoundingUp                    ->Write();
    histo_VVV_CMS_MVAMETBoundingDown                    ->Write();
    histo_WZ_CMS_MVAMETBoundingUp                    ->Write();
    histo_WZ_CMS_MVAMETBoundingDown                    ->Write();
    histo_ZZ_CMS_MVAMETBoundingUp                    ->Write();
    histo_ZZ_CMS_MVAMETBoundingDown                    ->Write();
    histo_ggZH_hinv_CMS_MVAMETBoundingUp            ->Write();
    histo_ggZH_hinv_CMS_MVAMETBoundingDown          ->Write();
    
    histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]            ->Write();
    histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]    ->Write(); 
    histo_VVV_CMS_MVAJESBoundingUp                    ->Write();
    histo_VVV_CMS_MVAJESBoundingDown                    ->Write();
    histo_WZ_CMS_MVAJESBoundingUp                     ->Write();
    histo_WZ_CMS_MVAJESBoundingDown                    ->Write();
    histo_ZZ_CMS_MVAJESBoundingUp                     ->Write();
    histo_ZZ_CMS_MVAJESBoundingDown                    ->Write();
    histo_ggZH_hinv_CMS_MVAJESBoundingUp            ->Write();
    histo_ggZH_hinv_CMS_MVAJESBoundingDown          ->Write();
    
    histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]            ->Write();
    histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]   ->Write(); 
    histo_VVV_CMS_MVABTAGBoundingUp                    ->Write();
    histo_VVV_CMS_MVABTAGBoundingDown                    ->Write();
    histo_WZ_CMS_MVABTAGBoundingUp                     ->Write();
    histo_WZ_CMS_MVABTAGBoundingDown                    ->Write();
    histo_ZZ_CMS_MVABTAGBoundingUp                     ->Write();
    histo_ZZ_CMS_MVABTAGBoundingDown                    ->Write();
    histo_ggZH_hinv_CMS_MVABTAGBoundingUp            ->Write();
    histo_ggZH_hinv_CMS_MVABTAGBoundingDown         ->Write();
    
    histo_ZH_hinv_CMS_PUBoundingUp[nModel]          ->Write();
    histo_ZH_hinv_CMS_PUBoundingDown[nModel]        ->Write();
    histo_VVV_CMS_PUBoundingUp                            ->Write();
    histo_VVV_CMS_PUBoundingDown                    ->Write();
    histo_WZ_CMS_PUBoundingUp                            ->Write();
    histo_WZ_CMS_PUBoundingDown                            ->Write();
    histo_ZZ_CMS_PUBoundingUp                            ->Write();
    histo_ZZ_CMS_PUBoundingDown                            ->Write();
    histo_ggZH_hinv_CMS_PUBoundingUp                ->Write();
    histo_ggZH_hinv_CMS_PUBoundingDown              ->Write();
    
    histo_WZ_CMS_EWKCorrUp                            ->Write();
    histo_WZ_CMS_EWKCorrDown                            ->Write();
    histo_ZZ_CMS_EWKCorrUp                            ->Write();
    histo_ZZ_CMS_EWKCorrDown                            ->Write();
    histo_ZZ_CMS_ggCorrUp                            ->Write();
    histo_ZZ_CMS_ggCorrDown                            ->Write();
    histo_ZH_hinv_CMS_EWKCorrUp[nModel]                    ->Write();
    histo_ZH_hinv_CMS_EWKCorrDown[nModel]            ->Write();
    histo_Zjets_CMS_ZjetsSystUp                            ->Write();
    histo_Zjets_CMS_ZjetsSystDown                    ->Write();
    if(useBDT) {
      bdt_syst_electronScaleUp_VVV       ->Write();
      bdt_syst_electronScaleUp_WZ        ->Write();
      bdt_syst_electronScaleUp_ZZ        ->Write();
      bdt_syst_electronScaleUp_ggZH_hinv ->Write();
      bdt_syst_electronScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_muonScaleUp_VVV           ->Write();
      bdt_syst_muonScaleUp_WZ            ->Write();
      bdt_syst_muonScaleUp_ZZ            ->Write();
      bdt_syst_muonScaleUp_ggZH_hinv     ->Write();
      bdt_syst_muonScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_METScaleUp_VVV            ->Write();
      bdt_syst_METScaleUp_WZ             ->Write();
      bdt_syst_METScaleUp_ZZ             ->Write();
      bdt_syst_METScaleUp_ggZH_hinv      ->Write();
      bdt_syst_METScaleUp_ZH_hinv[nModel] ->Write();
      bdt_syst_electronScaleDown_VVV       ->Write();
      bdt_syst_electronScaleDown_WZ        ->Write();
      bdt_syst_electronScaleDown_ZZ        ->Write();
      bdt_syst_electronScaleDown_ggZH_hinv ->Write();
      bdt_syst_electronScaleDown_ZH_hinv[nModel] ->Write();
      bdt_syst_muonScaleDown_VVV           ->Write();
      bdt_syst_muonScaleDown_WZ            ->Write();
      bdt_syst_muonScaleDown_ZZ            ->Write();
      bdt_syst_muonScaleDown_ggZH_hinv     ->Write();
      bdt_syst_muonScaleDown_ZH_hinv[nModel] ->Write();
      bdt_syst_METScaleDown_VVV            ->Write();
      bdt_syst_METScaleDown_WZ             ->Write();
      bdt_syst_METScaleDown_ZZ             ->Write();
      bdt_syst_METScaleDown_ggZH_hinv      ->Write();
      bdt_syst_METScaleDown_ZH_hinv[nModel] ->Write();   
    }
    outFileLimits->Close();
    
    // Record the BDT shape systs in a separate file for the caching (only save plotModel for now)
    if(useBDT && !useCachedBDTSystematics) assert(SaveBDTSystematics(plotModel));
    double lumiE = 1.025;
    double systLepResE[5] = {1.01,1.01,1.01,1.01,1.01};
    double systLepResM[5] = {1.01,1.01,1.01,1.01,1.01};
    double syst_btag = 1.02;
    double syst_WZl[2] = {1.010, 1.003};
    if(nJetsType > 0) syst_WZl[1] = 1.012;
    
    for(int nb=1; nb<=nBinMVA; nb++){
      double nggZHEvt = histo_ggZH_hinv->GetBinContent(nb);
      if(nModel != 0) nggZHEvt = 0.0;
      // QCD study
      double systQCDScale[5] = {TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nModel][0]->GetBinContent(nb)-histo_ZH_hinv[nModel]->GetBinContent(nb)),
                                TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]            ->GetBinContent(nb)-histo_VVV            ->GetBinContent(nb)),
                                TMath::Abs(histo_WZ_CMS_QCDScaleBounding[0]             ->GetBinContent(nb)-histo_WZ             ->GetBinContent(nb)),
                                TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0]             ->GetBinContent(nb)-histo_ZZ             ->GetBinContent(nb)),
                                TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[0]        ->GetBinContent(nb)-histo_ggZH_hinv      ->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nqcd]->GetBinContent(nb)  -histo_ZH_hinv[nModel]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_ZH_hinv_CMS_QCDScaleBounding[nModel][nqcd]->GetBinContent(nb)  -histo_ZH_hinv[nModel]  ->GetBinContent(nb));
        if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)              -histo_VVV    ->GetBinContent(nb))         > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)      -histo_VVV      ->GetBinContent(nb));
        if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)               -histo_WZ     ->GetBinContent(nb))         > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_WZ       ->GetBinContent(nb));
        if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)               -histo_ZZ     ->GetBinContent(nb))         > systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)       -histo_ZZ       ->GetBinContent(nb));
        if(TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)        -histo_ggZH_hinv->GetBinContent(nb))       > systQCDScale[4]) systQCDScale[4] = TMath::Abs(histo_ggZH_hinv_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggZH_hinv->GetBinContent(nb));
      }                 
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_ZH_hinv[nModel]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histo_VVV->GetBinContent(nb)             > 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_VVV            ->GetBinContent(nb); else systQCDScale[1] = 1;
      if(histo_WZ->GetBinContent(nb)              > 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_WZ             ->GetBinContent(nb); else systQCDScale[2] = 1;
      if(histo_ZZ->GetBinContent(nb)              > 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_ZZ             ->GetBinContent(nb); else systQCDScale[3] = 1;
      if(histo_ggZH_hinv->GetBinContent(nb)       > 0) systQCDScale[4] = 1 + systQCDScale[4]/histo_ggZH_hinv      ->GetBinContent(nb); else systQCDScale[4] = 1;
      for(int ntype=0; ntype<5; ntype++) if(systQCDScale[ntype] < 0) systQCDScale[ntype] = 1.0;
      if(verbose) printf("QCDScale(%d): %f %f %f %f %f - ",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[4]);

      if(systQCDScale[0]-1 > qcdScaleTotal[0]) systQCDScale[0] = 1 + sqrt((systQCDScale[0]-1.0)*(systQCDScale[0]-1.0)-qcdScaleTotal[0]*qcdScaleTotal[0]); else systQCDScale[0] = 1.0;
      if(systQCDScale[4]-1 > qcdScaleTotal[1]) systQCDScale[4] = 1 + sqrt((systQCDScale[4]-1.0)*(systQCDScale[4]-1.0)-qcdScaleTotal[1]*qcdScaleTotal[1]); else systQCDScale[4] = 1.0;
      if(verbose) printf("%f/%f %f/%f\n",systQCDScale[0],qcdScaleTotal[0],systQCDScale[4],qcdScaleTotal[1]);

      // PDF study
      double systPDFUp[5], systPDFDown[5];
      systPDFUp  [0] = (histo_ZH_hinv[nModel]  ->GetBinContent(nb)>0)?
        1.+sqrt(TMath::Max( pow((histo_ZH_hinv_CMS_PDFUp[nModel]->GetBinContent(nb) - histo_ZH_hinv[nModel]->GetBinContent(nb))/histo_ZH_hinv[nModel]->GetBinContent(nb),2) - pdfTotal[0]*pdfTotal[0], 0.)): 1.;
      systPDFUp  [0] = (histo_ZH_hinv[nModel]  ->GetBinContent(nb)>0)?
        1.-sqrt(TMath::Max( pow((histo_ZH_hinv[nModel]->GetBinContent(nb) - histo_ZH_hinv_CMS_PDFUp[nModel]->GetBinContent(nb) )/histo_ZH_hinv[nModel]->GetBinContent(nb),2) - pdfTotal[0]*pdfTotal[0], 0.)): 1.;
      systPDFUp  [1] = (histo_VVV      ->GetBinContent(nb)>0)? 1.+(histo_VVV_CMS_PDFUp      ->GetBinContent(nb) - histo_VVV    ->GetBinContent(nb))/histo_VVV    ->GetBinContent(nb):1.;
      systPDFDown[1] = (histo_VVV      ->GetBinContent(nb)>0)? 1.+(histo_VVV_CMS_PDFDown    ->GetBinContent(nb) - histo_VVV    ->GetBinContent(nb))/histo_VVV    ->GetBinContent(nb):1.;
      systPDFUp  [2] = (histo_WZ       ->GetBinContent(nb)>0)? 1.+(histo_WZ_CMS_PDFUp       ->GetBinContent(nb) - histo_WZ     ->GetBinContent(nb))/histo_WZ     ->GetBinContent(nb):1.;
      systPDFDown[2] = (histo_WZ       ->GetBinContent(nb)>0)? 1.+(histo_WZ_CMS_PDFDown     ->GetBinContent(nb) - histo_WZ     ->GetBinContent(nb))/histo_WZ     ->GetBinContent(nb):1.;
      systPDFUp  [3] = (histo_ZZ       ->GetBinContent(nb)>0)? 1.+(histo_ZZ_CMS_PDFUp       ->GetBinContent(nb) - histo_ZZ     ->GetBinContent(nb))/histo_ZZ     ->GetBinContent(nb):1.;
      systPDFDown[3] = (histo_ZZ       ->GetBinContent(nb)>0)? 1.+(histo_ZZ_CMS_PDFDown     ->GetBinContent(nb) - histo_ZZ     ->GetBinContent(nb))/histo_ZZ     ->GetBinContent(nb):1.;
      systPDFUp  [4] = (histo_ggZH_hinv  ->GetBinContent(nb)>0)?
        1.+sqrt(TMath::Max( pow((histo_ggZH_hinv_CMS_PDFUp->GetBinContent(nb) - histo_ggZH_hinv->GetBinContent(nb))/histo_ggZH_hinv->GetBinContent(nb),2) - pdfTotal[1]*pdfTotal[1], 0.)): 1.;
      systPDFUp  [4] = (histo_ggZH_hinv  ->GetBinContent(nb)>0)?
        1.-sqrt(TMath::Max( pow((histo_ggZH_hinv->GetBinContent(nb) - histo_ggZH_hinv_CMS_PDFUp->GetBinContent(nb) )/histo_ggZH_hinv->GetBinContent(nb),2) - pdfTotal[1]*pdfTotal[1], 0.)): 1.;
      if(verbose) printf("PDF(%d): %f/%f %f/%f %f/%f %f/%f %f/%f\n",nb,systPDFUp[0],systPDFDown[0],systPDFUp[1],systPDFDown[1],systPDFUp[2],systPDFDown[2],systPDFUp[3],systPDFDown[3],systPDFUp[4],systPDFDown[4]); // ZH, VVV, WZ, ZZ, ggZH
  
      double systLepEffM[5] = {1.0,1.0,1.0,1.0,1.0};
      if     (histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]     ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZH_hinv_CMS_MVALepEffMBoundingUp[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb);
      else if(histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_ZH_hinv_CMS_MVALepEffMBoundingAvg[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffMBoundingDown[nModel]->GetBinContent(nb);
      if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingUp         ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)        > 0 && histo_VVV_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffMBoundingUp          ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_WZ_CMS_MVALepEffMBoundingDown        ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp          ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)         > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown        ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[4] = histo_ggZH_hinv_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
      else if(histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ggZH_hinv_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[4] = histo_ggZH_hinv_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ggZH_hinv_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
  
      double systLepEffE[5] = {1.0,1.0,1.0,1.0,1.0};
      if     (histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]     ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZH_hinv_CMS_MVALepEffEBoundingUp[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb);
      else if(histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb)    > 0 && histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_ZH_hinv_CMS_MVALepEffEBoundingAvg[nModel]->GetBinContent(nb)/histo_ZH_hinv_CMS_MVALepEffEBoundingDown[nModel]->GetBinContent(nb);
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
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]     ->GetBinContent(nb) > 0) systMetUp  [0] = histo_ZH_hinv_CMS_MVAMETBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]   ->GetBinContent(nb) > 0) systMetDown[0] = histo_ZH_hinv_CMS_MVAMETBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
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
      for(int nmet=0; nmet<5; nmet++) if(systMetUp[nmet]   > 1.04) systMetUp[nmet]   = 1.04;
      for(int nmet=0; nmet<5; nmet++) if(systMetUp[nmet]   < 0.96) systMetUp[nmet]   = 0.96;
      for(int nmet=0; nmet<5; nmet++) if(systMetDown[nmet] > 1.04) systMetDown[nmet] = 1.04;
      for(int nmet=0; nmet<5; nmet++) if(systMetDown[nmet] < 0.96) systMetDown[nmet] = 0.96;
  
      double systJesUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systJesDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]        ->GetBinContent(nb) > 0) systJesUp  [0] = histo_ZH_hinv_CMS_MVAJESBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]        ->GetBinContent(nb) > 0) systJesDown[0] = histo_ZH_hinv_CMS_MVAJESBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)          > 0 && histo_VVV_CMS_MVAJESBoundingUp         ->GetBinContent(nb) > 0) systJesUp  [1] = histo_VVV_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)          > 0 && histo_VVV_CMS_MVAJESBoundingDown        ->GetBinContent(nb) > 0) systJesDown[1] = histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)          > 0 && histo_WZ_CMS_MVAJESBoundingUp          ->GetBinContent(nb) > 0) systJesUp  [2] = histo_WZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)          > 0 && histo_WZ_CMS_MVAJESBoundingDown        ->GetBinContent(nb) > 0) systJesDown[2] = histo_WZ_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)          > 0 && histo_ZZ_CMS_MVAJESBoundingUp          ->GetBinContent(nb) > 0) systJesUp  [3] = histo_ZZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)          > 0 && histo_ZZ_CMS_MVAJESBoundingDown        ->GetBinContent(nb) > 0) systJesDown[3] = histo_ZZ_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAJESBoundingUp        ->GetBinContent(nb) > 0) systJesUp  [4] = histo_ggZH_hinv_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJesDown[4] = histo_ggZH_hinv_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      for(int njes=0; njes<5; njes++) if(systJesUp[njes]   > 1.10) systJesUp[njes]   = 1.10;
      for(int njes=0; njes<5; njes++) if(systJesUp[njes]   < 0.90) systJesUp[njes]   = 0.90;
      for(int njes=0; njes<5; njes++) if(systJesDown[njes] > 1.10) systJesDown[njes] = 1.10;
      for(int njes=0; njes<5; njes++) if(systJesDown[njes] < 0.90) systJesDown[njes] = 0.90;

      double systBtagUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBtagDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]        ->GetBinContent(nb) > 0) systBtagUp  [0] = histo_ZH_hinv_CMS_MVABTAGBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]        ->GetBinContent(nb) > 0) systBtagDown[0] = histo_ZH_hinv_CMS_MVABTAGBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)          > 0 && histo_VVV_CMS_MVABTAGBoundingUp            ->GetBinContent(nb) > 0) systBtagUp  [1] = histo_VVV_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb)          > 0 && histo_VVV_CMS_MVABTAGBoundingDown           ->GetBinContent(nb) > 0) systBtagDown[1] = histo_VVV_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)          > 0 && histo_WZ_CMS_MVABTAGBoundingUp             ->GetBinContent(nb) > 0) systBtagUp  [2] = histo_WZ_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb)          > 0 && histo_WZ_CMS_MVABTAGBoundingDown           ->GetBinContent(nb) > 0) systBtagDown[2] = histo_WZ_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)          > 0 && histo_ZZ_CMS_MVABTAGBoundingUp             ->GetBinContent(nb) > 0) systBtagUp  [3] = histo_ZZ_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb)          > 0 && histo_ZZ_CMS_MVABTAGBoundingDown           ->GetBinContent(nb) > 0) systBtagDown[3] = histo_ZZ_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVABTAGBoundingUp   ->GetBinContent(nb) > 0) systBtagUp  [4] = histo_ggZH_hinv_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0 && histo_ggZH_hinv_CMS_MVABTAGBoundingDown ->GetBinContent(nb) > 0) systBtagDown[4] = histo_ggZH_hinv_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_ggZH_hinv->GetBinContent(nb);
      for(int njes=0; njes<5; njes++) if(systBtagUp[njes]   > 1.10) systBtagUp[njes]   = 1.10;
      for(int njes=0; njes<5; njes++) if(systBtagUp[njes]   < 0.90) systBtagUp[njes]   = 0.90;
      for(int njes=0; njes<5; njes++) if(systBtagDown[njes] > 1.10) systBtagDown[njes] = 1.10;
      for(int njes=0; njes<5; njes++) if(systBtagDown[njes] < 0.90) systBtagDown[njes] = 0.90;

      // BDT systematic study
      double systBDTMuonUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTMuonDown[5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTElectronUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTElectronDown[5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTMETUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systBDTMETDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(useBDT) { 
        systBDTElectronUp  [0] = bdt_syst_electronScaleUp_ZH_hinv[nModel]   ->GetBinContent(nb);
        systBDTElectronDown[0] = bdt_syst_electronScaleDown_ZH_hinv[nModel] ->GetBinContent(nb);
        systBDTMuonUp  [0]     = bdt_syst_muonScaleUp_ZH_hinv[nModel]       ->GetBinContent(nb);
        systBDTMuonDown[0]     = bdt_syst_muonScaleDown_ZH_hinv[nModel]     ->GetBinContent(nb);
        systBDTMETUp  [0]      = bdt_syst_METScaleUp_ZH_hinv[nModel]        ->GetBinContent(nb);
        systBDTMETDown[0]      = bdt_syst_METScaleDown_ZH_hinv[nModel]      ->GetBinContent(nb);
        systBDTElectronUp  [1] = bdt_syst_electronScaleUp_VVV   ->GetBinContent(nb);
        systBDTElectronDown[1] = bdt_syst_electronScaleDown_VVV ->GetBinContent(nb);
        systBDTMuonUp  [1]     = bdt_syst_muonScaleUp_VVV       ->GetBinContent(nb);
        systBDTMuonDown[1]     = bdt_syst_muonScaleDown_VVV     ->GetBinContent(nb);
        systBDTMETUp  [1]      = bdt_syst_METScaleUp_VVV        ->GetBinContent(nb);
        systBDTMETDown[1]      = bdt_syst_METScaleDown_VVV      ->GetBinContent(nb);
        systBDTElectronUp  [2] = bdt_syst_electronScaleUp_WZ   ->GetBinContent(nb);
        systBDTElectronDown[2] = bdt_syst_electronScaleDown_WZ ->GetBinContent(nb);
        systBDTMuonUp  [2]     = bdt_syst_muonScaleUp_WZ       ->GetBinContent(nb);
        systBDTMuonDown[2]     = bdt_syst_muonScaleDown_WZ     ->GetBinContent(nb);
        systBDTMETUp  [2]      = bdt_syst_METScaleUp_WZ        ->GetBinContent(nb);
        systBDTMETDown[2]      = bdt_syst_METScaleDown_WZ      ->GetBinContent(nb);
        systBDTElectronUp  [3] = bdt_syst_electronScaleUp_ZZ   ->GetBinContent(nb);
        systBDTElectronDown[3] = bdt_syst_electronScaleDown_ZZ ->GetBinContent(nb);
        systBDTMuonUp  [3]     = bdt_syst_muonScaleUp_ZZ       ->GetBinContent(nb);
        systBDTMuonDown[3]     = bdt_syst_muonScaleDown_ZZ     ->GetBinContent(nb);
        systBDTMETUp  [3]      = bdt_syst_METScaleUp_ZZ        ->GetBinContent(nb);
        systBDTMETDown[3]      = bdt_syst_METScaleDown_ZZ      ->GetBinContent(nb);
        systBDTElectronUp  [4] = bdt_syst_electronScaleUp_ggZH_hinv   ->GetBinContent(nb);
        systBDTElectronDown[4] = bdt_syst_electronScaleDown_ggZH_hinv ->GetBinContent(nb);
        systBDTMuonUp  [4]     = bdt_syst_muonScaleUp_ggZH_hinv       ->GetBinContent(nb);
        systBDTMuonDown[4]     = bdt_syst_muonScaleDown_ggZH_hinv     ->GetBinContent(nb);
        systBDTMETUp  [4]      = bdt_syst_METScaleUp_ggZH_hinv        ->GetBinContent(nb);
        systBDTMETDown[4]      = bdt_syst_METScaleDown_ggZH_hinv      ->GetBinContent(nb);
        
      }
   
      double systPUUp  [5] = {1.0,1.0,1.0,1.0,1.0};
      double systPUDown[5] = {1.0,1.0,1.0,1.0,1.0};
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_PUBoundingUp[nModel]     ->GetBinContent(nb) > 0) systPUUp  [0] = histo_ZH_hinv_CMS_PUBoundingUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb)   > 0 && histo_ZH_hinv_CMS_PUBoundingDown[nModel]   ->GetBinContent(nb) > 0) systPUDown[0] = histo_ZH_hinv_CMS_PUBoundingDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
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
  
      double syst_EWKCorrUp[4]   = {1.0,1.0,1.0,1.0}; // WZ, ZZ, ggZZ, ZH
      double syst_EWKCorrDown[4] = {1.0,1.0,1.0,1.0};
      if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [0] = histo_WZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_WZ->GetBinContent(nb) > 0 &&  histo_WZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[0] = histo_WZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [1] = histo_ZZ_CMS_EWKCorrUp  ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb) > 0) syst_EWKCorrDown[1] = histo_ZZ_CMS_EWKCorrDown->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb) > 0) syst_EWKCorrUp  [2] = histo_ZZ_CMS_ggCorrUp   ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZZ->GetBinContent(nb) > 0 &&  histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb) > 0) syst_EWKCorrDown[2] = histo_ZZ_CMS_ggCorrDown ->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0 &&  histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->GetBinContent(nb) > 0) syst_EWKCorrUp  [3] = histo_ZH_hinv_CMS_EWKCorrUp[nModel]  ->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0 &&  histo_ZH_hinv_CMS_EWKCorrDown[nModel]->GetBinContent(nb) > 0) syst_EWKCorrDown[3] = histo_ZH_hinv_CMS_EWKCorrDown[nModel]->GetBinContent(nb)/histo_ZH_hinv[nModel]->GetBinContent(nb);

      // If VV normalization is from data
      double systVV[2] = {0,0};
      if(useVVFromData){
        systVV[0] = sqrt(TMath::Power(syst_EWKCorrUp[0]-1,2)+TMath::Power(systQCDScale[2]-1,2)+TMath::Power((systPDFUp[2]-systPDFDown[2])/2.,2));
        systVV[1] = sqrt(TMath::Power(syst_EWKCorrUp[1]-1,2)+TMath::Power(systQCDScale[3]-1,2)+TMath::Power((systPDFUp[3]-systPDFDown[3])/2.,2));
        // this commented stuff is not right
        //printf("systVV_theo = %f(%f/%f/%f) %f(%f/%f/%f)\n",systVV[0],syst_EWKCorrUp[0]-1,systQCDScale[2]-1,systPDF[2]-1,
        //                                                   systVV[1],syst_EWKCorrUp[1]-1,systQCDScale[3]-1,systPDF[3]-1);
      }

      char outputLimitsShape[200];                                            
      sprintf(outputLimitsShape,"MitZHAnalysis/datacards%s/histo_limits_zll%shinv%s_%s_shape_%s_bin%d.txt",subdirectory.c_str(), addChan.Data(), finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb-1);
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
      newcardShape << Form("bin hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d hinv%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process ZH_hinv Zjets VVV WZ ZZ EM ggZH_hinv\n");
      newcardShape << Form("process 0 1 2 3 4 5 -1\n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",histo_ZH_hinv[nModel]->GetBinContent(nb),histo_Zjets->GetBinContent(nb),TMath::Max(histo_VVV->GetBinContent(nb),0.0),histo_WZ->GetBinContent(nb),histo_ZZ->GetBinContent(nb),histo_EM->GetBinContent(nb),nggZHEvt);
      newcardShape << Form("lumi_%4s                               lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE);                     
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effMName.Data(),systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3],systLepEffM[4]);
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",effEName.Data(),systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3],systLepEffE[4]);
      if(MVAVarType != 3) {
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momMName.Data(),systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3],systLepResM[4]);
      newcardShape << Form("%s                                     lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",momEName.Data(),systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3],systLepResE[4]);
      }
      newcardShape << Form("CMS_pu2016                             lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systPUUp[0],systPUDown[0],systPUUp[1],systPUDown[1],systPUUp[2],systPUDown[2],systPUUp[3],systPUDown[3],systPUUp[0],systPUDown[0]); // 0 --> 4
      newcardShape << Form("CMS_scale_met                          lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systMetUp[0],systMetDown[0],systMetUp[1],systMetDown[1],systMetUp[2],systMetDown[2],systMetUp[3],systMetDown[3],systMetUp[0],systMetDown[0]); // 0 --> 4
      newcardShape << Form("CMS_scale_j                            lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systJesUp[0],systJesDown[0],systJesUp[1],systJesDown[1],systJesUp[2],systJesDown[2],systJesUp[3],systJesDown[3],systJesUp[0],systJesDown[0]); // 0 --> 4                 
      newcardShape << Form("CMS_trigger2016                        lnN  %7.5f   -   %7.5f %7.5f %7.5f   -   %7.5f\n",1.01,1.01,1.01,1.01,1.01);
      newcardShape << Form("UEPS                                   lnN  1.030   -     -     -     -     -   1.030\n");
      newcardShape << Form("CMS_eff_b_2016                         lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -   %7.5f/%7.5f\n",systBtagUp[0],systBtagDown[0],systBtagUp[1],systBtagDown[1],systBtagUp[2],systBtagDown[2],systBtagUp[3],systBtagDown[3],systBtagUp[0],systBtagDown[0]); // 0 --> 4
      newcardShape << Form("pdf_qqbar_ACCEPT                       lnN  %7.5f/%7.5f   -   %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -     -  \n",systPDFUp[0],systPDFDown[0],systPDFUp[1],systPDFDown[1],systPDFUp[2],systPDFDown[2],systPDFUp[3],systPDFDown[3]);
      newcardShape << Form("pdf_gg_ACCEPT                          lnN    -     -     -     -     -     -   %7.5f/%7.5f\n",systPDFUp[4],systPDFDown[4]);
      newcardShape << Form("pdf_qqbar                              lnN  %7.5f   -     -     -     -     -     -  \n",1.0+pdfTotal[0]);
      if(histo_ggZH_hinv->GetBinContent(nb) > 0)
      newcardShape << Form("pdf_gg                                 lnN    -     -     -     -     -     -   %7.5f\n",1.0+pdfTotal[1]);
      if(systQCDScale[0] != 1.0)
      newcardShape << Form("QCDscale_VH_ACCEPT                         lnN  %7.5f   -     -     -     -     -     -  \n",systQCDScale[0]);  
      if(systQCDScale[4] != 1.0 && histo_ggZH_hinv->GetBinContent(nb) > 0)
      newcardShape << Form("QCDscale_ggVH_ACCEPT                         lnN    -     -     -     -     -     -   %7.5f\n",systQCDScale[4]);  
      newcardShape << Form("QCDscale_VH                                 lnN  %7.5f   -     -     -     -     -     -  \n",1.0+qcdScaleTotal[0]);  
      if(histo_ggZH_hinv->GetBinContent(nb) > 0)
      newcardShape << Form("QCDscale_ggVH                                 lnN    -     -     -     -     -     -   %7.5f\n",1.0+qcdScaleTotal[1]);  
      newcardShape << Form("QCDscale_VVV                                 lnN    -     -   %7.5f   -     -     -     -  \n",systQCDScale[1]);                
      newcardShape << Form("QCDscale_WZ                                   lnN    -     -     -   %7.5f      -    -      -  \n",systQCDScale[2]);
      newcardShape << Form("QCDscale_ZZ                                   lnN    -     -     -   -      %7.5f    -      -  \n",systQCDScale[3]);
      newcardShape << Form("CMS_zllhinv_ZH_EWKCorr                 lnN    %7.5f/%7.5f      -     -   -  -      -      -  \n",syst_EWKCorrUp[3],syst_EWKCorrDown[3]);                

      if(useVVFromData && nb != 1){
      if(useZZWZEWKUnc){
      newcardShape << Form("CMS_hinv_vvnorm_bin%d rateParam  * WZ 1 [0.1,10]\n",nb-1);                
      newcardShape << Form("CMS_hinv_vvnorm_bin%d rateParam  * ZZ 1 [0.1,10]\n",nb-1);        
      newcardShape << Form("CMS_zllhinv_ZZWZ_EWKCorr               lnN    -     -     -   -      %7.5f    -      -\n",1.+sqrt(0.02*0.02+(syst_EWKCorrUp[1]-1.0)*(syst_EWKCorrUp[1]-1.0)));                
      } else {
      newcardShape << Form("CMS_hinv_wznorm_bin%d rateParam  * WZ 1 [0.1,10]\n",nb-1);                
      newcardShape << Form("CMS_hinv_wznorm_bin%d param 1 %5.3f\n",nb-1, systVV[0]);                
      newcardShape << Form("CMS_hinv_zznorm_bin%d rateParam  * ZZ 1 [0.1,10]\n",nb-1);        
      newcardShape << Form("CMS_hinv_zznorm_bin%d param 1 %5.3f\n",nb-1, systVV[1]);        
      }
      }
      else if(!useVVFromData){
      newcardShape << Form("CMS_zllhinv_WZ_EWKCorr                 lnN    -     -     -   %7.5f/%7.5f   -      -      -  \n",syst_EWKCorrUp[0],syst_EWKCorrDown[0]);                
      newcardShape << Form("CMS_zllhinv_ZZ_EWKCorr                 lnN    -     -     -     -    %7.5f/%7.5f   -      -  \n",syst_EWKCorrUp[1],syst_EWKCorrDown[1]);                
      }

      newcardShape << Form("CMS_zllhinv_ggZZCorr                   lnN    -     -     -     -   %7.5f/%7.5f   -     -  \n",syst_EWKCorrUp[2],syst_EWKCorrDown[2]);                

      newcardShape << Form("CMS_zllhinv_WZ_lep2016                 lnN     -         -     -   %7.5f   -          -    -  \n",syst_WZl[0]);            
      newcardShape << Form("CMS_zllhinv_WZ_tau2016                 lnN     -         -     -   %7.5f   -          -    -  \n",syst_WZl[1]);            
      if(MVAVarType == 3 || MVAVarType==4) {
      newcardShape << Form("CMS_BDT_scale_m                   lnN        %7.5f/%7.5f -  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f -  %7.5f/%7.5f   \n", systBDTMuonUp[0], systBDTMuonDown[0], systBDTMuonUp[1], systBDTMuonDown[1], systBDTMuonUp[2], systBDTMuonDown[2], systBDTMuonUp[3], systBDTMuonDown[3], systBDTMuonUp[4], systBDTMuonDown[4]);            
      newcardShape << Form("CMS_BDT_scale_e                   lnN        %7.5f/%7.5f -  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f -  %7.5f/%7.5f   \n", systBDTElectronUp[0], systBDTElectronDown[0], systBDTElectronUp[1], systBDTElectronDown[1], systBDTElectronUp[2], systBDTElectronDown[2], systBDTElectronUp[3], systBDTElectronDown[3], systBDTElectronUp[4], systBDTElectronDown[4]);            
      newcardShape << Form("CMS_BDT_scale_MET                 lnN        %7.5f/%7.5f -  %7.5f/%7.5f   %7.5f/%7.5f   %7.5f/%7.5f -  %7.5f/%7.5f   \n", systBDTMETUp[0], systBDTMETDown[0], systBDTMETUp[1], systBDTMETDown[1], systBDTMETUp[2], systBDTMETDown[2], systBDTMETUp[3], systBDTMETDown[3], systBDTMETUp[4], systBDTMETDown[4]);            
      }
      if(nb>=3){
      newcardShape << Form("CMS_zllhinv_ZLLNorm2016_%s_%s              lnN        -   %7.5f   -          -     -     -     -  \n",finalStateName,ECMsb.Data(),2.0);            
      //newcardShape << Form("CMS_zllhinv_ZLLShape2016_%s_%s             lnN        -   %7.5f/%7.5f   -        -     -     -     -  \n",finalStateName,ECMsb.Data(),systZjetsUp[0],systZjetsDown[0]);            
      }
      if(MVAVarType == 1 || MVAVarType == 2 || MVAVarType == 4){
      newcardShape << Form("CMS_zllhinv_ZLLFit2016_%s_%s               lnU        -   %7.5f   -          -     -     -     -  \n",finalStateName,ECMsb.Data(),2.0);            
      }
      if(nb>=2)
      newcardShape << Form("CMS_zllhinv_EMSyst2016_%s_%s               lnN        -     -     -          -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),systEM[0]);               

      newcardShape << Form("CMS_zllhinv_EMNorm2016_%s_%s               lnU        -     -     -          -     -   %7.5f   -  \n",finalStateName,ECMsb.Data(),2.0);      
  
      if(histo_ZH_hinv[nModel]->GetBinContent(nb) > 0) newcardShape << Form("CMS_zllhinv%s_MVAZHStatBounding2016_%s_Bin%d    lnN    %7.5f -      -         -    -    -         -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ZH_hinv[nModel]->GetBinError(nb)  /histo_ZH_hinv[nModel]->GetBinContent(nb),0.999));
  
      if(histo_Zjets->GetBinContent(nb)           > 0) newcardShape << Form("CMS_zllhinv%s_MVAZjetsStatBounding2016_%s_Bin%d  lnN        -  %7.5f   -        -    -    -        -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_Zjets->GetBinError(nb)    /histo_Zjets->GetBinContent(nb),0.999));
  
      if(histo_VVV->GetBinContent(nb)             > 0) newcardShape << Form("CMS_zllhinv%s_MVAVVVStatBounding2016_%s_Bin%d    lnN        -   -  %7.5f   -    -         -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_VVV->GetBinError(nb)           /histo_VVV->GetBinContent(nb),0.999));
  
      if(histo_WZ->GetBinContent(nb)              > 0) newcardShape << Form("CMS_zllhinv%s_MVAWZStatBounding2016_%s_Bin%d     lnN      -  -     -  %7.5f  -    -     -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WZ->GetBinError(nb)        /histo_WZ->GetBinContent(nb),0.999));
  
      if(histo_ZZ->GetBinContent(nb)                  > 0) newcardShape << Form("CMS_zllhinv%s_MVAZZStatBounding2016_%s_Bin%d     lnN          -  -     -        -  %7.5f  -        -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ZZ->GetBinError(nb)            /histo_ZZ->GetBinContent(nb),0.999));
  
      if(histo_EM->GetBinContent(nb)              > 0) newcardShape << Form("CMS_zllhinv%s_MVAEMStatBounding2016_%s_Bin%d     lnN          -  -     -        -    -  %7.5f        -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_EM->GetBinError(nb)            /histo_EM->GetBinContent(nb),0.999));
  
      if(histo_ggZH_hinv->GetBinContent(nb)       > 0) newcardShape << Form("CMS_zllhinv%s_MVAggZHStatBounding2016_%s_Bin%d   lnN      -    -     -    -    -    -   %7.5f\n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ggZH_hinv->GetBinError(nb)/histo_ggZH_hinv->GetBinContent(nb),0.999));
    }

  } // end loop over models

  unsigned long int t2 = static_cast<unsigned long int>(time(NULL));
  printf("zhAnalysis::SaveHistos : Saved plots to \"%s\" and wrote the shapes (%lu seconds)\n", output, t2-t1);
  return true;
}
bool zhAnalysis::LoadCachedBDTSystematics() {
  printf("zhAnalysis::LoadCachedBDTSystematics : Opening cached systematics in \"%s\"\n", filenameBDTSysts);
  unsigned long int t1 = static_cast<unsigned long int>(time(NULL));
  cachedSystFile = new TFile(filenameBDTSysts,"read");
  if(!cachedSystFile || !cachedSystFile->IsOpen()) {
    printf("Error in zhAnalysis::LoadCachedBDTSystematics : the file does not exist\n");
    return false;
  }
  bdt_syst_electronScaleUp_VVV       = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleUp_VVV"      );
  bdt_syst_electronScaleUp_WZ        = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleUp_WZ"       );
  bdt_syst_electronScaleUp_ZZ        = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleUp_ZZ"       );
  bdt_syst_electronScaleUp_ggZH_hinv = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleUp_ggZH_hinv");
  bdt_syst_electronScaleUp_ZH_hinv[plotModel] = (TH1F*) cachedSystFile->Get(Form("bdt_syst_electronScaleUp_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_syst_muonScaleUp_VVV           = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleUp_VVV"          );
  bdt_syst_muonScaleUp_WZ            = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleUp_WZ"           );
  bdt_syst_muonScaleUp_ZZ            = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleUp_ZZ"           );
  bdt_syst_muonScaleUp_ggZH_hinv     = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleUp_ggZH_hinv"    );
  bdt_syst_muonScaleUp_ZH_hinv[plotModel] = (TH1F*) cachedSystFile->Get(Form("bdt_syst_muonScaleUp_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_syst_METScaleUp_VVV            = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleUp_VVV"           );
  bdt_syst_METScaleUp_WZ             = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleUp_WZ"            );
  bdt_syst_METScaleUp_ZZ             = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleUp_ZZ"            );
  bdt_syst_METScaleUp_ggZH_hinv      = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleUp_ggZH_hinv"     );
  bdt_syst_METScaleUp_ZH_hinv[plotModel] = (TH1F*) cachedSystFile->Get(Form("bdt_syst_METScaleUp_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_syst_electronScaleDown_VVV       = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleDown_VVV"      );
  bdt_syst_electronScaleDown_WZ        = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleDown_WZ"       );
  bdt_syst_electronScaleDown_ZZ        = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleDown_ZZ"       );
  bdt_syst_electronScaleDown_ggZH_hinv = (TH1F*) cachedSystFile->Get("bdt_syst_electronScaleDown_ggZH_hinv");
  bdt_syst_electronScaleDown_ZH_hinv[plotModel] = (TH1F*) cachedSystFile->Get(Form("bdt_syst_electronScaleDown_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_syst_muonScaleDown_VVV           = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleDown_VVV"          );
  bdt_syst_muonScaleDown_WZ            = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleDown_WZ"           );
  bdt_syst_muonScaleDown_ZZ            = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleDown_ZZ"           );
  bdt_syst_muonScaleDown_ggZH_hinv     = (TH1F*) cachedSystFile->Get("bdt_syst_muonScaleDown_ggZH_hinv"    );
  bdt_syst_muonScaleDown_ZH_hinv[plotModel] = (TH1F*) cachedSystFile->Get(Form("bdt_syst_muonScaleDown_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_syst_METScaleDown_VVV            = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleDown_VVV"           );
  bdt_syst_METScaleDown_WZ             = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleDown_WZ"            );
  bdt_syst_METScaleDown_ZZ             = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleDown_ZZ"            );
  bdt_syst_METScaleDown_ggZH_hinv      = (TH1F*) cachedSystFile->Get("bdt_syst_METScaleDown_ggZH_hinv"     );
  bdt_syst_METScaleDown_ZH_hinv[plotModel] = (TH1F*) cachedSystFile->Get(Form("bdt_syst_METScaleDown_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_syst_electronScaleUp_VVV       ->SetDirectory(0);
  bdt_syst_electronScaleUp_WZ        ->SetDirectory(0);
  bdt_syst_electronScaleUp_ZZ        ->SetDirectory(0);
  bdt_syst_electronScaleUp_ggZH_hinv ->SetDirectory(0);
  bdt_syst_electronScaleUp_ZH_hinv[plotModel] ->SetDirectory(0);
  bdt_syst_muonScaleUp_VVV           ->SetDirectory(0);
  bdt_syst_muonScaleUp_WZ            ->SetDirectory(0);
  bdt_syst_muonScaleUp_ZZ            ->SetDirectory(0);
  bdt_syst_muonScaleUp_ggZH_hinv     ->SetDirectory(0);
  bdt_syst_muonScaleUp_ZH_hinv[plotModel] ->SetDirectory(0);
  bdt_syst_METScaleUp_VVV            ->SetDirectory(0);
  bdt_syst_METScaleUp_WZ             ->SetDirectory(0);
  bdt_syst_METScaleUp_ZZ             ->SetDirectory(0);
  bdt_syst_METScaleUp_ggZH_hinv      ->SetDirectory(0);
  bdt_syst_METScaleUp_ZH_hinv[plotModel] ->SetDirectory(0);
  bdt_syst_electronScaleDown_VVV       ->SetDirectory(0);
  bdt_syst_electronScaleDown_WZ        ->SetDirectory(0);
  bdt_syst_electronScaleDown_ZZ        ->SetDirectory(0);
  bdt_syst_electronScaleDown_ggZH_hinv ->SetDirectory(0);
  bdt_syst_electronScaleDown_ZH_hinv[plotModel] ->SetDirectory(0);
  bdt_syst_muonScaleDown_VVV           ->SetDirectory(0);
  bdt_syst_muonScaleDown_WZ            ->SetDirectory(0);
  bdt_syst_muonScaleDown_ZZ            ->SetDirectory(0);
  bdt_syst_muonScaleDown_ggZH_hinv     ->SetDirectory(0);
  bdt_syst_muonScaleDown_ZH_hinv[plotModel] ->SetDirectory(0);
  bdt_syst_METScaleDown_VVV            ->SetDirectory(0);
  bdt_syst_METScaleDown_WZ             ->SetDirectory(0);
  bdt_syst_METScaleDown_ZZ             ->SetDirectory(0);
  bdt_syst_METScaleDown_ggZH_hinv      ->SetDirectory(0);
  bdt_syst_METScaleDown_ZH_hinv[plotModel] ->SetDirectory(0);
  cachedSystFile->Close();
  bdt_toy_envelope_electronScale_VVV       = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_electronScale_VVV"      );
  bdt_toy_envelope_electronScale_WZ        = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_electronScale_WZ"       );
  bdt_toy_envelope_electronScale_ZZ        = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_electronScale_ZZ"       );
  bdt_toy_envelope_electronScale_ggZH_hinv = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_electronScale_ggZH_hinv");
  bdt_toy_envelope_electronScale_ZH_hinv[plotModel] = (TH2F*) cachedSystFile->Get(Form("bdt_toy_envelope_electronScale_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_toy_envelope_muonScale_VVV           = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_muonScale_VVV"          );
  bdt_toy_envelope_muonScale_WZ            = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_muonScale_WZ"           );
  bdt_toy_envelope_muonScale_ZZ            = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_muonScale_ZZ"           );
  bdt_toy_envelope_muonScale_ggZH_hinv     = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_muonScale_ggZH_hinv"    );
  bdt_toy_envelope_muonScale_ZH_hinv[plotModel] = (TH2F*) cachedSystFile->Get(Form("bdt_toy_envelope_muonScale_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_toy_envelope_METScale_VVV            = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_METScale_VVV"           );
  bdt_toy_envelope_METScale_WZ             = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_METScale_WZ"            );
  bdt_toy_envelope_METScale_ZZ             = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_METScale_ZZ"            );
  bdt_toy_envelope_METScale_ggZH_hinv      = (TH2F*) cachedSystFile->Get("bdt_toy_envelope_METScale_ggZH_hinv"     );
  bdt_toy_envelope_METScale_ZH_hinv[plotModel] = (TH2F*) cachedSystFile->Get(Form("bdt_toy_envelope_METScale_ZH_hinv_%s",signalName_[plotModel].Data()));
  bdt_toy_envelope_electronScale_VVV       ->SetDirectory(0);
  bdt_toy_envelope_electronScale_WZ        ->SetDirectory(0);
  bdt_toy_envelope_electronScale_ZZ        ->SetDirectory(0);
  bdt_toy_envelope_electronScale_ggZH_hinv ->SetDirectory(0);
  bdt_toy_envelope_electronScale_ZH_hinv[plotModel] ->SetDirectory(0);
  bdt_toy_envelope_muonScale_VVV           ->SetDirectory(0);
  bdt_toy_envelope_muonScale_WZ            ->SetDirectory(0);
  bdt_toy_envelope_muonScale_ZZ            ->SetDirectory(0);
  bdt_toy_envelope_muonScale_ggZH_hinv     ->SetDirectory(0);
  bdt_toy_envelope_muonScale_ZH_hinv[plotModel] ->SetDirectory(0);
  bdt_toy_envelope_METScale_VVV            ->SetDirectory(0);
  bdt_toy_envelope_METScale_WZ             ->SetDirectory(0);
  bdt_toy_envelope_METScale_ZZ             ->SetDirectory(0);
  bdt_toy_envelope_METScale_ggZH_hinv      ->SetDirectory(0);
  bdt_toy_envelope_METScale_ZH_hinv[plotModel] ->SetDirectory(0);
  cachedSystFile->Close();
  unsigned long int t2 = static_cast<unsigned long int>(time(NULL));
  printf("zhAnalysis::LoadCachedBDTSystematics : Complete (%lu seconds)\n", t2-t1);
  return true;
}
bool zhAnalysis::SetupBDTSystematics() {
  printf("zhAnalysis::SetupBDTSystematics : Starting the setup\n");
  unsigned long int t1 = static_cast<unsigned long int>(time(NULL));
  cachedSystFile = new TFile(filenameBDTSysts,"recreate");
  if(!cachedSystFile || !cachedSystFile->IsOpen()) {
    printf("Error in zhAnalysis::SetupBDTSystematics : could not open file \"%s\" for writing\n", filenameBDTSysts);
    return false;
  }
  TRandom3 toy_machine(randomToySeed);
  bdt_toy_envelope_electronScale_VVV        = new TH2F("bdt_toy_envelope_electronScale_VVV"       , "bdt_toy_envelope_electronScale_VVV"       , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_electronScale_WZ         = new TH2F("bdt_toy_envelope_electronScale_WZ"        , "bdt_toy_envelope_electronScale_WZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_electronScale_ZZ         = new TH2F("bdt_toy_envelope_electronScale_ZZ"        , "bdt_toy_envelope_electronScale_ZZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_electronScale_ggZH_hinv  = new TH2F("bdt_toy_envelope_electronScale_ggZH_hinv" , "bdt_toy_envelope_electronScale_ggZH_hinv" , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_muonScale_VVV        = new TH2F("bdt_toy_envelope_muonScale_VVV"       , "bdt_toy_envelope_muonScale_VVV"       , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_muonScale_WZ         = new TH2F("bdt_toy_envelope_muonScale_WZ"        , "bdt_toy_envelope_muonScale_WZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_muonScale_ZZ         = new TH2F("bdt_toy_envelope_muonScale_ZZ"        , "bdt_toy_envelope_muonScale_ZZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_muonScale_ggZH_hinv  = new TH2F("bdt_toy_envelope_muonScale_ggZH_hinv" , "bdt_toy_envelope_muonScale_ggZH_hinv" , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_METScale_VVV        = new TH2F("bdt_toy_envelope_METScale_VVV"       , "bdt_toy_envelope_METScale_VVV"       , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_METScale_WZ         = new TH2F("bdt_toy_envelope_METScale_WZ"        , "bdt_toy_envelope_METScale_WZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_METScale_ZZ         = new TH2F("bdt_toy_envelope_METScale_ZZ"        , "bdt_toy_envelope_METScale_ZZ"        , nBinMVA, xbins, 1000 , 0, 4); 
  bdt_toy_envelope_METScale_ggZH_hinv  = new TH2F("bdt_toy_envelope_METScale_ggZH_hinv" , "bdt_toy_envelope_METScale_ggZH_hinv" , nBinMVA, xbins, 1000 , 0, 4); 
  histo_bdt_toys_electronScale_VVV       = new TH2F("histo_bdt_toys_electronScale_VVV"      , "histo_bdt_toys_electronScale_VVV"       , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_VVV->Sumw2();
  histo_bdt_toys_electronScale_WZ        = new TH2F("histo_bdt_toys_electronScale_WZ"       , "histo_bdt_toys_electronScale_WZ"        , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_WZ->Sumw2();
  histo_bdt_toys_electronScale_ZZ        = new TH2F("histo_bdt_toys_electronScale_ZZ"       , "histo_bdt_toys_electronScale_ZZ"        , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_ZZ->Sumw2();
  histo_bdt_toys_electronScale_ggZH_hinv = new TH2F("histo_bdt_toys_electronScale_ggZH_hinv", "histo_bdt_toys_electronScale_ggZH_hinv" , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_ggZH_hinv->Sumw2();
  histo_bdt_toys_muonScale_VVV           = new TH2F("histo_bdt_toys_muonScale_VVV"          , "histo_bdt_toys_muonScale_VVV"           , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_VVV->Sumw2();
  histo_bdt_toys_muonScale_WZ            = new TH2F("histo_bdt_toys_muonScale_WZ"           , "histo_bdt_toys_muonScale_WZ"            , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_WZ->Sumw2();
  histo_bdt_toys_muonScale_ZZ            = new TH2F("histo_bdt_toys_muonScale_ZZ"           , "histo_bdt_toys_muonScale_ZZ"            , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_ZZ->Sumw2();
  histo_bdt_toys_muonScale_ggZH_hinv     = new TH2F("histo_bdt_toys_muonScale_ggZH_hinv"    , "histo_bdt_toys_muonScale_ggZH_hinv"     , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_ggZH_hinv->Sumw2();
  histo_bdt_toys_METScale_VVV            = new TH2F("histo_bdt_toys_METScale_VVV"           , "histo_bdt_toys_METScale_VVV"            , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_VVV->Sumw2();
  histo_bdt_toys_METScale_WZ             = new TH2F("histo_bdt_toys_METScale_WZ"            , "histo_bdt_toys_METScale_WZ"             , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_WZ->Sumw2();
  histo_bdt_toys_METScale_ZZ             = new TH2F("histo_bdt_toys_METScale_ZZ"            , "histo_bdt_toys_METScale_ZZ"             , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_ZZ->Sumw2();
  histo_bdt_toys_METScale_ggZH_hinv      = new TH2F("histo_bdt_toys_METScale_ggZH_hinv"     , "histo_bdt_toys_METScale_ggZH_hinv"      , nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_ggZH_hinv->Sumw2();
  bdt_syst_electronScaleUp_VVV       = new TH1F("bdt_syst_electronScaleUp_VVV", "bdt_syst_electronScaleUp_VVV", nBinMVA, xbins);
  bdt_syst_electronScaleUp_WZ        = new TH1F("bdt_syst_electronScaleUp_WZ", "bdt_syst_electronScaleUp_WZ", nBinMVA, xbins);
  bdt_syst_electronScaleUp_ZZ        = new TH1F("bdt_syst_electronScaleUp_ZZ", "bdt_syst_electronScaleUp_ZZ", nBinMVA, xbins);
  bdt_syst_electronScaleUp_ggZH_hinv = new TH1F("bdt_syst_electronScaleUp_ggZH_hinv", "bdt_syst_electronScaleUp_ggZH_hinv", nBinMVA, xbins);
  bdt_syst_muonScaleUp_VVV           = new TH1F("bdt_syst_muonScaleUp_VVV", "bdt_syst_muonScaleUp_VVV", nBinMVA, xbins);
  bdt_syst_muonScaleUp_WZ            = new TH1F("bdt_syst_muonScaleUp_WZ", "bdt_syst_muonScaleUp_WZ", nBinMVA, xbins);
  bdt_syst_muonScaleUp_ZZ            = new TH1F("bdt_syst_muonScaleUp_ZZ", "bdt_syst_muonScaleUp_ZZ", nBinMVA, xbins);
  bdt_syst_muonScaleUp_ggZH_hinv     = new TH1F("bdt_syst_muonScaleUp_ggZH_hinv", "bdt_syst_muonScaleUp_ggZH_hinv", nBinMVA, xbins);
  bdt_syst_METScaleUp_VVV            = new TH1F("bdt_syst_METScaleUp_VVV", "bdt_syst_METScaleUp_VVV", nBinMVA, xbins);
  bdt_syst_METScaleUp_WZ             = new TH1F("bdt_syst_METScaleUp_WZ", "bdt_syst_METScaleUp_WZ", nBinMVA, xbins);
  bdt_syst_METScaleUp_ZZ             = new TH1F("bdt_syst_METScaleUp_ZZ", "bdt_syst_METScaleUp_ZZ", nBinMVA, xbins);
  bdt_syst_METScaleUp_ggZH_hinv      = new TH1F("bdt_syst_METScaleUp_ggZH_hinv", "bdt_syst_METScaleUp_ggZH_hinv", nBinMVA, xbins);
  bdt_syst_electronScaleDown_VVV       = new TH1F("bdt_syst_electronScaleDown_VVV", "bdt_syst_electronScaleDown_VVV", nBinMVA, xbins);
  bdt_syst_electronScaleDown_WZ        = new TH1F("bdt_syst_electronScaleDown_WZ", "bdt_syst_electronScaleDown_WZ", nBinMVA, xbins);
  bdt_syst_electronScaleDown_ZZ        = new TH1F("bdt_syst_electronScaleDown_ZZ", "bdt_syst_electronScaleDown_ZZ", nBinMVA, xbins);
  bdt_syst_electronScaleDown_ggZH_hinv = new TH1F("bdt_syst_electronScaleDown_ggZH_hinv", "bdt_syst_electronScaleDown_ggZH_hinv", nBinMVA, xbins);
  bdt_syst_muonScaleDown_VVV           = new TH1F("bdt_syst_muonScaleDown_VVV", "bdt_syst_muonScaleDown_VVV", nBinMVA, xbins);
  bdt_syst_muonScaleDown_WZ            = new TH1F("bdt_syst_muonScaleDown_WZ", "bdt_syst_muonScaleDown_WZ", nBinMVA, xbins);
  bdt_syst_muonScaleDown_ZZ            = new TH1F("bdt_syst_muonScaleDown_ZZ", "bdt_syst_muonScaleDown_ZZ", nBinMVA, xbins);
  bdt_syst_muonScaleDown_ggZH_hinv     = new TH1F("bdt_syst_muonScaleDown_ggZH_hinv", "bdt_syst_muonScaleDown_ggZH_hinv", nBinMVA, xbins);
  bdt_syst_METScaleDown_VVV            = new TH1F("bdt_syst_METScaleDown_VVV", "bdt_syst_METScaleDown_VVV", nBinMVA, xbins);
  bdt_syst_METScaleDown_WZ             = new TH1F("bdt_syst_METScaleDown_WZ", "bdt_syst_METScaleDown_WZ", nBinMVA, xbins);
  bdt_syst_METScaleDown_ZZ             = new TH1F("bdt_syst_METScaleDown_ZZ", "bdt_syst_METScaleDown_ZZ", nBinMVA, xbins);
  bdt_syst_METScaleDown_ggZH_hinv      = new TH1F("bdt_syst_METScaleDown_ggZH_hinv", "bdt_syst_METScaleDown_ggZH_hinv", nBinMVA, xbins);
  bdt_syst_electronScaleUp_VVV       ->SetDirectory(0);
  bdt_syst_electronScaleUp_WZ        ->SetDirectory(0);
  bdt_syst_electronScaleUp_ZZ        ->SetDirectory(0);
  bdt_syst_electronScaleUp_ggZH_hinv ->SetDirectory(0);
  bdt_syst_muonScaleUp_VVV           ->SetDirectory(0);
  bdt_syst_muonScaleUp_WZ            ->SetDirectory(0);
  bdt_syst_muonScaleUp_ZZ            ->SetDirectory(0);
  bdt_syst_muonScaleUp_ggZH_hinv     ->SetDirectory(0);
  bdt_syst_METScaleUp_VVV            ->SetDirectory(0);
  bdt_syst_METScaleUp_WZ             ->SetDirectory(0);
  bdt_syst_METScaleUp_ZZ             ->SetDirectory(0);
  bdt_syst_METScaleUp_ggZH_hinv      ->SetDirectory(0);
  bdt_syst_electronScaleDown_VVV       ->SetDirectory(0);
  bdt_syst_electronScaleDown_WZ        ->SetDirectory(0);
  bdt_syst_electronScaleDown_ZZ        ->SetDirectory(0);
  bdt_syst_electronScaleDown_ggZH_hinv ->SetDirectory(0);
  bdt_syst_muonScaleDown_VVV           ->SetDirectory(0);
  bdt_syst_muonScaleDown_WZ            ->SetDirectory(0);
  bdt_syst_muonScaleDown_ZZ            ->SetDirectory(0);
  bdt_syst_muonScaleDown_ggZH_hinv     ->SetDirectory(0);
  bdt_syst_METScaleDown_VVV            ->SetDirectory(0);
  bdt_syst_METScaleDown_WZ             ->SetDirectory(0);
  bdt_syst_METScaleDown_ZZ             ->SetDirectory(0);
  bdt_syst_METScaleDown_ggZH_hinv      ->SetDirectory(0);
  for(int nModel=0; nModel<nSigModels; nModel++) { 
    bdt_toy_envelope_electronScale_ZH_hinv[nModel] = new TH2F(Form("bdt_toy_envelope_electronScale_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_toy_envelope_electronScale_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins, 1000 , 0, 4); 
    bdt_toy_envelope_muonScale_ZH_hinv[nModel]     = new TH2F(Form("bdt_toy_envelope_muonScale_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_toy_envelope_muonScale_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins, 1000 , 0, 4); 
    bdt_toy_envelope_METScale_ZH_hinv[nModel]      = new TH2F(Form("bdt_toy_envelope_METScale_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_toy_envelope_METScale_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins, 1000 , 0, 4); 
    histo_bdt_toys_electronScale_ZH_hinv[nModel]   = new TH2F(Form("histo_bdt_toys_electronScale_ZH_hinv_%s",signalName_[nModel].Data()),   Form("histo_bdt_toys_electronScale_ZH_hinv_%s",signalName_[nModel].Data()),   nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_electronScale_ZH_hinv[nModel]->Sumw2();
    histo_bdt_toys_muonScale_ZH_hinv[nModel]   = new TH2F(Form("histo_bdt_toys_muonScale_ZH_hinv_%s",signalName_[nModel].Data()),   Form("histo_bdt_toys_muonScale_ZH_hinv_%s",signalName_[nModel].Data()),   nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_muonScale_ZH_hinv[nModel]->Sumw2();
    histo_bdt_toys_METScale_ZH_hinv[nModel]   = new TH2F(Form("histo_bdt_toys_METScale_ZH_hinv_%s",signalName_[nModel].Data()),   Form("histo_bdt_toys_METScale_ZH_hinv_%s",signalName_[nModel].Data()),   nBinMVA, xbins, num_bdt_toys, 0, num_bdt_toys); histo_bdt_toys_METScale_ZH_hinv[nModel]->Sumw2();
    bdt_syst_electronScaleUp_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_electronScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_electronScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
    bdt_syst_muonScaleUp_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_muonScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_muonScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
    bdt_syst_METScaleUp_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_METScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_METScaleUp_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
    bdt_syst_electronScaleDown_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_electronScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_electronScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
    bdt_syst_muonScaleDown_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_muonScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_muonScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
    bdt_syst_METScaleDown_ZH_hinv[nModel] = new TH1F(Form("bdt_syst_METScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), Form("bdt_syst_METScaleDown_ZH_hinv_%s",signalName_[nModel].Data()), nBinMVA, xbins);
    bdt_syst_electronScaleUp_ZH_hinv[nModel] ->SetDirectory(0);
    bdt_syst_muonScaleUp_ZH_hinv[nModel] ->SetDirectory(0);
    bdt_syst_METScaleUp_ZH_hinv[nModel] ->SetDirectory(0);
    bdt_syst_electronScaleDown_ZH_hinv[nModel] ->SetDirectory(0);
    bdt_syst_muonScaleDown_ZH_hinv[nModel] ->SetDirectory(0);
    bdt_syst_METScaleDown_ZH_hinv[nModel] ->SetDirectory(0);
  }
  for(int nb=1; nb<nBinMVA; nb++) {  
    bdt_syst_electronScaleUp_VVV       ->SetBinContent(nb,1.);
    bdt_syst_electronScaleUp_WZ        ->SetBinContent(nb,1.);
    bdt_syst_electronScaleUp_ZZ        ->SetBinContent(nb,1.);
    bdt_syst_electronScaleUp_ggZH_hinv ->SetBinContent(nb,1.);
    bdt_syst_muonScaleUp_VVV           ->SetBinContent(nb,1.);
    bdt_syst_muonScaleUp_WZ            ->SetBinContent(nb,1.);
    bdt_syst_muonScaleUp_ZZ            ->SetBinContent(nb,1.);
    bdt_syst_muonScaleUp_ggZH_hinv     ->SetBinContent(nb,1.);
    bdt_syst_METScaleUp_VVV            ->SetBinContent(nb,1.);
    bdt_syst_METScaleUp_WZ             ->SetBinContent(nb,1.);
    bdt_syst_METScaleUp_ZZ             ->SetBinContent(nb,1.);
    bdt_syst_METScaleUp_ggZH_hinv      ->SetBinContent(nb,1.);
    bdt_syst_electronScaleDown_VVV       ->SetBinContent(nb,1.);
    bdt_syst_electronScaleDown_WZ        ->SetBinContent(nb,1.);
    bdt_syst_electronScaleDown_ZZ        ->SetBinContent(nb,1.);
    bdt_syst_electronScaleDown_ggZH_hinv ->SetBinContent(nb,1.);
    bdt_syst_muonScaleDown_VVV           ->SetBinContent(nb,1.);
    bdt_syst_muonScaleDown_WZ            ->SetBinContent(nb,1.);
    bdt_syst_muonScaleDown_ZZ            ->SetBinContent(nb,1.);
    bdt_syst_muonScaleDown_ggZH_hinv     ->SetBinContent(nb,1.);
    bdt_syst_METScaleDown_VVV            ->SetBinContent(nb,1.);
    bdt_syst_METScaleDown_WZ             ->SetBinContent(nb,1.);
    bdt_syst_METScaleDown_ZZ             ->SetBinContent(nb,1.);
    bdt_syst_METScaleDown_ggZH_hinv      ->SetBinContent(nb,1.);
    for(int nModel=0; nModel<nSigModels; nModel++) {
      bdt_syst_electronScaleUp_ZH_hinv[nModel] ->SetBinContent(nb,1.);
      bdt_syst_muonScaleUp_ZH_hinv[nModel] ->SetBinContent(nb,1.);
      bdt_syst_METScaleUp_ZH_hinv[nModel] ->SetBinContent(nb,1.);
      bdt_syst_electronScaleDown_ZH_hinv[nModel] ->SetBinContent(nb,1.);
      bdt_syst_muonScaleDown_ZH_hinv[nModel] ->SetBinContent(nb,1.);
      bdt_syst_METScaleDown_ZH_hinv[nModel] ->SetBinContent(nb,1.);
    }
  }
  // Store a value from Gaussian distribution with mean=0, sigma=1
  // These values will get multiplied by the individual nuisance sizes
  for(unsigned int i_toy=0; i_toy<num_bdt_toys; i_toy++) {
    bdt_toy_scale[i_toy]=toy_machine.Gaus(0,1);
  }
  
  unsigned long int t2 = static_cast<unsigned long int>(time(NULL));
  printf("zhAnalysis::SetupBDTSystematics : Complete (%lu seconds)\n", t2-t1);
  return true;

}
bool zhAnalysis::SaveBDTSystematics(int nModel) {
  printf("zhAnalysis::SaveBDTSystematics : Saving the BDT systematics to file\n");
  unsigned long int t1 = static_cast<unsigned long int>(time(NULL));
  cachedSystFile->cd();
  histo_bdt_toys_electronScale_VVV       ->Write();
  histo_bdt_toys_electronScale_WZ        ->Write();
  histo_bdt_toys_electronScale_ZZ        ->Write();
  histo_bdt_toys_electronScale_ggZH_hinv ->Write();
  histo_bdt_toys_electronScale_ZH_hinv[nModel]->Write();
  histo_bdt_toys_muonScale_VVV           ->Write();
  histo_bdt_toys_muonScale_WZ            ->Write();
  histo_bdt_toys_muonScale_ZZ            ->Write();
  histo_bdt_toys_muonScale_ggZH_hinv     ->Write();
  histo_bdt_toys_muonScale_ZH_hinv[nModel]    ->Write();
  histo_bdt_toys_METScale_VVV            ->Write();
  histo_bdt_toys_METScale_WZ             ->Write();
  histo_bdt_toys_METScale_ZZ             ->Write();
  histo_bdt_toys_METScale_ggZH_hinv      ->Write();
  histo_bdt_toys_METScale_ZH_hinv[nModel]     ->Write();
  bdt_toy_envelope_electronScale_VVV->Write();
  bdt_toy_envelope_electronScale_WZ->Write();
  bdt_toy_envelope_electronScale_ZZ->Write();
  bdt_toy_envelope_electronScale_ZH_hinv[nModel]->Write();
  bdt_toy_envelope_electronScale_ggZH_hinv->Write();
  bdt_toy_envelope_muonScale_VVV->Write();
  bdt_toy_envelope_muonScale_WZ->Write();
  bdt_toy_envelope_muonScale_ZZ->Write();
  bdt_toy_envelope_muonScale_ZH_hinv[nModel]->Write();
  bdt_toy_envelope_muonScale_ggZH_hinv->Write();
  bdt_toy_envelope_METScale_VVV->Write();
  bdt_toy_envelope_METScale_WZ->Write();
  bdt_toy_envelope_METScale_ZZ->Write();
  bdt_toy_envelope_METScale_ZH_hinv[nModel]->Write();
  bdt_toy_envelope_METScale_ggZH_hinv->Write();
  bdt_syst_electronScaleUp_VVV       ->Write();
  bdt_syst_electronScaleUp_WZ        ->Write();
  bdt_syst_electronScaleUp_ZZ        ->Write();
  bdt_syst_electronScaleUp_ggZH_hinv ->Write();
  bdt_syst_electronScaleUp_ZH_hinv[nModel] ->Write();
  bdt_syst_muonScaleUp_VVV           ->Write();
  bdt_syst_muonScaleUp_WZ            ->Write();
  bdt_syst_muonScaleUp_ZZ            ->Write();
  bdt_syst_muonScaleUp_ggZH_hinv     ->Write();
  bdt_syst_muonScaleUp_ZH_hinv[nModel] ->Write();
  bdt_syst_METScaleUp_VVV            ->Write();
  bdt_syst_METScaleUp_WZ             ->Write();
  bdt_syst_METScaleUp_ZZ             ->Write();
  bdt_syst_METScaleUp_ggZH_hinv      ->Write();
  bdt_syst_METScaleUp_ZH_hinv[nModel] ->Write();
  bdt_syst_electronScaleDown_VVV       ->Write();
  bdt_syst_electronScaleDown_WZ        ->Write();
  bdt_syst_electronScaleDown_ZZ        ->Write();
  bdt_syst_electronScaleDown_ggZH_hinv ->Write();
  bdt_syst_electronScaleDown_ZH_hinv[nModel] ->Write();
  bdt_syst_muonScaleDown_VVV           ->Write();
  bdt_syst_muonScaleDown_WZ            ->Write();
  bdt_syst_muonScaleDown_ZZ            ->Write();
  bdt_syst_muonScaleDown_ggZH_hinv     ->Write();
  bdt_syst_muonScaleDown_ZH_hinv[nModel] ->Write();
  bdt_syst_METScaleDown_VVV            ->Write();
  bdt_syst_METScaleDown_WZ             ->Write();
  bdt_syst_METScaleDown_ZZ             ->Write();
  bdt_syst_METScaleDown_ggZH_hinv      ->Write();
  bdt_syst_METScaleDown_ZH_hinv[nModel] ->Write();
  cachedSystFile->Close();
  
  unsigned long int t2 = static_cast<unsigned long int>(time(NULL));
  printf("zhAnalysis::SaveBDTSystematics : Complete (%lu seconds)\n", t2-t1);
  return true;
}
void zhAnalysis::SetBranchAddresses(TTree *theInputTree, bool isData) {
  theInputTree->SetBranchAddress("runNumber", &gltEvent.runNumber);
  theInputTree->SetBranchAddress("lumiNumber", &gltEvent.lumiNumber);
  theInputTree->SetBranchAddress("eventNumber", &gltEvent.eventNumber);
  theInputTree->SetBranchAddress("npv", &gltEvent.npv);
  theInputTree->SetBranchAddress("pu", &gltEvent.pu);
  theInputTree->SetBranchAddress("metFilter", &gltEvent.metFilter);
  theInputTree->SetBranchAddress("egmFilter", &gltEvent.egmFilter);
  theInputTree->SetBranchAddress("nLooseLep", &gltEvent.nLooseLep);
  theInputTree->SetBranchAddress("looseLep1PdgId", &gltEvent.looseLep1PdgId);
  theInputTree->SetBranchAddress("looseLep2PdgId", &gltEvent.looseLep2PdgId);
  theInputTree->SetBranchAddress("looseLep3PdgId", &gltEvent.looseLep3PdgId);
  theInputTree->SetBranchAddress("looseLep4PdgId", &gltEvent.looseLep4PdgId);
  theInputTree->SetBranchAddress("looseLep1SelBit", &gltEvent.looseLep1SelBit);
  theInputTree->SetBranchAddress("looseLep2SelBit", &gltEvent.looseLep2SelBit);
  theInputTree->SetBranchAddress("looseLep3SelBit", &gltEvent.looseLep3SelBit);
  theInputTree->SetBranchAddress("looseLep4SelBit", &gltEvent.looseLep4SelBit);
  theInputTree->SetBranchAddress("looseLep1Pt", &gltEvent.looseLep1Pt);
  theInputTree->SetBranchAddress("looseLep2Pt", &gltEvent.looseLep2Pt);
  theInputTree->SetBranchAddress("looseLep3Pt", &gltEvent.looseLep3Pt);
  theInputTree->SetBranchAddress("looseLep4Pt", &gltEvent.looseLep4Pt);
  theInputTree->SetBranchAddress("looseLep1Eta", &gltEvent.looseLep1Eta);
  theInputTree->SetBranchAddress("looseLep2Eta", &gltEvent.looseLep2Eta);
  theInputTree->SetBranchAddress("looseLep3Eta", &gltEvent.looseLep3Eta);
  theInputTree->SetBranchAddress("looseLep4Eta", &gltEvent.looseLep4Eta);
  theInputTree->SetBranchAddress("looseLep1Phi", &gltEvent.looseLep1Phi);
  theInputTree->SetBranchAddress("looseLep2Phi", &gltEvent.looseLep2Phi);
  theInputTree->SetBranchAddress("looseLep3Phi", &gltEvent.looseLep3Phi);
  theInputTree->SetBranchAddress("looseLep4Phi", &gltEvent.looseLep4Phi);
  theInputTree->SetBranchAddress("nJet", &gltEvent.nJet);
  theInputTree->SetBranchAddress("jetNLBtags", &gltEvent.jetNLBtags);
  theInputTree->SetBranchAddress("jetNMBtags", &gltEvent.jetNMBtags);
  theInputTree->SetBranchAddress("jetNTBtags", &gltEvent.jetNTBtags);
  theInputTree->SetBranchAddress("jet1Pt", &gltEvent.jet1Pt);
  theInputTree->SetBranchAddress("jet2Pt", &gltEvent.jet2Pt);
  theInputTree->SetBranchAddress("jet3Pt", &gltEvent.jet3Pt);
  theInputTree->SetBranchAddress("jet4Pt", &gltEvent.jet4Pt);
  theInputTree->SetBranchAddress("jet1Eta", &gltEvent.jet1Eta);
  theInputTree->SetBranchAddress("jet2Eta", &gltEvent.jet2Eta);
  theInputTree->SetBranchAddress("jet3Eta", &gltEvent.jet3Eta);
  theInputTree->SetBranchAddress("jet4Eta", &gltEvent.jet4Eta);
  theInputTree->SetBranchAddress("jet1Phi", &gltEvent.jet1Phi);
  theInputTree->SetBranchAddress("jet2Phi", &gltEvent.jet2Phi);
  theInputTree->SetBranchAddress("jet3Phi", &gltEvent.jet3Phi);
  theInputTree->SetBranchAddress("jet4Phi", &gltEvent.jet4Phi);
  theInputTree->SetBranchAddress("jet1BTag", &gltEvent.jet1BTag);
  theInputTree->SetBranchAddress("jet2BTag", &gltEvent.jet2BTag);
  theInputTree->SetBranchAddress("jet3BTag", &gltEvent.jet3BTag);
  theInputTree->SetBranchAddress("jet4BTag", &gltEvent.jet4BTag);
  theInputTree->SetBranchAddress("jet1SelBit", &gltEvent.jet1SelBit);
  theInputTree->SetBranchAddress("jet2SelBit", &gltEvent.jet2SelBit);
  theInputTree->SetBranchAddress("jet3SelBit", &gltEvent.jet3SelBit);
  theInputTree->SetBranchAddress("jet4SelBit", &gltEvent.jet4SelBit);
  theInputTree->SetBranchAddress("jet1PtUp", &gltEvent.jet1PtUp);
  theInputTree->SetBranchAddress("jet2PtUp", &gltEvent.jet2PtUp);
  theInputTree->SetBranchAddress("jet3PtUp", &gltEvent.jet3PtUp);
  theInputTree->SetBranchAddress("jet4PtUp", &gltEvent.jet4PtUp);
  theInputTree->SetBranchAddress("jet1PtDown", &gltEvent.jet1PtDown);
  theInputTree->SetBranchAddress("jet2PtDown", &gltEvent.jet2PtDown);
  theInputTree->SetBranchAddress("jet3PtDown", &gltEvent.jet3PtDown);
  theInputTree->SetBranchAddress("jet4PtDown", &gltEvent.jet4PtDown);
  theInputTree->SetBranchAddress("jet1EtaUp", &gltEvent.jet1EtaUp);
  theInputTree->SetBranchAddress("jet2EtaUp", &gltEvent.jet2EtaUp);
  theInputTree->SetBranchAddress("jet3EtaUp", &gltEvent.jet3EtaUp);
  theInputTree->SetBranchAddress("jet4EtaUp", &gltEvent.jet4EtaUp);
  theInputTree->SetBranchAddress("jet1EtaDown", &gltEvent.jet1EtaDown);
  theInputTree->SetBranchAddress("jet2EtaDown", &gltEvent.jet2EtaDown);
  theInputTree->SetBranchAddress("jet3EtaDown", &gltEvent.jet3EtaDown);
  theInputTree->SetBranchAddress("jet4EtaDown", &gltEvent.jet4EtaDown);
  theInputTree->SetBranchAddress("pfmet", &gltEvent.pfmet);
  theInputTree->SetBranchAddress("pfmetphi", &gltEvent.pfmetphi);
  theInputTree->SetBranchAddress("pfmetRaw", &gltEvent.pfmetRaw);
  theInputTree->SetBranchAddress("pfmetUp", &gltEvent.pfmetUp);
  theInputTree->SetBranchAddress("pfmetDown", &gltEvent.pfmetDown);
  theInputTree->SetBranchAddress("pfmetnomu", &gltEvent.pfmetnomu);
  theInputTree->SetBranchAddress("puppimet", &gltEvent.puppimet);
  theInputTree->SetBranchAddress("puppimetphi", &gltEvent.puppimetphi);
  theInputTree->SetBranchAddress("calomet", &gltEvent.calomet);
  theInputTree->SetBranchAddress("calometphi", &gltEvent.calometphi);
  theInputTree->SetBranchAddress("trkmet", &gltEvent.trkmet);
  theInputTree->SetBranchAddress("trkmetphi", &gltEvent.trkmetphi);
  theInputTree->SetBranchAddress("dphipfmet", &gltEvent.dphipfmet);
  theInputTree->SetBranchAddress("dphipuppimet", &gltEvent.dphipuppimet);
  theInputTree->SetBranchAddress("nTau", &gltEvent.nTau);
  theInputTree->SetBranchAddress("nLoosePhoton", &gltEvent.nLoosePhoton);
  theInputTree->SetBranchAddress("loosePho1Pt", &gltEvent.loosePho1Pt);
  theInputTree->SetBranchAddress("loosePho1Eta", &gltEvent.loosePho1Eta);
  theInputTree->SetBranchAddress("loosePho1Phi", &gltEvent.loosePho1Phi);
  if(isData) theInputTree->SetBranchAddress("trigger", &gltEvent.trigger);
  else {
    theInputTree->SetBranchAddress("looseGenLep1PdgId", &gltEvent.looseGenLep1PdgId);
    theInputTree->SetBranchAddress("looseGenLep2PdgId", &gltEvent.looseGenLep2PdgId);
    theInputTree->SetBranchAddress("looseGenLep3PdgId", &gltEvent.looseGenLep3PdgId);
    theInputTree->SetBranchAddress("looseGenLep4PdgId", &gltEvent.looseGenLep4PdgId);
    theInputTree->SetBranchAddress("jet1GenPt", &gltEvent.jet1GenPt);
    theInputTree->SetBranchAddress("jet2GenPt", &gltEvent.jet2GenPt);
    theInputTree->SetBranchAddress("jet3GenPt", &gltEvent.jet3GenPt);
    theInputTree->SetBranchAddress("jet4GenPt", &gltEvent.jet4GenPt);
    theInputTree->SetBranchAddress("jet1Flav", &gltEvent.jet1Flav);
    theInputTree->SetBranchAddress("jet2Flav", &gltEvent.jet2Flav);
    theInputTree->SetBranchAddress("jet3Flav", &gltEvent.jet3Flav);
    theInputTree->SetBranchAddress("jet4Flav", &gltEvent.jet4Flav);
    theInputTree->SetBranchAddress("sf_pu"            , &gltEvent.sf_pu);
    theInputTree->SetBranchAddress("sf_puUp"          , &gltEvent.sf_puUp);
    theInputTree->SetBranchAddress("sf_puDown"        , &gltEvent.sf_puDown);
    theInputTree->SetBranchAddress("sf_zz"            , &gltEvent.sf_zz);
    theInputTree->SetBranchAddress("sf_zzUnc"         , &gltEvent.sf_zzUnc);
    theInputTree->SetBranchAddress("sf_wz"            , &gltEvent.sf_wz);
    theInputTree->SetBranchAddress("sf_zh"            , &gltEvent.sf_zh);
    theInputTree->SetBranchAddress("sf_zhUp"          , &gltEvent.sf_zhUp);
    theInputTree->SetBranchAddress("sf_zhDown"        , &gltEvent.sf_zhDown);
    theInputTree->SetBranchAddress("sf_tt"            , &gltEvent.sf_tt);
    theInputTree->SetBranchAddress("sf_trk1"          , &gltEvent.sf_trk1);
    theInputTree->SetBranchAddress("sf_trk2"          , &gltEvent.sf_trk2);
    theInputTree->SetBranchAddress("sf_trk3"          , &gltEvent.sf_trk3);
    theInputTree->SetBranchAddress("sf_trk4"          , &gltEvent.sf_trk4);
    theInputTree->SetBranchAddress("sf_loose1"        , &gltEvent.sf_loose1);
    theInputTree->SetBranchAddress("sf_loose2"        , &gltEvent.sf_loose2);
    theInputTree->SetBranchAddress("sf_loose3"        , &gltEvent.sf_loose3);
    theInputTree->SetBranchAddress("sf_loose4"        , &gltEvent.sf_loose4);
    theInputTree->SetBranchAddress("sf_medium1"       , &gltEvent.sf_medium1);
    theInputTree->SetBranchAddress("sf_medium2"       , &gltEvent.sf_medium2);
    theInputTree->SetBranchAddress("sf_medium3"       , &gltEvent.sf_medium3);
    theInputTree->SetBranchAddress("sf_medium4"       , &gltEvent.sf_medium4);
    theInputTree->SetBranchAddress("sf_tight1"        , &gltEvent.sf_tight1);
    theInputTree->SetBranchAddress("sf_tight2"        , &gltEvent.sf_tight2);
    theInputTree->SetBranchAddress("sf_tight3"        , &gltEvent.sf_tight3);
    theInputTree->SetBranchAddress("sf_tight4"        , &gltEvent.sf_tight4);
    theInputTree->SetBranchAddress("sf_unc1"          , &gltEvent.sf_unc1);
    theInputTree->SetBranchAddress("sf_unc2"          , &gltEvent.sf_unc2);
    theInputTree->SetBranchAddress("sf_unc3"          , &gltEvent.sf_unc3);
    theInputTree->SetBranchAddress("sf_unc4"          , &gltEvent.sf_unc4);
    theInputTree->SetBranchAddress("sf_btag0"         , &sf_btag0        );
    theInputTree->SetBranchAddress("sf_btag0BUp"      , &sf_btag0BUp     );
    theInputTree->SetBranchAddress("sf_btag0BDown"    , &sf_btag0BDown   );
    theInputTree->SetBranchAddress("sf_btag0MUp"      , &sf_btag0MUp     );
    theInputTree->SetBranchAddress("sf_btag0MDown"    , &sf_btag0MDown   );
    theInputTree->SetBranchAddress("normalizedWeight" , &normalizedWeight);
  }

}
