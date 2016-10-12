#include <TROOT.h>
#include <TLorentzVector.h>
#include <TMath.h>
//#include "TMVA/Reader.h"

// mvaNuisances:
//   Access the MVA variables
//   (optionally) vary relevant nuisances
void mvaNuisances(
  TLorentzVector lepton1,
  TLorentzVector lepton2,
  TLorentzVector MET,
  TLorentzVector jet1,
  Float_t  &mva_balance,
  Float_t  &mva_cos_theta_star_l1,
  Float_t  &mva_cos_theta_CS_l1,
  Float_t  &mva_delphi_ptll_MET,
  Float_t  &mva_delphi_ll,
  Float_t  &mva_delphi_jet_MET,
  Float_t  &mva_deltaR_ll,
  Float_t  &mva_etall,
  Float_t  &mva_etal1,
  Float_t  &mva_etal2,
  Float_t  &mva_MET,
  Float_t  &mva_mll_minus_mZ,
  Float_t  &mva_mTjetMET,
  Float_t  &mva_mTll,
  Float_t  &mva_mTl1MET,
  Float_t  &mva_mTl2MET,
  Float_t  &mva_ptll,
  Float_t  &mva_ptl1,
  Float_t  &mva_ptl2,
  Float_t  &mva_ptl1mptl2_over_ptll,
  Float_t  lepton1_scale_variation = 0, // Fractional lepton scale variation
  Float_t  lepton2_scale_variation = 0, // Fractional lepton scale variation
  Float_t  MET_scale_variation    = 0, // MET scale variation
  Float_t  jet_scale_variation    = 0  // jet scale variation
) {
  bool eventHasJet = jet1.E()>0;
  TLorentzVector d_lepton1(0,0,0,0), d_lepton2(0,0,0,0), d_jet1(0,0,0,0), d_MET; // Variations of the reconstructed objects
  TLorentzVector lepton1_new, lepton2_new, MET_new, jet1_new;
  if(lepton1_scale_variation!=0) {
    lepton1_new.SetPtEtaPhiE(
      lepton1.Pt() * (1.+lepton1_scale_variation),
      lepton1.Eta(),
      lepton1.Phi(),
      lepton1.E()  * (1.+lepton1_scale_variation)
    );
    d_lepton1 = lepton1_new - lepton1;
    lepton1 = lepton1_new;
  }
  if(lepton2_scale_variation!=0) {
    lepton2_new.SetPtEtaPhiE(
      lepton2.Pt() * (1.+lepton2_scale_variation),
      lepton2.Eta(),
      lepton2.Phi(),
      lepton2.E()  * (1.+lepton2_scale_variation)
    );
    d_lepton2 = lepton2_new - lepton2;
    lepton2 = lepton2_new;
  }
  if(jet_scale_variation!=0 && eventHasJet) {
    jet1_new.SetPtEtaPhiE(
      jet1.Pt() * (1.+jet_scale_variation),
      jet1.Eta(),
      jet1.Phi(),
      jet1.E()  * (1.+jet_scale_variation)
    );
    d_jet1 = jet1_new - jet1;
    jet1 = jet1_new;
  }
  if(MET_scale_variation!=0) {
    MET_new.SetPtEtaPhiE(
      MET.Pt() * (1.+MET_scale_variation),
      MET.Eta(),
      MET.Phi(),
      MET.E()  * (1.+MET_scale_variation)
    );
    d_MET = MET_new - MET; 
  }
  MET = MET + d_MET - d_lepton1 - d_lepton2 - d_jet1;
  TLorentzVector dilep = lepton1+lepton2;
  
  mva_balance             = TMath::Abs(dilep.Pt() - MET.Pt()) / dilep.Pt();
  mva_delphi_ptll_MET     = TMath::Abs(dilep.DeltaPhi(MET)); 
  mva_cos_theta_star_l1   = cos_theta_star( lepton1, lepton2, dilep+MET);
  mva_cos_theta_CS_l1     = cos_theta_collins_soper( lepton1, lepton2 );
  mva_deltaR_ll           = lepton1.DeltaR(lepton2);
  mva_delphi_ll           = TMath::Abs(lepton1.DeltaPhi(lepton2));
  mva_delphi_jet_MET      = eventHasJet ? TMath::Abs(jet1.DeltaPhi(MET)) : -1.;
  mva_etal1               = lepton1.Eta();
  mva_etal2               = lepton2.Eta();
  mva_etall               = dilep.Eta(); 
  mva_MET                 = MET.Pt();
  mva_mll_minus_mZ        = TMath::Abs(dilep.M() - 91.1876); 
  mva_mTjetMET            = eventHasJet ? TMath::Sqrt(2.0*jet1.Pt()*MET.Pt()*(1.0 - cos(TMath::Abs(jet1.DeltaPhi(MET))))) : -1.;
  mva_mTll                = TMath::Sqrt(2.0*dilep.Pt()*MET.Pt()*(1.0 - cos(TMath::Abs(dilep.DeltaPhi(MET)))));
  mva_mTl1MET             = TMath::Sqrt(2.0*lepton1.Pt()*MET.Pt()*(1.0 - cos(TMath::Abs(lepton1.DeltaPhi(MET))))); 
  mva_mTl2MET             = TMath::Sqrt(2.0*lepton2.Pt()*MET.Pt()*(1.0 - cos(TMath::Abs(lepton2.DeltaPhi(MET))))); 
  mva_ptll                = dilep.Pt(); 
  mva_ptl1                = lepton1.Pt();
  mva_ptl2                = lepton2.Pt();
  mva_ptl1mptl2_over_ptll = TMath::Abs(lepton1.Pt() - lepton2.Pt()) / dilep.Pt();

  if(mva_cos_theta_star_l1 != mva_cos_theta_star_l1) {
    printf("PROBLEM with mvaNuisances: mva_cos_theta_star_l1 = NaN (are we in an indian restaurant?)\n");
    printf("(pt, eta, phi, E) lepton1 (%f, %f, %f, %f) lepton2 (%f, %f, %f, %f)\n", lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), lepton1.E(), lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), lepton2.E());
    printf("(pt, eta, phi, E) d_lepton1 (%f, %f, %f, %f) d_lepton2 (%f, %f, %f, %f)\n", d_lepton1.Pt(), d_lepton1.Eta(), d_lepton1.Phi(), d_lepton1.E(), d_lepton2.Pt(), d_lepton2.Eta(), d_lepton2.Phi(), d_lepton2.E());
  }
  if(mva_cos_theta_CS_l1 != mva_cos_theta_CS_l1) {
    printf("PROBLEM with mvaNuisances: mva_cos_theta_star_l1 = NaN (are we in an indian restaurant?)\n");
    printf("(pt, eta, phi, E) lepton1 (%f, %f, %f, %f) lepton2 (%f, %f, %f, %f)\n", lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), lepton1.E(), lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), lepton2.E());
    printf("(pt, eta, phi, E) d_lepton1 (%f, %f, %f, %f) d_lepton2 (%f, %f, %f, %f)\n", d_lepton1.Pt(), d_lepton1.Eta(), d_lepton1.Phi(), d_lepton1.E(), d_lepton2.Pt(), d_lepton2.Eta(), d_lepton2.Phi(), d_lepton2.E());
  }
}
