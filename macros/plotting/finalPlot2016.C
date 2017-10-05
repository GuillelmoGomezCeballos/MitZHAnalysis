#include "TROOT.h"
#include "Math/ProbFuncMathCore.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "StandardPlot2016.C"
//#include "GoodStyle.C"
#include "TSystem.h"
#include "CMS_lumi.C"
#include "TRandom.h"

double scaling[8] = {1,1,1,1,1,1,1,1};
//double scaling[8] = {73.496/69.913,1.2,168.114/157.695,404.402/361.347,2.583/2.537,1,1,1}; // EM/Zjets/WZ/ZZ/VVV
//double scaling[8] = {0.94,1.04,1.04,1.01,1.00,1.02,1.01,1.05}; // EWK/QCD/VV/VVV/WS/WG/DPS/Wjets

void eraselabel(TPad *p,Double_t h){
  p->cd();
  TPad *pe = new TPad("pe","pe",0.02,0,p->GetLeftMargin(),h);	   
  pe->Draw();
  pe->SetFillColor(p->GetFillColor()); 
  pe->SetBorderMode(0);
}

void atributes(TH1D *histo, TString xtitle = "", TString ytitle = "Fraction", TString units = ""){

  histo->SetTitle("");
  //histo->SetMarkerStyle(20);
  //histo->SetMarkerSize(0.8);
  //histo->SetLineWidth(4);
  if(strcmp(units.Data(),"")==0){
    histo->GetXaxis()->SetTitle(xtitle.Data());
  } else {
    histo->GetXaxis()->SetTitle(Form("%s [%s]",xtitle.Data(),units.Data()));
  }
  histo->GetXaxis()->SetLabelFont  (   42);
  histo->GetXaxis()->SetLabelOffset(0.015);
  histo->GetXaxis()->SetLabelSize  (0.140);
  histo->GetXaxis()->SetNdivisions (  505);
  histo->GetXaxis()->SetTitleFont  (   42);
  histo->GetXaxis()->SetTitleOffset( 0.95);
  histo->GetXaxis()->SetTitleSize  (0.140);
  histo->GetXaxis()->SetTickLength (0.07 );

  histo->GetYaxis()->SetTitle(ytitle.Data());
  histo->GetYaxis()->SetLabelFont  (   42);
  histo->GetYaxis()->SetLabelOffset(0.015);
  histo->GetYaxis()->SetLabelSize  (0.120);
  histo->GetYaxis()->SetNdivisions (  505);
  histo->GetYaxis()->SetTitleFont  (   42);
  histo->GetYaxis()->SetTitleOffset(  0.5);
  histo->GetYaxis()->SetTitleSize  (0.120);
  histo->GetYaxis()->SetTickLength (0.03 );

  histo->SetLineColor  (kBlack);
  histo->SetMarkerStyle(kFullCircle);
}

void finalPlot2016(int nsel = 0, int ReBin = 1, TString XTitle = "N_{jets}", TString units = "", TString plotName = "histo_nice.root", TString outputName = "njets",
                bool isLogY = false, TString higgsLabel = "", double lumi = 1.0, bool isBlind = false, TString extraLabel = "",
		TString plotExtraName = "histo_nice2.root", TString higgs2Label = "", bool show2D = false) {

  gInterpreter->ExecuteMacro("GoodStyle.C");
  //gROOT->LoadMacro("StandardPlot2016.C");
  gStyle->SetOptStat(0);

  TFile* file = new TFile(plotName, "read");  if(!file) {printf("File %s does not exist\n",plotName.Data()); return;}
  TFile* fileExtra = new TFile(plotExtraName, "read"); bool hasFileExtra = true; if(!fileExtra->IsOpen()) hasFileExtra = false;

  StandardPlot2016 myPlot;
  myPlot.setLumi(lumi);
  myPlot.setLabel(XTitle);
  myPlot.addLabel(extraLabel.Data());
  myPlot.setHiggsLabel(higgsLabel.Data());
  myPlot.setHiggs2Label(higgs2Label.Data());
  myPlot.setUnits(units);

  TH1F* hOther;
  TH1F* hWWEWK;
  TH1F* hWWQCD;
  TH1F* hDPS  ;
  TH1F* hWS   ;
  TH1F* hWW   ;
  TH1F* hZJets;
  TH1F* hTop  ;
  TH1F* hVV   ;
  TH1F* hWJets;
  TH1F* hWG   ;
  TH1F* hHiggs;
  TH1F* hggZH ;
  TH1F* hData ;
  TH1F* hWZ   ;
  TH1F* hZZ   ;
  TH1F* hVVV  ;
  TH1F* hEM   ;
  TH1F* hBck  ;
  TH1F* hHiggs2;

  if     (nsel == 0){
    printf("ZLL style plots\n");
    hWW     = (TH1F*)file->Get("histo1");
    hZJets  = (TH1F*)file->Get("histo2");
    hTop    = (TH1F*)file->Get("histo3");
    hVV     = (TH1F*)file->Get("histo4");
    hWJets  = (TH1F*)file->Get("histo5");
    hWG     = (TH1F*)file->Get("histo6");
    hHiggs  = (TH1F*)file->Get("histo7");
    hData   = (TH1F*)file->Get("histo0");
    hZJets->Scale(lumi);

    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=1; i<=6; i++){
        hWW    ->SetBinContent(i,0);
        hZJets ->SetBinContent(i,0);
        hTop   ->SetBinContent(i,0);
        hVV    ->SetBinContent(i,0);
        hWJets ->SetBinContent(i,0);
        hWG    ->SetBinContent(i,0);
        hHiggs ->SetBinContent(i,0);
        hData  ->SetBinContent(i,0);
      }
    }

    hBck    = (TH1F*)hWW->Clone();
    hBck->Add(hZJets,1);
    hBck->Add(hTop  ,1);  
    hBck->Add(hVV   ,1); 
    hBck->Add(hWJets,1);
    hBck->Add(hWG   ,1);
    hBck->Add(hHiggs,1);
    
    bool doRescale = false;
    if(doRescale){
      double sf = hData->GetSumOfWeights()/hBck->GetSumOfWeights();
      hBck  ->Scale(sf);
      hWW   ->Scale(sf);
      hZJets->Scale(sf);
      hTop  ->Scale(sf);
      hVV   ->Scale(sf);
      hWJets->Scale(sf);
      hWG   ->Scale(sf);
      hHiggs->Scale(sf);
    }

  }
  else if(nsel == 1){
    printf("WZ style plots\n");
    hEM     = (TH1F*)file->Get("histo1");
    hZJets  = (TH1F*)file->Get("histo2");
    hWZ     = (TH1F*)file->Get("histo3");
    hZZ     = (TH1F*)file->Get("histo4");
    hVVV    = (TH1F*)file->Get("histo5");
    hData   = (TH1F*)file->Get("histo0");
    hHiggs  = (TH1F*)file->Get("histo6");
    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=hData->GetNbinsX()-185; i<=hData->GetNbinsX(); i++){
        hEM	->SetBinContent(i,0);
        hZJets  ->SetBinContent(i,0);
        hWZ	->SetBinContent(i,0);
        hZZ	->SetBinContent(i,0);
        hVVV	->SetBinContent(i,0);
        hData	->SetBinContent(i,0);
        hHiggs	->SetBinContent(i,0);
      }
    }

    hBck    = (TH1F*)hEM->Clone();
    hBck->Add(hZJets,1);
    hBck->Add(hWZ   ,1);  
    hBck->Add(hZZ   ,1);
    hBck->Add(hVVV  ,1);
  }
  else if(nsel == 2){
    printf("ZH style plots\n");
    hEM     = (TH1F*)file->Get("histo1");
    hZJets  = (TH1F*)file->Get("histo2");
    hWZ     = (TH1F*)file->Get("histo3");
    hZZ     = (TH1F*)file->Get("histo4");
    hVVV    = (TH1F*)file->Get("histo5");
    hHiggs  = (TH1F*)file->Get("histo6");
    hggZH   = (TH1F*)file->Get("histo7");
    hData   = (TH1F*)file->Get("histo0");

    hEM	  ->Scale(scaling[0]);
    hZJets->Scale(scaling[1]);
    hWZ	  ->Scale(scaling[2]);
    hZZ	  ->Scale(scaling[3]);
    hVVV  ->Scale(scaling[4]);

    if(hasFileExtra) hHiggs2 = (TH1F*)fileExtra->Get("histo6");
    hZJets->Scale(lumi);
    //hHiggs->Scale(lumi);
    hHiggs->Scale(0.1);
    hggZH ->Scale(1);
    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=hData->GetNbinsX()-100; i<=hData->GetNbinsX(); i++){
        hEM	->SetBinContent(i,0);
        hZJets  ->SetBinContent(i,0);
        hWZ	->SetBinContent(i,0);
        hZZ	->SetBinContent(i,0);
        hVVV	->SetBinContent(i,0);
        hHiggs  ->SetBinContent(i,0);
        hggZH   ->SetBinContent(i,0);
        hData	->SetBinContent(i,0);
      }
    }
    hHiggs->Add(hggZH);
    hBck    = (TH1F*)hEM->Clone();
    hBck->Add(hZJets,1);
    hBck->Add(hWZ   ,1);  
    hBck->Add(hZZ   ,1); 
    hBck->Add(hVVV  ,1);
  }
  else if(nsel == 3){
    printf("ZG style plots\n");
    hEM     = (TH1F*)file->Get("histo1");
    hVV     = (TH1F*)file->Get("histo2");
    hWG     = (TH1F*)file->Get("histo3");
    hZJets  = (TH1F*)file->Get("histo4");
    hData   = (TH1F*)file->Get("histo0");
    hVV->Scale(lumi);
    hWG   ->Scale(1.0);
    hZJets->Scale(1.0);

    hBck    = (TH1F*)hEM->Clone();
    hBck->Add(hVV   ,1);
    hBck->Add(hZJets,1);
    hBck->Add(hWG   ,1);  
  }
  else if(nsel == 4){
    printf("WWSS OLD style plots\n");
    hData   = (TH1F*)file->Get("histo0");
    hWWEWK  = (TH1F*)file->Get("histo1");
    hWWQCD  = (TH1F*)file->Get("histo2");
    hVV     = (TH1F*)file->Get("histo3");
    hZZ     = (TH1F*)file->Get("histo4");
    hVV->Add(hZZ);
    hVVV    = (TH1F*)file->Get("histo5");
    hWS     = (TH1F*)file->Get("histo6");
    hWG     = (TH1F*)file->Get("histo7");
    hDPS    = (TH1F*)file->Get("histo8");
    hWJets  = (TH1F*)file->Get("histo9");
    TH1F * hFakeE   = (TH1F*)file->Get("histo10");
    hWJets->Add(hFakeE);
    //hWS->Scale(lumi);

    hWWEWK ->Scale(scaling[0]);
    hWWQCD ->Scale(scaling[1]);
    hVV    ->Scale(scaling[2]);
    hVVV   ->Scale(scaling[3]);
    hWS    ->Scale(scaling[4]);
    hWG    ->Scale(scaling[5]);
    hDPS   ->Scale(scaling[6]);
    hWJets ->Scale(scaling[7]);

    if(hasFileExtra) hHiggs  = (TH1F*)file->Get("histo11");
    if(hasFileExtra) hHiggs2 = (TH1F*)fileExtra->Get("histo12");

    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=hData->GetNbinsX(); i>=25; i--){
        hData   ->SetBinContent(i,0);
        //hWWEWK  ->SetBinContent(i,0);
        //hWWQCD  ->SetBinContent(i,0);
        //hVV     ->SetBinContent(i,0);
        //hVVV    ->SetBinContent(i,0);
        //hWS     ->SetBinContent(i,0);
        //hWG     ->SetBinContent(i,0);
        //hDPS    ->SetBinContent(i,0);
        //hWJets  ->SetBinContent(i,0);
      }
    }
    hBck    = (TH1F*)hWWEWK->Clone();
    hBck->Add(hWWQCD,1);
    hBck->Add(hVV   ,1);  
    hBck->Add(hVVV  ,1); 
    hBck->Add(hWS   ,1);
    hBck->Add(hWG   ,1);
    hBck->Add(hDPS  ,1);  
    hBck->Add(hWJets,1); 
  }
  else if(nsel == 5){
    printf("WW style plots\n");
    hWW     = (TH1F*)file->Get("histo1"); TH1F *hWWAux = (TH1F*)file->Get("histo2"); hWW->Add(hWWAux);
    hTop    = (TH1F*)file->Get("histo3");
    hZJets  = (TH1F*)file->Get("histo4");
    hVV     = (TH1F*)file->Get("histo5"); TH1F *hVVAux = (TH1F*)file->Get("histo6"); hVV->Add(hVVAux);
    hWG     = (TH1F*)file->Get("histo7"); TH1F *hWGAux = (TH1F*)file->Get("histo8"); hWG->Add(hWGAux);
    hWJets  = (TH1F*)file->Get("histo9"); TH1F *hWJetsAux = (TH1F*)file->Get("histo10"); hWJets->Add(hWJetsAux);
    hHiggs  = (TH1F*)file->Get("histo11");
    hData   = (TH1F*)file->Get("histo0");

    //hWW->Scale(lumi);
    //hTop->Scale(0.91);

    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=1; i<=6; i++){
        hWW    ->SetBinContent(i,0);
        hTop   ->SetBinContent(i,0);
        hZJets ->SetBinContent(i,0);
        hVV    ->SetBinContent(i,0);
        hWG    ->SetBinContent(i,0);
        hWJets ->SetBinContent(i,0);
        hHiggs ->SetBinContent(i,0);
        hData  ->SetBinContent(i,0);
      }
    }

    hBck    = (TH1F*)hWW->Clone();
    hBck->Add(hTop  ,1);  
    hBck->Add(hZJets,1); 
    hBck->Add(hVV   ,1);
    hBck->Add(hWG   ,1);  
    hBck->Add(hWJets,1); 
    hBck->Add(hHiggs,1);
    
    bool doRescale = false;
    if(doRescale){
      double sf = hData->GetSumOfWeights()/hBck->GetSumOfWeights();
      hBck  ->Scale(sf);
      hWW   ->Scale(sf);
      hTop  ->Scale(sf);
      hZJets->Scale(sf);
      hVV   ->Scale(sf);
      hWG   ->Scale(sf);
      hWJets->Scale(sf);
      hHiggs->Scale(sf);
    }
  }
  else if(nsel == 6){
    printf("WWSS style plots\n");
    hData   = (TH1F*)file->Get("histo0");
    hWWEWK  = (TH1F*)file->Get("histo1");
    hWWQCD  = (TH1F*)file->Get("histo2");
    hWZ     = (TH1F*)file->Get("histo3");
    hZZ     = (TH1F*)file->Get("histo4");
    hVVV    = (TH1F*)file->Get("histo5");
    hWS     = (TH1F*)file->Get("histo6");
    hWG     = (TH1F*)file->Get("histo7");
    hDPS    = (TH1F*)file->Get("histo8");
    hWJets  = (TH1F*)file->Get("histo9");
    TH1F * hFakeE   = (TH1F*)file->Get("histo10");
    hWJets->Add(hFakeE);

    hOther  = (TH1F*)hWWQCD->Clone();
    hOther->Add(hZZ );
    hOther->Add(hVVV);
    hOther->Add(hWS );
    hOther->Add(hWG );
    hOther->Add(hDPS);

    hWWEWK ->Scale(scaling[0]);
    hOther ->Scale(scaling[1]);
    hWZ    ->Scale(scaling[2]);
    hWJets ->Scale(scaling[3]);

    if(hasFileExtra) hHiggs  = (TH1F*)file->Get("histo11");
    if(hasFileExtra) hHiggs2 = (TH1F*)fileExtra->Get("histo12");

    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=hData->GetNbinsX(); i>=25; i--){
        hData  ->SetBinContent(i,0);
        //hWWEWK->SetBinContent(i,0);
        //hOther->SetBinContent(i,0);
        //hWZ	->SetBinContent(i,0);
        //hWJets->SetBinContent(i,0);
      }
    }
    hBck    = (TH1F*)hWWEWK->Clone();
    hBck->Add(hOther,1);
    hBck->Add(hWZ   ,1);  
    hBck->Add(hWJets,1);
  }
  else {
    return;
  }

  if     (nsel == 0 || nsel == 5){
    if(hWW   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWW   , (TH1F*)hWW   ->Clone("hWW"));
    if(hZJets->GetSumOfWeights() > 0) myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));
    if(hTop  ->GetSumOfWeights() > 0) myPlot.setMCHist(iTop  , (TH1F*)hTop  ->Clone("hTop"));
    if(hVV   ->GetSumOfWeights() > 0) myPlot.setMCHist(iVV   , (TH1F*)hVV   ->Clone("hVV")); 
    if(hWJets->GetSumOfWeights() > 0) myPlot.setMCHist(iWJets, (TH1F*)hWJets->Clone("hWJets"));
    if(hWG   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWG   , (TH1F*)hWG   ->Clone("hWG"));
    if(hHiggs->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs, (TH1F*)hHiggs->Clone("hHiggs"));

    printf("%f + %f + %f + %f + %f + %f + %f = %f - %f\n",
          hWW->GetSumOfWeights(),hZJets->GetSumOfWeights(),hTop->GetSumOfWeights(),hVV->GetSumOfWeights(),
  	  hWJets->GetSumOfWeights(),hWG->GetSumOfWeights(),hHiggs->GetSumOfWeights(),
	  hBck->GetSumOfWeights(),
	  hData->GetSumOfWeights());
  }
  else if(nsel == 1){
    if(hEM   ->GetSumOfWeights() > 0) myPlot.setMCHist(iEM   , (TH1F*)hEM   ->Clone("hEM"));
    if(hZJets->GetSumOfWeights() > 0) myPlot.setMCHist(iWG   , (TH1F*)hZJets->Clone("hZJets"));
    if(hWZ   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWZ   , (TH1F*)hWZ   ->Clone("hWZ"));
    if(hZZ   ->GetSumOfWeights() > 0) myPlot.setMCHist(iZZ   , (TH1F*)hZZ   ->Clone("hZZ"));
    if(hVVV  ->GetSumOfWeights() > 0) myPlot.setMCHist(iVVV  , (TH1F*)hVVV  ->Clone("hVVV"));
    if(hHiggs->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs, (TH1F*)hHiggs->Clone("hHiggs"));
    myPlot.setOverlaid(false);
 
    printf("%f + %f + %f + %f + %f = %f - %f - %f\n",
          hEM->GetSumOfWeights(),hZJets->GetSumOfWeights(),hWZ->GetSumOfWeights(),hZZ->GetSumOfWeights(),hVVV->GetSumOfWeights(),
	  hBck->GetSumOfWeights(),hData->GetSumOfWeights(),hHiggs->GetSumOfWeights());
  }
  else if(nsel == 2){
    if(hEM    ->GetSumOfWeights() > 0) myPlot.setMCHist(iEM   , (TH1F*)hEM   ->Clone("hEM"));
    if(hZJets ->GetSumOfWeights() > 0) myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));
    if(hWZ    ->GetSumOfWeights() > 0) myPlot.setMCHist(iWZ   , (TH1F*)hWZ   ->Clone("hWZ"));
    if(hZZ    ->GetSumOfWeights() > 0) myPlot.setMCHist(iZZ   , (TH1F*)hZZ   ->Clone("hZZ")); 
    if(hVVV   ->GetSumOfWeights() > 0) myPlot.setMCHist(iVVV  , (TH1F*)hVVV  ->Clone("hVVV"));
    if(hHiggs ->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs, (TH1F*)hHiggs->Clone("hHiggs"));
    if(hasFileExtra && 
       hHiggs2->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs2,(TH1F*)hHiggs2->Clone("hHiggs2"));
    if(hasFileExtra && 
       hHiggs2->GetSumOfWeights() > 0) printf("Higgs2: %f\n",hHiggs2->GetSumOfWeights());
    
    myPlot.setOverlaid(false);
    myPlot.setLabelEM(" Nonresonant");
    myPlot.setLabelVVV("Other bkg.");

    printf("%f + %f + %f + %f + %f = %f - %f - %f\n",
          hEM->GetSumOfWeights(),hZJets->GetSumOfWeights(),hWZ->GetSumOfWeights(),hZZ->GetSumOfWeights(),hVVV->GetSumOfWeights(),
	  hBck->GetSumOfWeights(),
	  hData->GetSumOfWeights(),hHiggs->GetSumOfWeights());
  }
  else if(nsel == 3){
    if(hEM   ->GetSumOfWeights() > 0) myPlot.setMCHist(iEM   , (TH1F*)hEM   ->Clone("hEM"));
    if(hVV   ->GetSumOfWeights() > 0) myPlot.setMCHist(iVV   , (TH1F*)hVV   ->Clone("hVV"));
    if(hWG   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWG   , (TH1F*)hWG   ->Clone("hWG"));
    if(hZJets->GetSumOfWeights() > 0) myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));

    printf("%f + %f + %f + %f = %f - %f\n",
          hEM->GetSumOfWeights(),hVV->GetSumOfWeights(),hWG->GetSumOfWeights(),hZJets->GetSumOfWeights(),
	  hBck->GetSumOfWeights(),
	  hData->GetSumOfWeights());
  }
  else if(nsel == 4){
    if(hWWEWK->GetSumOfWeights() > 0) myPlot.setMCHist(iWWEWK, (TH1F*)hWWEWK->Clone("hWWEWK"));
    if(hWWQCD->GetSumOfWeights() > 0) myPlot.setMCHist(iWWQCD, (TH1F*)hWWQCD->Clone("hWWQCD"));
    if(hVV   ->GetSumOfWeights() > 0) myPlot.setMCHist(iVV   , (TH1F*)hVV   ->Clone("hVV"));
    if(hVVV  ->GetSumOfWeights() > 0) myPlot.setMCHist(iVVV  , (TH1F*)hVVV  ->Clone("hVVV"));
    if(hWS   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWS   , (TH1F*)hWS   ->Clone("hWS"));
    if(hWG   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWG   , (TH1F*)hWG   ->Clone("hWG"));
    if(hDPS  ->GetSumOfWeights() > 0) myPlot.setMCHist(iDPS  , (TH1F*)hDPS  ->Clone("hDPS"));
    if(hWJets->GetSumOfWeights() > 0) myPlot.setMCHist(iWJets, (TH1F*)hWJets->Clone("hWJets"));

    myPlot.setOverlaid(false);
    if(hasFileExtra && hHiggs ->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs,(TH1F*)hHiggs->Clone("hHiggs"));
    if(hasFileExtra && hHiggs ->GetSumOfWeights() > 0) printf("Higgs: %f\n",hHiggs->GetSumOfWeights());
    if(hasFileExtra && hHiggs2->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs2,(TH1F*)hHiggs2->Clone("hHiggs2"));
    if(hasFileExtra && hHiggs2->GetSumOfWeights() > 0) printf("Higgs2: %f\n",hHiggs2->GetSumOfWeights());

    printf("%f + %f + %f + %f + %f + %f + %f + %f = %f - %f\n",
          hWWEWK->GetSumOfWeights(),hWWQCD->GetSumOfWeights(),hVV->GetSumOfWeights(),hVVV->GetSumOfWeights(),
          hWS->GetSumOfWeights(),hWG->GetSumOfWeights(),hDPS->GetSumOfWeights(),hWJets->GetSumOfWeights(),
	  hBck->GetSumOfWeights(),
	  hData->GetSumOfWeights());
  }
  else if(nsel == 6){
    if(hWJets->GetSumOfWeights() > 0) myPlot.setMCHist(iWJets, (TH1F*)hWJets->Clone("hWJets"));
    if(hWWEWK->GetSumOfWeights() > 0) myPlot.setMCHist(iWWEWK, (TH1F*)hWWEWK->Clone("hWWEWK"));
    if(hOther->GetSumOfWeights() > 0) myPlot.setMCHist(iOther, (TH1F*)hOther->Clone("hOther"));
    if(hWZ   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWZ   , (TH1F*)hWZ   ->Clone("hWZ"));

    myPlot.setOverlaid(false);
    if(hasFileExtra && hHiggs ->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs,(TH1F*)hHiggs->Clone("hHiggs"));
    if(hasFileExtra && hHiggs ->GetSumOfWeights() > 0) printf("Higgs: %f\n",hHiggs->GetSumOfWeights());
    if(hasFileExtra && hHiggs2->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs2,(TH1F*)hHiggs2->Clone("hHiggs2"));
    if(hasFileExtra && hHiggs2->GetSumOfWeights() > 0) printf("Higgs2: %f\n",hHiggs2->GetSumOfWeights());

    printf("%f + %f + %f + %f = %f - %f\n",
          hWWEWK->GetSumOfWeights(),hOther->GetSumOfWeights(),hWZ->GetSumOfWeights(),hWJets->GetSumOfWeights(),
	  hBck->GetSumOfWeights(),
	  hData->GetSumOfWeights());
  }

  if(isBlind == false) myPlot.setDataHist(hData);

  TCanvas* c1 = new TCanvas("c1", "c1",5,5,500,500);

  if(show2D==false){
  if(isLogY == true) c1->SetLogy();
  myPlot.Draw(ReBin);  // Can pass a rebin 
  CMS_lumi( c1, 4, 12 );
  } else {
  c1->SetBottomMargin(0.1);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1",0.00,0.30,1.00,1.00);
  TPad *pad2 = new TPad("pad2", "pad2",0.00,0.00,1.00,0.30);
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  if(isLogY == true) c1->SetLogy();
  if(isLogY == true) pad1->SetLogy();
  myPlot.Draw(ReBin);
  CMS_lumi( pad1, 4, 33 );

  pad2->cd();
  //pad2->SetGridy();
  hBck ->Rebin(ReBin);
  TH1D* hTotalDivision = (TH1D*) hData->Clone();
  hTotalDivision->Reset();
  TH1D* hDataDivision = (TH1D*) hData->Clone();
  hDataDivision->Reset();
  TH1D* hRatio = (TH1D*) hData->Clone();
  hRatio->Reset();
  TH1D* hBand = (TH1D*) hData->Clone();
  hBand->Reset();

  hDataDivision ->Add(hData );
  hTotalDivision->Add(hBck  );
  hTotalDivision->Scale(hDataDivision->GetSumOfWeights()/hTotalDivision->GetSumOfWeights());
  TGraphAsymmErrors *g = new TGraphAsymmErrors(hDataDivision);
  bool showPulls = false;
  bool useGarwood = false;
  for(int i=1; i<=hDataDivision->GetNbinsX(); i++){
    if(showPulls){
      double pull = 0.0;
      if(useGarwood == false){
        if((hDataDivision->GetBinError(i) > 0 || hTotalDivision->GetBinError(i) > 0) && hDataDivision->GetBinContent(i) > 0){
          pull = (hDataDivision->GetBinContent(i)-hTotalDivision->GetBinContent(i))/sqrt(hDataDivision->GetBinError(i)*hDataDivision->GetBinError(i)+hTotalDivision->GetBinError(i)*hTotalDivision->GetBinError(i));
        }
      } 
      else {
        pull = TMath::Sqrt(2)*TMath::ErfInverse(-1+2*ROOT::Math::poisson_cdf(hDataDivision->GetBinContent(i),hTotalDivision->GetBinContent(i))); 
      }
      hRatio->SetBinContent(i,pull);
      hRatio->SetBinError(i,1.0);
      printf("pull(%3d): da: %f pred: %f --> %f\n",i,hDataDivision ->GetBinContent(i),hTotalDivision->GetBinContent(i),pull);
    } 
    else {
      double N = g->GetY()[i-1];
      double alpha=(1-0.6827);
      double L = (N==0) ? 0 : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      double U = ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
      //double diff = hDataDivision->GetBinError(i)
      double diff = (U-L)/2;
      if( N != hDataDivision ->GetBinContent(i)) cout << "PROBLEM" << endl;
      double pull = 1.0; double pullerr = 0.0;
      if(hDataDivision->GetBinContent(i) > 0 && hTotalDivision->GetBinContent(i) > 0){
        pull = (hDataDivision->GetBinContent(i)/hTotalDivision->GetBinContent(i));
	pullerr = pull*diff/hDataDivision->GetBinContent(i);
      }
      //if(pull>1.1 && i<10) pull = 1.1-gRandom->Uniform()*0.01;
      //if(pull<0.9 && i<10) pull = 0.9+gRandom->Uniform()*0.01;
      hRatio->SetBinContent(i,pull);
      hRatio->SetBinError(i,pullerr);
      hBand->SetBinContent(i,1);
      hBand->SetBinError(i,hTotalDivision->GetBinError(i)/hTotalDivision->GetBinContent(i));
      //printf("ratio(%3d): %f +/- %f --> da: %f +/- %f (%f) pred: %f +/- %f\n",i,pull,pullerr,hDataDivision ->GetBinContent(i),hDataDivision ->GetBinError(i),diff,hTotalDivision->GetBinContent(i),hTotalDivision->GetBinError(i));
    }
  }
  if(showPulls) atributes(hRatio,XTitle.Data(),"Pull",units.Data());
  else          atributes(hRatio,XTitle.Data(),"Data/Bkg.",units.Data());
  hRatio->Draw("e");
  hBand->SetFillColor(12);
  hBand->SetFillStyle(3002);
  hBand->SetMarkerSize(0);
  hBand->SetLineWidth(0);
  hBand->Draw("E2same");

  // Draw a line throgh y=0
  double theLines[2] = {1.0, 0.5};
  if(showPulls) {theLines[0] = 0.0; theLines[1] = 1.5;}
  TLine* baseline = new TLine(hRatio->GetXaxis()->GetXmin(), theLines[0],
                              hRatio->GetXaxis()->GetXmax(), theLines[0]);
  baseline->SetLineStyle(kDashed);
  baseline->Draw();
  // Set the y-axis range symmetric around y=0
  Double_t dy = TMath::Max(TMath::Abs(hRatio->GetMaximum()),
                           TMath::Abs(hRatio->GetMinimum())) + theLines[1];
  if(showPulls) hRatio->GetYaxis()->SetRangeUser(-dy, +dy);
  else          hRatio->GetYaxis()->SetRangeUser(0.2, +1.8);
  hRatio->GetYaxis()->CenterTitle();
  eraselabel(pad1,hData->GetXaxis()->GetLabelSize());
  }

  char CommandToExec[300];
  sprintf(CommandToExec,"mkdir -p plots");
  gSystem->Exec(CommandToExec);  

  if(strcmp(outputName.Data(),"") != 0){
    TString myOutputFile;
    myOutputFile = Form("plots/%s.eps",outputName.Data());
    //c1->SaveAs(myOutputFile.Data());
    myOutputFile = Form("plots/%s.png",outputName.Data());
    //c1->SaveAs(myOutputFile.Data());
    myOutputFile = Form("plots/%s.pdf",outputName.Data());
    c1->SaveAs(myOutputFile.Data());
    myOutputFile = Form("plots/%s.root",outputName.Data());
    c1->SaveAs(myOutputFile.Data());
  }

  bool computePU = false;
  if(computePU){
    hBck->Scale(hData->GetSumOfWeights()/hBck->GetSumOfWeights());
    TH1D * puWeights =  (TH1D*)hData->Clone("puWeights");
    puWeights->Sumw2();
    puWeights->Divide(hBck);
    TFile output("puWeights_13TeV.root","RECREATE");
    puWeights->Write();
    output.Close();
  }

}
