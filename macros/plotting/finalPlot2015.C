#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "StandardPlot2015.C"
#include "TSystem.h"
#include "CMS_lumi.C"

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
  histo->GetXaxis()->SetLabelSize  (0.100);
  histo->GetXaxis()->SetNdivisions (  505);
  histo->GetXaxis()->SetTitleFont  (   42);
  histo->GetXaxis()->SetTitleOffset(  1.1);
  histo->GetXaxis()->SetTitleSize  (0.100);
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

void finalPlot2015(int nsel = 0, int ReBin = 1, TString XTitle = "N_{jets}", TString units = "", TString plotName = "histo_nice.root", TString outputName = "njets",
                bool isLogY = false, TString higgsLabel = "", double lumi = 1.0, bool isBlind = false, TString extraLabel = "") {

  gInterpreter->ExecuteMacro("GoodStyle.C");
  gROOT->LoadMacro("StandardPlot2015.C");
  gStyle->SetOptStat(0);

  TFile* file = new TFile(plotName, "read");

  StandardPlot2015 myPlot;
  myPlot.setLumi(lumi);
  myPlot.setLabel(XTitle);
  myPlot.addLabel(extraLabel.Data());
  myPlot.setHiggsLabel(higgsLabel.Data());
  myPlot.setUnits(units);

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

  if     (nsel == 0){
    printf("WW style plots\n");
    hWW     = (TH1F*)file->Get("histo1");
    hZJets  = (TH1F*)file->Get("histo2");
    hTop    = (TH1F*)file->Get("histo3");
    hVV     = (TH1F*)file->Get("histo4");
    hWJets  = (TH1F*)file->Get("histo5");
    hWG     = (TH1F*)file->Get("histo6");
    hHiggs  = (TH1F*)file->Get("histo7");
    hData   = (TH1F*)file->Get("histo0");
    hZJets->Scale(lumi);

    hBck    = (TH1F*)hWW->Clone();
    hBck->Add(hZJets,1);
    hBck->Add(hTop  ,1);  
    hBck->Add(hVV   ,1); 
    hBck->Add(hWJets,1);
    hBck->Add(hWG   ,1);
    hBck->Add(hHiggs,1);
  }
  else if(nsel == 1){
    printf("WZ style plots\n");
    hEM     = (TH1F*)file->Get("histo1");
    hZJets  = (TH1F*)file->Get("histo2");
    hWZ     = (TH1F*)file->Get("histo3");
    hZZ     = (TH1F*)file->Get("histo4");
    hVVV    = (TH1F*)file->Get("histo5");
    hData   = (TH1F*)file->Get("histo0");
    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=hData->GetNbinsX()-185; i<=hData->GetNbinsX(); i++){
        hEM	->SetBinContent(i,0);
        hZJets  ->SetBinContent(i,0);
        hWZ	->SetBinContent(i,0);
        hZZ	->SetBinContent(i,0);
        hVVV	->SetBinContent(i,0);
        hData	->SetBinContent(i,0);
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
    hZJets->Scale(lumi);
    hHiggs->Scale(1);
    hggZH ->Scale(1);
    bool doRemoveBins = false;
    if(doRemoveBins){
      for(int i=hData->GetNbinsX()-11; i<=hData->GetNbinsX(); i++){
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
  else {
    return;
  }

  if     (nsel == 0){
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

    printf("%f + %f + %f + %f + %f = %f - %f\n",
          hEM->GetSumOfWeights(),hZJets->GetSumOfWeights(),hWZ->GetSumOfWeights(),hZZ->GetSumOfWeights(),hVVV->GetSumOfWeights(),
	  hBck->GetSumOfWeights(),
	  hData->GetSumOfWeights());
  }
  else if(nsel == 2){
    if(hEM   ->GetSumOfWeights() > 0) myPlot.setMCHist(iEM   , (TH1F*)hEM   ->Clone("hEM"));
    if(hZJets->GetSumOfWeights() > 0) myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));
    if(hWZ   ->GetSumOfWeights() > 0) myPlot.setMCHist(iWZ   , (TH1F*)hWZ   ->Clone("hWZ"));
    if(hZZ   ->GetSumOfWeights() > 0) myPlot.setMCHist(iZZ   , (TH1F*)hZZ   ->Clone("hZZ")); 
    if(hVVV  ->GetSumOfWeights() > 0) myPlot.setMCHist(iVVV  , (TH1F*)hVVV  ->Clone("hVVV"));
    if(hHiggs->GetSumOfWeights() > 0) myPlot.setMCHist(iHiggs, (TH1F*)hHiggs->Clone("hHiggs"));
    myPlot.setOverlaid(false);
    myPlot.setLabelEM("WW+top-quark");

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

  if(isBlind == false) myPlot.setDataHist(hData);

  TCanvas* c1 = new TCanvas("c1", "c1",5,5,500,500);
  bool show2D = false;

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
  CMS_lumi( c1, 4, 12 );

  pad2->cd();
  bool showPulls = false;
  //pad2->SetGridy();
  hBck ->Rebin(ReBin);
  TH1D* hTotalDivision = (TH1D*) hData->Clone();
  hTotalDivision->Reset();
  TH1D* hDataDivision = (TH1D*) hData->Clone();
  hDataDivision->Reset();
  TH1D* hRatio = (TH1D*) hData->Clone();
  hRatio->Reset();

  hDataDivision ->Add(hData );
  hTotalDivision->Add(hBck  );
  for(int i=1; i<=hDataDivision->GetNbinsX(); i++){
    if(showPulls){
      double pull = 0.0;
      if((hDataDivision->GetBinError(i) > 0 || hTotalDivision->GetBinError(i) > 0) && hDataDivision->GetBinContent(i) > 0){
        pull = (hDataDivision->GetBinContent(i)-hTotalDivision->GetBinContent(i))/sqrt(hDataDivision->GetBinError(i)*hDataDivision->GetBinError(i)+hTotalDivision->GetBinError(i)*hTotalDivision->GetBinError(i));
      }
      hRatio->SetBinContent(i,pull);
      hRatio->SetBinError(i,1.0);
    } else {
      double pull = 1.0; double pullerr = 0.0;
      if(hDataDivision->GetBinContent(i) > 0 && hTotalDivision->GetBinContent(i) > 0){
        pull = (hDataDivision->GetBinContent(i)/hTotalDivision->GetBinContent(i));
	pullerr = pull*sqrt(hDataDivision ->GetBinError(i)*hDataDivision ->GetBinError(i)/hDataDivision ->GetBinContent(i)/hDataDivision ->GetBinContent(i)+
	                    hTotalDivision->GetBinError(i)*hTotalDivision->GetBinError(i)/hTotalDivision->GetBinContent(i)/hTotalDivision->GetBinContent(i));
      }
      hRatio->SetBinContent(i,pull);
      hRatio->SetBinError(i,pullerr);
      //printf("ratio(%d): %f +/- %f\n",i,pull,pullerr);
    }
  }
  if(showPulls) atributes(hRatio,XTitle.Data(),"Pull",units.Data());
  else          atributes(hRatio,XTitle.Data(),"Data/Bkg.",units.Data());
  hRatio->Draw("e");

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
  else          hRatio->GetYaxis()->SetRangeUser(0.5, +1.5);
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
