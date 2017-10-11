#include<vector>

//#if !defined (__CINT__) || defined (__MAKECINT__)
#include "THStack.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TExec.h"
#include <iostream>
#include "TPaveText.h"
//#endif

Bool_t isHWWOverlaid = true;
enum samp { iWW, iOther, iVV, iWJets, iZJets, iEM, iVVV, iTop, iWWQCD, iHiggs, iHiggs2, iWZ, iZZ, iWG, iDPS,  iWS, iWWEWK, nSamples };

float xPos[nSamples+1] = {0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.30,0.30,0.30,0.30,0.30,0.30,0.30,0.30,0.30}; 
float yOff[nSamples+1] = {   0,   1,   2,   3,   4,   5,   6,   7,   8,   0,   1,   2,   3,   4,   5,   6,   7,   8};

const Float_t _tsize   = 0.04;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;

const bool doApplyBinWidth = true;
const double SFBinWidth = 1;

//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1F* h)
{
    Float_t maxWithErrors = 0;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if(doApplyBinWidth) binHeight = binHeight/h->GetBinWidth(i)*SFBinWidth;

        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
}

//------------------------------------------------------------------------------
// GetMinimumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMinimumIncludingErrors(TH1F* h)
{
    Float_t minWithErrors = 100000000000;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if(doApplyBinWidth) binHeight = binHeight/h->GetBinWidth(i)*SFBinWidth;

        if (binHeight < minWithErrors) minWithErrors = binHeight;
    }

    return minWithErrors;
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis*  axis,
        TString coordinate,
        TString title)
{
    axis->SetLabelFont  (   42);
    axis->SetLabelOffset(0.015);
    axis->SetLabelSize  (0.050);
    axis->SetNdivisions (  505);
    axis->SetTitleFont  (   42);
    axis->SetTitleOffset(  1.3);
    axis->SetTitleSize  (0.070);

    if      (coordinate == "y" && doApplyBinWidth == true) axis->SetTitleOffset(1.1);
    else if (coordinate == "y")                            axis->SetTitleOffset(1.1);

    axis->SetTitle(title);
}


//------------------------------------------------------------------------------
// THStackAxisFonts
//------------------------------------------------------------------------------
void THStackAxisFonts(THStack* h,
        TString  coordinate,
        TString  title)
{
    TAxis* axis = NULL;

    if (coordinate.Contains("x")) axis = h->GetHistogram()->GetXaxis();
    if (coordinate.Contains("y")) axis = h->GetHistogram()->GetYaxis();

    AxisFonts(axis, coordinate, title);
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend(Float_t x1,
        Float_t y1,
        TH1F*   hist,
        TString label,
        TString option)
{
    TLegend* legend = new TLegend(x1,
            y1,
            x1 + _xoffset,
            y1 + _yoffset);

    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (_tsize);

    legend->AddEntry(hist, label.Data(), option.Data());

    legend->Draw();
}

void DrawLegendTG(Float_t x1,
        Float_t y1,
        TGraph  *hist,
        TString label,
        TString option)
{
    TLegend* legend = new TLegend(x1,
            y1,
            x1 + _xoffset,
            y1 + _yoffset);

    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (_tsize);

    legend->AddEntry(hist, label.Data(), option.Data());

    legend->Draw();
}

class StandardPlot2016 {

    public: 
        StandardPlot2016() { _hist.resize(nSamples,0); _data = 0; _breakdown = false; _HiggsLabel = "";_Higgs2Label = "";_labelEM = " Nonprompt";_labelVVV = " VVV";}
        void setMCHist   (const samp &s, TH1F * h)  { _hist[s]       = h;  } 
        void setDataHist  (TH1F * h) { _data	      = h;  }
        void setWWHist    (TH1F * h) { setMCHist(iWW   ,h); }
        void setDSPHist   (TH1F * h) { setMCHist(iDPS  ,h); }
        void setWWEWKHist (TH1F * h) { setMCHist(iWWEWK,h); }
        void setWWQCDHist (TH1F * h) { setMCHist(iWWQCD,h); }
        void setWSHist    (TH1F * h) { setMCHist(iWS   ,h); }
        void setZJetsHist (TH1F * h) { setMCHist(iZJets,h); }
        void setTopHist   (TH1F * h) { setMCHist(iTop  ,h); }
        void setVVHist    (TH1F * h) { setMCHist(iVV   ,h); }
        void setWJetsHist (TH1F * h) { setMCHist(iWJets,h); }
        void setWGHist    (TH1F * h) { setMCHist(iWG   ,h); }
        void setHiggsHist (TH1F * h) { setMCHist(iHiggs,h); }
        void setHiggs2Hist(TH1F * h) { setMCHist(iHiggs2,h); }
        void setWZHist    (TH1F * h) { setMCHist(iWZ   ,h); }
        void setZZHist    (TH1F * h) { setMCHist(iZZ   ,h); }
        void setVVVHist   (TH1F * h) { setMCHist(iVVV  ,h); }
        void setEMHist    (TH1F * h) { setMCHist(iEM   ,h); }
        void setOtherHist (TH1F * h) { setMCHist(iOther   ,h); }
        void setOverlaid  (bool b)   { isHWWOverlaid = b;   }
        void setLabelEM   (TString s){ _labelEM  = s.Data();}
        void setLabelVVV  (TString s){ _labelVVV = s.Data();}

  TH1F* getDataHist() { return _data; }

        void setHiggsLabel (TString s) {_HiggsLabel=s;}
        void setHiggs2Label(TString s) {_Higgs2Label=s;}

        TH1* DrawAndRebinTo(const int &rebinTo) {

            if(rebinTo == 0) return Draw();
            int rebin = 0, nbins = 0;
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;

                nbins = _hist[i]->GetNbinsX();
            }
            if (nbins == 0) return Draw();

            rebin = nbins / rebinTo;
            while(nbins % rebin != 0) rebin--;
            return Draw(rebin);

        }

        TH1* Draw(const int &rebin=1, double xmin=0, double xmax=0) {

            Color_t _sampleColor[nSamples];
            _sampleColor[iWW    ] = kAzure-9;
            _sampleColor[iZJets ] = 901;//kGreen+2;
            _sampleColor[iTop   ] = kYellow;
            _sampleColor[iVV    ] = kAzure-2;
            _sampleColor[iWJets ] = TColor::GetColor(155,152,204);
            _sampleColor[iWG    ] = kViolet-9;
            _sampleColor[iHiggs ] = 809;//kOrange+7;
            _sampleColor[iHiggs2] = 419;//kOrange+7;
            _sampleColor[iWZ    ] = kAzure-4;//TColor::GetColor(222,90,106);
            _sampleColor[iZZ    ] = 842;//kBlue-3;
            _sampleColor[iVVV   ] = 832;//kRed+1;
            _sampleColor[iEM    ] = 798;//kYellow;
            _sampleColor[iWWEWK ] = TColor::GetColor(248,206,104);
            _sampleColor[iWWQCD ] = kYellow;
            _sampleColor[iDPS   ] = kGreen;;//kGreen+2;
            _sampleColor[iWS    ] = 419;
            _sampleColor[iOther ] = TColor::GetColor(250,202,255);

            //setUpStyle();
            //if(!gPad) new TCanvas();

            THStack* hstack = new THStack();
	    TH1D* hSum;
	    if     (_hist[iZJets]) hSum = (TH1D*)_hist[iZJets]->Clone();
	    else if(_hist[iWZ])    hSum = (TH1D*)_hist[iWZ]   ->Clone();
	    else if(_hist[iVV])    hSum = (TH1D*)_hist[iVV]   ->Clone();
	    else if(_data)         hSum = (TH1D*)_data        ->Clone();
	    else if(_hist[iZZ])    hSum = (TH1D*)_hist[iZZ]   ->Clone();
	    else if(_hist[iWJets]) hSum = (TH1D*)_hist[iWJets]->Clone();
	    else                   {printf("PROBLEM\n");assert(0);}
	    hSum->Rebin(rebin);
	    hSum->Scale(0.0);
  	    if(xmin!=xmax) hSum->GetXaxis()->SetRangeUser(xmin,xmax);
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;
  	    	bool modifyXAxis = false;
		if(modifyXAxis == true){
		  TAxis *xa =_hist[i]->GetXaxis();
  	    	  for(Int_t k=1;k<=_hist[i]->GetNbinsX();++k){
  	    	    xa->SetBinLabel(1 ,"2#mu+#mu");
  	    	    xa->SetBinLabel(2 ,"2#mu+e");
  	    	    xa->SetBinLabel(3 ,"2e+#mu");
  	    	    xa->SetBinLabel(4 ,"2e+e");
  	    	    xa->SetRangeUser(1,4);
  	    	  }
		}
                _hist[i]->Rebin(rebin);
                _hist[i]->SetLineColor(_sampleColor[i]);

                if(doApplyBinWidth == true){
                  for(int nb=1; nb<=_hist[i]->GetNbinsX(); nb++){
		    _hist[i]->SetBinContent(nb, _hist[i]->GetBinContent(nb)/_hist[i]->GetBinWidth(nb)*SFBinWidth);
		    _hist[i]->SetBinError  (nb, _hist[i]->GetBinError(nb)  /_hist[i]->GetBinWidth(nb)*SFBinWidth);
                  }
                }

                // signal gets overlaid
                if (i == iHiggs  && isHWWOverlaid == false) continue;
                if (i == iHiggs2) continue;

                _hist[i]->SetFillColor(_sampleColor[i]);
                _hist[i]->SetFillStyle(1001);
  	         if(xmin!=xmax) _hist[i]->GetXaxis()->SetRangeUser(xmin,xmax);

                hstack->Add(_hist[i]);
		hSum->Add(_hist[i]);
            }

            if(_hist[iHiggs] ) _hist[iHiggs ]->SetLineWidth(4);
            if(_hist[iHiggs2]) _hist[iHiggs2]->SetLineWidth(4);
            if(_hist[iHiggs2]) _hist[iHiggs2]->SetLineStyle(2);
            if(_data) _data->Rebin(rebin);
	    //_data->SetBinContent(1,_data->GetBinContent(1)+100);
	    //_data->SetBinContent(4,_data->GetBinContent(4)-50);
	    //_data->SetBinContent(5,_data->GetBinContent(5)-20);
	    //_data->SetBinContent(6,_data->GetBinContent(6)-50);
	    //_data->SetBinContent(9,_data->GetBinContent(9)+10);
	    //_data->SetBinContent(5,_data->GetBinContent(5)+2);
            if(_data) _data->SetLineColor  (kBlack);
            if(_data) _data->SetMarkerStyle(kFullCircle);
	    hstack->Draw("hist");
  	    if(xmin!=xmax) hstack->GetXaxis()->SetRangeUser(xmin,xmax);
            gPad->Update();

  	    TGraphAsymmErrors * gsyst = new TGraphAsymmErrors(hSum);
	    bool plotSystErrorBars = true;
	    if(plotSystErrorBars == true) {
              for (int i = 0; i < gsyst->GetN(); ++i) {
                double systBck = 0;
		if(_hist[iWW	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iWW    ]->GetBinContent(i+1),2);
		if(_hist[iWWEWK ]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iWWEWK ]->GetBinContent(i+1),2);
		if(_hist[iWWQCD ]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iWWQCD ]->GetBinContent(i+1),2);
		if(_hist[iDPS	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iDPS   ]->GetBinContent(i+1),2);
		if(_hist[iWS    ]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iWS    ]->GetBinContent(i+1),2);
		if(_hist[iZJets ]) systBck = systBck + 1.000*TMath::Power(1.000*_hist[iZJets ]->GetBinContent(i+1),2);
		if(_hist[iTop	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iTop   ]->GetBinContent(i+1),2);
		if(_hist[iVV	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iVV    ]->GetBinContent(i+1),2);
		if(_hist[iWJets ]) systBck = systBck + 1.000*TMath::Power(0.360*_hist[iWJets ]->GetBinContent(i+1),2);
		if(_hist[iWG	]) systBck = systBck + 1.000*TMath::Power(0.300*_hist[iWG    ]->GetBinContent(i+1),2);
		if(_hist[iWZ	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iWZ    ]->GetBinContent(i+1),2);
		if(_hist[iZZ	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iZZ    ]->GetBinContent(i+1),2);
		if(_hist[iVVV	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iVVV   ]->GetBinContent(i+1),2);
		if(_hist[iEM	]) systBck = systBck + 1.000*TMath::Power(0.400*_hist[iEM    ]->GetBinContent(i+1),2);
		if(_hist[iOther ]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iOther ]->GetBinContent(i+1),2);
		if(_hist[iHiggs	]) systBck = systBck + 1.000*TMath::Power(0.200*_hist[iHiggs ]->GetBinContent(i+1),2);
                double total = hSum->GetBinContent(i+1);
                if(total > 0) systBck = sqrt(systBck)/total; else systBck=0.0; 
		systBck=0.0;
                gsyst->SetPointEYlow (i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
                gsyst->SetPointEYhigh(i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
	      }
              gsyst->SetFillColor(12);
              gsyst->SetFillStyle(3345);
              gsyst->SetMarkerSize(0);
              gsyst->SetLineWidth(0);
              gsyst->SetLineColor(kWhite);
	      gsyst->Draw("E2same");
              //TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
              //setex1->Draw();
	    }

            if(_hist[iHiggs ] && isHWWOverlaid == false) _hist[iHiggs ]->Draw("hist,same");
            if(_hist[iHiggs2]) _hist[iHiggs2]->Draw("hist,same");

            if(_data) {
	      bool plotCorrectErrorBars = true;
	      if(plotCorrectErrorBars == true) {
  		TGraphAsymmErrors * g = new TGraphAsymmErrors(_data);
  		for (int i = 0; i < g->GetN(); ++i) {
                  double N = g->GetY()[i];
                  double alpha=(1-0.6827);
                  double L = (N==0) ? 0 : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
                  double U = ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
                  g->SetPointEYlow(i,double(N)-L);
                  g->SetPointEYhigh(i, U-double(N));
                  
                  g->SetPointEXlow (i, _data->GetBinWidth(i+1)/2.);
                  g->SetPointEXhigh(i, _data->GetBinWidth(i+1)/2.);
  		}
                if(doApplyBinWidth == true){
 		  for (int i = 0; i < g->GetN(); ++i) {
                    g->SetPointEYlow (i, g->GetErrorYlow(i) /_data->GetBinWidth(i+1)*SFBinWidth);
                    g->SetPointEYhigh(i, g->GetErrorYhigh(i)/_data->GetBinWidth(i+1)*SFBinWidth);
                    g->SetPoint      (i, g->GetX()[i], g->GetY()[i]/_data->GetBinWidth(i+1)*SFBinWidth);
  		  }
                }
  		g->Draw("ep");
	      }
	      else {
                if(doApplyBinWidth == true){
                  for(int nb=1; nb<=_data->GetNbinsX(); nb++){
		    _data->SetBinContent(nb, _data->GetBinContent(nb)/_data->GetBinWidth(nb)*SFBinWidth);
		    _data->SetBinError  (nb, _data->GetBinError(nb)  /_data->GetBinWidth(nb)*SFBinWidth);
                  }
                }
	        _data->Draw("ep,same");
	      }
            }
	    
            hstack->SetTitle("");
            //hstack->SetTitle("CMS");

   	    //TPaveText *pt = new TPaveText(0.61,0.8337762,0.9408059,0.8862238,"blNDC");
   	    //pt->SetName("title");
   	    //pt->SetBorderSize(0);
   	    //pt->SetFillColor(10);
   	    //pt->SetTextFont(42);
   	    //pt->SetTextSize(_tsize);
   	    //pt->AddText("CMS preliminary");
   	    //pt->Draw();
   
            Float_t theMax = hstack->GetMaximum();
            Float_t theMin = hstack->GetMinimum();

            if (_hist[iHiggs]) {
                if (_hist[iHiggs]->GetMaximum() > theMax) theMax = _hist[iHiggs]->GetMaximum();
                if (_hist[iHiggs]->GetMinimum() < theMin) theMin = _hist[iHiggs]->GetMinimum();
            }

            if (_data) {

                Float_t dataMax = GetMaximumIncludingErrors(_data);
                if (dataMax > theMax) theMax = dataMax;

                Float_t dataMin = GetMinimumIncludingErrors(_data);
                if (dataMin < theMin) theMin = dataMin;
            }

            if (gPad->GetLogy()) {
            	hstack->SetMaximum(250 * theMax);
            	hstack->SetMinimum(TMath::Max(0.9 * theMin,0.010));
            } else {
              hstack->SetMaximum(1.75 * theMax);
            }

            if(_breakdown) {
                THStackAxisFonts(hstack, "y", "Events");
                hstack->GetHistogram()->LabelsOption("v");
            } else {
                if(_units.Sizeof() == 1) {
                    THStackAxisFonts(hstack, "x", _xLabel.Data());
                    if     (doApplyBinWidth == true && SFBinWidth == 1) THStackAxisFonts(hstack, "y", "Events / GeV");
                    else if(doApplyBinWidth == true)                    THStackAxisFonts(hstack, "y", Form("Events / %.2f",SFBinWidth));
                    else                                                THStackAxisFonts(hstack, "y", "Events");
                } else {
		    if(_units.EndsWith("BIN") == false){
                      THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                      if(hSum->GetBinWidth(0) < 1) THStackAxisFonts(hstack, "y", TString::Format("Events / %.1f %s",hSum->GetBinWidth(0),_units.Data()));
		      else                         THStackAxisFonts(hstack, "y", TString::Format("Events / %.0f %s",hSum->GetBinWidth(0),_units.Data()));
                    }
		    else {
		      _units = _units.ReplaceAll("BIN","");
                      THStackAxisFonts(hstack, "x", TString::Format("%s (%s)",_xLabel.Data(),_units.Data()));
                      THStackAxisFonts(hstack, "y", TString::Format("Events / bin"));
		    }
                }
            }

            // total mess to get it nice, should be redone
            size_t j=0;
            TString higgsLabel  = Form("%s",_HiggsLabel.Data());
            TString higgs2Label = Form("%s",_Higgs2Label.Data());

            if(_data          && _data         ->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _data,	       "Data",    "lp"); j++; }
            if(_hist[iWWEWK ] && _hist[iWWEWK ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWWEWK ], "EW WW", "f" ); j++; }
            if(_hist[iZZ    ] && _hist[iZZ    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZZ    ], "ZZ",	  "f" ); j++; }
            if(_hist[iWZ    ] && _hist[iWZ    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWZ    ], "WZ",	  "f" ); j++; }
            if(_hist[iWW    ] && _hist[iWW    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWW    ], "WW",	  "f" ); j++; }
            if(_hist[iVVV   ] && _hist[iVVV   ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVVV   ], _labelVVV.Data(),	  "f" ); j++; }
            if(_hist[iEM    ] && _hist[iEM    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iEM    ], _labelEM.Data(), "f" ); j++; }
            if(_hist[iZJets ] && _hist[iZJets ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZJets ], "Drell-Yan", "f" ); j++; }
            if(_hist[iTop   ] && _hist[iTop   ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iTop   ], "t#bar{t}","f" ); j++; }
            if(_hist[iVV    ] && _hist[iVV    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVV    ], "VV","f" ); j++; }
            //if(_hist[iWJets ] && _hist[iWJets ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWJets ], "Nonprompt",  "f" ); j++; }
            if(_hist[iWJets ] && _hist[iWJets ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWJets ], "Nonresonant",  "f" ); j++; }
            if(_hist[iWG    ] && _hist[iWG    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWG    ], "V+#gamma", "f" ); j++; }
            if(_hist[iWWQCD ] && _hist[iWWQCD ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWWQCD ], "QCD WW", "f" ); j++; }
            if(_hist[iDPS   ] && _hist[iDPS   ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iDPS   ], "DPS WW",  "f" ); j++; }
            if(_hist[iWS    ] && _hist[iWS    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWS    ], "Wrong-sign",	  "f" ); j++; }
            if(_hist[iOther ] && _hist[iOther ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iOther ], "Others", "f" ); j++; }
            if     (_hist[iHiggs   ] && isHWWOverlaid) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHiggs   ], higgsLabel, "f" ); j++; }
            else if(_hist[iHiggs   ])                  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHiggs   ], higgsLabel, "l" ); j++; }
            if(_hist[iHiggs2]) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset*1.03, _hist[iHiggs2], higgs2Label, "l" ); j++; }

            if(plotSystErrorBars == true) {DrawLegendTG(xPos[j], 0.84 - yOff[j]*_yoffset*1.07,gsyst, "Bkg. unc.",  "f" ); j++;}

            //TLatex* luminosity = new TLatex(0.9, 0.8, TString::Format("L = %.1f fb^{-1}",_lumi));
            //luminosity->SetNDC();
            //luminosity->SetTextAlign(32);
            //luminosity->SetTextFont(42);
            //luminosity->SetTextSize(_tsize);
            //luminosity->Draw("same");
            if(_extraLabel) _extraLabel->Draw("same");

            return hstack->GetHistogram();
        }
        void setLumi(const float &l) { _lumi = l; }
        void setLabel(const TString &s) { _xLabel = s; }
        void setUnits(const TString &s) { _units = s; }
        void setBreakdown(const bool &b = true) { _breakdown = b; }
        void addLabel(const std::string &s) {
            _extraLabel = new TLatex(0.35, 0.70, TString(s));
            _extraLabel->SetNDC();
            _extraLabel->SetTextAlign(32);
            _extraLabel->SetTextFont(42);
            _extraLabel->SetTextSize(0.06);
        }

    private: 
        std::vector<TH1F*> _hist;
        TH1F* _data;

        //MWL
        float    _lumi;
        TString  _xLabel;
        TString  _units;
        TLatex * _extraLabel;
        bool     _breakdown;
        TString  _HiggsLabel;
        TString  _Higgs2Label;
        TString  _labelEM;
        TString  _labelVVV;

};
