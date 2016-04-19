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
enum samp { iWW, iZJets, iTop, iVV, iWJets, iHiggs, iWZ, iZZ, iWG, iVVV, iEM, nSamples };

float xPos[nSamples+1] = {0.19,0.19,0.19,0.19,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41}; 
float yOff[nSamples+1] = {   0,   1,   2,   3,   0,   1,   2,   3,   4,   5,   6,   7};

const Float_t _tsize   = 0.035;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;


//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1F* h)
{
    Float_t maxWithErrors = 0;

    for (Int_t i=1; i<=h->GetNbinsX(); i++) {

        Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

        if (binHeight > maxWithErrors) maxWithErrors = binHeight;
    }

    return maxWithErrors;
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
    axis->SetTitleOffset(  1.5);
    axis->SetTitleSize  (0.050);

    if (coordinate == "y") axis->SetTitleOffset(1.6);

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

class StandardPlot2015 {

    public: 
        StandardPlot2015() { _hist.resize(nSamples,0); _data = 0; _breakdown = false; _HiggsLabel = "";_labelEM = " Non-prompt";}
        void setMCHist   (const samp &s, TH1F * h)  { _hist[s]       = h;  } 
        void setDataHist  (TH1F * h) { _data	      = h;  }
        void setWWHist    (TH1F * h) { setMCHist(iWW   ,h); }
        void setZJetsHist (TH1F * h) { setMCHist(iZJets,h); }
        void setTopHist   (TH1F * h) { setMCHist(iTop  ,h); }
        void setVVHist    (TH1F * h) { setMCHist(iVV   ,h); }
        void setWJetsHist (TH1F * h) { setMCHist(iWJets,h); }
        void setWGHist    (TH1F * h) { setMCHist(iWG   ,h); }
        void setHiggsHist (TH1F * h) { setMCHist(iHiggs,h); }
        void setWZHist    (TH1F * h) { setMCHist(iWZ   ,h); }
        void setZZHist    (TH1F * h) { setMCHist(iZZ   ,h); }
        void setVVVHist   (TH1F * h) { setMCHist(iVVV  ,h); }
        void setEMHist    (TH1F * h) { setMCHist(iEM   ,h); }
        void setOverlaid  (bool b)   { isHWWOverlaid = b;   }
        void setLabelEM   (TString s){ _labelEM = s.Data();}

  TH1F* getDataHist() { return _data; }

        void setHiggsLabel(TString s) {_HiggsLabel=s;}

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

        TH1* Draw(const int &rebin=1) {

            Color_t _sampleColor[nSamples];
            _sampleColor[iWW    ] = kAzure-9;
            _sampleColor[iZJets ] = 901;//kGreen+2;
            _sampleColor[iTop   ] = kYellow;
            _sampleColor[iVV    ] = kAzure-2;
            _sampleColor[iWJets ] = kGray+1;
            _sampleColor[iWG    ] = kViolet-9;
            _sampleColor[iHiggs ] = 809;//kOrange+7;
            _sampleColor[iWZ    ] = 856;//kAzure-2;
            _sampleColor[iZZ    ] = 842;//kBlue-3;
            _sampleColor[iVVV   ] = 832;//kRed+1;
            _sampleColor[iEM    ] = 798;//kYellow;

            //setUpStyle();
            //if(!gPad) new TCanvas();

            THStack* hstack = new THStack();
	    TH1D* hSum;
	    if     (_hist[iZJets]) hSum = (TH1D*)_hist[iZJets]->Clone();
	    else if(_hist[iWZ])    hSum = (TH1D*)_hist[iWZ]   ->Clone();
	    else if(_hist[iVV])    hSum = (TH1D*)_hist[iVV]   ->Clone();
	    else if(_data)         hSum = (TH1D*)_data        ->Clone();
	    else                   {printf("PROBLEM\n");assert(0);}
	    hSum->Rebin(rebin);
	    hSum->Scale(0.0);
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

                // signal gets overlaid
                if (i == iHiggs && isHWWOverlaid == false) continue;

                _hist[i]->SetFillColor(_sampleColor[i]);
                _hist[i]->SetFillStyle(1001);

                hstack->Add(_hist[i]);
		hSum->Add(_hist[i]);
            }

            if(_hist[iHiggs]) _hist[iHiggs]->SetLineWidth(3);
            if(_data) _data->Rebin(rebin);
	    //_data->SetBinContent(1,_data->GetBinContent(1)-3);
	    //_data->SetBinContent(3,_data->GetBinContent(3)+1);
	    //_data->SetBinContent(5,_data->GetBinContent(5)+2);
            if(_data) _data->SetLineColor  (kBlack);
            if(_data) _data->SetMarkerStyle(kFullCircle);
	    hstack->Draw("hist");

	    bool plotSystErrorBars = true;
	    if(plotSystErrorBars == true) {
  	      TGraphAsymmErrors * gsyst = new TGraphAsymmErrors(hSum);
              for (int i = 0; i < gsyst->GetN(); ++i) {
                double systBck = 0;
		if(_hist[iWW	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iWW    ]->GetBinContent(i+1),2);
		if(_hist[iZJets ]) systBck = systBck + 1.000*TMath::Power(0.500*_hist[iZJets ]->GetBinContent(i+1),2);
		if(_hist[iTop	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iTop   ]->GetBinContent(i+1),2);
		if(_hist[iVV	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iVV    ]->GetBinContent(i+1),2);
		if(_hist[iWJets ]) systBck = systBck + 1.000*TMath::Power(0.360*_hist[iWJets ]->GetBinContent(i+1),2);
		if(_hist[iWG	]) systBck = systBck + 1.000*TMath::Power(0.300*_hist[iWG    ]->GetBinContent(i+1),2);
		if(_hist[iWZ	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iWZ    ]->GetBinContent(i+1),2);
		if(_hist[iZZ	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iZZ    ]->GetBinContent(i+1),2);
		if(_hist[iVVV	]) systBck = systBck + 1.000*TMath::Power(0.100*_hist[iVVV   ]->GetBinContent(i+1),2);
		if(_hist[iEM	]) systBck = systBck + 1.000*TMath::Power(0.400*_hist[iEM    ]->GetBinContent(i+1),2);
		if(_hist[iHiggs	]) systBck = systBck + 1.000*TMath::Power(0.200*_hist[iHiggs ]->GetBinContent(i+1),2);
                double total = hSum->GetBinContent(i+1);
                if(total > 0) systBck = sqrt(systBck)/total; else systBck=0.0; 
		systBck=0.0;
                gsyst->SetPointEYlow (i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
                gsyst->SetPointEYhigh(i, sqrt(hSum->GetBinError(i+1)*hSum->GetBinError(i+1)+hSum->GetBinContent(i+1)*hSum->GetBinContent(i+1)*systBck*systBck));
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

            if(_hist[iHiggs] && isHWWOverlaid == false) _hist[iHiggs]->Draw("hist,same");

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
                  
                  g->SetPointEXlow (i, 0);
                  g->SetPointEXhigh(i, 0);
  		}
  		g->Draw("P");
	      }
	      else {
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
            }

            if (gPad->GetLogy()) {
            	hstack->SetMaximum(2000 * theMax);
            	hstack->SetMinimum(TMath::Max(0.9 * theMin,0.50));
            } else {
              hstack->SetMaximum(1.75 * theMax);
            }

            if(_breakdown) {
                THStackAxisFonts(hstack, "y", "Events");
                hstack->GetHistogram()->LabelsOption("v");
            } else {
                THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                if(_units.Sizeof() == 1) {
                    THStackAxisFonts(hstack, "x", _xLabel.Data());
                    THStackAxisFonts(hstack, "y", "Events");
                } else {
                    THStackAxisFonts(hstack, "x", TString::Format("%s [%s]",_xLabel.Data(),_units.Data()));
                    if(hSum->GetBinWidth(0) < 1) THStackAxisFonts(hstack, "y", TString::Format("Events / %.1f %s",hSum->GetBinWidth(0),_units.Data()));
		    else                         THStackAxisFonts(hstack, "y", TString::Format("Events / %.0f %s",hSum->GetBinWidth(0),_units.Data()));
                }
            }

            // total mess to get it nice, should be redone
            size_t j=0;
            TString higgsLabel = Form("%s",_HiggsLabel.Data());

            if(_data          && _data         ->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _data,	       "data",    "lp"); j++; }
            if(_hist[iWW    ] && _hist[iWW    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWW    ], "WW",	  "f" ); j++; }
            if(_hist[iZJets ] && _hist[iZJets ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZJets ], "Z+jets/#gamma", "f" ); j++; }
            if(_hist[iTop   ] && _hist[iTop   ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iTop   ], "t#bar{t}","f" ); j++; }
            if(_hist[iVV    ] && _hist[iVV    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVV    ], "VV/VVV","f" ); j++; }
            if(_hist[iWJets ] && _hist[iWJets ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWJets ], "W+jets",  "f" ); j++; }
            if(_hist[iWG    ] && _hist[iWG    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWG    ], "V+#gamma^{(*)}", "f" ); j++; }
            if(_hist[iWZ    ] && _hist[iWZ    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iWZ    ], "WZ",	  "f" ); j++; }
            if(_hist[iZZ    ] && _hist[iZZ    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZZ    ], "ZZ",	  "f" ); j++; }
            if(_hist[iVVV   ] && _hist[iVVV   ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVVV   ], "VVV",	  "f" ); j++; }
            if(_hist[iEM    ] && _hist[iEM    ]->GetSumOfWeights() > 0) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iEM    ], _labelEM.Data(), "f" ); j++; }
            if     (_hist[iHiggs   ] && isHWWOverlaid) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHiggs   ], higgsLabel, "f" ); j++; }
            else if(_hist[iHiggs   ])                  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iHiggs   ], higgsLabel, "l" ); j++; }

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
            _extraLabel = new TLatex(0.9, 0.74, TString(s));
            _extraLabel->SetNDC();
            _extraLabel->SetTextAlign(32);
            _extraLabel->SetTextFont(42);
            _extraLabel->SetTextSize(_tsize);
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
        TString  _labelEM;

};
