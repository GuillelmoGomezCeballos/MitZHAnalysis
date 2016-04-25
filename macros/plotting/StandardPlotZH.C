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

const Bool_t isZHOverlaid = false;
enum samp { iEM, iVV, iZx, iZH, iVg, iZeg, iZjg, nSamples };

float xPos[nSamples+1] = {0.19,0.19,0.19,0.41,0.41}; 
float yOff[nSamples+1] = {0,1,2,0,1};

const Float_t _tsize   = 0.033;
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


class StandardPlotZH {

    public: 
        StandardPlotZH() { _hist.resize(nSamples,0); _data = 0; _breakdown = false; _massH = 0; _massChi = 0;}
        void setMCHist   (const samp &s, TH1F * h)  { _hist[s]       = h;  } 
        void setDataHist  (TH1F * h)                 { _data          = h;  } 
        void setEMHist    (TH1F * h)                 { setMCHist(iEM   ,h); } 
        void setVVHist    (TH1F * h)                 { setMCHist(iVV   ,h); } 
        void setZxHist    (TH1F * h)                 { setMCHist(iZx   ,h); } 
        void setZegHist   (TH1F * h)		     { setMCHist(iZeg  ,h); } 
        void setZjgHist   (TH1F * h)		     { setMCHist(iZjg  ,h); } 
        void setZHHist    (TH1F * h)                 { setMCHist(iZH   ,h); } 
        void setVgHist    (TH1F * h)                 { setMCHist(iVg   ,h); } 

  TH1F* getDataHist() { return _data; }

        void setMassH  (const int &m) {_massH  =m;}
        void setMassChi(const int &m) {_massChi=m;}

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
            _sampleColor[iEM    ] = kGreen+2;
            _sampleColor[iVV    ] = kAzure-9;
            _sampleColor[iZx    ] = kAzure-2;
            _sampleColor[iZH    ] = kRed+1;
            _sampleColor[iZeg	] = kCyan;
            _sampleColor[iZjg	] = kAzure-9;
            _sampleColor[iVg    ] = kAzure-2;

            //setUpStyle();
            //if(!gPad) new TCanvas();

            THStack* hstack = new THStack();
	    TH1D* hSum = (TH1D*)_data->Clone();
	    hSum->Rebin(rebin);
	    hSum->Scale(0.0);
            for (int i=0; i<nSamples; i++) {

                // in case the user doesn't set it
                if( !_hist[i] ) continue;
                _hist[i]->Rebin(rebin);
                _hist[i]->SetLineColor(_sampleColor[i]);

                // signal gets overlaid
                if (i == iZH && isZHOverlaid == false) continue;

                _hist[i]->SetFillColor(_sampleColor[i]);
                _hist[i]->SetFillStyle(1001);

                hstack->Add(_hist[i]);
		hSum->Add(_hist[i]);
            }

            if(_hist[iZH]) _hist[iZH]->SetLineWidth(3);
            if(_data) _data->Rebin(rebin);
            if(_data) _data->SetLineColor  (kBlack);
            if(_data) _data->SetMarkerStyle(kFullCircle);
	    hstack->Draw("hist");

	    bool plotSystErrorBars = true;
	    if(plotSystErrorBars == true) {
  	      TGraphAsymmErrors * gsyst = new TGraphAsymmErrors(hSum);
              for (int i = 0; i < gsyst->GetN(); ++i) {
                double systBck = 0;
		if(_hist[iZx	]) systBck = systBck + 0.500*0.500*_hist[iZx    ]->GetBinContent(i+1)*_hist[iZx	   ]->GetBinContent(i+1);
		if(_hist[iVV	]) systBck = systBck + 0.200*0.200*_hist[iVV    ]->GetBinContent(i+1)*_hist[iVV	   ]->GetBinContent(i+1);
		if(_hist[iEM    ]) systBck = systBck + 0.750*0.750*_hist[iEM    ]->GetBinContent(i+1)*_hist[iEM    ]->GetBinContent(i+1);
                double total = hSum->GetBinContent(i+1);
                if(total > 0) systBck = sqrt(systBck)/total;
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

            if(_hist[iZH] && isZHOverlaid == false) _hist[iZH]->Draw("hist,same");

            if(_data) {
	      bool plotCorrectErrorBars = true;
	      if(plotCorrectErrorBars == true) {
  		TGraphAsymmErrors * g = new TGraphAsymmErrors(_data);
                double alpha=(1-0.6827);
                for (int i = 0; i < g->GetN(); ++i) {
                  int N = g->GetY()[i];
                  double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
                  double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
                  g->SetPointEYlow(i, double(N)-L);
                  g->SetPointEYhigh(i, U-double(N));
                }
  		g->Draw("P");
	      }
	      else {
	        _data->Draw("ep,same");
	      }
            }
	    
            hstack->SetTitle("");
   
            Float_t theMax = hstack->GetMaximum();
            Float_t theMin = hstack->GetMinimum();

            if (_hist[iZH]) {
                if (_hist[iZH]->GetMaximum() > theMax) theMax = _hist[iZH]->GetMaximum();
                if (_hist[iZH]->GetMinimum() < theMin) theMin = _hist[iZH]->GetMinimum();
            }

            if (_data) {

                Float_t dataMax = GetMaximumIncludingErrors(_data);

                if (dataMax > theMax) theMax = dataMax;
            }

            if (gPad->GetLogy()) {
            	hstack->SetMaximum(18 * theMax);
            	hstack->SetMinimum(TMath::Max(0.9 * theMin,0.50));
            } else {
              hstack->SetMaximum(1.55 * theMax);
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
                    if(_data->GetBinWidth(0) < 1) THStackAxisFonts(hstack, "y", TString::Format("Events / %.1f %s", _data->GetBinWidth(0),_units.Data()));
		    else                          THStackAxisFonts(hstack, "y", TString::Format("Events / %.0f %s", _data->GetBinWidth(0),_units.Data()));
                }
            }

            // total mess to get it nice, should be redone
            size_t j=0;
            TString higgsLabel = " ZH";
            if(_massH != 0 && _massChi != 0) higgsLabel.Form(" ZH(%d,%d)",_massH,_massChi);

            if(_data         ) { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _data,          " Data",    "lp"); j++; }
            if     (_hist[iZH   ] && isZHOverlaid)  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZH  ], higgsLabel,	    "f" ); j++; }
            else if(_hist[iZH   ])                  { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZH  ], higgsLabel,	    "l" ); j++; }
            if(_hist[iEM    ])                      { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iEM  ], " WW + top-quark",   "f" ); j++; }
            if(_hist[iZx    ])                      { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZx  ], " Z#gamma + Zjets", "f" ); j++; }
            if(_hist[iVV    ])                      { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVV  ], " WZ + ZZ",	    "f" ); j++; }
            if(_hist[iZeg   ])  		    { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZeg ], " e #rightarrow #gamma", "f" ); j++; }
            if(_hist[iZjg   ])  		    { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iZjg ], " j #rightarrow #gamma", "f" ); j++; }
            if(_hist[iVg    ])  		    { DrawLegend(xPos[j], 0.84 - yOff[j]*_yoffset, _hist[iVg  ], " V#gamma", "f" ); j++; }

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
        int      _massH;
        int      _massChi;


};
