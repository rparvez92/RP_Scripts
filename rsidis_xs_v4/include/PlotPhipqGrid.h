#ifndef RSIDIS_XS_V4_PLOT_PHIPQ_GRID_H
#define RSIDIS_XS_V4_PLOT_PHIPQ_GRID_H

#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLatex.h"
#include "TStyle.h"

namespace rsidis_xs_v4 {

struct ZSeries {
  double zLo;
  double zHi;
  std::unique_ptr<TGraphErrors> gr;
  std::string label;
};

inline std::unique_ptr<TGraphErrors> HistToGraphErrors(const TH1* h, bool useBinWidthAsXErr=false) {
  if(!h) return nullptr;
  const int n = h->GetNbinsX();
  auto gr = std::make_unique<TGraphErrors>(n);
  for(int i=1;i<=n;i++){
    const double x  = h->GetBinCenter(i);
    const double y  = h->GetBinContent(i);
    const double ex = useBinWidthAsXErr ? 0.5*h->GetBinWidth(i) : 0.0;
    const double ey = h->GetBinError(i);
    gr->SetPoint(i-1, x, y);
    gr->SetPointError(i-1, ex, ey);
  }
  return gr;
}

inline std::unique_ptr<TCanvas> MakePhipqGridCanvas(
    const std::vector<double>& ptEdges,
    const std::vector< std::vector<ZSeries> >& seriesByPt,
    const std::string& canvasName,
    const std::string& yTitle,
    const std::string& mainTitle)
{
  const int nPt = (int)ptEdges.size()-1;
  const int nCols = 2;
  const int nRows = (nPt + nCols - 1)/nCols;

  auto c = std::make_unique<TCanvas>(canvasName.c_str(), mainTitle.c_str(), 1200, 900);
  c->Divide(nCols, nRows, 0.001, 0.001);

  gStyle->SetOptStat(0);

  for(int ipt=0; ipt<nPt; ipt++){
    c->cd(ipt+1);
    gPad->SetTicks(1,1);
    gPad->SetMargin(0.12, 0.04, 0.12, 0.08);

    // Determine y range
    double yMax = 0.0;
    for(const auto& s : seriesByPt[ipt]){
      if(!s.gr) continue;
      for(int i=0;i<s.gr->GetN();i++){
        double x,y;
        s.gr->GetPoint(i,x,y);
        yMax = std::max(yMax, y + 1.2*s.gr->GetErrorY(i));
      }
    }
    if(yMax<=0) yMax = 1.0;

    auto frame = gPad->DrawFrame(0.0, 0.0, 2.0*M_PI, 1.15*yMax);
    frame->GetXaxis()->SetTitle("#phi_{pq} (rad)");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.05);

    // Legend
    auto leg = new TLegend(0.15, 0.72, 0.55, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    int colorIdx = 1;
    for(const auto& s : seriesByPt[ipt]){
      if(!s.gr) continue;
      s.gr->SetMarkerStyle(20 + (colorIdx%5));
      s.gr->SetMarkerSize(1.0);
      s.gr->SetLineWidth(2);
      s.gr->SetLineColor(colorIdx+1);
      s.gr->SetMarkerColor(colorIdx+1);
      s.gr->Draw("P SAME");
      leg->AddEntry(s.gr.get(), s.label.c_str(), "P");
      colorIdx++;
    }
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.05);
    lat.DrawLatex(0.15, 0.64, Form("p_{T} [%.2f, %.2f] GeV", ptEdges[ipt], ptEdges[ipt+1]));
  }

  // Title on pad 1
  c->cd(1);
  TLatex title;
  title.SetNDC();
  title.SetTextSize(0.05);
  title.DrawLatex(0.12, 0.96, mainTitle.c_str());

  return c;
}

} // namespace rsidis_xs_v4

#endif
