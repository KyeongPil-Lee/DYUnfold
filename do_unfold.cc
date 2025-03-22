#include "CommonTool.h"
#include "TUnfold.h"
#include "SimplePlotTools.h" // -- for plotting

void do_unfold() {
  TFile* f_input = TFile::Open("respM_DYMC.root");
  TH1D* h_mll_gen = (TH1D *)f_input->Get("h_mll_gen");
  TH1D* h_mll_rec = (TH1D *)f_input->Get("h_mll_rec");
  TH2D* h_respM = (TH2D *)f_input->Get("h_respM");

  // -- define TUnfold object
  TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;

  TUnfoldDensity unfold(h_respM, TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);

  // -- DY "fake" events: true event doesn't exist, but reco. event exist
  // ---- e.g. true event is out of acceptance but reco. event is within acceptance
  // -- this type of events are considered as background: need to be subtracted in the unfolding procedure
  TH1D* h_DYFake = CommonTool::MakeHist_Underflow("true", h_respM);

  double backScale = 1.0;      // -- background scale factor (default: 1.0: no changes)
  double backScaleError = 0.0; // -- error on the background scale
  unfold.SubtractBackground(h_DYFake, "DYFake", backScale, backScaleError);

  // -- do unfolding
  double tau = 0; // -- regularization strength (0: no reg.)
  unfold.DoUnfold(tau, h_mll_rec); // -- unfolding reco mll
  TH1D* h_unfolded = (TH1D *)unfold.GetOutput("h_unfolded");
  // -- get covariance matrix
  TH2 *h_covM = unfold.GetEmatrixTotal("h_covM");
  // -- assign error to each bin using the covariance matrix
  for (Int_t i = 0; i < h_unfolded->GetNbinsX(); i++) {
    Int_t i_bin = i + 1;
    Double_t error = std::sqrt(h_covM->GetBinContent(i_bin, i_bin)); // -- diagonal term: unc^2
    h_unfolded->SetBinError(i_bin, error);
  }

  // -- plotting
  TH1D* h_mll_rec_rebin = (TH1D*)h_mll_gen->Clone(); // -- to have same binning with gen for comparison plots
  h_mll_rec_rebin->Reset("ICES");
  for(Int_t i=0; i<h_mll_rec_rebin->GetNbinsX(); ++i) {
    Int_t i_bin = i + 1;
    Double_t lowerEdge = h_mll_rec_rebin->GetBinLowEdge(i_bin);
    Double_t upperEdge = h_mll_rec_rebin->GetBinLowEdge(i_bin+1);

    for(Int_t j=0; j<h_mll_rec->GetNbinsX(); ++j) {
      Int_t j_bin = j + 1;
      Double_t lowerEdge_rec = h_mll_rec->GetBinLowEdge(j_bin);
      Double_t upperEdge_rec = h_mll_rec->GetBinLowEdge(j_bin+1);
      if (lowerEdge == lowerEdge_rec && upperEdge == upperEdge_rec) {
        h_mll_rec_rebin->SetBinContent(i_bin, h_mll_rec->GetBinContent(j_bin));
        h_mll_rec_rebin->SetBinError(i_bin, h_mll_rec->GetBinError(j_bin));
        break;
      }
    } // -- end of j loop
  } // -- end of i loop

  PlotTool::HistCanvaswRatio *canvas = new PlotTool::HistCanvaswRatio("c_closureTest", 1, 1);
  canvas->SetTitle("m [GeV]", "# events", "ratio over true");

  canvas->Register(h_mll_gen,       "Gen. distribution", kBlack);
  canvas->Register(h_mll_rec_rebin, "Rec. distribution", kBlue);
  canvas->Register(h_unfolded,      "Unfolded distribution", kRed);

  canvas->SetLegendPosition(0.50, 0.73, 0.95, 0.95);

  canvas->SetRangeX(40, 1000);
  // canvas->SetRangeY(minY, maxY);
  // canvas->SetRangeRatio(minRatio, maxRatio);
  canvas->SetAutoRangeY();
  canvas->SetAutoRangeRatio();

  canvas->Latex_CMSInternal();
  // canvas->RegisterLatex(0.16, 0.91, "#font[42]{#scale[0.6]{Gaussian distributions}}"); // same with above

  canvas->Draw();

  canvas->SetCanvasName("c_closureTest_zoomIn");
  canvas->SetRangeRatio(0.93, 1.07);
  canvas->Draw();

  // cout << "h_mll_gen" << endl;
  // PlotTool::Print_Histogram(h_mll_gen);

  // cout << "h_mll_rec" << endl;
  // PlotTool::Print_Histogram(h_mll_rec);
}