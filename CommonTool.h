namespace CommonTool {
  Double_t mass_mu = 105.6583745*1e-3; // -- PDG value

  Bool_t Pass_Acceptance(std::vector<TLorentzVector>& vec_vecP_gen) {
    // -- check if the two muons are in the acceptance
    Double_t pt_lead = std::max(vec_vecP_gen[0].Pt(), vec_vecP_gen[1].Pt());
    Double_t pt_sub = std::min(vec_vecP_gen[0].Pt(), vec_vecP_gen[1].Pt());
    if( pt_lead < 20 || pt_sub < 15 ) return kFALSE;
    if( std::abs(vec_vecP_gen[0].Eta()) > 2.4 || std::abs(vec_vecP_gen[1].Eta()) > 2.4 ) return kFALSE;

    return kTRUE;
  }

  TH1D *MakeHist_Underflow(TString type, TH2D *h2D)
  {
    TH1D *h_return = nullptr;
    if (type == "reco")
      h_return = h2D->ProjectionY();
    else if (type == "true")
      h_return = h2D->ProjectionX();
    else
      throw std::invalid_argument("[MakeHist_Underflow] type = " + type + " is not supported");

    h_return->Reset("ICES");

    Int_t nBin = (type == "reco") ? h2D->GetNbinsY() : h2D->GetNbinsX();

    for (Int_t i = 0; i < nBin; ++i) {
      Int_t i_bin = i + 1;

      Double_t value = (type == "reco") ? h2D->GetBinContent(0, i_bin) : h2D->GetBinContent(i_bin, 0);
      Double_t error = (type == "reco") ? h2D->GetBinError(0, i_bin) : h2D->GetBinError(i_bin, 0);

      h_return->SetBinContent(i_bin, value);
      h_return->SetBinError(i_bin, error);
    }

    return h_return;
  }
}