namespace CommonTool {
  Double_t mass_mu = 105.6583745*1e-3; // -- PDG value

  Bool_t Pass_Acceptance(std::vector<TLorentzVector>& vec_vecP_gen) {
    // -- check if the two muons are in the acceptance
    if( vec_vecP_gen[0].Pt() < 20 || vec_vecP_gen[1].Pt() < 15 ) return kFALSE;
    if( std::abs(vec_vecP_gen[0].Eta()) > 2.4 || std::abs(vec_vecP_gen[1].Eta()) > 2.4 ) return kFALSE;

    return kTRUE;
  }
}