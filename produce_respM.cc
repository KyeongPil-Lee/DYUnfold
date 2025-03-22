#include "CommonTool.h"
#include "TH2D.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"

void produce_respM() {
  // -- define response matrix
  const Int_t nMllBin_gen = 36;
  Double_t arr_mllBinEdge_gen[nMllBin_gen+1] = {
    40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,
    126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000};

  // -- number of rec. bin should be (equal) or larger than numbe of gen. bins for TUnfold
  // -- add one more bin in the under & overflow region w.r.t. gen bin edges
  const Int_t nMllBin_rec = 38;
  Double_t arr_mllBinEdge_rec[nMllBin_rec+1] = {
      35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120,
      126, 133, 141, 150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 510, 600, 700, 830, 1000, 1200};

  TH1::SetDefaultSumw2();
  TH1D* h_mll_gen = new TH1D("h_mll_gen", "h_mll_gen", nMllBin_gen, arr_mllBinEdge_gen); // -- gen-level invariant mass
  TH1D* h_mll_rec = new TH1D("h_mll_rec", "h_mll_rec", nMllBin_rec, arr_mllBinEdge_rec); // -- rec-level invariant mass
  TH2D* h_respM = new TH2D("h_respM", "h_respM", nMllBin_rec, arr_mllBinEdge_rec, nMllBin_gen, arr_mllBinEdge_gen); // -- response matrix

  cout << "Reading TTree..." << endl;

  // -- reading TTree
  TChain* chain = new TChain("Events");
  // -- under KNU
  // TString path_DYMC = "/pnfs/knu.ac.kr/data/cms/store/user/sungwon/DY_Run2_UL_NanoAOD/2018/DYJetsToMuMu_M-50/DYJetsToMuMu_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/DYJetsToMuMu_M-50/241017_133720/0000";
  // // chain->Add(path_DYMC+"/*.root"); // -- use all files
  chain->Add(path_DYMC+"/tree_52.root"); // -- largest file
  chain->Add(path_DYMC+"/tree_108.root"); // -- 2nd largest file
  // chain->Add("/Users/kplee/Research/Analysis/DYUnfoldExample/tree_1.root");
  // chain->Add("/Users/kplee/Research/Analysis/DYUnfoldExample/tree_67.root");
  TTreeReader reader(chain);

  // -- generator level variables -- //
  TTreeReaderValue<Float_t> genWeight(reader, "genWeight");

  TTreeReaderArray<Float_t> GenDressedLepton_pt(reader, "GenDressedLepton_pt");
  TTreeReaderArray<Float_t> GenDressedLepton_eta(reader, "GenDressedLepton_eta");
  TTreeReaderArray<Float_t> GenDressedLepton_phi(reader, "GenDressedLepton_phi");
  TTreeReaderArray<Int_t>   GenDressedLepton_pdgId(reader, "GenDressedLepton_pdgId");
 
  // -- reconstruction level variables -- //
  // -- trigger
  TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");

  // -- reconstructed muon
  TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
  TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
  TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");

  TTreeReaderArray<Bool_t> Muon_tightId(reader, "Muon_tightId");
  TTreeReaderArray<Float_t> Muon_pfRelIso04_all(reader, "Muon_pfRelIso04_all");

  // -- correction only for reco-level events
  TTreeReaderValue<Float_t> L1PreFiringWeight_Nom(reader, "L1PreFiringWeight_Nom");

  Int_t nEvent_tot = reader.GetEntries(kTRUE);
  cout << "nEvent_tot = " << nEvent_tot << endl;

  // -- for testing; remove this line for processing all events
  // if( nEvent_tot > 1e6 ) nEvent_tot = 1e6;

  for (Int_t i_ev = 0; i_ev < nEvent_tot; ++i_ev) {
    reader.SetEntry(i_ev);

    // -- progress bar
    // if( i_ev % (nEvent_tot/100) == 0 ) {
    //   cout << "Processing " << i_ev << " / " << nEvent_tot << " events" << endl;
    // }

    Double_t eventWeight = (*genWeight) * (*L1PreFiringWeight_Nom);

    // -- check gen-level muons (dressed leptons)
    std::vector<TLorentzVector> vec_vecP_gen;
    for (Int_t i_lep = 0; i_lep < GenDressedLepton_pt.GetSize(); ++i_lep) {
      if (std::abs(GenDressedLepton_pdgId[i_lep]) == 13) {
        TLorentzVector vecP_gen;
        vecP_gen.SetPtEtaPhiM(GenDressedLepton_pt[i_lep], GenDressedLepton_eta[i_lep], GenDressedLepton_phi[i_lep], CommonTool::mass_mu);
        vec_vecP_gen.push_back(vecP_gen);
      }
    }

    Bool_t hasDY_gen = kFALSE;
    Double_t mll_gen = -1.0;
    if( vec_vecP_gen.size() == 2 && 
        CommonTool::Pass_Acceptance(vec_vecP_gen) ) {
      hasDY_gen = kTRUE;
      mll_gen = (vec_vecP_gen[0] + vec_vecP_gen[1]).M();
    }

    // -- check reconstruction level -- //

    // -- check muon pair -- //
    std::vector<TLorentzVector> vec_vecP_rec;
    for (Int_t i_lep = 0; i_lep < Muon_pt.GetSize(); ++i_lep) {
      if (Muon_tightId[i_lep] && Muon_pfRelIso04_all[i_lep] < 0.15 && 
          Muon_pt[i_lep] > 15.0 && std::abs(Muon_eta[i_lep]) < 2.4) {
        TLorentzVector vecP_rec;
        vecP_rec.SetPtEtaPhiM(Muon_pt[i_lep], Muon_eta[i_lep], Muon_phi[i_lep], CommonTool::mass_mu);
        vec_vecP_rec.push_back(vecP_rec);
      }
    }

    Bool_t hasDY_rec = kFALSE;
    Double_t mll_rec = -1.0;
    if( *HLT_IsoMu24 &&
        vec_vecP_rec.size() == 2 && // -- could be improved e.g. select the best muon pair
        CommonTool::Pass_Acceptance(vec_vecP_rec) ) {
      hasDY_rec = kTRUE;
      mll_rec = (vec_vecP_rec[0] + vec_vecP_rec[1]).M();
    }

    // -- fill mass histograms
    if( hasDY_gen ) h_mll_gen->Fill(mll_gen, *genWeight); // -- only apply gen-level weight
    if( hasDY_rec ) h_mll_rec->Fill(mll_rec, eventWeight);

    // -- fill response matrix
    if( hasDY_gen && hasDY_rec ) { // -- if both gen- and rec-level pairs are found
      h_respM->Fill(mll_rec, mll_gen, eventWeight);
      // -- one more fill rec-underflow bin with (gen-weight - tot-weight)
      // -- to have correct sum(weight) when 2D histogram is projected to gen-level axis
      // https://www.desy.de/~sschmitt/TUnfold/tunfold_manual_v17.9.pdf, page 10
      h_respM->Fill(-1.0,    mll_gen, *genWeight-eventWeight);
    }
    else {
      if( hasDY_gen )  { // -- if only gen-level pair is found
        h_respM->Fill(-1.0,    mll_gen, *genWeight); // -- rec: fill the underflow bin
      }
      if( hasDY_rec )  { // -- if only rec-level pair is found
        h_respM->Fill(mll_rec, -1.0, eventWeight); // -- gen: fill the underflow bin
        h_respM->Fill(-1.0,    -1.0, *genWeight-eventWeight);
      }
    }    
  } // -- end of event loop

  TFile* f_out = TFile::Open("respM_DYMC.root", "RECREATE");
  f_out->cd();
  h_mll_gen->Write();
  h_mll_rec->Write();
  h_respM->Write();
  f_out->Close();

  cout << "Done!" << endl;
}