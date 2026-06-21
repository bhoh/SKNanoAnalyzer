#include "AtobbMLTree.h"

AtobbMLTree::AtobbMLTree() :
  const_top_mass(172.5),
  const_top_width(1.5),
  const_w_mass(80.4),
  const_w_width(2.085)
{
}

AtobbMLTree::~AtobbMLTree() {}

void AtobbMLTree::initializeAnalyzer() {

  MuonIDs = {
    Muon::MuonID::POG_TIGHT,
    Muon::MuonID::POG_MEDIUM_PROMPT,
    Muon::MuonID::POG_MVA_MU_TIGHT
  };

  MuonIDISOSFKeys = {
    "NUM_TightID_DEN_TrackerMuons",
    "NUM_TightPFIso_DEN_TightID",
    "NUM_MediumPromptID_DEN_TrackerMuons",
    "NUM_TightPFIso_DEN_MediumPromptID",
    "NUM_TightMvaMuID_DEN_TrackerMuons"
  };

  if (DataEra == "2016preVFP" || DataEra == "2016postVFP" || DataEra == "2018") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  }
  else if (DataEra == "2017") {
    IsoMuTriggerName = "HLT_IsoMu27";
    TriggerSafePtCut = 29.;
  }
  else if (DataEra == "2022" || DataEra == "2022EE") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  }
  else if (DataEra == "2023") {
    IsoMuTriggerName = "";
    TriggerSafePtCut = 26.;
  }
  else if (DataEra == "2024") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  }
  else {
    cerr << "[AtobbMLTree::initializeAnalyzer] DataEra is not set properly" << endl;
    exit(EXIT_FAILURE);
  }

  cout << "[AtobbMLTree::initializeAnalyzer] IsoMuTriggerName = "
       << IsoMuTriggerName << endl;
  cout << "[AtobbMLTree::initializeAnalyzer] TriggerSafePtCut = "
       << TriggerSafePtCut << endl;

  myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
  myCorr->SetTaggingParam(
      JetTagging::JetFlavTagger::ParT,
      JetTagging::JetFlavTaggerWP::Medium
  );

  string SKNANO_HOME = getenv("SKNANO_HOME");
  if (IsDATA) {
    systHelper = std::make_unique<SystematicHelper>(
        SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra
    );
  }
  else {
    systHelper = std::make_unique<SystematicHelper>(
        SKNANO_HOME + "/docs/noSyst.yaml", MCSample, DataEra
    );
  }

  NewTree("Training_Tree", {""}, {"*"});
}

void AtobbMLTree::executeEvent() {

  AllMuonViews = GetAllMuonViews();
  AllElectronViews = GetAllElectronViews();
  AllJetViews = GetAllJetViews();
  AllGenViews = GetAllGenViews();
  ev = GetEvent();

  if (!PassMetFilter(AllJetViews, ev)) return;
  // trigger cut is not necessary for MC only study
  //if (IsoMuTriggerName != "" && !(ev.PassTrigger(IsoMuTriggerName))) return;

  Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI);

  bool use_pog_tight_muon_id = true;
  bool use_pog_mva_tight_muon_id = false;
  bool correct_b_jet_pt = true;
  bool use_UParT_JEC = false;
  bool apply_pu_id = false;

  Muon::MuonID this_muon_id = MuonIDs[0];
  TString this_muon_id_sf_key = MuonIDISOSFKeys[0];
  TString this_muon_iso_sf_key = MuonIDISOSFKeys[1];
  TString this_muon_trig_sf_key = "";

  if (use_pog_mva_tight_muon_id) {
    this_muon_id = MuonIDs[2];
    this_muon_id_sf_key = MuonIDISOSFKeys[4];
  }
  else if (use_pog_tight_muon_id) {
    this_muon_id = MuonIDs[0];
    this_muon_id_sf_key = MuonIDISOSFKeys[0];
    this_muon_iso_sf_key = MuonIDISOSFKeys[1];
    this_muon_trig_sf_key = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";
  }
  else {
    this_muon_id = MuonIDs[1];
    this_muon_id_sf_key = MuonIDISOSFKeys[2];
    this_muon_iso_sf_key = MuonIDISOSFKeys[3];
  }

  std::vector<size_t> SelectedMuonIndices_id_only =
      SelectMuonIndices(AllMuonViews, this_muon_id, 15., 2.4);

  std::vector<size_t> SelectedMuonIndices = {};
  if (use_pog_mva_tight_muon_id) {
    SelectedMuonIndices = SelectedMuonIndices_id_only;
  }
  else {
    SelectedMuonIndices = SelectMuonIndices(
        AllMuonViews,
        SelectedMuonIndices_id_only,
        Muon::MuonID::POG_PFISO_TIGHT,
        15.,
        2.4
    );
  }

  std::vector<size_t> SelectedElectronIndices =
      SelectElectronIndices(AllElectronViews, Electron::ElectronID::POG_MVAISO_WP80, 15., 2.5);

  if (SelectedMuonIndices.size() + SelectedElectronIndices.size() != 1) return;

  RVec<Muon> muons = MaterializeMuons(AllMuonViews, SelectedMuonIndices);
  RVec<Electron> electrons = MaterializeElectrons(AllElectronViews, SelectedElectronIndices);

  Lepton lepton;
  if (muons.size() == 1) {
    lepton = Lepton(muons.at(0));
    if (muons.at(0).Pt() <= TriggerSafePtCut) return;
  } else if (electrons.size() == 1) {
    lepton = Lepton(electrons.at(0));
    if (electrons.at(0).Pt() <= 32) return; //TODO: set electron trigger safe cut properly
  } else {
    return;
  }


  MyCorrection::variation jes_variation = MyCorrection::variation::nom;
  MyCorrection::variation btag_jes_variation = MyCorrection::variation::nom;
  TString btag_source = "total";

  auto jet_id = apply_pu_id ? Jet::JetID::PUID_LOOSE : Jet::JetID::TIGHT;

  std::vector<size_t> SelectedJetIndices =
      SelectJetIndices(AllJetViews, jet_id, 0., 5.191, jes_variation, MyCorrection::variation::nom);

  RVec<Jet> jets =
      MaterializeJets(AllJetViews, SelectedJetIndices, jes_variation, MyCorrection::variation::nom);

  jets = JetsVetoLeptonInside(jets, electrons, muons, 0.3);

  std::sort(jets.begin(), jets.end(), PtComparing);

  if (jets.size() < 6) return;
  if (!PassJetVetoMap(AllJetViews, AllMuonViews, "jetvetomap_fpix")) return;

  float btag_wp_cut = myCorr->GetBTaggingWP();

  for (auto& jet : jets) {
    double bscore = jet.GetTaggerResult(
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerScoreType::B
    );

    if (bscore > btag_wp_cut) {
      double current_pt  = jet.Pt();
      double raw_pt      = jet.GetRawPt();
      double current_eta = jet.Eta();
      double current_phi = jet.Phi();
      double current_m   = jet.M();

      double L1L2L3Res = myCorr->GetJESSF(
          jet.GetArea(),
          current_eta,
          raw_pt,
          current_phi,
          Rho_fixedGridRhoFastjetAll,
          ev.run()
      );

      double L2L3Residual = myCorr->GetJESSF(
          jet.GetArea(),
          current_eta,
          raw_pt,
          current_phi,
          Rho_fixedGridRhoFastjetAll,
          ev.run(),
          "L2L3Residual"
      );

      double UParTAK4RegPtRawCorr = jet.UParTAK4RegPtRawCorr();
      double UParTAK4RegPtRawCorrNeutrino = jet.UParTAK4RegPtRawCorrNeutrino();
      double UParT_ratio = UParTAK4RegPtRawCorrNeutrino / UParTAK4RegPtRawCorr;

      double modified_pt = current_pt;
      double modified_m  = current_m;

      if (use_UParT_JEC) {
        modified_pt = current_pt * (1. / L1L2L3Res)
                    * UParTAK4RegPtRawCorrNeutrino
                    * L2L3Residual;
        modified_m = current_m * modified_pt / current_pt;
      }
      else {
        modified_pt = current_pt * UParT_ratio;
        modified_m  = current_m  * UParT_ratio;
      }

      if (correct_b_jet_pt) {
        jet.SetPtEtaPhiM(modified_pt, current_eta, current_phi, modified_m);
      }
    }
  }

  std::vector<bool> btag_vector;
  int n_bjets = 0;
  int n_jets_for_count = 0;

  for (unsigned int ij = 0; ij < jets.size(); ij++) {
    double bscore = jets.at(ij).GetTaggerResult(
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerScoreType::B
    );

    bool is_btag = bscore > btag_wp_cut;
    btag_vector.push_back(is_btag);

    if (is_btag) {
      n_bjets++;
      if (jets.at(ij).Pt() > 25.) n_jets_for_count++;
    }
    else {
      if (jets.at(ij).Pt() > 30.) n_jets_for_count++;
    }
  }

  if (n_bjets < 3) return;
  if (n_jets_for_count < 6) return;

  float weight_mc = 1.f;
  float weight_trigger_lumi = 1.f;
  float weight_pu = 1.f;
  float weight_btag = 1.f;
  float weight_muon_id = 1.f;
  float weight_muon_iso = 1.f;
  float weight_muon_trig = 1.f;

  if (!IsDATA) {
    weight_mc = MCweight();
    weight_trigger_lumi = ev.GetTriggerLumi("Full");

    weight_muon_id = myCorr->GetMuonIDSF(
        this_muon_id_sf_key,
        muons,
        MyCorrection::variation::nom
    );

    if (!use_pog_mva_tight_muon_id) {
      weight_muon_iso = myCorr->GetMuonIDSF(
          this_muon_iso_sf_key,
          muons,
          MyCorrection::variation::nom
      );
    }

    if (use_pog_tight_muon_id) {
      weight_muon_trig = myCorr->GetMuonTriggerSF(
          this_muon_trig_sf_key,
          muons,
          MyCorrection::variation::nom
      );
    }

    weight_pu = myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    RVec<Jet> jets_2p5 = RVec<Jet>();
    for (auto jet : jets) {
      if (abs(jet.Eta()) < 2.499 && jet.Pt() > 20.) {
        jets_2p5.push_back(jet);
      }
    }
    weight_btag = myCorr->GetBTaggingSF(
        jets_2p5,
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerWP::Medium,
        JetTagging::JetTaggingSFMethod::comb,
        btag_jes_variation,
        btag_source
    );
  }

  float weight_train = 1.f;
  if (!IsDATA) {
    weight_train = weight_mc * weight_trigger_lumi * weight_pu * weight_btag;
  }

  float weight_nominal = 1.f;
  if (!IsDATA) {
    weight_nominal = weight_mc
                   * weight_trigger_lumi
                   * weight_pu
                   * weight_btag
                   * weight_muon_id
                   * weight_muon_iso
                   * weight_muon_trig;
  }

  float best_chi2 = -999.f;
  float had_W_mass = -999.f;
  float had_top_mass = -999.f;
  float lep_W_mass = -999.f;
  float lep_top_mass = -999.f;

  std::vector<unsigned int> top_b_jet_candidates;
  std::vector<unsigned int> had_W_candidates;

  unsigned int n_check = std::min<unsigned int>(6, jets.size());

  for (unsigned int ij = 0; ij < n_check; ij++) {
    if (btag_vector.at(ij)) {
      if (top_b_jet_candidates.size() < 2) {
        top_b_jet_candidates.push_back(ij);
      }
    }
    else {
      if (had_W_candidates.size() < 3) {
        had_W_candidates.push_back(ij);
      }
    }
  }

  if (top_b_jet_candidates.size() >= 2 && had_W_candidates.size() >= 2) {
    std::vector<ttCombinatoric> combinatorics;

    for (unsigned int w1 = 0; w1 < had_W_candidates.size() - 1; w1++) {
      for (unsigned int w2 = w1 + 1; w2 < had_W_candidates.size(); w2++) {
        for (unsigned int b_idx = 0; b_idx < 2; b_idx++) {
          ttCombinatoric comb;
          comb.lepton = &lepton;
          comb.jets = &jets;
          comb.met = &METv;
          comb.had_W_jet_idx_1 = had_W_candidates.at(w1);
          comb.had_W_jet_idx_2 = had_W_candidates.at(w2);
          comb.had_top_b_jet_idx = top_b_jet_candidates.at(b_idx);
          comb.lep_top_b_jet_idx = top_b_jet_candidates.at(1 - b_idx);

          this->EvalChi2(comb);
          combinatorics.push_back(comb);
        }
      }
    }

    if (!combinatorics.empty()) {
      ttCombinatoric* best = &combinatorics.at(0);
      for (unsigned int i = 1; i < combinatorics.size(); i++) {
        if (combinatorics.at(i).best_chi2 < best->best_chi2) {
          best = &combinatorics.at(i);
        }
      }

      best_chi2 = best->best_chi2;
      had_W_mass = best->had_W_mass;
      had_top_mass = best->had_top_mass;
      lep_W_mass = best->best_lep_W_mass;
      lep_top_mass = best->best_lep_top_mass;
    }
  }

  FillTreeBranches(
      lepton,
      jets,
      METv,
      n_bjets,
      n_jets_for_count,
      btag_wp_cut,
      weight_train,
      weight_nominal,
      weight_mc * weight_trigger_lumi,
      weight_pu,
      weight_btag,
      weight_muon_id,
      weight_muon_iso,
      weight_muon_trig,
      best_chi2,
      had_W_mass,
      had_top_mass,
      lep_W_mass,
      lep_top_mass
  );
  
}

void AtobbMLTree::FillJetBranches(const TString& prefix, const Jet& jet, float btag_wp_cut) {
  float bscore = jet.GetTaggerResult(
      JetTagging::JetFlavTagger::ParT,
      JetTagging::JetFlavTaggerScoreType::B
  );

  SetBranch("Training_Tree", prefix + "_Pt",     float(jet.Pt()));
  SetBranch("Training_Tree", prefix + "_Eta",    float(jet.Eta()));
  SetBranch("Training_Tree", prefix + "_CosPhi", float(std::cos(jet.Phi())));
  SetBranch("Training_Tree", prefix + "_SinPhi", float(std::sin(jet.Phi())));
  SetBranch("Training_Tree", prefix + "_M",      float(jet.M()));
  SetBranch("Training_Tree", prefix + "_Bscore", bscore);
  SetBranch("Training_Tree", prefix + "_isBtag", int(bscore > btag_wp_cut));
  SetBranch("Training_Tree", prefix + "_isJet", 1);
}

void AtobbMLTree::FillDummyJetBranches(const TString& prefix) {
  SetBranch("Training_Tree", prefix + "_Pt",     -999.f);
  SetBranch("Training_Tree", prefix + "_Eta",    -999.f);
  SetBranch("Training_Tree", prefix + "_CosPhi", -999.f);
  SetBranch("Training_Tree", prefix + "_SinPhi", -999.f);
  SetBranch("Training_Tree", prefix + "_M",      -999.f);
  SetBranch("Training_Tree", prefix + "_Bscore", -999.f);
  SetBranch("Training_Tree", prefix + "_isBtag", 0);
  SetBranch("Training_Tree", prefix + "_isJet", 0);
}

void AtobbMLTree::FillTreeBranches(
    Lepton& lepton,
    RVec<Jet>& jets,
    Particle& METv,
    int n_bjets,
    int n_jets_for_count,
    float btag_wp_cut,
    float weight_train,
    float weight_nominal,
    float weight_mc,
    float weight_pu,
    float weight_btag,
    float weight_muon_id,
    float weight_muon_iso,
    float weight_muon_trig,
    float best_chi2,
    float had_W_mass,
    float had_top_mass,
    float lep_W_mass,
    float lep_top_mass
) {
  int label = 0;
  if (MCSample.Contains("AToBB") || MCSample.Contains("HcToWA")) {
    label = 1;
  }

  SetBranch("Training_Tree", "label", label);
  SetBranch("Training_Tree", "isData", int(IsDATA));

  SetBranch("Training_Tree", "weight_train", weight_train);
  SetBranch("Training_Tree", "weight_nominal", weight_nominal);
  SetBranch("Training_Tree", "weight_mc", weight_mc);
  SetBranch("Training_Tree", "weight_pu", weight_pu);
  SetBranch("Training_Tree", "weight_btag", weight_btag);
  SetBranch("Training_Tree", "weight_muon_id", weight_muon_id);
  SetBranch("Training_Tree", "weight_muon_iso", weight_muon_iso);
  SetBranch("Training_Tree", "weight_muon_trig", weight_muon_trig);

  SetBranch("Training_Tree", "nJets", int(jets.size()));
  SetBranch("Training_Tree", "nJets_for_count", n_jets_for_count);
  SetBranch("Training_Tree", "nBJets", n_bjets);

  SetBranch("Training_Tree", "Lepton0_Pt",     float(lepton.Pt()));
  SetBranch("Training_Tree", "Lepton0_Eta",    float(lepton.Eta()));
  SetBranch("Training_Tree", "Lepton0_CosPhi", float(std::cos(lepton.Phi())));
  SetBranch("Training_Tree", "Lepton0_SinPhi", float(std::sin(lepton.Phi())));
  SetBranch("Training_Tree", "Lepton0_M",      float(lepton.M()));
  SetBranch("Training_Tree", "Lepton0_isJet",  0);

  SetBranch("Training_Tree", "MET_Pt",  float(METv.Pt()));
  SetBranch("Training_Tree", "MET_Phi", float(METv.Phi()));

  float HT = 0.f;
  for (const auto& jet : jets) {
    HT += float(jet.Pt());
  }
  SetBranch("Training_Tree", "HT", HT);
  SetBranch("Training_Tree", "ST", HT + float(lepton.Pt()) + float(METv.Pt()));

  for (unsigned int i = 0; i < 6; i++) {
    TString prefix = "Jet" + TString(std::to_string(i));
    if (i < jets.size()) FillJetBranches(prefix, jets.at(i), btag_wp_cut);
    else FillDummyJetBranches(prefix);
  }

  std::vector<unsigned int> bjet_indices;
  for (unsigned int i = 0; i < jets.size(); i++) {
    float bscore = jets.at(i).GetTaggerResult(
        JetTagging::JetFlavTagger::ParT,
        JetTagging::JetFlavTaggerScoreType::B
    );
    if (bscore > btag_wp_cut) bjet_indices.push_back(i);
  }

  float lnu_mass = float((static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
  SetBranch("Training_Tree", "lnu_mass", lnu_mass);

  // Initialize with a default value (e.g., -999.0 or 0.0) in case there are fewer than 2 non-b-jets
  float jj_mass = 0.f; 

  // Vector to store indices of non-b-tagged jets
  std::vector<int> non_bjet_indices;
  // Iterate through all jets to find non-b-jets
  // Since 'jets' is already pT-ordered, the first ones found have the highest pT
  for (size_t i = 0; i < jets.size(); ++i) {
    // Check if the current jet index is NOT in the bjet_indices list
    if (std::find(bjet_indices.begin(), bjet_indices.end(), i) == bjet_indices.end()) {
      non_bjet_indices.push_back(i);
    }
    
    // Stop the loop once we find the two highest pT non-b-jets
    if (non_bjet_indices.size() >= 2) {
      break;
    }
  }
  if (non_bjet_indices.size() >= 2) {
    TLorentzVector j0 = static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0)));
    TLorentzVector j1 = static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)));
    jj_mass = float((j0 + j1).M());
    SetBranch("Training_Tree", "jj_mass", jj_mass);
  }
  else {
    SetBranch("Training_Tree", "jj_mass", jj_mass);
  }

  float bb_mass_01 = -999.f;
  float bb_mass_02 = -999.f;
  float bb_mass_12 = -999.f;
  float bb_mass_03 = 0.f;
  float bb_mass_13 = 0.f;
  float bb_mass_23 = 0.f;
  float bb_dr_01 = -999.f;
  float bb_dr_02 = -999.f;
  float bb_dr_12 = -999.f;
  float bb_dr_03 = 0.f;
  float bb_dr_13 = 0.f;
  float bb_dr_23 = 0.f;

  float bblnu_mass_01 = -999.f;
  float bblnu_mass_02 = -999.f;
  float bblnu_mass_12 = -999.f;
  float bblnu_mass_03 = 0.f;
  float bblnu_mass_13 = 0.f;
  float bblnu_mass_23 = 0.f;

  float bbjj_mass_01 = -999.f;
  float bbjj_mass_02 = -999.f;
  float bbjj_mass_12 = -999.f;
  float bbjj_mass_03 = 0.f;
  float bbjj_mass_13 = 0.f;
  float bbjj_mass_23 = 0.f;

  float bbbW_mass_012 = 0.f;
  float bbbW_mass_013 = 0.f;
  float bbbW_mass_023 = 0.f;
  float bbbW_mass_123 = 0.f;

  float bW_mass_013 = 0.f;
  float bW_mass_023 = 0.f;
  float bW_mass_123 = 0.f;


// Calculate kinematics for the first two b-jets
  if (bjet_indices.size() >= 2) {
    TLorentzVector b0 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(0)));
    TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(1)));
    
    bb_mass_01 = float((b0 + b1).M());
    bb_dr_01 = float(b0.DeltaR(b1));
    bblnu_mass_01 = float((b0 + b1 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
    if (non_bjet_indices.size() >= 2) {
      bbjj_mass_01 = float((b0 + b1 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
    }
  }

  // Calculate kinematics involving the third b-jet
  if (bjet_indices.size() >= 3) {
    TLorentzVector b0 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(0)));
    TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(1)));
    TLorentzVector b2 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(2)));

    bb_mass_02 = float((b0 + b2).M());
    bb_mass_12 = float((b1 + b2).M());

    bb_dr_02 = float(b0.DeltaR(b2));
    bb_dr_12 = float(b1.DeltaR(b2));

    bblnu_mass_02 = float((b0 + b2 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
    bblnu_mass_12 = float((b1 + b2 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
    if (non_bjet_indices.size() >= 2) {
      bbjj_mass_02 = float((b0 + b2 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      bbjj_mass_12 = float((b1 + b2 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());

      float bbbW_mass_012_cand1 = float((b0 + b1 + b2 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
      float bbbW_mass_012_cand2 = float((b0 + b1 + b2 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      if (std::abs(bbbW_mass_012_cand1 - const_top_mass) < std::abs(bbbW_mass_012_cand2 - const_top_mass)) {
        bbbW_mass_012 = bbbW_mass_012_cand1;
      }
      else {
        bbbW_mass_012 = bbbW_mass_012_cand2;
      }
    }
  }

  // Calculate kinematics involving the fourth b-jet
  if (bjet_indices.size() >= 4) {
    TLorentzVector b0 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(0)));
    TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(1)));
    TLorentzVector b2 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(2)));
    TLorentzVector b3 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(3)));

    bb_mass_03 = float((b0 + b3).M());
    bb_mass_13 = float((b1 + b3).M());
    bb_mass_23 = float((b2 + b3).M());

    bb_dr_03 = float(b0.DeltaR(b3));
    bb_dr_13 = float(b1.DeltaR(b3));
    bb_dr_23 = float(b2.DeltaR(b3));

    bblnu_mass_03 = float((b0 + b3 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
    bblnu_mass_13 = float((b1 + b3 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
    bblnu_mass_23 = float((b2 + b3 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());

    if (non_bjet_indices.size() >= 2) {
      bbjj_mass_03 = float((b0 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      bbjj_mass_13 = float((b1 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      bbjj_mass_23 = float((b2 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());

      float bbbW_mass_013_cand1 = float((b0 + b1 + b3 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
      float bbbW_mass_013_cand2 = float((b0 + b1 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      float bW_mass_013_cand1 = float((b2 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      float bW_mass_013_cand2 = float((b2 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
      if (std::abs(bbbW_mass_013_cand1 - const_top_mass) < std::abs(bbbW_mass_013_cand2 - const_top_mass)) {
        bbbW_mass_013 = bbbW_mass_013_cand1;
        bW_mass_013 = bW_mass_013_cand1;
      }
      else {
        bbbW_mass_013 = bbbW_mass_013_cand2;
        bW_mass_013 = bW_mass_013_cand2;
      }
      float bbbW_mass_023_cand1 = float((b0 + b2 + b3 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
      float bbbW_mass_023_cand2 = float((b0 + b2 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      float bW_mass_023_cand1 = float((b1 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      float bW_mass_023_cand2 = float((b1 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
      if (std::abs(bbbW_mass_023_cand1 - const_top_mass) < std::abs(bbbW_mass_023_cand2 - const_top_mass)) {
        bbbW_mass_023 = bbbW_mass_023_cand1;
        bW_mass_023 = bW_mass_023_cand1;
      }
      else {
        bbbW_mass_023 = bbbW_mass_023_cand2;
        bW_mass_023 = bW_mass_023_cand2;
      }
      float bbbW_mass_123_cand1 = float((b1 + b2 + b3 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
      float bbbW_mass_123_cand2 = float((b1 + b2 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      float bW_mass_123_cand1 = float((b0 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      float bW_mass_123_cand2 = float((b0 + static_cast<TLorentzVector>(lepton) + static_cast<TLorentzVector>(METv)).M());
      if (std::abs(bbbW_mass_123_cand1 - const_top_mass) < std::abs(bbbW_mass_123_cand2 - const_top_mass)) {
        bbbW_mass_123 = bbbW_mass_123_cand1;
        bW_mass_123 = bW_mass_123_cand1;
      }
      else {
        bbbW_mass_123 = bbbW_mass_123_cand2;
        bW_mass_123 = bW_mass_123_cand2;
      }
    }
  }

  // Set branches for invariant masses
  SetBranch("Training_Tree", "bb_mass_01", bb_mass_01);
  SetBranch("Training_Tree", "bb_mass_02", bb_mass_02);
  SetBranch("Training_Tree", "bb_mass_12", bb_mass_12);
  SetBranch("Training_Tree", "bb_mass_03", bb_mass_03);
  SetBranch("Training_Tree", "bb_mass_13", bb_mass_13);
  SetBranch("Training_Tree", "bb_mass_23", bb_mass_23);

  // Set branches for Delta R
  SetBranch("Training_Tree", "bb_dr_01", bb_dr_01);
  SetBranch("Training_Tree", "bb_dr_02", bb_dr_02);
  SetBranch("Training_Tree", "bb_dr_12", bb_dr_12);
  SetBranch("Training_Tree", "bb_dr_03", bb_dr_03);
  SetBranch("Training_Tree", "bb_dr_13", bb_dr_13);
  SetBranch("Training_Tree", "bb_dr_23", bb_dr_23);

  SetBranch("Training_Tree", "bblnu_mass_01", bblnu_mass_01);
  SetBranch("Training_Tree", "bblnu_mass_02", bblnu_mass_02);
  SetBranch("Training_Tree", "bblnu_mass_12", bblnu_mass_12);
  SetBranch("Training_Tree", "bblnu_mass_03", bblnu_mass_03);
  SetBranch("Training_Tree", "bblnu_mass_13", bblnu_mass_13);
  SetBranch("Training_Tree", "bblnu_mass_23", bblnu_mass_23);
  SetBranch("Training_Tree", "bbjj_mass_01", bbjj_mass_01);
  SetBranch("Training_Tree", "bbjj_mass_02", bbjj_mass_02);
  SetBranch("Training_Tree", "bbjj_mass_12", bbjj_mass_12);
  SetBranch("Training_Tree", "bbjj_mass_03", bbjj_mass_03);
  SetBranch("Training_Tree", "bbjj_mass_13", bbjj_mass_13);
  SetBranch("Training_Tree", "bbjj_mass_23", bbjj_mass_23);

  SetBranch("Training_Tree", "bbbW_mass_012", bbbW_mass_012);
  SetBranch("Training_Tree", "bbbW_mass_013", bbbW_mass_013);
  SetBranch("Training_Tree", "bbbW_mass_023", bbbW_mass_023);
  SetBranch("Training_Tree", "bbbW_mass_123", bbbW_mass_123);
  SetBranch("Training_Tree", "bW_mass_013", bW_mass_013);
  SetBranch("Training_Tree", "bW_mass_023", bW_mass_023);
  SetBranch("Training_Tree", "bW_mass_123", bW_mass_123);

  // 1. Construct all valid combinations of (b_cand_idx, w_cand_idx)
  struct JetCombination {
    std::vector<unsigned int> b_cands;
    std::vector<unsigned int> w_cands;
  };
  std::vector<JetCombination> jet_combinations;

  std::vector<unsigned int> b_indices = bjet_indices;
  std::vector<unsigned int> non_b_indices;
  for (unsigned int i = 0; i < std::min<unsigned int>(6, jets.size()); i++) {
    if (std::find(b_indices.begin(), b_indices.end(), i) == b_indices.end()) {
      non_b_indices.push_back(i);
    }
  }

  if (b_indices.size() >= 4) {
    JetCombination job;
    job.b_cands = {b_indices[0], b_indices[1], b_indices[2], b_indices[3]};
    for (unsigned int i = 0; i < std::min<unsigned int>(6, jets.size()); i++) {
      if (std::find(job.b_cands.begin(), job.b_cands.end(), i) == job.b_cands.end()) {
        job.w_cands.push_back(i);
      }
    }
    jet_combinations.push_back(job);
  } else if (b_indices.size() == 3) {
    for (unsigned int extra_b : non_b_indices) {
      JetCombination job;
      job.b_cands = {b_indices[0], b_indices[1], b_indices[2], extra_b};
      for (unsigned int i = 0; i < std::min<unsigned int>(6, jets.size()); i++) {
        if (std::find(job.b_cands.begin(), job.b_cands.end(), i) == job.b_cands.end()) {
          job.w_cands.push_back(i);
        }
      }
      if (job.w_cands.size() >= 2) {
        jet_combinations.push_back(job);
      }
    }
  }

  std::vector<float> vec_AWb_mass;
  std::vector<float> vec_AW_mass;
  std::vector<float> vec_A_mass;
  std::vector<float> vec_W_mass0;
  std::vector<float> vec_Wb_mass;
  std::vector<float> vec_W_mass1;
  std::vector<float> vec_A_dr;
  std::vector<float> vec_AW_dr;
  std::vector<float> vec_AWb_dr;

  // New constituent index vectors
  std::vector<int> A_cands;
  std::vector<int> W_signal_cands;
  std::vector<int> top_b_signal_cands;
  std::vector<int> W_cands;
  std::vector<int> top_b_cands;

  TLorentzVector lep = static_cast<TLorentzVector>(lepton);
  TLorentzVector met = static_cast<TLorentzVector>(METv);
  TLorentzVector w_lep = lep + met;

  for (const auto& comb : jet_combinations) {
    const auto& b_cand_idx = comb.b_cands;
    const auto& w_cand_idx = comb.w_cands;

    TLorentzVector b_vec[4];
    for (int idx = 0; idx < 4; ++idx) {
      b_vec[idx] = static_cast<TLorentzVector>(jets.at(b_cand_idx[idx]));
    }

    int comb_pairs[6][2] = {
      {0, 1}, {0, 2}, {0, 3},
      {1, 2}, {1, 3},
      {2, 3}
    };

    for (int p = 0; p < 6; ++p) {
      int i = comb_pairs[p][0];
      int j = comb_pairs[p][1];

      int k = -1, l = -1;
      for (int idx = 0; idx < 4; ++idx) {
        if (idx != i && idx != j) {
          if (k == -1) k = idx;
          else l = idx;
        }
      }

      int rem_perms[2][2] = {
        {k, l},
        {l, k}
      };

      for (int perm = 0; perm < 2; ++perm) {
        int idx_top_H = rem_perms[perm][0];
        int idx_top_other = rem_perms[perm][1];

        TLorentzVector b_A0 = b_vec[i];
        TLorentzVector b_A1 = b_vec[j];
        TLorentzVector b_top_H = b_vec[idx_top_H];
        TLorentzVector b_top_other = b_vec[idx_top_other];

        TLorentzVector A = b_A0 + b_A1;
        float A_mass = float(A.M());
        float A_dr = float(b_A0.DeltaR(b_A1));

        for (size_t m = 0; m < w_cand_idx.size() - 1; ++m) {
          for (size_t n = m + 1; n < w_cand_idx.size(); ++n) {
            TLorentzVector j0 = static_cast<TLorentzVector>(jets.at(w_cand_idx[m]));
            TLorentzVector j1 = static_cast<TLorentzVector>(jets.at(w_cand_idx[n]));
            TLorentzVector w_had = j0 + j1;

            // Mode 1: W_signal ( leptonic ) / W ( hadronic )
            {
              TLorentzVector W_signal = w_lep;
              TLorentzVector W = w_had;

              TLorentzVector AW = A + W_signal;
              TLorentzVector AWb = AW + b_top_H;
              TLorentzVector Wb = W + b_top_other;

              vec_A_mass.push_back(A_mass);
              vec_AW_mass.push_back(float(AW.M()));
              vec_AWb_mass.push_back(float(AWb.M()));
              vec_W_mass0.push_back(float(W_signal.M()));
              vec_W_mass1.push_back(float(W.M()));
              vec_Wb_mass.push_back(float(Wb.M()));

              vec_A_dr.push_back(A_dr);
              vec_AW_dr.push_back(float(A.DeltaR(W_signal)));
              vec_AWb_dr.push_back(float(AW.DeltaR(b_top_H)));

              // Index vectors
              A_cands.push_back(int(b_cand_idx[i]) + 2);
              A_cands.push_back(int(b_cand_idx[j]) + 2);

              W_signal_cands.push_back(0); // Lepton
              W_signal_cands.push_back(1); // MET

              top_b_signal_cands.push_back(int(b_cand_idx[idx_top_H]) + 2);

              W_cands.push_back(int(w_cand_idx[m]) + 2);
              W_cands.push_back(int(w_cand_idx[n]) + 2);

              top_b_cands.push_back(int(b_cand_idx[idx_top_other]) + 2);
            }

            // Mode 2: W_signal ( hadronic ) / W ( leptonic )
            {
              TLorentzVector W_signal = w_had;
              TLorentzVector W = w_lep;

              TLorentzVector AW = A + W_signal;
              TLorentzVector AWb = AW + b_top_H;
              TLorentzVector Wb = W + b_top_other;

              vec_A_mass.push_back(A_mass);
              vec_AW_mass.push_back(float(AW.M()));
              vec_AWb_mass.push_back(float(AWb.M()));
              vec_W_mass0.push_back(float(W_signal.M()));
              vec_W_mass1.push_back(float(W.M()));
              vec_Wb_mass.push_back(float(Wb.M()));

              vec_A_dr.push_back(A_dr);
              vec_AW_dr.push_back(float(A.DeltaR(W_signal)));
              vec_AWb_dr.push_back(float(AW.DeltaR(b_top_H)));

              // Index vectors
              A_cands.push_back(int(b_cand_idx[i]) + 2);
              A_cands.push_back(int(b_cand_idx[j]) + 2);

              W_signal_cands.push_back(int(w_cand_idx[m]) + 2);
              W_signal_cands.push_back(int(w_cand_idx[n]) + 2);

              top_b_signal_cands.push_back(int(b_cand_idx[idx_top_H]) + 2);

              W_cands.push_back(0); // Lepton
              W_cands.push_back(1); // MET

              top_b_cands.push_back(int(b_cand_idx[idx_top_other]) + 2);
            }
          }
        }
      }
    }
  }

  // Set vector branches for combinations
  SetBranch("Training_Tree", "vec_AWb_mass", vec_AWb_mass);
  SetBranch("Training_Tree", "vec_AW_mass",  vec_AW_mass);
  SetBranch("Training_Tree", "vec_A_mass",   vec_A_mass);
  SetBranch("Training_Tree", "vec_W_mass0",  vec_W_mass0);
  SetBranch("Training_Tree", "vec_Wb_mass",  vec_Wb_mass);
  SetBranch("Training_Tree", "vec_W_mass1",  vec_W_mass1);
  SetBranch("Training_Tree", "vec_A_dr",     vec_A_dr);
  SetBranch("Training_Tree", "vec_AW_dr",    vec_AW_dr);
  SetBranch("Training_Tree", "vec_AWb_dr",   vec_AWb_dr);

  SetBranch("Training_Tree", "A_cands",            A_cands);
  SetBranch("Training_Tree", "W_signal_cands",     W_signal_cands);
  SetBranch("Training_Tree", "top_b_signal_cands", top_b_signal_cands);
  SetBranch("Training_Tree", "W_cands",            W_cands);
  SetBranch("Training_Tree", "top_b_cands",        top_b_cands);

  std::vector<float> vec_pt;
  std::vector<float> vec_eta;
  std::vector<float> vec_cos_phi;
  std::vector<float> vec_sin_phi;
  std::vector<float> vec_M;
  std::vector<bool> vec_btagged;

  // Entry 0: Lepton
  vec_pt.push_back(float(lepton.Pt()));
  vec_eta.push_back(float(lepton.Eta()));
  vec_cos_phi.push_back(float(std::cos(lepton.Phi())));
  vec_sin_phi.push_back(float(std::sin(lepton.Phi())));
  vec_M.push_back(float(lepton.M()));
  vec_btagged.push_back(false);

  // Entry 1: MET
  vec_pt.push_back(float(METv.Pt()));
  vec_eta.push_back(float(METv.Eta()));
  vec_cos_phi.push_back(float(std::cos(METv.Phi())));
  vec_sin_phi.push_back(float(std::sin(METv.Phi())));
  vec_M.push_back(float(METv.M()));
  vec_btagged.push_back(false);

  // Remaining entries: Jets
  for (unsigned int i = 0; i < jets.size(); i++) {
    const auto& jet = jets.at(i);
    vec_pt.push_back(float(jet.Pt()));
    vec_eta.push_back(float(jet.Eta()));
    vec_cos_phi.push_back(float(std::cos(jet.Phi())));
    vec_sin_phi.push_back(float(std::sin(jet.Phi())));
    vec_M.push_back(float(jet.M()));

    bool is_btag = std::find(bjet_indices.begin(), bjet_indices.end(), i) != bjet_indices.end();
    vec_btagged.push_back(is_btag);
  }

  // Set raw vector branches
  SetBranch("Training_Tree", "vec_pt",      vec_pt);
  SetBranch("Training_Tree", "vec_eta",     vec_eta);
  SetBranch("Training_Tree", "vec_cos_phi", vec_cos_phi);
  SetBranch("Training_Tree", "vec_sin_phi", vec_sin_phi);
  SetBranch("Training_Tree", "vec_M",       vec_M);
  SetBranch("Training_Tree", "vec_btagged", vec_btagged);

  SetBranch("Training_Tree", "best_chi2", best_chi2);
  SetBranch("Training_Tree", "had_W_mass", had_W_mass);
  SetBranch("Training_Tree", "had_top_mass", had_top_mass);
  SetBranch("Training_Tree", "lep_W_mass", lep_W_mass);
  SetBranch("Training_Tree", "lep_top_mass", lep_top_mass);

  FillTrees("Training_Tree");
}

void AtobbMLTree::EvalChi2(ttCombinatoric& tt_combinatoric) {
  tt_combinatoric.EvalHadronicPart();
  tt_combinatoric.EvalLeptonicPart();

  tt_combinatoric.best_chi2 = 1e9;

  for (unsigned int i = 0; i < tt_combinatoric.neu_pz.size(); i++) {
    double chi2 = this->Chi2Function(
        tt_combinatoric.had_top_mass,
        tt_combinatoric.had_W_mass,
        tt_combinatoric.lep_top_mass.at(i),
        tt_combinatoric.lep_W_mass.at(i)
    );

    tt_combinatoric.chi2.push_back(chi2);

    if (chi2 < tt_combinatoric.best_chi2) {
      tt_combinatoric.best_lep_top_mass = tt_combinatoric.lep_top_mass.at(i);
      tt_combinatoric.best_lep_W_mass   = tt_combinatoric.lep_W_mass.at(i);
      tt_combinatoric.best_neu_pz       = tt_combinatoric.neu_pz.at(i);
      tt_combinatoric.best_chi2         = chi2;
    }
  }
}

double AtobbMLTree::Chi2Function(
    double had_top_mass,
    double had_W_mass,
    double lep_top_mass,
    double lep_W_mass
) {
  double chi2 = 0.;
  chi2 += TMath::Power((had_top_mass - const_top_mass) / const_top_width, 2);
  chi2 += TMath::Power((had_W_mass   - const_w_mass)   / const_w_width,   2);
  chi2 += TMath::Power((lep_top_mass - const_top_mass) / const_top_width, 2);
  chi2 += TMath::Power((lep_W_mass   - const_w_mass)   / const_w_width,   2);
  return chi2;
}

void AtobbMLTree::ttCombinatoric::EvalHadronicPart() {
  TLorentzVector had_W_vector =
      static_cast<TLorentzVector>(jets->at(had_W_jet_idx_1))
    + static_cast<TLorentzVector>(jets->at(had_W_jet_idx_2));

  TLorentzVector had_top_vector =
      had_W_vector
    + static_cast<TLorentzVector>(jets->at(had_top_b_jet_idx));

  had_W_mass = had_W_vector.M();
  had_top_mass = had_top_vector.M();
}

void AtobbMLTree::ttCombinatoric::EvalLeptonicPart() {
  lep_top_mass.clear();
  lep_W_mass.clear();
  neu_pz.clear();
  chi2.clear();

  double step = 5.;
  double pz = -700.;

  while (pz <= 700.) {
    neu_pz.push_back(pz);
    pz += step;
  }

  for (auto& a_neu_pz : neu_pz) {
    double E_neu = TMath::Sqrt(met->E() * met->E() + a_neu_pz * a_neu_pz);
    TLorentzVector neutrino_vector(met->Px(), met->Py(), a_neu_pz, E_neu);

    TLorentzVector lep_W_vector =
        neutrino_vector
      + static_cast<TLorentzVector>(*lepton);

    TLorentzVector lep_top_vector =
        lep_W_vector
      + static_cast<TLorentzVector>(jets->at(lep_top_b_jet_idx));

    lep_W_mass.push_back(lep_W_vector.M());
    lep_top_mass.push_back(lep_top_vector.M());
  }
}
