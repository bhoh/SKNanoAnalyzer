#include "Tutorial_reco_tt.h"
#include "lester_mt2_bisect.h"

Tutorial_reco_tt::Tutorial_reco_tt() :
  const_top_mass(172.5),
  const_top_width(1.5),
  const_w_mass(80.4),
  const_w_width(2.085)
{
}

Tutorial_reco_tt::~Tutorial_reco_tt() {}

void Tutorial_reco_tt::initializeAnalyzer() {

  MuonIDs = { Muon::MuonID::POG_TIGHT, Muon::MuonID::POG_MEDIUM_PROMPT, Muon::MuonID::POG_MVA_MU_TIGHT };
  MuonIDISOSFKeys = { "NUM_TightID_DEN_TrackerMuons", "NUM_TightPFIso_DEN_TightID", "NUM_MediumPromptID_DEN_TrackerMuons", "NUM_TightPFIso_DEN_MediumPromptID", "NUM_TightMvaMuID_DEN_TrackerMuons" };

  if (DataEra == "2016preVFP" || DataEra == "2016postVFP" ||
      DataEra == "2018") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  } else if (DataEra == "2017") {
    IsoMuTriggerName = "HLT_IsoMu27";
    TriggerSafePtCut = 29.;
  } else if (DataEra == "2022") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  } else if (DataEra == "2022EE") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  } else if (DataEra == "2023") {
    IsoMuTriggerName = "";
    TriggerSafePtCut = 26.;
  } else if (DataEra == "2024") {
    IsoMuTriggerName = "HLT_IsoMu24";
    TriggerSafePtCut = 26.;
  } else {
    cerr << "[ExampleRun::initializeAnalyzer] DataEra is not set properly"
         << endl;
    exit(EXIT_FAILURE);
  }

  // init B-Tagging (DeepJet Medium WP example)
  myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
  myCorr->SetTaggingParam(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);

  // init SystematicHelper
  string SKNANO_HOME = getenv("SKNANO_HOME");
  if (IsDATA) {
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
  } else {
    //systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/JesTotal.yaml", MCSample, DataEra);
  }
}

void Tutorial_reco_tt::executeEvent() {

  AllMuonViews = GetAllMuonViews();
  AllElectronViews = GetAllElectronViews();
  AllJetViews = GetAllJetViews();
  AllGenViews = GetAllGenViews();
  AllGenJetViews = GetAllGenJetViews();

  ev = GetEvent();

  // Check this for Run3
  weight_Prefire = 1.; 

  // Systematic sources from YAML

  for (const auto &syst_dummy : *systHelper) {
    executeEventFromParameter();
  }
}

void Tutorial_reco_tt::executeEventFromParameter() {

  bool draw_include_pu_jets = false;
  bool correct_b_jet_pt = true;
  bool apply_pu_id = false;
  bool use_pog_tight_muon_id = true; // if not use medium prompt ID
  bool use_pog_mva_tight_muon_id = false;
  bool eval_top_pt_reweight_normalization = false;

  const TString this_syst = systHelper->getCurrentSysName();
  if (IsDATA && this_syst != "Central") return;

  if(eval_top_pt_reweight_normalization && MCSample.Contains("TT") && this_syst == "Central") {
    auto [firstTopIdx, firstAntiTopIdx, lastTopIdx, lastAntiTopIdx] =
    GetTopAndAntiTopIndices(AllGenViews);

    const TLorentzVector top = AllGenViews[firstTopIdx].P4();
    const TLorentzVector antiTop = AllGenViews[firstAntiTopIdx].P4();
    float w_toppt = myCorr->GetTopPtReweight(top, antiTop);
    FillHist(this_syst + "/count/w_toppt" + this_syst, 0.5, 1., 2, 0., 2.);
    FillHist(this_syst + "/count/w_toppt" + this_syst, 1.5, w_toppt, 2, 0., 2.);
  }


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

  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 0.5, 1., 6, 0., 6.);

  //==== MET Filter & Trigger
  if (!PassMetFilter(AllJetViews, ev)) return;
  if (!(ev.PassTrigger(IsoMuTriggerName))) return;

  Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI); 

  //==== Lepton Selection
  std::vector<size_t> SelectedMuonIndices_id_only = SelectMuonIndices(AllMuonViews, this_muon_id, 15., 2.4);
  std::vector<size_t> SelectedMuonIndices = {};
  if (use_pog_mva_tight_muon_id){
    SelectedMuonIndices = SelectedMuonIndices_id_only;
  }
  else{
    SelectedMuonIndices = SelectMuonIndices(AllMuonViews, SelectedMuonIndices_id_only, Muon::MuonID::POG_PFISO_TIGHT, 15., 2.4);
  }
  std::vector<size_t> SelectedElectronIndices = SelectElectronIndices(AllElectronViews, Electron::ElectronID::POG_LOOSE, 15., 2.5);

  if (SelectedMuonIndices.size() + SelectedElectronIndices.size() != 2) return;

  RVec<Muon> muons = MaterializeMuons(AllMuonViews, SelectedMuonIndices);
  RVec<Electron> electrons = MaterializeElectrons(AllElectronViews, SelectedElectronIndices);

  //==== Jet Selection
  MyCorrection::variation jes_variation = MyCorrection::variation::nom;
  MyCorrection::variation btag_jes_variation = MyCorrection::variation::nom;
  TString btag_source = "total";
  if (this_syst.Contains("JESTotal")) {
    ApplyJetScaleVariation(AllJetViews, "total");
    if (this_syst.Contains("Up")) {
      jes_variation = MyCorrection::variation::up;
      btag_jes_variation = MyCorrection::variation::nom;
      btag_source = "total";
    } else if (this_syst.Contains("Down")) {
      jes_variation = MyCorrection::variation::down;
      btag_jes_variation = MyCorrection::variation::nom;
      btag_source = "total";
    }
  }
  auto jet_id = apply_pu_id ? Jet::JetID::PUID_LOOSE : Jet::JetID::TIGHT;
  std::vector<size_t> SelectedJetIndices = SelectJetIndices(AllJetViews, jet_id, 25., 2.4, jes_variation, MyCorrection::variation::nom);
  RVec<Jet> jets = MaterializeJets(AllJetViews, SelectedJetIndices, jes_variation, MyCorrection::variation::nom);
  jets = JetsVetoLeptonInside(jets, electrons, muons, 0.3);

  //==== Sorting
  sort(muons.begin(), muons.end(), PtComparing);
  sort(jets.begin(), jets.end(), PtComparing);

  //==== Event selections
  if (muons.size() != 2) return;
  if (electrons.size() != 0) return;
  if (muons.at(1).Pt() <= TriggerSafePtCut) return;
  if (jets.size() < 4) return;
  //if (METv.Pt() <= 20) return;
  //require OS pair
  if (muons.at(0).Charge() * muons.at(1).Charge() >= 0) return;
  float mll = (muons.at(0) + muons.at(1)).M();
  // Z mass window veto
  if (mll > 76. && mll < 106.) return;
  // meason resonance veto
  if (mll < 15.) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 1.5, 1., 6, 0., 6.);

  if (!PassJetVetoMap(AllJetViews, AllMuonViews, "jetvetomap_fpix")) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 2.5, 1., 6, 0., 6.);

  //==== B-Tagging (DeepJet Medium WP example)
  int NBJets = 0;
  int njets_pt30_non_btagged = 0;
  float btag_wp_cut = myCorr->GetBTaggingWP(); 
  std::vector<bool> btag_vector;

  for (unsigned int ij = 0; ij < jets.size(); ij++) {
    double this_discr = jets.at(ij).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    if (this_discr > btag_wp_cut) {
      NBJets++;
      btag_vector.push_back(true);
      njets_pt30_non_btagged++;
    } else {
      btag_vector.push_back(false);
      if (jets.at(ij).Pt() > 30){
        njets_pt30_non_btagged++;
      }
    }
  }





  if (NBJets != 3) return;
  if (njets_pt30_non_btagged < 4) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 3.5, 1., 6, 0., 6.);

  //==== Event Weight
  float weight = 1.;
  if (!IsDATA) {
    weight *= MCweight();
    weight *= ev.GetTriggerLumi("Full");
    // muon official trigger SF is not available yet.
    float muon_id_sf = myCorr->GetMuonIDSF(this_muon_id_sf_key, muons, MyCorrection::variation::nom);
    weight *= muon_id_sf;
    float muon_iso_sf = 1;
    if (use_pog_mva_tight_muon_id) {
      // isolation SF is already included in the ID SF for the MVA tight WP, so we don't apply it separately
    }
    else {
      muon_iso_sf = myCorr->GetMuonIDSF(this_muon_iso_sf_key, muons, MyCorrection::variation::nom);
      weight *= muon_iso_sf;
    }
    if (use_pog_tight_muon_id) {
        float muon_trig_sf = myCorr->GetMuonTriggerSF(this_muon_trig_sf_key, muons, MyCorrection::variation::nom);
        weight *= muon_trig_sf;
    }
    float pu_weight = myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    weight *= pu_weight;
    float btag_sf = myCorr->GetBTaggingSF(jets, 
            JetTagging::JetFlavTagger::ParT, 
            JetTagging::JetFlavTaggerWP::Medium,
            JetTagging::JetTaggingSFMethod::comb,
            btag_jes_variation, btag_source
        );
    weight *= btag_sf;

    if (MCSample.Contains("TT")) {
      auto [firstTopIdx, firstAntiTopIdx, lastTopIdx, lastAntiTopIdx] =
          GetTopAndAntiTopIndices(AllGenViews);

      const TLorentzVector top = AllGenViews[firstTopIdx].P4();
      const TLorentzVector antiTop = AllGenViews[firstAntiTopIdx].P4();
      float w_toppt = myCorr->GetTopPtReweight(top, antiTop);
      weight *= w_toppt * 1.2360; // 1.2360 for top_pt_reweight normalization correction

      auto [topIdx, WTopIdx, BHadTopIdx, antiTopIdx, WAntiTopIdx,
          BHadAntiTopIdx] = myCorr->GetGenIdxofTopDecayProducts(AllGenViews);
      float weight_bfrag = 1.f;
      float weight_bfrag_up = 1.f;
      float xb = -1.f;
      float xb_anti = -1.f;
      if ((BHadTopIdx == std::numeric_limits<std::size_t>::max()) ||
          (BHadAntiTopIdx == std::numeric_limits<std::size_t>::max())) {
        weight_bfrag = -1.f;
        weight_bfrag_up = -1.f;
      } else {
        auto LastCopyTop = AllGenViews[topIdx].P4();
        auto LastCopyAntiTop = AllGenViews[antiTopIdx].P4();
        auto LastCopyWPlus = AllGenViews[WTopIdx].P4();
        auto LastCopyWMinus = AllGenViews[WAntiTopIdx].P4();
        auto FirstCopyAntiTopBHad = AllGenViews[BHadAntiTopIdx].P4();
        auto FirstCopyTopBHad = AllGenViews[BHadTopIdx].P4();

        const float x_e_top =
            2 * FirstCopyTopBHad * LastCopyTop / LastCopyTop.M2();
        const float x_e_antitop =
            2 * FirstCopyAntiTopBHad * LastCopyAntiTop / LastCopyAntiTop.M2();
        const float w_top = LastCopyWPlus.M2() / LastCopyTop.M2();
        const float w_antitop = LastCopyWMinus.M2() / LastCopyAntiTop.M2();
        const float clip_value = 1.2f;
        const float x_b_top = std::min(x_e_top / (1 - w_top), clip_value);
        const float x_b_antitop =
            std::min(x_e_antitop / (1 - w_antitop), clip_value);
        xb = x_b_top;
        xb_anti = x_b_antitop;

        weight_bfrag = myCorr->GetBFragReweight(
            LastCopyTop, LastCopyAntiTop, LastCopyWPlus, LastCopyWMinus,
            FirstCopyTopBHad, FirstCopyAntiTopBHad, MyCorrection::variation::nom);
        weight_bfrag_up = myCorr->GetBFragReweight(
            LastCopyTop, LastCopyAntiTop, LastCopyWPlus, LastCopyWMinus,
            FirstCopyTopBHad, FirstCopyAntiTopBHad, MyCorrection::variation::up);
      }
      //weight *= weight_bfrag;

      if(genTtbarId%100>=51 && genTtbarId%100<=55){
        weight *= 1.36;
      }
      else if(genTtbarId%100>=41 && genTtbarId%100<=45){
        weight *= 1.11;
      }
    }

  }

  unordered_map<int, int> matched_genjet_idx = GenJetMatching(jets, MaterializeGenJets(AllGenJetViews), Rho_fixedGridRhoFastjetAll);
  bool isPileupJet = false;
  for (auto &[reco_idx, gen_idx] : matched_genjet_idx) {
    if (gen_idx == -999) {
      isPileupJet = true;
      break;
    }
  }

  for (auto& jet : jets) {
    // Get the b-tagging discriminator score
    double this_discr = jet.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    
    // Check if the jet is b-tagged
    if (this_discr > btag_wp_cut) {
      FillHist(this_syst + "/baseLineCut/btagged_jet_pt0_" + this_syst, float(jet.Pt()), weight, 80, 0., 400.);
      break;
    }
  }

  // Loop through the jets and apply corrections to b-tagged ones
  for (auto& jet : jets) {
    // Get the b-tagging discriminator score
    double this_discr = jet.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    
    // Check if the jet is b-tagged
    if (this_discr > btag_wp_cut) {
        
        double current_pt  = jet.Pt();
        double current_eta = jet.Eta();
        double current_phi = jet.Phi();
        double current_m   = jet.M();

        double UParTAK4RegPtRawCorr = jet.UParTAK4RegPtRawCorr();
        double UParTAK4RegPtRawCorrNeutrino = jet.UParTAK4RegPtRawCorrNeutrino();
        double UParT_ratio = UParTAK4RegPtRawCorrNeutrino / UParTAK4RegPtRawCorr;
        FillHist(this_syst + "/corrections/UParTAK4RegPtRawCorr_" + this_syst, UParTAK4RegPtRawCorr, weight, 80, 0., 2);
        FillHist(this_syst + "/corrections/UParTAK4RegPtRawCorrNeutrino_" + this_syst, UParTAK4RegPtRawCorrNeutrino, weight, 80, 0., 2);
        FillHist(this_syst + "/corrections/UParT_ratio_" + this_syst, UParT_ratio, weight, 80, 0., 2);
        // Calculate your modified pT
        double modified_pt = current_pt * UParT_ratio;
        double modified_m  = current_m  * UParT_ratio;

        // Update the LorentzVector with the new pT
        if(correct_b_jet_pt){
          jet.SetPtEtaPhiM(modified_pt, current_eta, current_phi, modified_m);
        }
        
    }
  }

  for (auto& jet : jets) {
    // Get the b-tagging discriminator score
    double this_discr = jet.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    
    // Check if the jet is b-tagged
    if (this_discr > btag_wp_cut) {
      FillHist(this_syst + "/baseLineCut/btagged_RegCorr_jet_pt0_" + this_syst, float(jet.Pt()), weight, 80, 0., 400.);
      break;
    }
  }

  float muon_pt0 = muons.at(0).Pt();
  float muon_eta0 = muons.at(0).Eta();
  float muon_pt1 = muons.at(1).Pt();
  float muon_eta1 = muons.at(1).Eta();
  
  float jet_pt0 = jets.at(0).Pt();
  float jet_eta0 = jets.at(0).Eta();
  float njets = njets_pt30_non_btagged;
  float MET_pt = METv.Pt();
  float MET_phi = METv.Phi();
  FillHist(this_syst + "/baseLineCut/muon_pt0_" + this_syst, muon_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/baseLineCut/muon_eta0_" + this_syst, muon_eta0, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/baseLineCut/muon_pt1_" + this_syst, muon_pt1, weight, 80, 0., 400.);
  FillHist(this_syst + "/baseLineCut/muon_eta1_" + this_syst, muon_eta1, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/baseLineCut/njets_" + this_syst, njets, weight, 10, 0., 10.);
  FillHist(this_syst + "/baseLineCut/mll_" + this_syst, mll, weight, 80, 0., 400.);


  std::vector<size_t> SelectedJetIndices2 = SelectJetIndices(AllJetViews, jet_id, 40., 2.4, jes_variation, MyCorrection::variation::nom);
  RVec<Jet> jets2 = MaterializeJets(AllJetViews, SelectedJetIndices2, jes_variation, MyCorrection::variation::nom);
  jets2 = JetsVetoLeptonInside(jets2, electrons, muons, 0.3);
  FillHist(this_syst + "/baseLineCut/njets2_" + this_syst, float(jets2.size()), weight, 10, 0., 10.);
  FillHist(this_syst + "/baseLineCut/jet_pt0_" + this_syst, jet_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/baseLineCut/jet_eta0_" + this_syst, jet_eta0, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/baseLineCut/MET_pt_" + this_syst, MET_pt, weight, 40, 0., 200.);
  FillHist(this_syst + "/baseLineCut/MET_phi_" + this_syst, MET_phi, weight, 40, -3.14, 3.14);

  // MT2 calculation
  
  // 1. Define the hypothesised mass of the invisible particles.
  //    (Set to 0.0 for neutrinos).
  double mInvis = 0.0; 
  
  // ==========================================
  // PAIRING 1: (muon 0 + jet 0) and (muon 1 + jet 1)
  // ==========================================
  TLorentzVector visA_1 = muons.at(0) + jets.at(0);
  TLorentzVector visB_1 = muons.at(1) + jets.at(1);
  
  double mt2_pairing1 = asymm_mt2_lester_bisect::get_mT2(
      visA_1.M(), visA_1.Px(), visA_1.Py(),
      visB_1.M(), visB_1.Px(), visB_1.Py(),
      METv.Px(), METv.Py(),
      mInvis, mInvis
  );
  
  // ==========================================
  // PAIRING 2: (muon 0 + jet 1) and (muon 1 + jet 0)
  // ==========================================
  TLorentzVector visA_2 = muons.at(0) + jets.at(1);
  TLorentzVector visB_2 = muons.at(1) + jets.at(0);
  
  double mt2_pairing2 = asymm_mt2_lester_bisect::get_mT2(
      visA_2.M(), visA_2.Px(), visA_2.Py(),
      visB_2.M(), visB_2.Px(), visB_2.Py(),
      METv.Px(), METv.Py(),
      mInvis, mInvis
  );
  
  // ==========================================
  // FINAL RESULT: Resolve ambiguity
  // ==========================================
  // The correct physical MT2 is the minimum of the possible valid groupings.
  double mt2_0 = 0;
  double mt2_1 = 0;
  double mbl_0 = 0;
  double mbl_1 = 0;
  if(mt2_pairing1 < mt2_pairing2){
    mt2_0 = mt2_pairing1;
    mt2_1 = mt2_pairing2;
    mbl_0 = (visA_1.M() + visB_1.M())/2.;
    mbl_1 = (visA_2.M() + visB_2.M())/2.;
  }
  else{
    mt2_0 = mt2_pairing2;
    mt2_1 = mt2_pairing1;
    mbl_0 = (visA_2.M() + visB_2.M())/2.;
    mbl_1 = (visA_1.M() + visB_1.M())/2.;
  }


  double final_mt2 = std::min(mt2_pairing1, mt2_pairing2);
  FillHist(this_syst + "/baseLineCut/MT2_0_" + this_syst, mt2_0, weight, 80, 0., 400.);
  FillHist(this_syst + "/baseLineCut/MT2_1_" + this_syst, mt2_1, weight, 80, 0., 800.);
  FillHist(this_syst + "/baseLineCut/Mbl_0_" + this_syst, mbl_0, weight, 40, 0., 200.);
  FillHist(this_syst + "/baseLineCut/Mbl_1_" + this_syst, mbl_1, weight, 40, 0., 400.);

  if (!IsDATA && draw_include_pu_jets) {
    if(isPileupJet){
      FillHist("Pileup/" + this_syst + "/baseLineCut/muon_pt0_" + this_syst, muon_pt0, weight, 80, 0., 400.);
      FillHist("Pileup/" + this_syst + "/baseLineCut/muon_eta0_" + this_syst, muon_eta0, weight, 40, -2.4, 2.4);
      FillHist("Pileup/" + this_syst + "/baseLineCut/njets_" + this_syst, njets, weight, 10, 0., 10.);
      FillHist("Pileup/" + this_syst + "/baseLineCut/njets2_" + this_syst, float(jets2.size()), weight, 10, 0., 10.);
      FillHist("Pileup/" + this_syst + "/baseLineCut/jet_pt0_" + this_syst, jet_pt0, weight, 80, 0., 400.);
      FillHist("Pileup/" + this_syst + "/baseLineCut/jet_eta0_" + this_syst, jet_eta0, weight, 40, -2.4, 2.4);
      FillHist("Pileup/" + this_syst + "/baseLineCut/MET_pt_" + this_syst, MET_pt, weight, 40, 0., 200.);
      FillHist("Pileup/" + this_syst + "/baseLineCut/MET_phi_" + this_syst, MET_phi, weight, 40, -3.14, 3.14); 
    }
    else{
      FillHist("noPileup/" + this_syst + "/baseLineCut/muon_pt0_" + this_syst, muon_pt0, weight, 80, 0., 400.);
      FillHist("noPileup/" + this_syst + "/baseLineCut/muon_eta0_" + this_syst, muon_eta0, weight, 40, -2.4, 2.4);
      FillHist("noPileup/" + this_syst + "/baseLineCut/njets_" + this_syst, njets, weight, 10, 0., 10.);
      FillHist("noPileup/" + this_syst + "/baseLineCut/njets2_" + this_syst, float(jets2.size()), weight, 10, 0., 10.);
      FillHist("noPileup/" + this_syst + "/baseLineCut/jet_pt0_" + this_syst, jet_pt0, weight, 80, 0., 400.);
      FillHist("noPileup/" + this_syst + "/baseLineCut/jet_eta0_" + this_syst, jet_eta0, weight, 40, -2.4, 2.4);
      FillHist("noPileup/" + this_syst + "/baseLineCut/MET_pt_" + this_syst, MET_pt, weight, 40, 0., 200.);
      FillHist("noPileup/" + this_syst + "/baseLineCut/MET_phi_" + this_syst, MET_phi, weight, 40, -3.14, 3.14); 
     
    }
  }
// MT2 cut
if (mt2_0 > 180) return;

  FillHist(this_syst + "/MT2Cut/muon_pt0_" + this_syst, muon_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/MT2Cut/muon_eta0_" + this_syst, muon_eta0, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/MT2Cut/muon_pt1_" + this_syst, muon_pt1, weight, 80, 0., 400.);
  FillHist(this_syst + "/MT2Cut/muon_eta1_" + this_syst, muon_eta1, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/MT2Cut/njets_" + this_syst, njets, weight, 10, 0., 10.);
  FillHist(this_syst + "/MT2Cut/mll_" + this_syst, mll, weight, 80, 0., 400.);


  FillHist(this_syst + "/MT2Cut/njets2_" + this_syst, float(jets2.size()), weight, 10, 0., 10.);
  FillHist(this_syst + "/MT2Cut/jet_pt0_" + this_syst, jet_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/MT2Cut/jet_eta0_" + this_syst, jet_eta0, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/MT2Cut/MET_pt_" + this_syst, MET_pt, weight, 40, 0., 200.);
  FillHist(this_syst + "/MT2Cut/MET_phi_" + this_syst, MET_phi, weight, 40, -3.14, 3.14);

  FillHist(this_syst + "/MT2Cut/MT2_0_" + this_syst, mt2_0, weight, 80, 0., 400.);
  FillHist(this_syst + "/MT2Cut/MT2_1_" + this_syst, mt2_1, weight, 80, 0., 800.);
  FillHist(this_syst + "/MT2Cut/Mbl_0_" + this_syst, mbl_0, weight, 40, 0., 200.);
  FillHist(this_syst + "/MT2Cut/Mbl_1_" + this_syst, mbl_1, weight, 40, 0., 400.);



}

void Tutorial_reco_tt::EvalChi2(ttCombinatoric& tt_combinatoric) {
  tt_combinatoric.EvalHadronicPart();
  tt_combinatoric.EvalLeptonicPart();

  tt_combinatoric.best_chi2 = 1e9;
  for(unsigned int i(0); i<tt_combinatoric.neu_pz.size(); i++){
    double chi2 = this->Chi2Function(
            tt_combinatoric.had_top_mass,
            tt_combinatoric.had_W_mass,
            tt_combinatoric.lep_top_mass.at(i),
            tt_combinatoric.lep_W_mass.at(i)
            );

    tt_combinatoric.chi2.push_back(chi2);

    if(chi2 < tt_combinatoric.best_chi2){
      tt_combinatoric.best_lep_top_mass = tt_combinatoric.lep_top_mass.at(i);
      tt_combinatoric.best_lep_W_mass   = tt_combinatoric.lep_W_mass  .at(i);
      tt_combinatoric.best_neu_pz       = tt_combinatoric.neu_pz      .at(i);
      tt_combinatoric.best_chi2         = chi2;
    }
  }
}

double Tutorial_reco_tt::Chi2Function(double had_top_mass, double had_W_mass, double lep_top_mass, double lep_W_mass) {
  double chi2 = 0.;
  chi2 += TMath::Power( (had_top_mass - const_top_mass )/const_top_width,  2);
  chi2 += TMath::Power( (had_W_mass   - const_w_mass   )/const_w_width,    2);
  chi2 += TMath::Power( (lep_top_mass - const_top_mass )/const_top_width,  2);
  chi2 += TMath::Power( (lep_W_mass   - const_w_mass   )/const_w_width,    2);
  return chi2;
}

//=== define Tutorial_reco_tt::ttCombinatoric
void Tutorial_reco_tt::ttCombinatoric::EvalHadronicPart() {
  TLorentzVector had_W_vector   = static_cast<TLorentzVector>(jets->at( had_W_jet_idx_1 ))
                                + static_cast<TLorentzVector>(jets->at( had_W_jet_idx_2 ));

  TLorentzVector had_top_vector = had_W_vector
                                + static_cast<TLorentzVector>(jets->at( had_top_b_jet_idx ));

  had_W_mass   = had_W_vector.M();
  had_top_mass = had_top_vector.M();
}

void Tutorial_reco_tt::ttCombinatoric::EvalLeptonicPart() {
  lep_top_mass.clear();
  lep_W_mass  .clear();
  neu_pz      .clear();
  chi2        .clear();

  double step = 5.;
  double pz   = -700.;
  while(pz<=700){
    neu_pz.push_back(pz);
    pz += step;
  }

  for(auto& a_neu_pz : neu_pz){
    double E_neu = TMath::Sqrt( met->E() * met->E() + a_neu_pz * a_neu_pz );
    TLorentzVector neutrino_vector(met->Px(), met->Py(), a_neu_pz, E_neu);

    TLorentzVector lep_W_vector   = neutrino_vector
                                  + static_cast<TLorentzVector>(*lepton);
    TLorentzVector lep_top_vector = lep_W_vector
                                  + static_cast<TLorentzVector>(jets->at( lep_top_b_jet_idx ));

    lep_W_mass  .push_back(lep_W_vector.M());
    lep_top_mass.push_back(lep_top_vector.M());
  }
}

array<size_t, 4> Tutorial_reco_tt::GetTopAndAntiTopIndices(const GenViewCollection &gens) {
  constexpr size_t npos = std::numeric_limits<size_t>::max();

  size_t FirstCopyTopIndex = npos;
  size_t FirstCopyAntiTopIndex = npos;
  size_t LastCopyTopIndex = npos;
  size_t LastCopyAntiTopIndex = npos;

  const size_t n = gens.size();

  constexpr unsigned long FIRST_COPY_BIT = 1UL << 12;
  constexpr unsigned long LAST_COPY_BIT = 1UL << 13;

  for (size_t idx = 0; idx < n; ++idx) {
    const GenView &gen = gens[idx];

    const int pdg = gen.PdgId();
    const auto flags = gen.StatusFlags();

    const bool isFirstCopy = (flags & FIRST_COPY_BIT) != 0;
    const bool isLastCopy = (flags & LAST_COPY_BIT) != 0;

    if (pdg == 6) { // top
      if (isFirstCopy) {
        assert(FirstCopyTopIndex == npos &&
               "Multiple first-copy tops found in event");
        FirstCopyTopIndex = idx;
      }
      if (isLastCopy) {
        assert(LastCopyTopIndex == npos &&
               "Multiple last-copy tops found in event");
        LastCopyTopIndex = idx;
      }
    } else if (pdg == -6) { // anti-top
      if (isFirstCopy) {
        assert(FirstCopyAntiTopIndex == npos &&
               "Multiple first-copy antitops found in event");
        FirstCopyAntiTopIndex = idx;
      }
      if (isLastCopy) {
        assert(LastCopyAntiTopIndex == npos &&
               "Multiple last-copy antitops found in event");
        LastCopyAntiTopIndex = idx;
      }
    }
  }

  assert(FirstCopyTopIndex != npos && "No first-copy top found in event");
  assert(FirstCopyAntiTopIndex != npos &&
         "No first-copy antitop found in event");
  assert(LastCopyTopIndex != npos && "No last-copy top found in event");
  assert(LastCopyAntiTopIndex != npos && "No last-copy antitop found in event");

  return {FirstCopyTopIndex, FirstCopyAntiTopIndex, LastCopyTopIndex,
          LastCopyAntiTopIndex};
}