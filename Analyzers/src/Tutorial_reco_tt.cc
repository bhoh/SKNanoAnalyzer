#include "Tutorial_reco_tt.h"

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
  bool use_UParT_JEC = false;

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

  if (SelectedMuonIndices.size() + SelectedElectronIndices.size() != 1) return;

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
  std::vector<size_t> SelectedJetIndices = SelectJetIndices(AllJetViews, jet_id, 0., 5.191, jes_variation, MyCorrection::variation::nom);
  RVec<Jet> jets = MaterializeJets(AllJetViews, SelectedJetIndices, jes_variation, MyCorrection::variation::nom);
  jets = JetsVetoLeptonInside(jets, electrons, muons, 0.3);
  //==== Sorting
  sort(muons.begin(), muons.end(), PtComparing);
  sort(jets.begin(), jets.end(), PtComparing);

  //==== Event selections
  if (muons.size() != 1) return;
  if (electrons.size() != 0) return;
  if (muons.at(0).Pt() <= TriggerSafePtCut) return;
  if (jets.size() < 6) return;
  //if (METv.Pt() <= 20) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 1.5, 1., 6, 0., 6.);

  if (!PassJetVetoMap(AllJetViews, AllMuonViews, "jetvetomap_fpix")) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 2.5, 1., 6, 0., 6.);
  //==== B-Tagging (DeepJet Medium WP example)
  int NBJets = 0;
  int njets_pt30_non_btagged = 0;
  float btag_wp_cut = myCorr->GetBTaggingWP();


  // Loop through the jets and apply corrections to b-tagged ones
  for (auto& jet : jets) {
    // Get the b-tagging discriminator score
    double this_discr = jet.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    
    // Check if the jet is b-tagged
    if (this_discr > btag_wp_cut) {
        
        double current_pt  = jet.Pt();
        double raw_pt = jet.GetRawPt();
        double unsmeared_pt = jet.GetUnsmearedP4().Pt();
        double current_eta = jet.Eta();
        double current_phi = jet.Phi();
        double current_m   = jet.M();

        double L1L2L3Res = myCorr->GetJESSF(jet.GetArea(), current_eta, raw_pt, current_phi, Rho_fixedGridRhoFastjetAll, ev.run());
        double L2L3Residual = myCorr->GetJESSF(jet.GetArea(), current_eta, raw_pt, current_phi, Rho_fixedGridRhoFastjetAll, ev.run(), "L2L3Residual");

        double UParTAK4RegPtRawCorr = jet.UParTAK4RegPtRawCorr();
        double UParTAK4RegPtRawCorrNeutrino = jet.UParTAK4RegPtRawCorrNeutrino();
        double UParT_ratio = UParTAK4RegPtRawCorrNeutrino / UParTAK4RegPtRawCorr;
        FillHist(this_syst + "/corrections/UParTAK4RegPtRawCorr_" + this_syst, UParTAK4RegPtRawCorr, 1., 80, 0., 2);
        FillHist(this_syst + "/corrections/UParTAK4RegPtRawCorrNeutrino_" + this_syst, UParTAK4RegPtRawCorrNeutrino, 1., 80, 0., 2);
        FillHist(this_syst + "/corrections/UParT_ratio_" + this_syst, UParT_ratio, 1., 80, 0., 2);
        // Calculate your modified pT
        double modified_pt = current_pt;
        double modified_m  = current_m;
        if (use_UParT_JEC) {
          modified_pt = current_pt * (1/L1L2L3Res) * UParTAK4RegPtRawCorrNeutrino * L2L3Residual;
          modified_m  = current_m * modified_pt / current_pt;
        }
        else{
          modified_pt = current_pt * UParT_ratio;
          modified_m  = current_m  * UParT_ratio;
        }

        // Update the LorentzVector with the new pT
        if(correct_b_jet_pt){
          jet.SetPtEtaPhiM(modified_pt, current_eta, current_phi, modified_m);
        }
        
    }
    else{
      if(use_UParT_JEC){
        double current_pt  = jet.Pt();
        double raw_pt = jet.GetRawPt();
        double unsmeared_pt = jet.GetUnsmearedP4().Pt();
        double current_eta = jet.Eta();
        double current_phi = jet.Phi();
        double current_m   = jet.M();

        double L1L2L3Res = myCorr->GetJESSF(jet.GetArea(), current_eta, raw_pt, current_phi, Rho_fixedGridRhoFastjetAll, ev.run());
        double L2L3Residual = myCorr->GetJESSF(jet.GetArea(), current_eta, raw_pt, current_phi, Rho_fixedGridRhoFastjetAll, ev.run(), "L2L3Residual");

        double UParTAK4RegPtRawCorr = jet.UParTAK4RegPtRawCorr();
        double modified_pt = current_pt * (1/L1L2L3Res) * UParTAK4RegPtRawCorr * L2L3Residual;
        double modified_m  = current_m * modified_pt / current_pt;
        jet.SetPtEtaPhiM(modified_pt, current_eta, current_phi, modified_m);
      }
    }
  }
  std::vector<bool> btag_vector;

  for (unsigned int ij = 0; ij < jets.size(); ij++) {
    double this_discr = jets.at(ij).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    if (this_discr > btag_wp_cut) {
      NBJets++;
      btag_vector.push_back(true);
      if (jets.at(ij).Pt() > 25){
        njets_pt30_non_btagged++;
      }

    } else {
      btag_vector.push_back(false);
      if (jets.at(ij).Pt() > 30 && abs(jets.at(ij).Eta()) < 2.5){
        njets_pt30_non_btagged++;
      }
      else if (jets.at(ij).Pt() > 30 && abs(jets.at(ij).Eta()) >= 2.5){
        njets_pt30_non_btagged++;
      }
    }
  }





  if (NBJets != 3) return;
  if (njets_pt30_non_btagged < 6) return;
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
    RVec<Jet> jets_2p5 = RVec<Jet>();
    for (auto& jet : jets) {
      if (abs(jet.Eta()) < 2.499 && jet.Pt() > 20.) {
        jets_2p5.push_back(jet);
      }
    }
    float btag_sf = myCorr->GetBTaggingSF(jets_2p5, 
            JetTagging::JetFlavTagger::ParT, 
            JetTagging::JetFlavTaggerWP::Medium,
            JetTagging::JetTaggingSFMethod::comb,
            btag_jes_variation, btag_source
        );
    weight *= btag_sf;

    if (MCSample.Contains("powheg") && MCSample.Contains("TT")) {
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
  float jet_pt0 = jets.at(0).Pt();
  float jet_eta0 = jets.at(0).Eta();
  float njets = njets_pt30_non_btagged;
  float MET_pt = METv.Pt();
  float MET_phi = METv.Phi();
  FillHist(this_syst + "/baseLineCut/muon_pt0_" + this_syst, muon_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/baseLineCut/muon_eta0_" + this_syst, muon_eta0, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/baseLineCut/njets_" + this_syst, njets, weight, 10, 0., 10.);


  std::vector<size_t> SelectedJetIndices2 = SelectJetIndices(AllJetViews, jet_id, 40., 2.4, jes_variation, MyCorrection::variation::nom);
  RVec<Jet> jets2 = MaterializeJets(AllJetViews, SelectedJetIndices2, jes_variation, MyCorrection::variation::nom);
  jets2 = JetsVetoLeptonInside(jets2, electrons, muons, 0.3);
  FillHist(this_syst + "/baseLineCut/njets2_" + this_syst, float(jets2.size()), weight, 10, 0., 10.);
  FillHist(this_syst + "/baseLineCut/jet_pt0_" + this_syst, jet_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/baseLineCut/jet_eta0_" + this_syst, jet_eta0, weight, 40, -5.2, 5.2);
  FillHist(this_syst + "/baseLineCut/MET_pt_" + this_syst, MET_pt, weight, 40, 0., 200.);
  FillHist(this_syst + "/baseLineCut/MET_phi_" + this_syst, MET_phi, weight, 40, -3.14, 3.14);

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


//==== Take leading five jets in pT
  std::vector<unsigned int> top_b_jet_candidates;
  std::vector<unsigned int> had_W_candidates;

  // Check up to the leading 5 jets
  for(unsigned int ij(0); ij < 5; ij++){
    if(btag_vector.at(ij)){
      // We still need exactly 2 b-jets for the ttbar system
      if(top_b_jet_candidates.size() < 2) top_b_jet_candidates.push_back(ij);
      // allow soft b-tagged jet for mistag c (W -> cs)
      else if(top_b_jet_candidates.size() >= 2) had_W_candidates.push_back(ij);
    }
    else{
      // Expanded to keep up to 3 hadronic W jet candidates
      if(had_W_candidates.size() < 3) had_W_candidates.push_back(ij);
    }
  }

  // Require at least 3 had_W candidates and 2 top_b candidates
  if(had_W_candidates.size() < 2 || top_b_jet_candidates.size() < 2) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 4.5, 1., 6, 0., 6.);

  //==== Combinatorics
  // Array to hold the 6 combinations: 3C2 (W jets) * 2P2 (b jets) = 6
  Tutorial_reco_tt::ttCombinatoric combinatorics[6];
  int comb_idx = 0;

  // Loop over the 3 ways to choose 2 W jets out of the 3 candidates: (0,1), (0,2), (1,2)
  for(unsigned int w1 = 0; w1 < had_W_candidates.size()-1; w1++){
    for(unsigned int w2 = w1 + 1; w2 < had_W_candidates.size(); w2++){
      // Loop over the 2 ways to assign the b-jets (hadronic vs leptonic)
      for(unsigned int b_idx = 0; b_idx < 2; b_idx++){
        combinatorics[comb_idx].lepton            = &(muons.at(0));
        combinatorics[comb_idx].jets              = &jets;
        combinatorics[comb_idx].met               = &METv;
        combinatorics[comb_idx].had_W_jet_idx_1   = had_W_candidates.at(w1);
        combinatorics[comb_idx].had_W_jet_idx_2   = had_W_candidates.at(w2);
        combinatorics[comb_idx].had_top_b_jet_idx = top_b_jet_candidates.at(b_idx);
        combinatorics[comb_idx].lep_top_b_jet_idx = top_b_jet_candidates.at(1 - b_idx);

        // Evaluate Chi2 for this specific combination
        this->EvalChi2(combinatorics[comb_idx]);
        comb_idx++;
      }
    }
  }

  //==== Find the best combinatoric based on the minimum chi2
  Tutorial_reco_tt::ttCombinatoric* best_combinatoric = &combinatorics[0];
  for(unsigned int i = 1; i < 6; i++){
    if(combinatorics[i].best_chi2 < best_combinatoric->best_chi2){
      best_combinatoric = &combinatorics[i];
    }
  }

  FillHist(this_syst + "/noChi2Cut/had_W_mass_" + this_syst, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
  FillHist(this_syst + "/noChi2Cut/had_top_mass_" + this_syst, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
  FillHist(this_syst + "/noChi2Cut/lep_W_mass_" + this_syst, best_combinatoric->best_lep_W_mass, weight, 40, 0., 200.);
  FillHist(this_syst + "/noChi2Cut/lep_top_mass_" + this_syst, best_combinatoric->best_lep_top_mass, weight, 80, 0., 400.);
  FillHist(this_syst + "/noChi2Cut/chi2_" + this_syst, best_combinatoric->best_chi2, weight, 50, 0., 100000.);
  FillHist(this_syst + "/noChi2Cut/njets_" + this_syst, njets, weight, 10, 0., 10.);

  if(best_combinatoric->best_chi2 >= 2e3) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 5.5, 1., 6, 0., 6.);
  
  FillHist(this_syst + "/Chi2Cut/had_W_mass_" + this_syst, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
  FillHist(this_syst + "/Chi2Cut/had_top_mass_" + this_syst, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
  FillHist(this_syst + "/Chi2Cut/lep_W_mass_" + this_syst, best_combinatoric->best_lep_W_mass, weight, 40, 0., 200.);
  FillHist(this_syst + "/Chi2Cut/lep_top_mass_" + this_syst, best_combinatoric->best_lep_top_mass, weight, 80, 0., 400.);
  FillHist(this_syst + "/Chi2Cut/chi2_" + this_syst, best_combinatoric->best_chi2, weight, 50, 0., 2000.);

  FillHist(this_syst + "/Chi2Cut/muon_pt0_" + this_syst, muon_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/Chi2Cut/muon_eta0_" + this_syst, muon_eta0, weight, 40, -2.4, 2.4);
  FillHist(this_syst + "/Chi2Cut/njets_" + this_syst, njets, weight, 10, 0., 10.);
  FillHist(this_syst + "/Chi2Cut/njets2_" + this_syst, float(jets2.size()), weight, 10, 0., 10.);
  FillHist(this_syst + "/Chi2Cut/jet_pt0_" + this_syst, jet_pt0, weight, 80, 0., 400.);
  FillHist(this_syst + "/Chi2Cut/jet_eta0_" + this_syst, jet_eta0, weight, 40, -5.2, 5.2);
  FillHist(this_syst + "/Chi2Cut/MET_pt_" + this_syst, MET_pt, weight, 40, 0., 200.);
  FillHist(this_syst + "/Chi2Cut/MET_phi_" + this_syst, MET_phi, weight, 40, -3.14, 3.14);

  for (auto& jet : jets) {
    // Get the b-tagging discriminator score
    double this_discr = jet.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    
    // Check if the jet is b-tagged
    if (this_discr > btag_wp_cut) {
      FillHist(this_syst + "/Chi2Cut/btagged_RegCorr_jet_pt0_" + this_syst, float(jet.Pt()), weight, 80, 0., 400.);
      break;
    }
  }

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