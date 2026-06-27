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
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/TutorialSystematic.yaml", MCSample, DataEra);
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
    const auto top_indices = GetTopAndAntiTopIndices(AllGenViews);
    constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
    if (top_indices[0] != npos && top_indices[1] != npos) {
      const TLorentzVector top = AllGenViews[top_indices[0]].P4();
      const TLorentzVector antiTop = AllGenViews[top_indices[1]].P4();
      float w_toppt = myCorr->GetTopPtReweight(top, antiTop);
      FillHist(this_syst + "/count/w_toppt" + this_syst, 0.5, 1., 2, 0., 2.);
      FillHist(this_syst + "/count/w_toppt" + this_syst, 1.5, w_toppt, 2, 0., 2.);
    }
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
  MyCorrection::variation jer_variation = MyCorrection::variation::nom;
  if (this_syst.Contains("JESTotal")) {
    ApplyJetScaleVariation(AllJetViews, "total");
    if (this_syst.Contains("Up")) {
      jes_variation = MyCorrection::variation::up;
    } else if (this_syst.Contains("Down")) {
      jes_variation = MyCorrection::variation::down;
    }
  } else if (this_syst.Contains("JER")) {
    if (this_syst.Contains("Up")) {
      jer_variation = MyCorrection::variation::up;
    } else if (this_syst.Contains("Down")) {
      jer_variation = MyCorrection::variation::down;
    }
  }
  auto jet_id = apply_pu_id ? Jet::JetID::PUID_LOOSE : Jet::JetID::TIGHT;
  std::vector<size_t> SelectedJetIndices = SelectJetIndices(AllJetViews, jet_id, 0., 5.191, jes_variation, jer_variation);
  RVec<Jet> jets = MaterializeJets(AllJetViews, SelectedJetIndices, jes_variation, jer_variation);
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
  float leading_btagged_jet_pt_before_reg = -1.f;
  float leading_btagged_jet_pt_after_reg = -1.f;


  // Loop through the jets and apply corrections to b-tagged ones
  for (auto& jet : jets) {
    // Get the b-tagging discriminator score
    double this_discr = jet.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    
    // Check if the jet is b-tagged
    if (this_discr > btag_wp_cut) {

        if (leading_btagged_jet_pt_before_reg < 0.f) {
          leading_btagged_jet_pt_before_reg = static_cast<float>(jet.Pt());
        }
        double current_pt  = jet.Pt();
        double raw_pt = jet.GetRawPt();
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
        if (this_syst.Contains("BFragmentation")) {
          double bfrag_scale = this_syst.Contains("Up") ? 1.01 : 0.99;
          modified_pt *= bfrag_scale;
          modified_m *= bfrag_scale;
          current_pt *= bfrag_scale;
          current_m *= bfrag_scale;
        }

        if(correct_b_jet_pt){
          jet.SetPtEtaPhiM(modified_pt, current_eta, current_phi, modified_m);
        }
        else if (this_syst.Contains("BFragmentation")) {
          jet.SetPtEtaPhiM(current_pt, current_eta, current_phi, current_m);
        }
        if (leading_btagged_jet_pt_after_reg < 0.f) {
          leading_btagged_jet_pt_after_reg = static_cast<float>(jet.Pt());
        }

    }
    else{
      if(use_UParT_JEC){
        double current_pt  = jet.Pt();
        double raw_pt = jet.GetRawPt();
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

  //==== Event Weight and Systematic Weight Map
  float base_weight = 1.;
  unordered_map<std::string, float> weight_map;

  if (IsDATA) {
    weight_map["Central"] = 1.f;
  } else {
    base_weight *= MCweight();
    base_weight *= ev.GetTriggerLumi("Full");

    if (MCSample.Contains("powheg") && MCSample.Contains("TT")) {
      base_weight *= 1.2360; // top_pt_reweight normalization factor
      if (genTtbarId % 100 >= 51 && genTtbarId % 100 <= 55) {
        base_weight *= 1.36;
      } else if (genTtbarId % 100 >= 41 && genTtbarId % 100 <= 45) {
        base_weight *= 1.11;
      }
    }

    // Assign weight functions to SystematicHelper
    unordered_map<std::string, std::variant<std::function<float(MyCorrection::variation, TString)>, std::function<float()>>> weight_function_map;

    std::function<float(MyCorrection::variation, TString)> central_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)syst;
          (void)source;
          return 1.f;
        };
    weight_function_map["Central"] = central_lambda;

    std::function<float(MyCorrection::variation, TString)> muon_id_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          return myCorr->GetMuonIDSF(this_muon_id_sf_key, muons, syst);
        };
    weight_function_map["Muon_ID"] = muon_id_lambda;

    std::function<float(MyCorrection::variation, TString)> muon_iso_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          if (use_pog_mva_tight_muon_id) {
            return 1.f;
          } else {
            return myCorr->GetMuonIDSF(this_muon_iso_sf_key, muons, syst);
          }
        };
    weight_function_map["Muon_Iso"] = muon_iso_lambda;

    std::function<float(MyCorrection::variation, TString)> muon_trig_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          if (use_pog_tight_muon_id) {
            return myCorr->GetMuonTriggerSF(this_muon_trig_sf_key, muons, syst);
          }
          return 1.f;
        };
    weight_function_map["Muon_Trig"] = muon_trig_lambda;

    std::function<float(MyCorrection::variation, TString)> pileup_lambda =
        [&](MyCorrection::variation syst, TString source) {
          return myCorr->GetPUWeight(ev.nTrueInt(), syst, source);
        };
    weight_function_map["Pileup"] = pileup_lambda;

    std::function<float(MyCorrection::variation, TString)> btag_lambda =
        [&](MyCorrection::variation syst, TString source) {
          RVec<Jet> jets_2p5 = RVec<Jet>();
          for (auto& jet : jets) {
            if (abs(jet.Eta()) < 2.499 && jet.Pt() > 20.) {
              jets_2p5.push_back(jet);
            }
          }
          return myCorr->GetBTaggingSF(jets_2p5,
              JetTagging::JetFlavTagger::ParT,
              JetTagging::JetFlavTaggerWP::Medium,
              JetTagging::JetTaggingSFMethod::comb,
              syst, source
          );
        };
    weight_function_map["BTagSF"] = btag_lambda;

    std::function<float()> toppt_lambda = [&]() {
      if (MCSample.Contains("powheg") && MCSample.Contains("TT")) {
        const auto top_indices = GetTopAndAntiTopIndices(AllGenViews);
        constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
        if (top_indices[0] == npos || top_indices[1] == npos) return 1.f;
        const TLorentzVector top = AllGenViews[top_indices[0]].P4();
        const TLorentzVector antiTop = AllGenViews[top_indices[1]].P4();
        return myCorr->GetTopPtReweight(top, antiTop);
      }
      return 1.f;
    };
    weight_function_map["Top_Pt_Reweight"] = toppt_lambda;

    systHelper->assignWeightFunctionMap(weight_function_map);
    weight_map = systHelper->calculateWeight();
  }

  unordered_map<int, int> matched_genjet_idx = GenJetMatching(jets, MaterializeGenJets(AllGenJetViews), Rho_fixedGridRhoFastjetAll);
  bool isPileupJet = false;
  for (auto &[reco_idx, gen_idx] : matched_genjet_idx) {
    if (gen_idx == -999) {
      isPileupJet = true;
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

  std::vector<size_t> SelectedJetIndices2 = SelectJetIndices(AllJetViews, jet_id, 40., 2.4, jes_variation, jer_variation);
  RVec<Jet> jets2 = MaterializeJets(AllJetViews, SelectedJetIndices2, jes_variation, jer_variation);
  jets2 = JetsVetoLeptonInside(jets2, electrons, muons, 0.3);

  for (const auto &w : weight_map) {
    TString systName = w.first;
    float weight = base_weight * w.second;

    if (leading_btagged_jet_pt_before_reg >= 0.f) {
      FillHist(this_syst + "/baseLineCut/btagged_jet_pt0_" + systName, leading_btagged_jet_pt_before_reg, weight, 80, 0., 400.);
    }
    if (leading_btagged_jet_pt_after_reg >= 0.f) {
      FillHist(this_syst + "/baseLineCut/btagged_RegCorr_jet_pt0_" + systName, leading_btagged_jet_pt_after_reg, weight, 80, 0., 400.);
    }

    FillHist(this_syst + "/baseLineCut/muon_pt0_" + systName, muon_pt0, weight, 80, 0., 400.);
    FillHist(this_syst + "/baseLineCut/muon_eta0_" + systName, muon_eta0, weight, 40, -2.4, 2.4);
    FillHist(this_syst + "/baseLineCut/njets_" + systName, njets, weight, 10, 0., 10.);
    FillHist(this_syst + "/baseLineCut/njets2_" + systName, float(jets2.size()), weight, 10, 0., 10.);
    FillHist(this_syst + "/baseLineCut/jet_pt0_" + systName, jet_pt0, weight, 80, 0., 400.);
    FillHist(this_syst + "/baseLineCut/jet_eta0_" + systName, jet_eta0, weight, 40, -5.2, 5.2);
    FillHist(this_syst + "/baseLineCut/MET_pt_" + systName, MET_pt, weight, 40, 0., 200.);
    FillHist(this_syst + "/baseLineCut/MET_phi_" + systName, MET_phi, weight, 40, -3.14, 3.14);

    if (!IsDATA && draw_include_pu_jets) {
      if(isPileupJet){
        FillHist("Pileup/" + this_syst + "/baseLineCut/muon_pt0_" + systName, muon_pt0, weight, 80, 0., 400.);
        FillHist("Pileup/" + this_syst + "/baseLineCut/muon_eta0_" + systName, muon_eta0, weight, 40, -2.4, 2.4);
        FillHist("Pileup/" + this_syst + "/baseLineCut/njets_" + systName, njets, weight, 10, 0., 10.);
        FillHist("Pileup/" + this_syst + "/baseLineCut/njets2_" + systName, float(jets2.size()), weight, 10, 0., 10.);
        FillHist("Pileup/" + this_syst + "/baseLineCut/jet_pt0_" + systName, jet_pt0, weight, 80, 0., 400.);
        FillHist("Pileup/" + this_syst + "/baseLineCut/jet_eta0_" + systName, jet_eta0, weight, 40, -2.4, 2.4);
        FillHist("Pileup/" + this_syst + "/baseLineCut/MET_pt_" + systName, MET_pt, weight, 40, 0., 200.);
        FillHist("Pileup/" + this_syst + "/baseLineCut/MET_phi_" + systName, MET_phi, weight, 40, -3.14, 3.14);
      }
      else{
        FillHist("noPileup/" + this_syst + "/baseLineCut/muon_pt0_" + systName, muon_pt0, weight, 80, 0., 400.);
        FillHist("noPileup/" + this_syst + "/baseLineCut/muon_eta0_" + systName, muon_eta0, weight, 40, -2.4, 2.4);
        FillHist("noPileup/" + this_syst + "/baseLineCut/njets_" + systName, njets, weight, 10, 0., 10.);
        FillHist("noPileup/" + this_syst + "/baseLineCut/njets2_" + systName, float(jets2.size()), weight, 10, 0., 10.);
        FillHist("noPileup/" + this_syst + "/baseLineCut/jet_pt0_" + systName, jet_pt0, weight, 80, 0., 400.);
        FillHist("noPileup/" + this_syst + "/baseLineCut/jet_eta0_" + systName, jet_eta0, weight, 40, -2.4, 2.4);
        FillHist("noPileup/" + this_syst + "/baseLineCut/MET_pt_" + systName, MET_pt, weight, 40, 0., 200.);
        FillHist("noPileup/" + this_syst + "/baseLineCut/MET_phi_" + systName, MET_phi, weight, 40, -3.14, 3.14);
      }
    }
  }

  //==== Take W and top-b candidates from the leading five jets in pT
  std::vector<unsigned int> top_b_jet_candidates;
  std::vector<unsigned int> had_W_candidates;

  // Check up to the leading 5 jets.
  for(unsigned int ij(0); ij < 5; ij++){
    if(btag_vector.at(ij)){
      // We still need exactly 2 b-jets for the ttbar system
      if(top_b_jet_candidates.size() < 2) top_b_jet_candidates.push_back(ij);
      // allow soft b-tagged jet for mistag c (W -> cs)
      else if(top_b_jet_candidates.size() >= 2) had_W_candidates.push_back(ij);
    }
    else{
      // Keep up to three hadronic W jet candidates.
      if(had_W_candidates.size() < 3) had_W_candidates.push_back(ij);
    }
  }

  // Require at least two W candidates and two top-b candidates.
  if(had_W_candidates.size() < 2 || top_b_jet_candidates.size() < 2) return;
  FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 4.5, 1., 6, 0., 6.);

  //==== Combinatorics
  // Up to six combinations: C(nW, 2) W-jet pairs times two b-jet assignments, with nW <= 3.
  Tutorial_reco_tt::ttCombinatoric combinatorics[6];
  int comb_idx = 0;

  // Loop over all ways to choose 2 W jets from the available candidates.
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
  for(int i = 1; i < comb_idx; i++){
    if(combinatorics[i].best_chi2 < best_combinatoric->best_chi2){
      best_combinatoric = &combinatorics[i];
    }
  }

  for (const auto &w : weight_map) {
    TString systName = w.first;
    float weight = base_weight * w.second;

    FillHist(this_syst + "/noChi2Cut/had_W_mass_" + systName, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
    FillHist(this_syst + "/noChi2Cut/had_top_mass_" + systName, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
    FillHist(this_syst + "/noChi2Cut/lep_W_mass_" + systName, best_combinatoric->best_lep_W_mass, weight, 40, 0., 200.);
    FillHist(this_syst + "/noChi2Cut/lep_top_mass_" + systName, best_combinatoric->best_lep_top_mass, weight, 80, 0., 400.);
    FillHist(this_syst + "/noChi2Cut/chi2_" + systName, best_combinatoric->best_chi2, weight, 50, 0., 100000.);
    FillHist(this_syst + "/noChi2Cut/njets_" + systName, njets, weight, 10, 0., 10.);

    if(best_combinatoric->best_chi2 < 2e3) {
      FillHist(this_syst + "/Chi2Cut/had_W_mass_" + systName, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
      FillHist(this_syst + "/Chi2Cut/had_top_mass_" + systName, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
      FillHist(this_syst + "/Chi2Cut/lep_W_mass_" + systName, best_combinatoric->best_lep_W_mass, weight, 40, 0., 200.);
      FillHist(this_syst + "/Chi2Cut/lep_top_mass_" + systName, best_combinatoric->best_lep_top_mass, weight, 80, 0., 400.);
      FillHist(this_syst + "/Chi2Cut/chi2_" + systName, best_combinatoric->best_chi2, weight, 50, 0., 2000.);

      FillHist(this_syst + "/Chi2Cut/muon_pt0_" + systName, muon_pt0, weight, 80, 0., 400.);
      FillHist(this_syst + "/Chi2Cut/muon_eta0_" + systName, muon_eta0, weight, 40, -2.4, 2.4);
      FillHist(this_syst + "/Chi2Cut/njets_" + systName, njets, weight, 10, 0., 10.);
      FillHist(this_syst + "/Chi2Cut/njets2_" + systName, float(jets2.size()), weight, 10, 0., 10.);
      FillHist(this_syst + "/Chi2Cut/jet_pt0_" + systName, jet_pt0, weight, 80, 0., 400.);
      FillHist(this_syst + "/Chi2Cut/jet_eta0_" + systName, jet_eta0, weight, 40, -5.2, 5.2);
      FillHist(this_syst + "/Chi2Cut/MET_pt_" + systName, MET_pt, weight, 40, 0., 200.);
      FillHist(this_syst + "/Chi2Cut/MET_phi_" + systName, MET_phi, weight, 40, -3.14, 3.14);

      if (leading_btagged_jet_pt_before_reg >= 0.f) {
        FillHist(this_syst + "/Chi2Cut/btagged_jet_pt0_" + systName, leading_btagged_jet_pt_before_reg, weight, 80, 0., 400.);
      }
      if (leading_btagged_jet_pt_after_reg >= 0.f) {
        FillHist(this_syst + "/Chi2Cut/btagged_RegCorr_jet_pt0_" + systName, leading_btagged_jet_pt_after_reg, weight, 80, 0., 400.);
      }
    }
  }

  if(best_combinatoric->best_chi2 < 2e3) {
    FillHist(this_syst + "/cutflow/cutflow_" + this_syst, 5.5, 1., 6, 0., 6.);
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

  auto warn_duplicate = [](const char *message, size_t current_idx,
                           size_t duplicate_idx) {
    cerr << "[Tutorial_reco_tt::GetTopAndAntiTopIndices] " << message
         << ": keeping index " << current_idx << ", ignoring index "
         << duplicate_idx << endl;
  };

  auto warn_missing = [](const char *message) {
    cerr << "[Tutorial_reco_tt::GetTopAndAntiTopIndices] " << message << endl;
  };

  for (size_t idx = 0; idx < n; ++idx) {
    const GenView &gen = gens[idx];

    const int pdg = gen.PdgId();
    const auto flags = gen.StatusFlags();

    const bool isFirstCopy = (flags & FIRST_COPY_BIT) != 0;
    const bool isLastCopy = (flags & LAST_COPY_BIT) != 0;

    if (pdg == 6) { // top
      if (isFirstCopy) {
        if (FirstCopyTopIndex == npos) {
          FirstCopyTopIndex = idx;
        } else {
          warn_duplicate("Multiple first-copy tops found in event",
                         FirstCopyTopIndex, idx);
          assert(FirstCopyTopIndex == npos &&
                 "Multiple first-copy tops found in event");
        }
      }
      if (isLastCopy) {
        if (LastCopyTopIndex == npos) {
          LastCopyTopIndex = idx;
        } else {
          warn_duplicate("Multiple last-copy tops found in event",
                         LastCopyTopIndex, idx);
          assert(LastCopyTopIndex == npos &&
                 "Multiple last-copy tops found in event");
        }
      }
    } else if (pdg == -6) { // anti-top
      if (isFirstCopy) {
        if (FirstCopyAntiTopIndex == npos) {
          FirstCopyAntiTopIndex = idx;
        } else {
          warn_duplicate("Multiple first-copy antitops found in event",
                         FirstCopyAntiTopIndex, idx);
          assert(FirstCopyAntiTopIndex == npos &&
                 "Multiple first-copy antitops found in event");
        }
      }
      if (isLastCopy) {
        if (LastCopyAntiTopIndex == npos) {
          LastCopyAntiTopIndex = idx;
        } else {
          warn_duplicate("Multiple last-copy antitops found in event",
                         LastCopyAntiTopIndex, idx);
          assert(LastCopyAntiTopIndex == npos &&
                 "Multiple last-copy antitops found in event");
        }
      }
    }
  }

  if (FirstCopyTopIndex == npos) {
    warn_missing("No first-copy top found in event");
    assert(FirstCopyTopIndex != npos && "No first-copy top found in event");
  }
  if (FirstCopyAntiTopIndex == npos) {
    warn_missing("No first-copy antitop found in event");
    assert(FirstCopyAntiTopIndex != npos &&
           "No first-copy antitop found in event");
  }
  if (LastCopyTopIndex == npos) {
    warn_missing("No last-copy top found in event");
    assert(LastCopyTopIndex != npos && "No last-copy top found in event");
  }
  if (LastCopyAntiTopIndex == npos) {
    warn_missing("No last-copy antitop found in event");
    assert(LastCopyAntiTopIndex != npos &&
           "No last-copy antitop found in event");
  }

  return {FirstCopyTopIndex, FirstCopyAntiTopIndex, LastCopyTopIndex,
          LastCopyAntiTopIndex};
}
