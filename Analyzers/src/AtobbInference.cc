#include "AtobbInference.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <limits>
#include <functional>
#include <numeric>

AtobbInference::AtobbInference() :
  const_top_mass(172.5),
  const_top_width(1.5),
  const_w_mass(80.4),
  const_w_width(2.085)
{
}

AtobbInference::~AtobbInference() {}

void AtobbInference::initializeAnalyzer() {

  MuonIDs = { Muon::MuonID::POG_TIGHT, Muon::MuonID::POG_MEDIUM_PROMPT, Muon::MuonID::POG_MVA_MU_TIGHT };
  MuonIDISOSFKeys = { "NUM_TightID_DEN_TrackerMuons", "NUM_TightPFIso_DEN_TightID", "NUM_MediumPromptID_DEN_TrackerMuons", "NUM_TightPFIso_DEN_MediumPromptID", "NUM_TightMvaMuID_DEN_TrackerMuons" };
  ElectronIDs = { Electron::ElectronID::POG_MVAISO_WP80 };
  ElectronIDSFKeys = { "wp80iso" };
  ElectronTriggerSFKeys = { "HLT_SF_Ele30_MVAiso80ID" };

  EleTriggerName = "HLT_Ele30_WPTight_Gsf";
  EleTriggerSafePtCut = 32.;

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
    cerr << "[AtobbInference::initializeAnalyzer] DataEra is not set properly"
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

  // Load ONNX BDT models for each mA mass point
  ma_points = {"MA18", "MA25", "MA36", "MA50", "MA70", "MA90", "MA110", "MA125"};
  std::string model_dir = "/home/bhoh/kkimeun/MachineLearning/models/";

  for (const auto& ma : ma_points) {
    std::string model_path = model_dir + "xgboost_atobb_MHc130_" + std::string(ma.Data()) + "_ttlj.onnx";
    try {
      ml_helpers[ma] = std::make_unique<MLHelper>(model_path, MLHelper::ModelType::ONNX);
      cout << "[AtobbInference::initializeAnalyzer] Loaded ONNX model for " << ma << " from " << model_path << endl;
    } catch (const std::exception& e) {
      cerr << "[AtobbInference::initializeAnalyzer] Failed to load ONNX model for " << ma << ": " << e.what() << endl;
    }
  }

}

void AtobbInference::executeEvent() {

  AllMuonViews = GetAllMuonViews();
  AllElectronViews = GetAllElectronViews();
  AllJetViews = GetAllJetViews();
  AllGenViews = GetAllGenViews();
  AllGenJetViews = GetAllGenJetViews();

  ev = GetEvent();

  weight_Prefire = 1.; 

  for (const auto &syst_dummy : *systHelper) {
    executeEventFromParameter();
  }
}

void AtobbInference::executeEventFromParameter() {

  bool correct_b_jet_pt = true;
  bool apply_pu_id = false;
  bool use_pog_tight_muon_id = true;
  bool use_pog_mva_tight_muon_id = false;
  bool eval_top_pt_reweight_normalization = false;
  bool use_UParT_JEC = false;

  const TString this_syst = systHelper->getCurrentSysName();
  if (IsDATA && this_syst != "Central") return;

  Muon::MuonID this_muon_id = MuonIDs[0];
  TString this_muon_id_sf_key = MuonIDISOSFKeys[0];
  TString this_muon_iso_sf_key = MuonIDISOSFKeys[1];
  TString this_muon_trig_sf_key = "";
  Electron::ElectronID this_electron_id = ElectronIDs[0];
  TString this_electron_id_sf_key = ElectronIDSFKeys[0];
  TString this_electron_trig_sf_key = ElectronTriggerSFKeys[0];

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

  //==== MET Filter & Trigger
  if (!PassMetFilter(AllJetViews, ev)) return;

  const bool use_met_unclustered = this_syst.Contains("UnclusteredEnergy");
  MyCorrection::variation met_variation = MyCorrection::variation::nom;
  if (use_met_unclustered) {
    if (this_syst.Contains("Up")) {
      met_variation = MyCorrection::variation::up;
    } else if (this_syst.Contains("Down")) {
      met_variation = MyCorrection::variation::down;
    }
  }
  Particle METv = use_met_unclustered
                      ? ev.GetMETVector(Event::MET_Type::PUPPI, met_variation,
                                        Event::MET_Syst::UE)
                      : ev.GetMETVector(Event::MET_Type::PUPPI);

  //==== Lepton Selection
  std::vector<size_t> SelectedMuonIndices_id_only = SelectMuonIndices(AllMuonViews, this_muon_id, 15., 2.4);
  std::vector<size_t> SelectedMuonIndices = {};
  if (use_pog_mva_tight_muon_id){
    SelectedMuonIndices = SelectedMuonIndices_id_only;
  }
  else{
    SelectedMuonIndices = SelectMuonIndices(AllMuonViews, SelectedMuonIndices_id_only, Muon::MuonID::POG_PFISO_TIGHT, 15., 2.4);
  }
  std::vector<size_t> SelectedElectronVetoIndices = SelectElectronIndices(AllElectronViews, Electron::ElectronID::POG_MVAISO_WP80, 15., 2.5);
  std::vector<size_t> SelectedElectronIndices = SelectElectronIndices(AllElectronViews, SelectedElectronVetoIndices, this_electron_id, 15., 2.5);

  if (SelectedMuonIndices.size() + SelectedElectronVetoIndices.size() != 1) return;
  if (SelectedMuonIndices.size() == 1 && SelectedElectronVetoIndices.size() != 0) return;
  if (SelectedElectronIndices.size() == 1 && SelectedElectronVetoIndices.size() != 1) return;

  RVec<Muon> muons = MaterializeMuons(AllMuonViews, SelectedMuonIndices);
  RVec<Electron> electrons = MaterializeElectrons(AllElectronViews, SelectedElectronIndices);

  if (this_syst.Contains("MuonEn")) {
    MyCorrection::variation syst = MyCorrection::variation::nom;
    if (this_syst.Contains("Up")) syst = MyCorrection::variation::up;
    else if (this_syst.Contains("Down")) syst = MyCorrection::variation::down;
    muons = ScaleMuons(muons, syst);
  }

  if (this_syst.Contains("ElectronEn")) {
    MyCorrection::variation syst = MyCorrection::variation::nom;
    if (this_syst.Contains("Up")) syst = MyCorrection::variation::up;
    else if (this_syst.Contains("Down")) syst = MyCorrection::variation::down;
    electrons = ScaleElectrons(ev, electrons, syst);
  } else if (this_syst.Contains("ElectronRes")) {
    MyCorrection::variation syst = MyCorrection::variation::nom;
    if (this_syst.Contains("Up")) syst = MyCorrection::variation::up;
    else if (this_syst.Contains("Down")) syst = MyCorrection::variation::down;
    electrons = SmearElectrons(electrons, syst);
  }

  sort(muons.begin(), muons.end(), PtComparing);
  sort(electrons.begin(), electrons.end(), PtComparing);

  bool is_muon_channel = muons.size() == 1;
  bool is_electron_channel = electrons.size() == 1;
  Lepton* selected_lepton = nullptr;
  TString channel_name = "";

  if (is_muon_channel) {
    if (!(ev.PassTrigger(IsoMuTriggerName))) return;
    if (muons.at(0).Pt() <= TriggerSafePtCut) return;
    selected_lepton = &(muons.at(0));
    channel_name = "Muon";
  } else if (is_electron_channel) {
    if (!(ev.PassTrigger(EleTriggerName))) return;
    if (electrons.at(0).Pt() <= EleTriggerSafePtCut) return;
    selected_lepton = &(electrons.at(0));
    channel_name = "Electron";
  } else {
    return;
  }

  TString sample_name = IsDATA ? "Data" : MCSample;
  if (!IsDATA && (MCSample == "TTLJ_powheg" || MCSample == "TTLL_powheg")) {
    int ttbar_flavor_id = genTtbarId % 100;
    if (ttbar_flavor_id < 0) ttbar_flavor_id += 100;

    if (ttbar_flavor_id >= 53 && ttbar_flavor_id <= 55) {
      sample_name = MCSample + "_bb";
    } else if (ttbar_flavor_id >= 51 && ttbar_flavor_id <= 52) {
      sample_name = MCSample + "_bj";
    } else if (ttbar_flavor_id >= 43 && ttbar_flavor_id <= 45) {
      sample_name = MCSample + "_cc";
    } else if (ttbar_flavor_id >= 41 && ttbar_flavor_id <= 42) {
      sample_name = MCSample + "_cj";
    } else {
      sample_name = MCSample + "_jj";
    }
  }
  const TString hist_sample_prefix = "histo_" + sample_name;
  TString event_hist_name = hist_sample_prefix;
  if (this_syst != "Central") event_hist_name += "_" + this_syst;

  FillHist(channel_name + "/cutflow/" + event_hist_name, 0.5, 1., 6, 0., 6.);

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
  sort(jets.begin(), jets.end(), PtComparing);

  //==== Event selections
  if (jets.size() < 6) return;
  FillHist(channel_name + "/cutflow/" + event_hist_name, 1.5, 1., 6, 0., 6.);

  if (!PassJetVetoMap(AllJetViews, AllMuonViews, "jetvetomap_fpix")) return;
  FillHist(channel_name + "/cutflow/" + event_hist_name, 2.5, 1., 6, 0., 6.);

  //==== B-Tagging (DeepJet Medium WP example)
  int NBJets = 0;
  int njets_pt30_non_btagged = 0;
  float btag_wp_cut = myCorr->GetBTaggingWP();

  for (auto& jet : jets) {
    double this_discr = jet.GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    if (this_discr > btag_wp_cut) {
        double current_pt  = jet.Pt();
        double current_eta = jet.Eta();
        double current_phi = jet.Phi();
        double current_m   = jet.M();

        double UParTAK4RegPtRawCorr = jet.UParTAK4RegPtRawCorr();
        double UParTAK4RegPtRawCorrNeutrino = jet.UParTAK4RegPtRawCorrNeutrino();
        double UParT_ratio = UParTAK4RegPtRawCorrNeutrino / UParTAK4RegPtRawCorr;

        double modified_pt = current_pt * UParT_ratio;
        double modified_m  = current_m  * UParT_ratio;

        if (this_syst.Contains("BFragmentation")) {
          double bfrag_scale = this_syst.Contains("Up") ? 1.01 : 0.99;
          modified_pt *= bfrag_scale;
          modified_m *= bfrag_scale;
          current_pt *= bfrag_scale;
          current_m *= bfrag_scale;
        }

        if (correct_b_jet_pt) {
          jet.SetPtEtaPhiM(modified_pt, current_eta, current_phi, modified_m);
        } else if (this_syst.Contains("BFragmentation")) {
          jet.SetPtEtaPhiM(current_pt, current_eta, current_phi, current_m);
        }
    }
  }

  std::vector<unsigned int> bjet_indices;
  for (unsigned int ij = 0; ij < jets.size(); ij++) {
    double this_discr = jets.at(ij).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    if (this_discr > btag_wp_cut) {
      NBJets++;
      bjet_indices.push_back(ij);
      if (jets.at(ij).Pt() > 25) {
        njets_pt30_non_btagged++;
      }
    } else {
      if (jets.at(ij).Pt() > 30) {
        njets_pt30_non_btagged++;
      }
    }
  }

  // Require >= 3 b-tagged jets per AtobbMLTree specification (instead of NBJets != 3)
  if (NBJets < 3) return;
  if (njets_pt30_non_btagged < 6) return;
  FillHist(channel_name + "/cutflow/" + event_hist_name, 3.5, 1., 6, 0., 6.);

  //==== Event Weight and Systematic Weight Map
  float base_weight = 1.;
  unordered_map<std::string, float> weight_map;

  if (IsDATA) {
    weight_map["Central"] = 1.f;
  } else {
    base_weight *= MCweight();
    base_weight *= ev.GetTriggerLumi("Full");

    if (MCSample.Contains("powheg") && MCSample.Contains("TT")) {
      base_weight *= 1.2360;
      if (genTtbarId % 100 >= 51 && genTtbarId % 100 <= 55) {
        base_weight *= 1.36;
      } else if (genTtbarId % 100 >= 41 && genTtbarId % 100 <= 45) {
        base_weight *= 1.11;
      }
    }

    unordered_map<std::string, std::variant<std::function<float(MyCorrection::variation, TString)>, std::function<float()>>> weight_function_map;

    weight_function_map["Central"] = [&](MyCorrection::variation syst, TString source) { (void)syst; (void)source; return 1.f; };

    weight_function_map["Muon_ID"] = [&](MyCorrection::variation syst, TString source) {
      (void)source;
      if (!is_muon_channel) return 1.f;
      return myCorr->GetMuonIDSF(this_muon_id_sf_key, muons, syst);
    };

    weight_function_map["Muon_Iso"] = [&](MyCorrection::variation syst, TString source) {
      (void)source;
      if (!is_muon_channel || use_pog_mva_tight_muon_id) return 1.f;
      return myCorr->GetMuonIDSF(this_muon_iso_sf_key, muons, syst);
    };

    weight_function_map["Muon_Trig"] = [&](MyCorrection::variation syst, TString source) {
      (void)source;
      if (is_muon_channel && use_pog_tight_muon_id) {
        return myCorr->GetMuonTriggerSF(this_muon_trig_sf_key, muons, syst);
      }
      return 1.f;
    };

    weight_function_map["Electron_ID"] = [&](MyCorrection::variation syst, TString source) {
      (void)source;
      if (!is_electron_channel) return 1.f;
      return myCorr->GetElectronIDSF(this_electron_id_sf_key, electrons, syst);
    };

    weight_function_map["Electron_Reco"] = [&](MyCorrection::variation syst, TString source) {
      (void)source;
      if (!is_electron_channel) return 1.f;
      return myCorr->GetElectronRECOSF(electrons, syst);
    };

    weight_function_map["Electron_Trig"] = [&](MyCorrection::variation syst, TString source) {
      (void)source;
      if (!is_electron_channel) return 1.f;
      const Electron &electron = electrons.at(0);
      return myCorr->GetElectronTriggerSF(this_electron_trig_sf_key, electron.scEta(), electron.Pt(), electron.Phi(), syst);
    };

    weight_function_map["Pileup"] = [&](MyCorrection::variation syst, TString source) {
      return myCorr->GetPUWeight(ev.nTrueInt(), syst, source);
    };

    weight_function_map["BTagSF"] = [&](MyCorrection::variation syst, TString source) {
      RVec<Jet> jets_2p5 = RVec<Jet>();
      for (auto& jet : jets) {
        if (abs(jet.Eta()) < 2.499 && jet.Pt() > 20.) {
          jets_2p5.push_back(jet);
        }
      }
      return myCorr->GetBTaggingSF(jets_2p5, JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium, JetTagging::JetTaggingSFMethod::comb, syst, source);
    };

    weight_function_map["Top_Pt_Reweight"] = [&]() {
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

    systHelper->assignWeightFunctionMap(weight_function_map);
    weight_map = systHelper->calculateWeight();

    if (!use_met_unclustered && this_syst == "Central") {
      auto safe_lhe_pdf_weight = [&](int index) {
        if (!LHEPdfWeight.valid() || nLHEPdfWeight <= index) return 1.f;
        const float weight = LHEPdfWeight[index];
        return std::isfinite(weight) ? weight : 1.f;
      };
      auto safe_lhe_scale_weight =
          [&](MyCorrection::variation muF_syst,
              MyCorrection::variation muR_syst) {
            if (!LHEScaleWeight.valid() || nLHEScaleWeight < 8) return 1.f;

            if (nLHEScaleWeight >= 9) {
              int index = -1;
              if (muF_syst == MyCorrection::variation::down &&
                  muR_syst == MyCorrection::variation::down) {
                index = 0;
              } else if (muF_syst == MyCorrection::variation::nom &&
                         muR_syst == MyCorrection::variation::down) {
                index = 1;
              } else if (muF_syst == MyCorrection::variation::up &&
                         muR_syst == MyCorrection::variation::down) {
                index = 2;
              } else if (muF_syst == MyCorrection::variation::down &&
                         muR_syst == MyCorrection::variation::nom) {
                index = 3;
              } else if (muF_syst == MyCorrection::variation::nom &&
                         muR_syst == MyCorrection::variation::nom) {
                return 1.f;
              } else if (muF_syst == MyCorrection::variation::up &&
                         muR_syst == MyCorrection::variation::nom) {
                index = 5;
              } else if (muF_syst == MyCorrection::variation::down &&
                         muR_syst == MyCorrection::variation::up) {
                index = 6;
              } else if (muF_syst == MyCorrection::variation::nom &&
                         muR_syst == MyCorrection::variation::up) {
                index = 7;
              } else if (muF_syst == MyCorrection::variation::up &&
                         muR_syst == MyCorrection::variation::up) {
                index = 8;
              }
              if (index < 0 || nLHEScaleWeight <= index) return 1.f;
              const float weight = LHEScaleWeight[index];
              return std::isfinite(weight) ? weight : 1.f;
            }

            const float weight = GetScaleVariation(muF_syst, muR_syst);
            return std::isfinite(weight) ? weight : 1.f;
          };

      const float nominal_systematic_weight =
          weight_map.count("Central") ? weight_map["Central"] : 1.f;
      for (int member = 1; member <= 100; ++member) {
        std::ostringstream pdf_name;
        pdf_name << "PDF_" << std::setw(3) << std::setfill('0') << member;
        weight_map[pdf_name.str()] =
            nominal_systematic_weight * safe_lhe_pdf_weight(member);
      }
      weight_map["AlphaS_Down"] =
          nominal_systematic_weight * safe_lhe_pdf_weight(101);
      weight_map["AlphaS_Up"] =
          nominal_systematic_weight * safe_lhe_pdf_weight(102);

      if (LHEScaleWeight.valid() && nLHEScaleWeight >= 8) {
        weight_map["LHEScale_MuF_Up"] =
            nominal_systematic_weight *
            safe_lhe_scale_weight(MyCorrection::variation::up,
                                  MyCorrection::variation::nom);
        weight_map["LHEScale_MuF_Down"] =
            nominal_systematic_weight *
            safe_lhe_scale_weight(MyCorrection::variation::down,
                                  MyCorrection::variation::nom);
        weight_map["LHEScale_MuR_Up"] =
            nominal_systematic_weight *
            safe_lhe_scale_weight(MyCorrection::variation::nom,
                                  MyCorrection::variation::up);
        weight_map["LHEScale_MuR_Down"] =
            nominal_systematic_weight *
            safe_lhe_scale_weight(MyCorrection::variation::nom,
                                  MyCorrection::variation::down);
        weight_map["LHEScale_MuFMuR_Up"] =
            nominal_systematic_weight *
            safe_lhe_scale_weight(MyCorrection::variation::up,
                                  MyCorrection::variation::up);
        weight_map["LHEScale_MuFMuR_Down"] =
            nominal_systematic_weight *
            safe_lhe_scale_weight(MyCorrection::variation::down,
                                  MyCorrection::variation::down);
      }

      if (PSWeight.valid() && nPSWeight >= 4) {
        weight_map["FSR_Up"] =
            nominal_systematic_weight *
            GetPSWeight(MyCorrection::variation::nom,
                        MyCorrection::variation::up);
        weight_map["FSR_Down"] =
            nominal_systematic_weight *
            GetPSWeight(MyCorrection::variation::nom,
                        MyCorrection::variation::down);
        weight_map["ISR_Up"] =
            nominal_systematic_weight *
            GetPSWeight(MyCorrection::variation::up,
                        MyCorrection::variation::nom);
        weight_map["ISR_Down"] =
            nominal_systematic_weight *
            GetPSWeight(MyCorrection::variation::down,
                        MyCorrection::variation::nom);
      }

      if (MCSample.Contains("powheg") && MCSample.Contains("TT")) {
        const auto top_indices = GetTopAndAntiTopIndices(AllGenViews);
        constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
        if (top_indices[0] != npos && top_indices[1] != npos) {
          const TLorentzVector first_copy_top =
              AllGenViews[top_indices[0]].P4();
          const TLorentzVector first_copy_antitop =
              AllGenViews[top_indices[1]].P4();
          weight_map["hdamp_Up"] =
              nominal_systematic_weight *
              myCorr->GethDampReweight(first_copy_top, first_copy_antitop,
                                       MyCorrection::variation::up);
          weight_map["hdamp_Down"] =
              nominal_systematic_weight *
              myCorr->GethDampReweight(first_copy_top, first_copy_antitop,
                                       MyCorrection::variation::down);
        }
      }
    }
  }

  //==== Kinematic Fitter / Combinatorics
  std::vector<unsigned int> top_b_jet_candidates;
  std::vector<unsigned int> had_W_candidates;

  for (unsigned int ij = 0; ij < 5; ij++) {
    if (std::find(bjet_indices.begin(), bjet_indices.end(), ij) != bjet_indices.end()) {
      if (top_b_jet_candidates.size() < 2) top_b_jet_candidates.push_back(ij);
      else if (top_b_jet_candidates.size() >= 2) had_W_candidates.push_back(ij);
    } else {
      if (had_W_candidates.size() < 3) had_W_candidates.push_back(ij);
    }
  }

  if (had_W_candidates.size() < 2 || top_b_jet_candidates.size() < 2) return;
  FillHist(channel_name + "/cutflow/" + event_hist_name, 4.5, 1., 6, 0., 6.);

  AtobbInference::ttCombinatoric combinatorics[6];
  int comb_idx = 0;

  for (unsigned int w1 = 0; w1 < had_W_candidates.size() - 1; w1++) {
    for (unsigned int w2 = w1 + 1; w2 < had_W_candidates.size(); w2++) {
      for (unsigned int b_idx = 0; b_idx < 2; b_idx++) {
        combinatorics[comb_idx].lepton            = selected_lepton;
        combinatorics[comb_idx].jets              = &jets;
        combinatorics[comb_idx].met               = &METv;
        combinatorics[comb_idx].had_W_jet_idx_1   = had_W_candidates.at(w1);
        combinatorics[comb_idx].had_W_jet_idx_2   = had_W_candidates.at(w2);
        combinatorics[comb_idx].had_top_b_jet_idx = top_b_jet_candidates.at(b_idx);
        combinatorics[comb_idx].lep_top_b_jet_idx = top_b_jet_candidates.at(1 - b_idx);

        this->EvalChi2(combinatorics[comb_idx]);
        comb_idx++;
      }
    }
  }

  AtobbInference::ttCombinatoric* best_combinatoric = &combinatorics[0];
  for (int i = 1; i < comb_idx; i++) {
    if (combinatorics[i].best_chi2 < best_combinatoric->best_chi2) {
      best_combinatoric = &combinatorics[i];
    }
  }

  // NOTE: Chi2Cut filtering is intentionally NOT applied per user directive.
  FillHist(channel_name + "/cutflow/" + event_hist_name, 5.5, 1., 6, 0., 6.);

  //==== 31 Feature Extraction for XGBoost ONNX Models
  std::vector<int> non_bjet_indices;
  for (size_t i = 0; i < jets.size(); ++i) {
    if (std::find(bjet_indices.begin(), bjet_indices.end(), i) == bjet_indices.end()) {
      non_bjet_indices.push_back(i);
    }
    if (non_bjet_indices.size() >= 2) break;
  }

  float HT = 0.f;
  for (const auto& jet : jets) HT += float(jet.Pt());
  float ST = HT + float(selected_lepton->Pt()) + float(METv.Pt());
  float lnu_mass = float((static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());

  float jj_mass = 0.f;
  if (non_bjet_indices.size() >= 2) {
    TLorentzVector j0 = static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0)));
    TLorentzVector j1 = static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)));
    jj_mass = float((j0 + j1).M());
  }

  float bb_mass_01 = -999.f, bb_mass_02 = -999.f, bb_mass_12 = -999.f;
  float bb_mass_03 = 0.f, bb_mass_13 = 0.f, bb_mass_23 = 0.f;

  float bb_dr_01 = -999.f, bb_dr_02 = -999.f, bb_dr_12 = -999.f;
  float bb_dr_03 = 0.f, bb_dr_13 = 0.f, bb_dr_23 = 0.f;

  float bblnu_mass_01 = -999.f, bblnu_mass_02 = -999.f, bblnu_mass_12 = -999.f;
  float bblnu_mass_03 = 0.f, bblnu_mass_13 = 0.f, bblnu_mass_23 = 0.f;

  float bbjj_mass_01 = -999.f, bbjj_mass_02 = -999.f, bbjj_mass_12 = -999.f;
  float bbjj_mass_03 = 0.f, bbjj_mass_13 = 0.f, bbjj_mass_23 = 0.f;

  float bbbW_mass_012 = 0.f;

  if (bjet_indices.size() >= 2) {
    TLorentzVector b0 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(0)));
    TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(1)));
    bb_mass_01 = float((b0 + b1).M());
    bb_dr_01 = float(b0.DeltaR(b1));
    bblnu_mass_01 = float((b0 + b1 + static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());
    if (non_bjet_indices.size() >= 2) {
      bbjj_mass_01 = float((b0 + b1 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
    }
  }

  if (bjet_indices.size() >= 3) {
    TLorentzVector b0 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(0)));
    TLorentzVector b1 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(1)));
    TLorentzVector b2 = static_cast<TLorentzVector>(jets.at(bjet_indices.at(2)));
    bb_mass_02 = float((b0 + b2).M());
    bb_mass_12 = float((b1 + b2).M());
    bb_dr_02 = float(b0.DeltaR(b2));
    bb_dr_12 = float(b1.DeltaR(b2));

    bblnu_mass_02 = float((b0 + b2 + static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());
    bblnu_mass_12 = float((b1 + b2 + static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());

    if (non_bjet_indices.size() >= 2) {
      bbjj_mass_02 = float((b0 + b2 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      bbjj_mass_12 = float((b1 + b2 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());

      float bbbW_mass_012_cand1 = float((b0 + b1 + b2 + static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());
      float bbbW_mass_012_cand2 = float((b0 + b1 + b2 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      if (std::abs(bbbW_mass_012_cand1 - const_top_mass) < std::abs(bbbW_mass_012_cand2 - const_top_mass)) {
        bbbW_mass_012 = bbbW_mass_012_cand1;
      } else {
        bbbW_mass_012 = bbbW_mass_012_cand2;
      }
    }
  }

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

    bblnu_mass_03 = float((b0 + b3 + static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());
    bblnu_mass_13 = float((b1 + b3 + static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());
    bblnu_mass_23 = float((b2 + b3 + static_cast<TLorentzVector>(*selected_lepton) + static_cast<TLorentzVector>(METv)).M());

    if (non_bjet_indices.size() >= 2) {
      bbjj_mass_03 = float((b0 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      bbjj_mass_13 = float((b1 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
      bbjj_mass_23 = float((b2 + b3 + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(0))) + static_cast<TLorentzVector>(jets.at(non_bjet_indices.at(1)))).M());
    }
  }

  // Derive bb_summary_features
  float pair_masses[6] = {bb_mass_01, bb_mass_02, bb_mass_12, bb_mass_03, bb_mass_13, bb_mass_23};
  float pair_drs[6]    = {bb_dr_01,   bb_dr_02,   bb_dr_12,   bb_dr_03,   bb_dr_13,   bb_dr_23};

  float bb_mass_min_all = 1e9f;
  float bb_mass_max_all = -1e9f;
  float bb_dr_min_all = 1e9f;
  float bb_mass_at_min_dr = -999.f;
  bool has_valid_pair = false;

  for (int p = 0; p < 6; p++) {
    float m = pair_masses[p];
    float dr = pair_drs[p];
    if (std::isfinite(m) && std::isfinite(dr) && m > 0.f && dr > 0.f) {
      has_valid_pair = true;
      if (m < bb_mass_min_all) bb_mass_min_all = m;
      if (m > bb_mass_max_all) bb_mass_max_all = m;
      if (dr < bb_dr_min_all) {
        bb_dr_min_all = dr;
        bb_mass_at_min_dr = m;
      }
    }
  }
  if (!has_valid_pair) {
    bb_mass_min_all = std::numeric_limits<float>::quiet_NaN();
    bb_mass_max_all = std::numeric_limits<float>::quiet_NaN();
    bb_dr_min_all = std::numeric_limits<float>::quiet_NaN();
    bb_mass_at_min_dr = std::numeric_limits<float>::quiet_NaN();
  }

  // Construct 31 features array matching feature_indices in xgboost_atobb_MHc130_*_features.json
  FloatArray features(31);
  features[0]  = HT;
  features[1]  = ST;
  features[2]  = float(METv.Pt());
  features[3]  = float(jets.size());
  features[4]  = float(NBJets);
  features[5]  = float(selected_lepton->Pt());
  features[6]  = float(selected_lepton->Eta());
  features[7]  = bb_mass_min_all;
  features[8]  = bb_mass_max_all;
  features[9]  = bb_mass_at_min_dr;
  features[10] = bb_dr_min_all;
  features[11] = bb_mass_01;
  features[12] = bb_mass_02;
  features[13] = bb_mass_12;
  features[14] = bb_dr_01;
  features[15] = bb_dr_02;
  features[16] = bb_dr_12;
  features[17] = bblnu_mass_01;
  features[18] = bblnu_mass_02;
  features[19] = bblnu_mass_12;
  features[20] = bbjj_mass_01;
  features[21] = bbjj_mass_02;
  features[22] = bbjj_mass_12;
  features[23] = bbbW_mass_012;
  features[24] = lnu_mass;
  features[25] = jj_mass;
  features[26] = float(best_combinatoric->best_chi2);
  features[27] = float(best_combinatoric->had_W_mass);
  features[28] = float(best_combinatoric->had_top_mass);
  features[29] = float(best_combinatoric->best_lep_W_mass);
  features[30] = float(best_combinatoric->best_lep_top_mass);

  std::unordered_map<std::string, VariousArray> inputDataMap;
  inputDataMap["input"] = features;

  static const std::unordered_map<std::string, IntArray> inputDataShapeMap = {
      {"input", IntArray{1, 31}}
  };

  //==== BDT Inference for each mA mass point
  std::map<TString, float> ma_scores;
  std::vector<float> all_scores;
  all_scores.reserve(ma_points.size());

  for (const auto& ma : ma_points) {
    if (ml_helpers.count(ma) && ml_helpers[ma]) {
      try {
        auto outputMap = ml_helpers[ma]->Run_ONNX_Model(inputDataMap, inputDataShapeMap);
        if (outputMap.count("probabilities") && outputMap["probabilities"].size() >= 2) {
          float bdt_score = outputMap["probabilities"][1]; // signal class score
          ma_scores[ma] = bdt_score;
          all_scores.push_back(bdt_score);
        }
      } catch (const std::exception& e) {
        static bool print_once = true;
        if (print_once) {
          cerr << "[AtobbInference::executeEventFromParameter] ONNX inference error for " << ma << ": " << e.what() << endl;
          print_once = false;
        }
      }
    }
  }

  //==== Compute aggregated score metrics: topk1, topk2, topk3, avg
  std::vector<float> sorted_scores = all_scores;
  std::sort(sorted_scores.begin(), sorted_scores.end(), std::greater<float>());

  float score_topk1 = !sorted_scores.empty() ? sorted_scores[0] : 0.f;
  float score_topk2 = 0.f;
  if (sorted_scores.size() >= 2) {
    score_topk2 = (sorted_scores[0] + sorted_scores[1]) / 2.f;
  } else if (!sorted_scores.empty()) {
    score_topk2 = sorted_scores[0];
  }

  float score_topk3 = 0.f;
  if (sorted_scores.size() >= 3) {
    score_topk3 = (sorted_scores[0] + sorted_scores[1] + sorted_scores[2]) / 3.f;
  } else if (!sorted_scores.empty()) {
    float sum_k = 0.f;
    for (float s : sorted_scores) sum_k += s;
    score_topk3 = sum_k / float(sorted_scores.size());
  }

  float score_avg = 0.f;
  if (!all_scores.empty()) {
    float sum_all = 0.f;
    for (float s : all_scores) sum_all += s;
    score_avg = sum_all / float(all_scores.size());
  }

  //==== Filling Score Histograms for each Systematic Weight
  for (const auto &w : weight_map) {
    TString systName = w.first;
    float weight = base_weight * w.second;
    TString hist_name = hist_sample_prefix;
    if (systName != "Central") hist_name += "_" + systName;

    for (const auto& ma_pair : ma_scores) {
      FillHist(channel_name + "/bdt_score_" + ma_pair.first + "/" + hist_name, ma_pair.second, weight, 100, 0., 1.);
    }

    if (!all_scores.empty()) {
      FillHist(channel_name + "/bdt_score_topk1/" + hist_name, score_topk1, weight, 100, 0., 1.);
      FillHist(channel_name + "/bdt_score_topk2/" + hist_name, score_topk2, weight, 100, 0., 1.);
      FillHist(channel_name + "/bdt_score_topk3/" + hist_name, score_topk3, weight, 100, 0., 1.);
      FillHist(channel_name + "/bdt_score_avg/" + hist_name, score_avg, weight, 100, 0., 1.);
    }
  }
}

void AtobbInference::EvalChi2(ttCombinatoric& tt_combinatoric) {
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

double AtobbInference::Chi2Function(double had_top_mass, double had_W_mass, double lep_top_mass, double lep_W_mass) {
  double chi2 = 0.;
  chi2 += TMath::Power( (had_top_mass - const_top_mass )/const_top_width,  2);
  chi2 += TMath::Power( (had_W_mass   - const_w_mass   )/const_w_width,    2);
  chi2 += TMath::Power( (lep_top_mass - const_top_mass )/const_top_width,  2);
  chi2 += TMath::Power( (lep_W_mass   - const_w_mass   )/const_w_width,    2);
  return chi2;
}

void AtobbInference::ttCombinatoric::EvalHadronicPart() {
  TLorentzVector had_W_vector   = static_cast<TLorentzVector>(jets->at( had_W_jet_idx_1 ))
                                + static_cast<TLorentzVector>(jets->at( had_W_jet_idx_2 ));

  TLorentzVector had_top_vector = had_W_vector
                                + static_cast<TLorentzVector>(jets->at( had_top_b_jet_idx ));

  had_W_mass   = had_W_vector.M();
  had_top_mass = had_top_vector.M();
}

void AtobbInference::ttCombinatoric::EvalLeptonicPart() {
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

array<size_t, 4> AtobbInference::GetTopAndAntiTopIndices(const GenViewCollection &gens) {
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
    cerr << "[AtobbInference::GetTopAndAntiTopIndices] " << message
         << ": keeping index " << current_idx << ", ignoring index "
         << duplicate_idx << endl;
  };

  auto warn_missing = [](const char *message) {
    cerr << "[AtobbInference::GetTopAndAntiTopIndices] " << message << endl;
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
