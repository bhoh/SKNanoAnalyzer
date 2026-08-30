#include "Tutorial_dilep.h"

Tutorial_dilep::Tutorial_dilep() {}

Tutorial_dilep::~Tutorial_dilep() {}

void Tutorial_dilep::initializeAnalyzer() {

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
    cerr << "[Tutorial_dilep::initializeAnalyzer] DataEra is not set properly"
         << endl;
    exit(EXIT_FAILURE);
  }

  myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
  myCorr->SetTaggingParam(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);

  string SKNANO_HOME = getenv("SKNANO_HOME");
  if (IsDATA) {
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
  } else {
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/TutorialDilepSystematic.yaml", MCSample, DataEra);
  }
}

void Tutorial_dilep::executeEvent() {
  AllMuonViews = GetAllMuonViews();
  AllElectronViews = GetAllElectronViews();
  AllGenViews = GetAllGenViews();

  ev = GetEvent();

  for (const auto &syst_dummy : *systHelper) {
    executeEventFromParameter();
  }
}

void Tutorial_dilep::executeEventFromParameter() {
  const TString this_syst = systHelper->getCurrentSysName();
  if (IsDATA && this_syst != "Central") return;

  Muon::MuonID this_muon_id = MuonIDs[0];
  TString this_muon_id_sf_key = MuonIDISOSFKeys[0];
  TString this_muon_iso_sf_key = MuonIDISOSFKeys[1];
  TString this_muon_trig_sf_key = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";

  Electron::ElectronID this_electron_id = ElectronIDs[0];
  TString this_electron_id_sf_key = ElectronIDSFKeys[0];
  TString this_electron_trig_sf_key = ElectronTriggerSFKeys[0];

  std::vector<size_t> SelectedMuonIndices_id_only = SelectMuonIndices(AllMuonViews, this_muon_id, 15., 2.4);
  std::vector<size_t> SelectedMuonIndices = SelectMuonIndices(AllMuonViews, SelectedMuonIndices_id_only, Muon::MuonID::POG_PFISO_TIGHT, 15., 2.4);

  std::vector<size_t> SelectedElectronVetoIndices = SelectElectronIndices(AllElectronViews, Electron::ElectronID::POG_MVAISO_WP80, 15., 2.5);
  std::vector<size_t> SelectedElectronIndices = SelectElectronIndices(AllElectronViews, SelectedElectronVetoIndices, this_electron_id, 15., 2.5);

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

  bool is_muon_channel = (muons.size() == 2 && SelectedElectronVetoIndices.size() == 0);
  bool is_electron_channel = (electrons.size() == 2 && SelectedMuonIndices.size() == 0);

  TString channel_name = "";
  double dilepton_mass = 0.;

  if (is_muon_channel) {
    if (muons[0].Charge() * muons[1].Charge() >= 0) return;
    if (!(ev.PassTrigger(IsoMuTriggerName))) return;
    if (muons.at(0).Pt() <= TriggerSafePtCut) return;
    channel_name = "Muon";
    dilepton_mass = (muons[0] + muons[1]).M();
  } else if (is_electron_channel) {
    if (electrons[0].Charge() * electrons[1].Charge() >= 0) return;
    if (!(ev.PassTrigger(EleTriggerName))) return;
    if (electrons.at(0).Pt() <= EleTriggerSafePtCut) return;
    channel_name = "Electron";
    dilepton_mass = (electrons[0] + electrons[1]).M();
  } else {
    return;
  }

  TString sample_name = IsDATA ? "Data" : MCSample;

  // Setup Scale Factors
  unordered_map<std::string, std::variant<std::function<float(MyCorrection::variation, TString)>, std::function<float()>>> weight_function_map;

  if (!IsDATA) {
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
          if (!is_muon_channel) return 1.f;
          return myCorr->GetMuonIDSF(this_muon_id_sf_key, muons, syst);
        };
    weight_function_map["Muon_ID"] = muon_id_lambda;

    std::function<float(MyCorrection::variation, TString)> muon_iso_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          if (!is_muon_channel) return 1.f;
          return myCorr->GetMuonISOSF(this_muon_iso_sf_key, muons, syst);
        };
    weight_function_map["Muon_Iso"] = muon_iso_lambda;

    std::function<float(MyCorrection::variation, TString)> muon_trig_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          if (!is_muon_channel) return 1.f;
          return myCorr->GetMuonTriggerSF(this_muon_trig_sf_key, muons.at(0), syst);
        };
    weight_function_map["Muon_Trig"] = muon_trig_lambda;

    std::function<float(MyCorrection::variation, TString)> electron_id_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          if (!is_electron_channel) return 1.f;
          return myCorr->GetElectronIDSF(this_electron_id_sf_key, electrons, syst);
        };
    weight_function_map["Electron_ID"] = electron_id_lambda;

    std::function<float(MyCorrection::variation, TString)> electron_reco_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          if (!is_electron_channel) return 1.f;
          return myCorr->GetElectronRECOSF(electrons, syst);
        };
    weight_function_map["Electron_Reco"] = electron_reco_lambda;

    std::function<float(MyCorrection::variation, TString)> electron_trig_lambda =
        [&](MyCorrection::variation syst, TString source) {
          (void)source;
          if (!is_electron_channel) return 1.f;
          const Electron &electron = electrons.at(0);
          return myCorr->GetElectronTriggerSF(this_electron_trig_sf_key,
                                             electron.scEta(), electron.Pt(),
                                             electron.Phi(), syst);
        };
    weight_function_map["Electron_Trig"] = electron_trig_lambda;

    std::function<float(MyCorrection::variation, TString)> pileup_lambda =
        [&](MyCorrection::variation syst, TString source) {
          return myCorr->GetPUWeight(ev.nTrueInt(), syst, source);
        };
    weight_function_map["Pileup"] = pileup_lambda;
  }

  systHelper->assignWeightFunctionMap(weight_function_map);

  float base_weight = 1.;
  if (!IsDATA) {
    base_weight *= MCweight();
    base_weight *= ev.GetTriggerLumi("Full");
  }

  unordered_map<std::string, float> weight_map;
  if (IsDATA) {
    weight_map["Central"] = 1.f;
  } else {
    weight_map = systHelper->calculateWeight();

    // Fill PDF/Scale/Parton Shower weights if needed
    if (this_syst == "Central") {
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

      const float central_weight = weight_map.count("Central") ? weight_map["Central"] : 1.f;

      for (int member = 1; member <= 100; ++member) {
        std::ostringstream pdf_name;
        pdf_name << "PDF_" << std::setw(3) << std::setfill('0') << member;
        weight_map[pdf_name.str()] = central_weight * safe_lhe_pdf_weight(member);
      }
      weight_map["AlphaS_Down"] = central_weight * safe_lhe_pdf_weight(101);
      weight_map["AlphaS_Up"] = central_weight * safe_lhe_pdf_weight(102);

      if (LHEScaleWeight.valid() && nLHEScaleWeight >= 8) {
        weight_map["LHEScale_MuF_Up"] = central_weight * safe_lhe_scale_weight(MyCorrection::variation::up, MyCorrection::variation::nom);
        weight_map["LHEScale_MuF_Down"] = central_weight * safe_lhe_scale_weight(MyCorrection::variation::down, MyCorrection::variation::nom);
        weight_map["LHEScale_MuR_Up"] = central_weight * safe_lhe_scale_weight(MyCorrection::variation::nom, MyCorrection::variation::up);
        weight_map["LHEScale_MuR_Down"] = central_weight * safe_lhe_scale_weight(MyCorrection::variation::nom, MyCorrection::variation::down);
        weight_map["LHEScale_MuFMuR_Up"] = central_weight * safe_lhe_scale_weight(MyCorrection::variation::up, MyCorrection::variation::up);
        weight_map["LHEScale_MuFMuR_Down"] = central_weight * safe_lhe_scale_weight(MyCorrection::variation::down, MyCorrection::variation::down);
      }

      if (PSWeight.valid() && nPSWeight >= 4) {
        weight_map["FSR_Up"] = central_weight * GetPSWeight(MyCorrection::variation::nom, MyCorrection::variation::up);
        weight_map["FSR_Down"] = central_weight * GetPSWeight(MyCorrection::variation::nom, MyCorrection::variation::down);
        weight_map["ISR_Up"] = central_weight * GetPSWeight(MyCorrection::variation::up, MyCorrection::variation::nom);
        weight_map["ISR_Down"] = central_weight * GetPSWeight(MyCorrection::variation::down, MyCorrection::variation::nom);
      }
    }
  }

  for (const auto &w : weight_map) {
    TString systName = w.first;
    float weight = base_weight * w.second;
    TString hist_name = "histo_" + sample_name;
    if (systName != "Central") hist_name += "_" + systName;

    FillHist(channel_name + "/dilepton_mass/" + hist_name, dilepton_mass, weight, 100, 0., 200.);
    FillHist(channel_name + "/dilepton_mass_Zpeak/" + hist_name, dilepton_mass, weight, 60, 60., 120.);
  }
}
