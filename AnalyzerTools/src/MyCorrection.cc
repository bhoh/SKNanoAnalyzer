#include "MyCorrection.h"
#include "MLHelper.h"

#include <TLorentzVector.h>
#include <execution>
#include <iostream>
#include <numeric>
#include <vector>

MyCorrection::MyCorrection() {}
MyCorrection::MyCorrection(const TString &era, const TString &period,
                           const TString &sample, const bool IsData,
                           const string &btagging_eff_file,
                           const string &ctagging_eff_file,
                           const string &btagging_R_file,
                           const string &ctagging_R_file) {
  cout << "[MyCorrection::MyCorrection] MyCorrection created for " << era
       << endl;
  SetEra(era);
  SetPeriod(period);
  SetSample(sample);
  setIsData(IsData);

  EraConfig config = GetEraConfig(era, btagging_eff_file, ctagging_eff_file,
                                  btagging_R_file, ctagging_R_file);
  struct CorrectionInfo {
    string name;
    string path;
    unique_ptr<CorrectionSet> &cset;
    bool isOptional;
  };

  // load correction sets
  // only jetid and jet, jerc are available for 2024 for now, I change optional
  // flag to true for other sets
  loadCorrectionSet("muon SF", config.json_muon, cset_muon, true);
  loadCorrectionSet("puWeights", config.json_puWeights, cset_puWeights, true);
  loadCorrectionSet("btagging", config.json_btagging, cset_btagging, true);
  loadCorrectionSet("ctagging", config.json_ctagging, cset_ctagging, true);
  loadCorrectionSet("btagging eff", config.json_btagging_eff, cset_btagging_eff,
                    true);
  loadCorrectionSet("ctagging eff", config.json_ctagging_eff, cset_ctagging_eff,
                    true);
  loadCorrectionSet("electron", config.json_electron, cset_electron, true);
  loadCorrectionSet("electron variation", config.json_electron_variation,
                    cset_electron_variation, true);
  loadCorrectionSet("photon", config.json_photon, cset_photon, true);
  loadCorrectionSet("jetid", config.json_jetid, cset_jetid, false);
  loadCorrectionSet("jerc", config.json_jerc, cset_jerc, false);
  loadCorrectionSet("jerc_fatjet", config.json_jerc_fatjet, cset_jerc_fatjet,
                    true);
  loadCorrectionSet("jetvetomap", config.json_jetvetomap, cset_jetvetomap,
                    false);
  // Optional files
  loadRoccoR(config.txt_roccor, true);
  loadGoldenJson(config.golden_json, true);
  loadCorrectionSet("jmar", config.json_jmar, cset_jmar, true);
  loadCorrectionSet("muon trig eff", config.json_muon_trig_eff,
                    cset_muon_trig_eff, true);
  loadCorrectionSet("muon trig sf", config.json_muon_trig_sf, cset_muon_trig_sf,
                    true); // temporary due to no mu trig sf for 2024
  loadCorrectionSet("electron hlt", config.json_electron_hlt, cset_electron_hlt,
                    true);
  loadCorrectionSet("met", config.json_met, cset_met, true);
  loadCorrectionSet("btagging R", config.json_btagging_R, cset_btagging_R,
                    true);
  loadCorrectionSet("ctagging R", config.json_ctagging_R, cset_ctagging_R,
                    true);
  loadCorrectionSet("muon TopHNT idsf", config.json_muon_TopHNT_idsf,
                    cset_muon_TopHNT_idsf, true);
  loadCorrectionSet("muon TopHNT dblmu leg1 eff",
                    config.json_muon_TopHNT_dblmu_leg1_eff,
                    cset_muon_TopHNT_dblmu_leg1_eff, true);
  loadCorrectionSet("muon TopHNT dblmu leg2 eff",
                    config.json_muon_TopHNT_dblmu_leg2_eff,
                    cset_muon_TopHNT_dblmu_leg2_eff, true);
  loadCorrectionSet("muon TopHNT emu leg1 eff",
                    config.json_muon_TopHNT_emu_leg1_eff,
                    cset_muon_TopHNT_emu_leg1_eff, true);
  loadCorrectionSet("muon TopHNT emu leg2 eff",
                    config.json_muon_TopHNT_emu_leg2_eff,
                    cset_muon_TopHNT_emu_leg2_eff, true);
  loadCorrectionSet("electron TopHNT idsf", config.json_electron_TopHNT_idsf,
                    cset_electron_TopHNT_idsf, true);
  loadCorrectionSet("electron TopHNT emu leg1 eff",
                    config.json_electron_TopHNT_emu_leg1_eff,
                    cset_electron_TopHNT_emu_leg1_eff, true);
  loadCorrectionSet("electron TopHNT emu leg2 eff",
                    config.json_electron_TopHNT_emu_leg2_eff,
                    cset_electron_TopHNT_emu_leg2_eff, true);

  MLHelper_hDampUp =
      make_unique<MLHelper>(config.onnx_hDampUp, MLHelper::ModelType::ONNX);
  MLHelper_hDampDown =
      make_unique<MLHelper>(config.onnx_hDampDown, MLHelper::ModelType::ONNX);
  MLHelper_TopPtReweight = make_unique<MLHelper>(config.onnx_toppt_reweight,
                                                 MLHelper::ModelType::ONNX);
  MLHelper_rBnom =
      make_unique<MLHelper>(config.onnx_rBnom, MLHelper::ModelType::ONNX);
  MLHelper_rBUp =
      make_unique<MLHelper>(config.onnx_rBUp, MLHelper::ModelType::ONNX);

  LUM_keys["2024"] = "Collisions24_BCDEFGHI_goldenJSON";
  LUM_keys["2023BPix"] = "Collisions2023_369803_370790_eraD_GoldenJson";
  LUM_keys["2023"] = "Collisions2023_366403_369802_eraBC_GoldenJson";
  LUM_keys["2022EE"] = "Collisions2022_359022_362760_eraEFG_GoldenJson";
  LUM_keys["2022"] = "Collisions2022_355100_357900_eraBCD_GoldenJson";
  LUM_keys["2018"] = "Collisions18_UltraLegacy_goldenJSON";
  LUM_keys["2017"] = "Collisions17_UltraLegacy_goldenJSON";
  LUM_keys["2016postVFP"] = "Collisions16_UltraLegacy_goldenJSON";
  LUM_keys["2016preVFP"] = "Collisions16_UltraLegacy_goldenJSON";

  EGM_keys["2024"] = "2024Prompt";
  EGM_keys["2023BPix"] = "2023PromptD";
  EGM_keys["2023"] = "2023PromptC";
  EGM_keys["2022EE"] = "2022Re-recoE+PromptFG";
  EGM_keys["2022"] = "2022Re-recoBCD";
  EGM_keys["2016preVFP"] = "2016preVFP";
  EGM_keys["2016postVFP"] = "2016postVFP";
  EGM_keys["2017"] = "2017";
  EGM_keys["2018"] = "2018";

  // Please use ####### as placeholder
  if (!IsData) {
    JME_JER_GT["2024"] =
        "Summer24Prompt24_JRV1_MC_######_AK4PFPuppi";
                                                               // because real
                                                               // content of
                                                               // file is this
    JME_JES_GT["2024"] = "Summer24Prompt24_V3_MC_######_AK4PFPuppi";
  } else {
    // JME_JER_GT["2024"] =
    // "Summer23BPixPrompt23_RunD_JRV1_DATA_######_AK4PFPuppi"; // this is
    // because real content of file is this
    JME_JES_GT["2024"] = "Summer24Prompt24_V3_DATA_######_AK4PFPuppi";
  }

  JME_vetomap_keys["2024"] = "Summer24Prompt24_RunBCDEFGHI_V1";

  JME_PILEUP_keys["2016preVFP"] = "PUJetID_eff";
  JME_PILEUP_keys["2016postVFP"] = "PUJetID_eff";
  JME_PILEUP_keys["2017"] = "PUJetID_eff";
  JME_PILEUP_keys["2018"] = "PUJetID_eff";
}

MyCorrection::~MyCorrection() {}

MyCorrection::EraConfig
MyCorrection::GetEraConfig(TString era, const string &btagging_eff_file,
                           const string &ctagging_eff_file,
                           const string &btagging_R_file,
                           const string &ctagging_R_file) const {
  EraConfig config;

  const char *json_pog_path = getenv("JSONPOG_REPO_PATH");
  const char *sknano_data = getenv("SKNANO_DATA");
  const char *external_roccor = getenv("ROCCOR_PATH");
  cout << "[MyCorrection::GetEraConfig] json_pog_path: "
       << (json_pog_path ? json_pog_path : "NULL") << endl;
  cout << "[MyCorrection::GetEraConfig] sknano_data: "
       << (sknano_data ? sknano_data : "NULL") << endl;
  cout << "[MyCorrection::GetEraConfig] external_roccor: "
       << (external_roccor ? external_roccor : "NULL") << endl;
  if (!json_pog_path || !sknano_data || !external_roccor) {
    throw runtime_error(
        "JSONPOG_REPO_PATH or SKNANO_DATA or ROCCOR_PATH is not set");
  }

  const string json_pog_path_str(json_pog_path);
  const string sknano_data_str(sknano_data);
  const string external_roccor_str(external_roccor);

  config.json_muon = json_pog_path_str + "/MUO";
  config.json_muon_trig_eff = sknano_data_str;
  config.json_muon_trig_sf =
      json_pog_path_str + "/MUO"; // temporary due to no mu trig sf for 2024
  config.json_puWeights = json_pog_path_str + "/LUM";
  config.json_btagging = json_pog_path_str + "/BTV";
  config.json_ctagging = json_pog_path_str + "/BTV";
  config.json_btagging_eff = sknano_data_str;
  config.json_ctagging_eff = sknano_data_str;
  config.json_btagging_R = sknano_data_str;
  config.json_ctagging_R = sknano_data_str;
  config.json_electron = json_pog_path_str + "/EGM";
  config.json_electron_hlt = config.json_electron;
  config.json_electron_variation = json_pog_path_str + "/EGM";
  config.json_photon = json_pog_path_str + "/EGM";
  config.json_jetid = json_pog_path_str + "/JME";
  config.json_jerc = json_pog_path_str + "/JME";
  config.json_jerc_fatjet = json_pog_path_str + "/JME";
  config.json_jetvetomap = json_pog_path_str + "/JME";
  config.json_jmar = json_pog_path_str + "/JME";
  config.json_met = json_pog_path_str + "/JME";
  config.txt_roccor = external_roccor_str;
  config.golden_json = sknano_data_str;

  config.onnx_hDampDown = sknano_data_str;
  config.onnx_hDampUp = sknano_data_str;
  config.onnx_toppt_reweight = sknano_data_str;
  config.onnx_rBnom = sknano_data_str;
  config.onnx_rBUp = sknano_data_str;

  // config.json_muon_custom_TopHNT_idsf = sknano_data_str;
  // config.json_muon_custom_dblmu_leg1_eff = sknano_data_str;
  // config.json_muon_custom_dblmu_leg2_eff = sknano_data_str;
  // config.json_muon_custom_emu_leg1_eff = sknano_data_str;
  // config.json_muon_custom_emu_leg2_eff = sknano_data_str;
  // config.json_electron_custom_TopHNT_idsf = sknano_data_str;
  // config.json_electron_custom_emu_leg1_eff = sknano_data_str;
  // config.json_electron_custom_emu_leg2_eff = sknano_data_str;

  if (era == "2024") {
    const string tag = "/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/";
    const string tag_jme_temp = "/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2025-12-02/";
    const string tag_temp = "/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/";
    config.json_muon += tag + "muon_Z.json.gz";
    config.json_muon_trig_eff += "/2024/MUO/muon_trig.json";
    config.json_muon_trig_sf += tag + "muon_Z.json.gz";
    config.json_puWeights += tag + "puWeights_BCDEFGHI.json.gz";
    config.json_btagging += tag + "btagging.json.gz";
    // config.json_ctagging += "/2023_Summer23BPix/ctagging.json.gz";
    config.json_btagging_eff += "/2024/BTV/" + btagging_eff_file;
    // config.json_ctagging_eff += "/2023BPix/BTV/" + ctagging_eff_file;
    // config.json_btagging_R += "/2023BPix/BTV/" + btagging_R_file;
    // config.json_ctagging_R += "/2023BPix/BTV/" + ctagging_R_file;
    config.json_electron += tag + "electron.json.gz";
    config.json_electron_variation =
        tag + "electronSS_EtDependent.json.gz";
    config.json_electron_hlt += tag + "electronHlt.json.gz";
    // config.json_photon += "/2023_Summer23BPix/photon.json.gz";
    config.json_jetid += tag_jme_temp + "jetid.json.gz";
    config.json_jerc += tag + "jet_jerc.json.gz";
    // config.json_jerc_fatjet += "/2023_Summer23BPix/fatJet_jerc.json.gz";
    config.json_jetvetomap += tag_jme_temp + "jetvetomaps.json.gz";
    // config.json_met += "/2023_Summer23BPix/met.json.gz";
    config.txt_roccor += "/RoccoR2023BPix.txt";
    config.golden_json +=
        "/2024/LUM/Cert_Collisions2024_378981_386951_Golden.json";

    config.onnx_hDampDown += "/2024/ONNX/mymodel12_hdamp_down_13.6TeV.onnx";
    config.onnx_hDampUp += "/2024/ONNX/mymodel12_hdamp_up_13.6TeV.onnx";
    config.onnx_toppt_reweight +=
        "/2024/ONNX/mymodel12_13TeV_MiNNLO_afterShower.onnx";
    config.onnx_rBnom += "/2024/ONNX/mymodel12_rB_nom_CP5_2M.onnx";
    config.onnx_rBUp += "/2024/ONNX/mymodel12_rB_up_CP5_2M.onnx";

    // print in red
    cout << "\033[1;31m[MyCorrection::GetEraConfig] Warning: ONNX models for "
            "TopPt reweight is for 13TeV! Please update the models for "
            "13.6TeV!\033[0m"
         << endl;

  } else {
    throw invalid_argument("[MyCorrection::GetEraConfig] Invalid era: " + era);
  }

  return config;
}

// GoldenLumi
bool MyCorrection::IsGoldenLumi(const unsigned int runNumber,
                                const unsigned int lumiSection) const {
  if(!IsDATA) return true;  
  return golden_json_parser->isGood(runNumber, lumiSection);
}

// Muon
// Rochestor correction
float MyCorrection::GetMuonScaleSF(const Muon &muon, const variation syst,
                                   const float matched_pt) const {
  float roccor = 1.;
  float roccor_err = 0.;

  // few GeVs of muon shuold not use this correction, because the authors did
  // not consider the radiations of low pt muons inside detectors still true for
  // Run3?
  if (muon.Pt() < 10.)
    return 1.;

  if (IsDATA) {
    roccor =
        rc.kScaleDT(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(), 0, 0);
    roccor_err =
        rc.kScaleDTerror(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi());
  } else {
    // Random seed is initialized in SKNanoLoader::Init()
    gRandom->SetSeed(int(muon.Pt() / muon.Eta()));
    float u = gRandom->Rndm();
    if (matched_pt > 0) { // matched
      roccor = rc.kSpreadMC(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(),
                            matched_pt, 0, 0);
      roccor_err = rc.kSpreadMCerror(muon.Charge(), muon.Pt(), muon.Eta(),
                                     muon.Phi(), matched_pt);
    } else {
      // roccor = rc.kScaleMC(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(),
      // 0, 0); roccor_err = 0.; roccor_err = rc.kScaleMCerror(muon.Charge(),
      // muon.Pt(), muon.Eta(), muon.Phi());
      roccor = rc.kSmearMC(muon.Charge(), muon.Pt(), muon.Eta(), muon.Phi(),
                           muon.nTrackerLayers(), u, 0, 0);
      roccor_err = rc.kSmearMCerror(muon.Charge(), muon.Pt(), muon.Eta(),
                                    muon.Phi(), muon.nTrackerLayers(), u);
    }
  }

  if (syst == variation::nom)
    return roccor;
  else if (syst == variation::up)
    return roccor + roccor_err;
  else if (syst == variation::down)
    return roccor - roccor_err;
  else {
    throw runtime_error("[MyCorrection::GetMuonScaleSF] Invalid syst value");
  }
}

float MyCorrection::GetMuonRECOSF(const Muon &muon,
                                  const variation syst) const {
  // No correction for Run3, see
  // https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_5_2
  if (Run == 3)
    return 1.;

  // For RECO efficiency, used 40-60 GeV muons due to the large background in
  // Z-peak. Plaetue already reached in a few GeV, okay to use for [10, 200] GeV
  // muons.
  auto cset = cset_muon->at("NUM_TrackerMuons_DEN_genTracks");
  return safeEvaluate(cset, "GetMuonRECOSF",
                      {muon.Eta(),
                       (muon.MiniAODPt() < 40. ? 40. : muon.MiniAODPt()),
                       getSystString_MUO(syst)});
}

float MyCorrection::GetMuonRECOSF(const RVec<Muon> &muons,
                                  const variation syst) const {
  float weight = 1.;
  for (const auto &muon : muons) {
    weight *= GetMuonRECOSF(muon, syst);
  }
  return weight;
}

float MyCorrection::GetMuonIDSF(const TString &Muon_ID_SF_Key, const Muon &muon,
                                const variation syst) const {
  if (Muon_ID_SF_Key == "TopHNT") {
    auto cset = cset_muon_TopHNT_idsf->at("sf");
    if (syst == variation::nom) {
      return safeEvaluate(cset, "GetMuonIDSF",
                          {fabs(muon.Eta()), muon.MiniAODPt(), "nom"});
    } else if (syst == variation::up) {
      return safeEvaluate(cset, "GetMuonIDSF",
                          {fabs(muon.Eta()), muon.MiniAODPt(), "up"});
    } else if (syst == variation::down) {
      return safeEvaluate(cset, "GetMuonIDSF",
                          {fabs(muon.Eta()), muon.MiniAODPt(), "down"});
    } else {
      throw runtime_error("[MyCorrection::GetElectronIDSF] Invalid syst value");
    }
  } else {
    auto cset = cset_muon->at(string(Muon_ID_SF_Key));
    return safeEvaluate(
        cset, "GetMuonIDSF",
        {fabs(muon.Eta()), muon.MiniAODPt(), getSystString_MUO(syst)});
  }
}

float MyCorrection::GetMuonIDSF(const TString &Muon_ID_SF_Key,
                                const RVec<Muon> &muons,
                                const variation syst) const {
  float weight = 1.;
  for (const auto &muon : muons) {
    weight *= GetMuonIDSF(Muon_ID_SF_Key, muon, syst);
  }
  return weight;
}

// Electron
// For Run2, scale uncertainty is not stored in the NanoAODv9.
// Should patch from https://github.com/cms-egamma/ScaleFactorsJSON
// https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3
float MyCorrection::GetElectronScaleUnc(const float scEta,
                                        const unsigned char seedGain,
                                        const unsigned int runNumber,
                                        const float r9, const float pt,
                                        const variation syst) const {
  if (IsDATA)
    return 1.0;

  switch (Run) {
  case 2: {
    if (syst == variation::nom)
      return 1.;
    auto cset = cset_electron_variation->at("UL-EGM_ScaleUnc");
    vector<correction::Variable::Type> args = {
        GetEra().Data(), getSystString_EGMScale(syst), scEta,
        static_cast<int>(seedGain)};
    return safeEvaluate(cset, "GetElectronScaleSF", args);
  }
  case 3: {
    if (syst == variation::nom)
      return 1.;
    const string key = (GetEra().Contains("2022"))
                           ? "Scale"
                           : EGM_keys.at(GetEra().Data()) + "_ScaleJSON";
    auto cset = cset_electron_variation->at(key);
    vector<correction::Variable::Type> args = {"total_uncertainty",
                                               static_cast<int>(seedGain),
                                               static_cast<float>(runNumber),
                                               scEta,
                                               r9,
                                               pt};
    const float unc = safeEvaluate(cset, "GetElectronScaleSF", args);
    if (syst == variation::up)
      return 1. + unc;
    else if (syst == variation::down)
      return 1. - unc;
    else
      throw runtime_error(
          "[MyCorrection::GetElectronScaleUnc] Invalid syst value");
  }
  default:
    throw runtime_error(
        "[MyCorrection::GetElectronScaleUnc] Invalid run number");
  }

  // This should never be reached, but added to avoid compiler warning
  return 1.0;
}

float MyCorrection::GetElectronSmearUnc(const Electron &electron,
                                        const variation syst,
                                        const unsigned int seed) const {
  if (IsDATA)
    return 1.0; // No smearing for data, only applied to MC
  if (Run == 2)
    throw runtime_error("[MyCorrection::GetElectronSmearUnc] Run2 is not "
                        "supported by NanoAODv9");

  const string key = (GetEra().Contains("2022"))
                         ? "Smearing"
                         : EGM_keys.at(GetEra().Data()) + "_SmearingJSON";
  auto cset = cset_electron_variation->at(key);
  vector<correction::Variable::Type> args = {"rho", electron.scEta(),
                                             electron.r9()};
  const float rho = safeEvaluate(cset, "GetElectronScaleSF", args);

  TRandom3 rng(seed);

  // Handle different variation cases
  if (syst == variation::nom) {
    // For nominal case, apply normal smearing
    return rng.Gaus(1.0, rho);
  } else if (syst == variation::up) {
    // For up variation, increase the width of the Gaussian
    return rng.Gaus(
        1.0, rho + safeEvaluate(cset, "GetElectronScaleSF",
                                {"err_rho", electron.scEta(), electron.r9()}));
  } else if (syst == variation::down) {
    // For down variation, decrease the width of the Gaussian
    return rng.Gaus(
        1.0, rho - safeEvaluate(cset, "GetElectronScaleSF",
                                {"err_rho", electron.scEta(), electron.r9()}));
  } else {
    throw runtime_error(
        "[MyCorrection::GetElectronSmearUnc] Invalid syst value");
  }
}

float MyCorrection::GetElectronRECOSF(const float eta, const float pt,
                                      const float phi,
                                      const variation syst) const {
  auto cset = Run == 3 ? cset_electron->at("Electron-ID-SF")
                       : cset_electron->at("UL-Electron-ID-SF");
  switch (Run) {
  case 2:
    if (pt < 20.)
      return GetElectronIDSF("RecoBelow20", eta, pt, phi, syst);
    else
      return GetElectronIDSF("RecoAbove20", eta, pt, phi, syst);
    break;
  case 3:
    if (pt < 20.)
      return safeEvaluate(cset, "GetElectronRECOSF",
                          {GetEra().Data() + std::string("Prompt"),
                           getSystString_EGM(syst), "RecoBelow20", eta, pt});
    else if (pt < 75.)
      return safeEvaluate(cset, "GetElectronRECOSF",
                          {GetEra().Data() + std::string("Prompt"),
                           getSystString_EGM(syst), "Reco20to75", eta, pt});
    else
      return safeEvaluate(cset, "GetElectronRECOSF",
                          {GetEra().Data() + std::string("Prompt"),
                           getSystString_EGM(syst), "RecoAbove75", eta, pt});
    break;
  default:
    throw runtime_error("[MyCorrection::GetElectronRECOSF] Invalid run number");
  }
}

float MyCorrection::GetElectronRECOSF(const RVec<Electron> &electrons,
                                      const variation syst) const {
  float weight = 1.;
  for (const auto &electron : electrons) {
    weight *= GetElectronRECOSF(fabs(electron.Eta()), electron.Pt(),
                                electron.Phi(), syst);
  }
  return weight;
}

float MyCorrection::GetElectronIDSF(const TString &Electron_ID_SF_Key,
                                    const float eta, const float pt,
                                    const float phi,
                                    const variation syst) const {
  if (Electron_ID_SF_Key == "TopHNT") {
    auto cset = cset_electron_TopHNT_idsf->at("sf");
    if (syst == variation::nom) {
      return safeEvaluate(cset, "GetElectronRECOSF", {eta, pt, "nom"});
    } else if (syst == variation::up) {
      return safeEvaluate(cset, "GetElectronRECOSF", {eta, pt, "up"});
    } else if (syst == variation::down) {
      return safeEvaluate(cset, "GetElectronRECOSF", {eta, pt, "down"});
    } else {
      throw runtime_error("[MyCorrection::GetElectronIDSF] Invalid syst value");
    }
  } else {
    // POG IDs
    string key;
    if (Run == 2)
      key = "UL-Electron-ID-SF";
    else if (Run == 3)
      key = "Electron-ID-SF";
    else
      throw runtime_error("[MyCorrection::GetElectronIDSF] Invalid run number");

    auto cset = cset_electron->at(key);
    return safeEvaluate(cset, "GetElectronIDSF",
                        {DataEra.Data() + std::string("Prompt"), getSystString_EGM(syst),
                         string(Electron_ID_SF_Key), eta,
                         pt < 999.9f ? pt : 999.9f});
  }
}

float MyCorrection::GetElectronIDSF(const TString &Electron_ID_SF_Key,
                                    const RVec<Electron> &electrons,
                                    const variation syst) const {
  float weight = 1.;
  for (const auto &electron : electrons) {
    if (Electron_ID_SF_Key == "TopHNT") {
      weight *= GetElectronIDSF(Electron_ID_SF_Key, electron.scEta(),
                                electron.Pt(), 0., syst);
    } else {
      weight *= GetElectronIDSF(Electron_ID_SF_Key, fabs(electron.Eta()),
                                electron.Pt(), electron.Phi(), syst);
    }
  }
  return weight;
}

// Trigger
float MyCorrection::GetMuonTriggerEff(const TString &Muon_Trigger_Eff_Key,
                                      const float eta, const float pt,
                                      const bool isData,
                                      const variation syst) const {
  static bool warned_missing_trig_eff = false;
  if (!cset_muon_trig_eff) {
    if (!warned_missing_trig_eff) {
      cerr << "[MyCorrection::GetMuonTriggerEff] Warning: trigger efficiency "
              "correction set is not loaded, returning 1."
           << endl;
      warned_missing_trig_eff = true;
    }
    return 1.;
  }

  correction::Correction::Ref cset;
  try {
    cset = cset_muon_trig_eff->at(string(Muon_Trigger_Eff_Key));
  } catch (const std::out_of_range &e) {
    if (!warned_missing_trig_eff) {
      cerr << "[MyCorrection::GetMuonTriggerEff] Warning: key "
           << Muon_Trigger_Eff_Key
           << " not found in trigger efficiency set, returning 1. (" << e.what()
           << ")" << endl;
      warned_missing_trig_eff = true;
    }
    return 1.;
  }
  if (isData)
    return safeEvaluate(cset, "GetTriggerEff",
                        {"data", getSystString_MUO(syst), eta, pt});
  else
    return safeEvaluate(cset, "GetTriggerEff",
                        {"mc", getSystString_MUO(syst), eta, pt});
}

float MyCorrection::GetMuonTriggerSF(const TString &Muon_Trigger_SF_Key,
                                     const RVec<Muon> &muons,
                                     const variation syst) const {
  if (!cset_muon_trig_eff) {
    float weight = 1.;
    for (const auto &muon : muons) {
      weight *= GetMuonTriggerSF(Muon_Trigger_SF_Key, muon, syst);
    }
    return weight;
  }

  float eff_data = 1.;
  float eff_mc = 1.;
  float weight = 1.;
  for (const auto &muon : muons) {
    eff_data *= 1. - GetMuonTriggerEff(Muon_Trigger_SF_Key, fabs(muon.Eta()),
                                       muon.Pt(), true, syst);
    eff_mc *= 1. - GetMuonTriggerEff(Muon_Trigger_SF_Key, fabs(muon.Eta()),
                                     muon.Pt(), false, syst);
  }
  weight = (1. - eff_data) / (1. - eff_mc);
  return weight;
}

float MyCorrection::GetMuonTriggerSF(const TString &Muon_Trigger_SF_Key,
                                     const Muon muon,
                                     const variation syst) const {
  if (IsDATA)
    return 1.0;

  float weight = 1.f;
  auto try_eval = [&](const unique_ptr<CorrectionSet> &set) -> bool {
    if (!set)
      return false;
    try {
      auto cset = set->at(Muon_Trigger_SF_Key.Data());
      weight = safeEvaluate(cset, "GetMuonTriggerSF",
                            {muon.Eta(),
                             muon.MiniAODPt() > 26.f ? muon.MiniAODPt() : 26.f,
                             getSystString_MUO(syst)});
      return true;
    } catch (const std::out_of_range &) {
      return false;
    }
  };

  if (try_eval(cset_muon_trig_sf))
    return weight;
  if (try_eval(cset_muon))
    return weight;

  static bool warned_missing_trig_sf = false;
  if (!warned_missing_trig_sf) {
    cerr << "[MyCorrection::GetMuonTriggerSF] Warning: trigger SF key "
         << Muon_Trigger_SF_Key << " not found for era " << DataEra
         << ", returning 1." << endl;
    warned_missing_trig_sf = true;
  }
  return 1.0;
}

float MyCorrection::GetElectronTriggerEff(
    const TString &Electron_Trigger_SF_Key, const float eta, const float pt,
    const float phi, const bool isDATA, const variation syst) const {
  string key = isDATA ? "Electron-HLT-DataEff" : "Electron-HLT-McEff";
  auto cset = cset_electron_hlt->at(key);
  // hardcoded replacemet
  string ValType = getSystString_EGM(syst);
  if (ValType == "sf")
    ValType = "nom";
  else if (ValType == "sfup")
    ValType = "up";
  else if (ValType == "sfdown")
    ValType = "down";
  else
    throw runtime_error(
        "[MyCorrection::GetElectronTriggerEff] Invalid syst value");
  try {
    if (!isInputInCorrection("phi", cset)) {
      return safeEvaluate(cset, "GetTriggerEff",
                          {EGM_keys.at(GetEra().Data()), ValType,
                           string(Electron_Trigger_SF_Key), eta, pt});
    } else {
      return safeEvaluate(cset, "GetTriggerEff",
                          {EGM_keys.at(GetEra().Data()), ValType,
                           string(Electron_Trigger_SF_Key), eta, pt, phi});
    }
  } catch (exception &e) {
    cerr << "[MyCorrection::GetElectronTriggerEff] " << e.what() << endl;
    throw e;
  }
}

float MyCorrection::GetElectronTriggerSF(const TString &Electron_Trigger_SF_Key,
                                         const float eta, const float pt,
                                         const float phi,
                                         const variation syst) const {
  auto cset = cset_electron_hlt->at("Electron-HLT-SF");
  try {
    if (!isInputInCorrection("phi", cset)) {
      return safeEvaluate(cset, "GetTriggerEff",
                          {EGM_keys.at(GetEra().Data()),
                           getSystString_EGM(syst),
                           string(Electron_Trigger_SF_Key), eta, pt});
    } else {
      return safeEvaluate(cset, "GetTriggerEff",
                          {EGM_keys.at(GetEra().Data()),
                           getSystString_EGM(syst),
                           string(Electron_Trigger_SF_Key), eta, pt, phi});
    }
  } catch (exception &e) {
    cerr << "[MyCorrection::GetElectronTriggerSF] " << e.what() << endl;
    throw e;
  }
}

// double lepton triggers
// This function is used for leptons passing TopHNT ID
float MyCorrection::GetTriggerEff(const Muon &muon,
                                  const TString &trigger_leg_key,
                                  const bool isData,
                                  const variation syst) const {
  if (trigger_leg_key == "DblMu_Mu17Leg") {
    const string jsonkey = isData ? "data" : "sim";
    auto cset = cset_muon_TopHNT_dblmu_leg1_eff->at(jsonkey);
    float eff = safeEvaluate(
        cset, "GetTriggerEff",
        {fabs(muon.Eta()), muon.MiniAODPt(), getSystString_CUSTOM(syst)});
    return eff < 1. ? eff : 1.;
  } else if (trigger_leg_key == "DblMu_Mu8Leg") {
    const string jsonkey = isData ? "data" : "sim";
    auto cset = cset_muon_TopHNT_dblmu_leg2_eff->at(jsonkey);
    float eff = safeEvaluate(
        cset, "GetTriggerEff",
        {fabs(muon.Eta()), muon.MiniAODPt(), getSystString_CUSTOM(syst)});
    return eff < 1. ? eff : 1.;
  } else if (trigger_leg_key == "EMu_Mu23Leg") {
    const string jsonkey = isData ? "Mu23El12_Data" : "Mu23El12_MC";
    auto cset = cset_muon_TopHNT_emu_leg1_eff->at(jsonkey);
    float eff = safeEvaluate(
        cset, "GetTriggerEff",
        {fabs(muon.Eta()), muon.MiniAODPt(), getSystString_CUSTOM(syst)});
    return eff < 1. ? eff : 1.;
  } else if (trigger_leg_key == "EMu_Mu8Leg") {
    const string jsonkey = isData ? "Mu8El23_Data" : "Mu8El23_MC";
    auto cset = cset_muon_TopHNT_emu_leg2_eff->at(jsonkey);
    float eff = safeEvaluate(
        cset, "GetTriggerEff",
        {fabs(muon.Eta()), muon.MiniAODPt(), getSystString_CUSTOM(syst)});
    return eff < 1. ? eff : 1.;
  } else {
    throw runtime_error(
        "[MyCorrection::GetTriggerEff] Invalid trigger leg key");
  }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetTriggerEff(const Electron &electron,
                                  const TString &trigger_leg_key,
                                  const bool isData,
                                  const variation syst) const {
  if (trigger_leg_key == "EMu_El23Leg") {
    const string jsonkey = isData ? "Mu8El23_Data" : "Mu8El23_MC";
    auto cset = cset_electron_TopHNT_emu_leg1_eff->at(jsonkey);
    float eff = safeEvaluate(
        cset, "GetTriggerEff",
        {fabs(electron.scEta()), electron.Pt(), getSystString_CUSTOM(syst)});
    return eff < 1. ? eff : 1.;
  } else if (trigger_leg_key == "EMu_El12Leg") {
    const string jsonkey = isData ? "Mu23El12_Data" : "Mu23El12_MC";
    auto cset = cset_electron_TopHNT_emu_leg2_eff->at(jsonkey);
    float eff = safeEvaluate(
        cset, "GetTriggerEff",
        {fabs(electron.scEta()), electron.Pt(), getSystString_CUSTOM(syst)});
    return eff < 1. ? eff : 1.;
  } else {
    throw runtime_error(
        "[MyCorrection::GetTriggerEff] Invalid trigger leg key");
  }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetDblMuTriggerEff(const RVec<Muon> &muons,
                                       const bool isData,
                                       const variation syst) const {
  const float filter_eff = GetPairwiseFilterEff("DblMu", isData);
  if (muons.size() == 2) {
    float leg1_eff = GetTriggerEff(muons.at(0), "DblMu_Mu17Leg", isData, syst);
    float leg2_eff = GetTriggerEff(muons.at(1), "DblMu_Mu8Leg", isData, syst);
    return leg1_eff * leg2_eff * filter_eff;
  } else if (muons.size() == 3) {
    const auto &mu1 = muons.at(0);
    const auto &mu2 = muons.at(1);
    const auto &mu3 = muons.at(2);

    float case1 = GetTriggerEff(mu1, "DblMu_Mu17Leg", isData, syst);
    case1 *= GetTriggerEff(mu2, "DblMu_Mu8Leg", isData, syst);
    case1 *= filter_eff;
    float case2 = 1. - GetTriggerEff(mu1, "DblMu_Mu17Leg", isData, syst);
    case2 *= GetTriggerEff(mu2, "DblMu_Mu17Leg", isData, syst);
    case2 *= GetTriggerEff(mu3, "DblMu_Mu8Leg", isData, syst);
    case2 *= filter_eff;
    float case3 = GetTriggerEff(mu1, "DblMu_Mu17Leg", isData, syst);
    case3 *=
        (1. - GetTriggerEff(mu2, "DblMu_Mu8Leg", isData, syst)) * filter_eff;
    case3 *= GetTriggerEff(mu3, "DblMu_Mu8Leg", isData, syst) * filter_eff;
    return case1 + case2 + case3;
  } else {
    throw runtime_error("[MyCorrection::GetDblMuTriggerEff] Invalid muon and "
                        "electron size configuration");
  }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetDblMuTriggerSF(const RVec<Muon> &muons,
                                      const variation syst) const {
  float eff_data = -999.;
  float eff_mc = -999.;
  if (syst == variation::nom) {
    eff_data = GetDblMuTriggerEff(muons, true, syst);
    eff_mc = GetDblMuTriggerEff(muons, false, syst);
  } else if (syst == variation::up) {
    eff_data = GetDblMuTriggerEff(muons, true, variation::up);
    eff_mc = GetDblMuTriggerEff(muons, false, variation::down);
  } else if (syst == variation::down) {
    eff_data = GetDblMuTriggerEff(muons, true, variation::down);
    eff_mc = GetDblMuTriggerEff(muons, false, variation::up);
  } else {
    throw runtime_error("[MyCorrection::GetDblMuTriggerSF] Invalid syst value");
  }
  return eff_mc > 0. ? eff_data / eff_mc : 1.;
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetEMuTriggerEff(const RVec<Electron> &electrons,
                                     const RVec<Muon> &muons, const bool isData,
                                     const variation syst) const {
  const float filter_eff = GetPairwiseFilterEff("EMu", isData);
  if (electrons.size() == 1 && muons.size() == 1) {
    const auto &el = electrons.at(0);
    const auto &mu = muons.at(0);
    float eff_el = mu.Pt() > 25.
                       ? GetTriggerEff(el, "EMu_El12Leg", isData, syst)
                       : GetTriggerEff(el, "EMu_El23Leg", isData, syst);
    float eff_mu = el.Pt() > 25.
                       ? GetTriggerEff(mu, "EMu_Mu8Leg", isData, syst)
                       : GetTriggerEff(mu, "EMu_Mu23Leg", isData, syst);
    return eff_el * eff_mu * filter_eff;
  } else if (electrons.size() == 1 && muons.size() == 2) {
    const auto &el = electrons.at(0);
    const auto &mu1 = muons.at(0);
    const auto &mu2 = muons.at(1);

    float case1 = GetTriggerEff(mu1, "EMu_Mu8Leg", isData, syst) +
                  (1. - GetTriggerEff(mu1, "EMu_Mu8Leg", isData, syst)) *
                      GetTriggerEff(mu2, "EMu_Mu8Leg", isData, syst) *
                      filter_eff;
    float case2 =
        mu2.Pt() > 25.
            ? GetTriggerEff(mu1, "EMu_Mu23Leg", isData, syst) +
                  (1. - GetTriggerEff(mu1, "EMu_Mu23Leg", isData, syst)) *
                      GetTriggerEff(mu2, "EMu_Mu23Leg", isData, syst) *
                      filter_eff
            : GetTriggerEff(mu1, "EMu_Mu23Leg", isData, syst);
    float eff_el = (mu1.Pt() > 25. || mu2.Pt() > 25.)
                       ? GetTriggerEff(el, "EMu_El12Leg", isData, syst)
                       : GetTriggerEff(el, "EMu_El23Leg", isData, syst);
    float eff_mu = el.Pt() > 25. ? case1 : case2;
    return eff_el * eff_mu * filter_eff;
  } else {
    throw runtime_error("[MyCorrection::GetEMuTriggerEff] Invalid electron and "
                        "muon size configuration");
  }
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetEMuTriggerSF(const RVec<Electron> &electrons,
                                    const RVec<Muon> &muons,
                                    const variation syst) const {
  float eff_data = -999.;
  float eff_mc = -999.;
  if (syst == variation::nom) {
    eff_data = GetEMuTriggerEff(electrons, muons, true, syst);
    eff_mc = GetEMuTriggerEff(electrons, muons, false, syst);
  } else if (syst == variation::up) {
    eff_data = GetEMuTriggerEff(electrons, muons, true, variation::up);
    eff_mc = GetEMuTriggerEff(electrons, muons, false, variation::down);
  } else if (syst == variation::down) {
    eff_data = GetEMuTriggerEff(electrons, muons, true, variation::down);
    eff_mc = GetEMuTriggerEff(electrons, muons, false, variation::up);
  } else {
    throw runtime_error("[MyCorrection::GetEMuTriggerSF] Invalid syst value");
  }
  return eff_mc > 0. ? eff_data / eff_mc : 1.;
}

// This function is used for leptons passing TopHNT ID
float MyCorrection::GetPairwiseFilterEff(const TString &filter_name,
                                         const bool isData) const {
  if (filter_name.Contains("DblMu")) {
    if (GetEra() == "2016preVFP") {
      return 1.;
    } else if (GetEra() == "2016postVFP") {
      return isData ? 0.9798 : 0.9968;
    } else if (GetEra() == "2017") {
      return isData ? 0.9961 : 0.9958;
    } else if (GetEra() == "2018") {
      return isData ? 0.9988 : 0.9998;
    } else if (GetEra() == "2023") {
      return isData ? 0.9993 : 0.9998;
    } else {
      cerr << "[MyCorrection::GetPairwiseFilterEff] " << filter_name
           << " is not implemented for " << GetEra() << endl;
      return 1.;
    }
  } else if (filter_name.Contains("EMu")) {
    if (GetEra() == "2016preVFP") {
      return 1.;
    } else if (GetEra() == "2016postVFP") {
      return isData ? 0.9638 : 0.9878;
    } else if (GetEra() == "2017") {
      return isData ? 0.9989 : 0.9955;
    } else if (GetEra() == "2018") {
      return isData ? 0.9946 : 0.9981;
    } else if (GetEra() == "2023") {
      return isData ? 0.9944 : 0.9976;
    } else {
      cerr << "[MyCorrection::GetPairwiseFilterEff] " << filter_name
           << " is not implemented for " << GetEra() << endl;
      return 1.;
    }
  } else {
    throw runtime_error(
        "[MyCorrection::GetPairwiseFilterEff] Invalid filter name");
  }
}

// Pileup
float MyCorrection::GetPUWeight(const float nTrueInt, const variation syst,
                                const TString &source) const {
  // nota bene: Input should be nTrueInt, not nPileUp
  correction::Correction::Ref cset = nullptr;
  cset = cset_puWeights->at(LUM_keys.at(GetEra().Data()));
  try {
    return safeEvaluate(cset, "GetPUWeight",
                        {nTrueInt, getSystString_LUM(syst)});
  } catch (exception &e) {
    cerr << "[MyCorrection::GetPUWeight] " << e.what() << endl;
    return 1.;
  }
}

// Heavy flavor tagging
void MyCorrection::SetTaggingParam(JetTagging::JetFlavTagger tagger,
                                   JetTagging::JetFlavTaggerWP wp) {
  global_tagger = tagger;
  global_wp = wp;
  global_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger);
  global_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp);
}

float MyCorrection::GetBTaggingWP() const {
  try {
    correction::Correction::Ref cset =
        cset_btagging->at(global_taggerStr + "_wp_values");
    return safeEvaluate(cset, "GetBTaggingWP", {global_wpStr});
  } catch (const exception &e) {
    cerr << "[Correction::GetBTaggingWP] Warning: Failed to evaluate WP '"
         << global_wpStr << "' for tagger '" << global_taggerStr << endl;
    throw runtime_error(e.what());
    return 1.f;
  }
}

float MyCorrection::GetBTaggingWP(JetTagging::JetFlavTagger tagger,
                                  JetTagging::JetFlavTaggerWP wp) const {
  // Convert enumerations to strings
  string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();

  try {
    correction::Correction::Ref cset =
        cset_btagging->at(this_taggerStr + "_wp_values");
    return safeEvaluate(cset, "GetBTaggingWP", {this_wpStr});
  } catch (const exception &e) {
    cerr << "[Correction::GetBTaggingWP] Warning: Failed to evaluate WP '"
         << this_wpStr << "' for tagger '" << this_taggerStr << endl;
    throw runtime_error(e.what());
    return 1.f;
  }
}

float MyCorrection::GetBTaggingEff(const float eta, const float pt,
                                   const int flav,
                                   JetTagging::JetFlavTagger tagger,
                                   JetTagging::JetFlavTaggerWP wp,
                                   const variation syst) {
  string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
  auto cset = cset_btagging_eff->at(this_taggerStr);
  return safeEvaluate(cset, "GetBTaggingSF",
                      {"central", this_wpStr, flav, fabs(eta), pt});
}

float MyCorrection::GetBTaggingSF(const RVec<Jet> &jets,
                                  const JetTagging::JetFlavTagger tagger,
                                  const JetTagging::JetFlavTaggerWP wp,
                                  const JetTagging::JetTaggingSFMethod method,
                                  const variation syst, const TString &source) {
  if (Run == 2 && tagger != JetTagging::JetFlavTagger::DeepJet) {
    cerr << "[MyCorrection::GetBTaggingSF] DeepJet is the only supported "
            "tagger for 2016preVFP, 2016postVFP, 2017, 2018, and 2018UL"
         << endl;
    return 1.;
  }

  const std::size_t nJets = jets.size();
  if (nJets == 0)
    return 1.f;

  const string this_taggerStr =
      JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  const string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
  string syst_str = getSystString_BTV(syst);
  if (source != "total")
    syst_str = getSystString_BTV(syst) + "_" + source;
  const string nominal_str = getSystString_BTV(MyCorrection::variation::nom);
  const unordered_set<string> c_flav_source = {"cferr1", "cferr2"};
  const bool source_targets_cferr =
      c_flav_source.find(source.Data()) != c_flav_source.end();

  std::vector<std::size_t> indices(nJets);
  std::iota(indices.begin(), indices.end(), 0);
  std::vector<int> flavours(nJets, 0);
  std::vector<float> absEtas(nJets, 0.f);
  std::vector<float> pts(nJets, 0.f);
  std::vector<float> scores(nJets, -1.f);

  std::for_each(std::execution::unseq, indices.begin(), indices.end(),
                [&](std::size_t idx) {
                  const auto &jet = jets[idx];
                  const int hadFlavour = std::abs(jet.hadronFlavour());
                  if (hadFlavour == 5 || hadFlavour == 4)
                    flavours[idx] = hadFlavour;
                  absEtas[idx] = std::fabs(jet.Eta());
                  pts[idx] = jet.Pt();
                  scores[idx] = jet.GetTaggerResult(
                      tagger, JetTagging::JetFlavTaggerScoreType::B);
                });

  float weight = 1.f;
  if (method == JetTagging::JetTaggingSFMethod::shape) {
    auto cset =
        cset_btagging->at(this_taggerStr + "_" +
                          JetTagging::GetJetTaggingSFMethodStr(method).Data());
    const string syst_for_cjets = source_targets_cferr ? syst_str : nominal_str;
    const string syst_for_non_cjets =
        source_targets_cferr ? nominal_str : syst_str;

    for (std::size_t idx = 0; idx < nJets; ++idx) {
      if (scores[idx] < 0.f)
        continue;
      const int flav = flavours[idx];
      const string &requested_syst =
          (flav == 4) ? syst_for_cjets : syst_for_non_cjets;
      weight *= safeEvaluate(
          cset, "GetBTaggingSF",
          {requested_syst, flav, absEtas[idx], pts[idx], scores[idx]});
    }
    return weight;
  } else if (method == JetTagging::JetTaggingSFMethod::comb ||
             method == JetTagging::JetTaggingSFMethod::mujets) {
    auto cset =
        cset_btagging->at(this_taggerStr + "_" +
                          JetTagging::GetJetTaggingSFMethodStr(method).Data());
    string light_str;
    if (Run == 2)
      light_str = this_taggerStr + "_incl";
    else if (Run == 3)
      light_str = this_taggerStr + "_light";
    auto cset_light = cset_btagging->at(light_str);
    const float this_cut = GetBTaggingWP(tagger, wp);

    std::vector<char> passesWP(nJets, 0);
    std::vector<float> efficiencies(nJets, 0.f);

    std::for_each(std::execution::unseq, indices.begin(), indices.end(),
                  [&](std::size_t idx) {
                    passesWP[idx] = static_cast<char>(scores[idx] > this_cut);
                  });

    std::transform(indices.begin(), indices.end(), efficiencies.begin(),
                   [&](std::size_t idx) {
                     return GetBTaggingEff(absEtas[idx], pts[idx],
                                           flavours[idx], tagger, wp, syst);
                   });

    for (std::size_t idx = 0; idx < nJets; ++idx) {
      const bool is_heavy = flavours[idx] == 5 || flavours[idx] == 4;
      auto this_cset = is_heavy ? cset : cset_light;
      const int flav = is_heavy ? flavours[idx] : 0;
      //XXX up_jes, down_jes not supported for light jet
      const float sf =
          safeEvaluate(this_cset, "GetBTaggingSF",
                       {syst_str, this_wpStr, flav, absEtas[idx], pts[idx]});
      if (passesWP[idx]) {
        weight *= sf;
      } else {
        const float eff = efficiencies[idx];
        weight *= (1.f - eff * sf) / (1.f - eff);
      }
    }
    return weight;
  } else {
    cout << "[MyCorrection::GetBTaggingSF] method "
         << JetTagging::GetJetTaggingSFMethodStr(method)
         << " is not implemented" << endl;
    exit(ENODATA);
  }
}

float MyCorrection::GetBTaggingR(const RVec<Jet> &jets,
                                 const JetTagging::JetFlavTagger tagger,
                                 string &processName, const variation syst,
                                 const TString &source) const {
  string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  string syst_str = getSystString_BTV(syst);

  if (source != "total")
    syst_str = getSystString_BTV(syst) + "_" + source;
  float weight = 1.;

  this_taggerStr += (string("_") + processName);

  // First, check cset_btagging_R is loaded or nullptr
  if (!cset_btagging_R)
    throw std::runtime_error("[MyCorrection::GetBTaggingR] cset_btagging_R is "
                             "not loaded or nullptr");
  auto cset = cset_btagging_R->at(this_taggerStr);

  for (const auto &jet : jets) {
    int this_flav = 0;
    if (abs(jet.hadronFlavour()) == 5)
      this_flav = 5;
    else if (abs(jet.hadronFlavour()) == 4)
      this_flav = 4;

    weight *= safeEvaluate(cset, "GetBTaggingEff",
                           {syst_str, this_flav, jet.Pt(), fabs(jet.Eta())});
  }
  return weight;
}

pair<float, float> MyCorrection::GetCTaggingWP() const {
  try {
    correction::Correction::Ref cset =
        cset_ctagging->at(global_taggerStr + "_wp_values");
    float valCvB = safeEvaluate(cset, "GetCTaggingWP", {global_wpStr, "CvB"});
    float valCvL = safeEvaluate(cset, "GetCTaggingWP", {global_wpStr, "CvL"});
    return make_pair(valCvB, valCvL);
  } catch (const exception &e) {
    // If the requested WP is not found or any other error occurs,
    // log a warning and return (1.f, 1.f) as a fallback.
    cerr << "[Correction::GetCTaggingWP] Warning: Failed to evaluate WP '"
         << global_wpStr << "' for tagger '" << global_taggerStr << endl;
    throw runtime_error(e.what());
    return make_pair(1.f, 1.f);
  }
}

pair<float, float>
MyCorrection::GetCTaggingWP(JetTagging::JetFlavTagger tagger,
                            JetTagging::JetFlavTaggerWP wp) const {
  // Convert enumerations to strings using your existing utility functions
  string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();

  // Retrieve the relevant correction set
  correction::Correction::Ref cset =
      cset_ctagging->at(this_taggerStr + "_wp_values");

  try {
    // Evaluate the corrections. If the WP does not exist, an exception might be
    // thrown
    float valCvB = safeEvaluate(cset, "GetCTaggingWP", {this_wpStr, "CvB"});
    float valCvL = safeEvaluate(cset, "GetCTaggingWP", {this_wpStr, "CvL"});

    // If everything is fine, return the pair
    return make_pair(valCvB, valCvL);
  } catch (const exception &e) {
    // In case the WP is not found (or any other error occurs),
    // print a warning (optional) and return default values
    cerr << "[Correction::GetCTaggingWP] Warning: WP '" << this_wpStr
         << "' not found for tagger '" << this_taggerStr << endl;
    throw runtime_error(e.what());
    return make_pair(1.f, 1.f);
  }
}

float MyCorrection::GetCTaggingEff(const float eta, const float pt,
                                   const int flav,
                                   JetTagging::JetFlavTagger tagger,
                                   JetTagging::JetFlavTaggerWP wp,
                                   const variation syst) {
  return 1.;
  string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
  correction::Correction::Ref cset = cset_btagging_eff->at(this_taggerStr);
  return safeEvaluate(cset, "GetBTaggingEff",
                      {
                          getSystString_BTV(syst),
                          this_wpStr,
                      });
}

float MyCorrection::GetCTaggingSF(const RVec<Jet> &jets,
                                  const JetTagging::JetFlavTagger tagger,
                                  const JetTagging::JetFlavTaggerWP wp,
                                  const JetTagging::JetTaggingSFMethod method,
                                  const variation syst, const TString &source) {
  if (Run == 2 && tagger != JetTagging::JetFlavTagger::DeepJet) {
    cerr << "[MyCorrection::GetCTaggingSF] Run2 only supports DeepJet tagger"
         << endl;
    return 1.;
  }
  if (Run == 3) {
    // cerr << "[MyCorrection::GetCTaggingSF] Run3 C-Tagger SF not provided by
    // POG yet" << endl;
    return 1.;
  }
  string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  string this_wpStr = JetTagging::GetTaggerCorrectionWPStr(wp).Data();
  string syst_str = getSystString_BTV(syst);
  if (source != "total")
    syst_str = getSystString_BTV(syst) + "_" + source;
  float weight = 1.;
  // method is comb, mujets, or shape
  if (method == JetTagging::JetTaggingSFMethod::shape) {
    auto cset =
        cset_ctagging->at(this_taggerStr + "_" +
                          JetTagging::GetJetTaggingSFMethodStr(method).Data());
    for (const auto &jet : jets) {
      int this_flav = 0;
      if (abs(jet.hadronFlavour()) == 5)
        this_flav = 5;
      else if (abs(jet.hadronFlavour()) == 4)
        this_flav = 4;
      pair<float, float> this_ctag = make_pair(
          jet.GetTaggerResult(tagger, JetTagging::JetFlavTaggerScoreType::CvB),
          jet.GetTaggerResult(tagger, JetTagging::JetFlavTaggerScoreType::CvL));
      if (this_ctag.first < 0.f || this_ctag.second < 0.f)
        continue; // defaulted jet
      weight *= safeEvaluate(
          cset, "GetCTaggingSF",
          {syst_str, this_flav, this_ctag.second, this_ctag.first});
    }
    return weight;
  } else if (method == JetTagging::JetTaggingSFMethod::wp) {
    auto cset =
        cset_btagging->at(this_taggerStr + "_" +
                          JetTagging::GetJetTaggingSFMethodStr(method).Data());
    auto cset_light = cset_btagging->at(this_taggerStr + "_incl");
    pair<float, float> this_cut = GetCTaggingWP(tagger, wp);
    pair<float, float> this_score;
    for (const auto &jet : jets) {
      this_score = make_pair(
          jet.GetTaggerResult(tagger, JetTagging::JetFlavTaggerScoreType::CvB),
          jet.GetTaggerResult(tagger, JetTagging::JetFlavTaggerScoreType::CvL));
      int this_flav = 0;
      if (abs(jet.hadronFlavour()) == 5)
        this_flav = 5;
      else if (abs(jet.hadronFlavour()) == 4)
        this_flav = 4;
      float eff =
          GetCTaggingEff(jet.Eta(), jet.Pt(), this_flav, tagger, wp, syst);
      float sf;
      if (this_flav == 0)
        sf = cset_light->evaluate(
            {syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
      else
        sf = safeEvaluate(
            cset, "GetCTaggingSF",
            {syst_str, this_wpStr, this_flav, fabs(jet.Eta()), jet.Pt()});
      if (this_score.first > this_cut.first &&
          this_score.second > this_cut.second) {
        weight *= sf;
      } else {
        weight *= (1. - eff * sf) / (1. - eff);
      }
    }
    return weight;
  } else {
    cout << "[MyCorrection::GetCTaggingSF] method "
         << JetTagging::GetJetTaggingSFMethodStr(method)
         << " is not implemented" << endl;
    exit(ENODATA);
  }
}

float MyCorrection::GetCTaggingR(const float nTrueInt, const float HT,
                                 const JetTagging::JetFlavTagger tagger,
                                 const TString &processName,
                                 const TString &ttBarCategory,
                                 const TString &syst_str) const {
  string this_taggerStr = JetTagging::GetTaggerCorrectionLibStr(tagger).Data();
  if (processName != "")
    this_taggerStr += "_" + processName;
  else
    this_taggerStr += (string("_") + Sample.Data());
  auto cset = cset_ctagging_R->at(this_taggerStr);
  return safeEvaluate(cset, "GetTopPtReweight",
                      {syst_str.Data(), ttBarCategory.Data(), nTrueInt, HT});
}

// Pileup Jet ID
float MyCorrection::GetPileupJetIDSF(const RVec<Jet> &jets,
                                     const unordered_map<int, int> &matched_idx,
                                     const TString &wp, const variation syst) {
  // Should pass jets after PUID, no mistag rate correction
  if (Run == 3)
    return 1.;

  float weight = 1.;
  auto cset = cset_jmar->at(JME_PILEUP_keys.at(GetEra().Data()));
  string wp_str;
  if (wp == "tight")
    wp_str = "T";
  else if (wp == "medium")
    wp_str = "M";
  else if (wp == "loose")
    wp_str = "L";
  else
    throw std::invalid_argument("[MCCorrection::GetPileupJetIDSF] Invalid WP");

  for (int i = 0; i < jets.size(); i++) {
    if (jets.at(i).Pt() > 50.)
      continue;
    if (matched_idx.find(i) == matched_idx.end() || matched_idx.at(i) < 0)
      continue; // not matched
    float this_eff = safeEvaluate(
        cset, "GetJetVetoMapEff",
        {jets.at(i).Eta(), jets.at(i).Pt(), getSystString_JME(syst), wp_str});
    weight *= this_eff;
  }
  return weight;
}

// JetID

bool MyCorrection::PassJetID(const Jet &jet, const Jet::JetID &id) const {
  correction::Correction::Ref cset = nullptr;
  float out;
  float puid_score = 1;
  float puid_WP = -1;
  switch (id) {
  case Jet::JetID::TIGHT:
    cset = cset_jetid->at("AK4PUPPI_Tight");
    out = cset->evaluate(
        {fabs(jet.Eta()), jet.chHEF(), jet.neHEF(), jet.chEmEF(), jet.neEmEF(),
         jet.muEF(), static_cast<int>(jet.chMultiplicity()),
         static_cast<int>(jet.neMultiplicity()),
         static_cast<int>(jet.chMultiplicity() + jet.neMultiplicity())});
    return out > 0.5; // return is real
    break;
  case Jet::JetID::TIGHTLEPVETO:
    cset = cset_jetid->at("AK4PUPPI_TightLeptonVeto");
    out = cset->evaluate(
        {fabs(jet.Eta()), jet.chHEF(), jet.neHEF(), jet.chEmEF(), jet.neEmEF(),
         jet.muEF(), static_cast<int>(jet.chMultiplicity()),
         static_cast<int>(jet.neMultiplicity()),
         static_cast<int>(jet.chMultiplicity() + jet.neMultiplicity())});
    return out > 0.5; // return is real
    break;
  case Jet::JetID::PUID_LOOSE:
    puid_score = jet.PuIdDisc();
    // WP definition for |eta| < 2.5
    if (jet.Pt() < 10){
      puid_WP = -1;
    }
    else if (jet.Pt() < 20){
      puid_WP = -0.51;
    }
    else if (jet.Pt() < 30){
      puid_WP = -0.72;
    }
    else if (jet.Pt() < 40){
      puid_WP = -0.60;
    }
    else if (jet.Pt() < 50){
      puid_WP = -0.40;
    }
    else{
      puid_WP = -1;
    }
    cset = cset_jetid->at("AK4PUPPI_Tight");
    out = cset->evaluate(
        {fabs(jet.Eta()), jet.chHEF(), jet.neHEF(), jet.chEmEF(), jet.neEmEF(),
         jet.muEF(), static_cast<int>(jet.chMultiplicity()),
         static_cast<int>(jet.neMultiplicity()),
         static_cast<int>(jet.chMultiplicity() + jet.neMultiplicity())});
    return ((out > 0.5) && (puid_score > puid_WP)); // return is real
    break;
  case Jet::JetID::PUID_MEDIUM:
    return true;
  case Jet::JetID::PUID_TIGHT:
    return true;
  case Jet::JetID::NOCUT:
    // No cut, always return true
    return true;
    break;
  default:
    throw runtime_error("[MyCorrection::PassJetID] Invalid JetID type");
  }
}

bool MyCorrection::PassJetID(const JetView &jet, const Jet::JetID &id) const {
  correction::Correction::Ref cset = nullptr;
  float out = 0.f;
  float puid_score = 1;
  float puid_WP = -1;
  switch (id) {
  case Jet::JetID::TIGHT:
    cset = cset_jetid->at("AK4PUPPI_Tight");
    out = cset->evaluate(
        {fabs(jet.Eta()), jet.ChHEF(), jet.NeHEF(), jet.ChEmEF(), jet.NeEmEF(),
         jet.MuEF(), static_cast<int>(jet.ChMultiplicity()),
         static_cast<int>(jet.NeMultiplicity()),
         static_cast<int>(jet.ChMultiplicity() + jet.NeMultiplicity())});
    return out > 0.5;
  case Jet::JetID::TIGHTLEPVETO:
    cset = cset_jetid->at("AK4PUPPI_TightLeptonVeto");
    out = cset->evaluate(
        {fabs(jet.Eta()), jet.ChHEF(), jet.NeHEF(), jet.ChEmEF(), jet.NeEmEF(),
         jet.MuEF(), static_cast<int>(jet.ChMultiplicity()),
         static_cast<int>(jet.NeMultiplicity()),
         static_cast<int>(jet.ChMultiplicity() + jet.NeMultiplicity())});
    return out > 0.5;
  case Jet::JetID::PUID_LOOSE:
    puid_score = jet.PuIdDisc();
    // WP definition for |eta| < 2.5
    if (jet.Pt() < 10){
      puid_WP = -1;
    }
    else if (jet.Pt() < 20){
      puid_WP = -0.51;
    }
    else if (jet.Pt() < 30){
      puid_WP = -0.72;
    }
    else if (jet.Pt() < 40){
      puid_WP = -0.60;
    }
    else if (jet.Pt() < 50){
      puid_WP = -0.40;
    }
    else{
      puid_WP = -1;
    }
    cset = cset_jetid->at("AK4PUPPI_Tight");
    out = cset->evaluate(
        {fabs(jet.Eta()), jet.ChHEF(), jet.NeHEF(), jet.ChEmEF(), jet.NeEmEF(),
         jet.MuEF(), static_cast<int>(jet.ChMultiplicity()),
         static_cast<int>(jet.NeMultiplicity()),
         static_cast<int>(jet.ChMultiplicity() + jet.NeMultiplicity())});
    return ((out > 0.5) && (puid_score > puid_WP)); // return is real
    break;

  case Jet::JetID::NOCUT:
    return true;
  default:
    throw runtime_error("[MyCorrection::PassJetID] Invalid JetID type");
  }
}

bool MyCorrection::PassFatJetID(const FatJet &fatjet,
                                const FatJet::FatJetID &id) const {
  correction::Correction::Ref cset = nullptr;
  float out;
  switch (id) {
  case FatJet::FatJetID::TIGHT:
    cset = cset_jetid->at("AK8PUPPI_Tight");
    out = cset->evaluate({fabs(fatjet.Eta()), fatjet.chHEF(), fatjet.neHEF(),
                          fatjet.chEmEF(), fatjet.neEmEF(), fatjet.muEF(),
                          fatjet.chMultiplicity(), fatjet.neMultiplicity(),
                          fatjet.chMultiplicity() + fatjet.neMultiplicity()});
    return out > 0.5; // return is real
    break;
  case FatJet::FatJetID::TIGHTLEPVETO:
    cset = cset_jetid->at("AK8PUPPI_TightLeptonVeto");
    out = cset->evaluate({fabs(fatjet.Eta()), fatjet.chHEF(), fatjet.neHEF(),
                          fatjet.chEmEF(), fatjet.neEmEF(), fatjet.muEF(),
                          fatjet.chMultiplicity(), fatjet.neMultiplicity(),
                          fatjet.chMultiplicity() + fatjet.neMultiplicity()});
    return out > 0.5; // return is real
    break;
  case FatJet::FatJetID::NOCUT:
    // No cut, always return true
    return true;
    break;
  default:
    throw runtime_error("[MyCorrection::PassFatJetID] Invalid JetID type");
  }
}

// JERC
float MyCorrection::GetJER(const float eta, const float pt,
                           const float rho) const {
  correction::Correction::Ref cset = nullptr;
  string cset_string = JME_JER_GT.at("2024");
  if (std::string(GetEra().Data()) != "2024") {
    throw std::runtime_error(
        std::string("[MyCorrection::GetJER] Error: Expected era 2024, but got ") + 
        GetEra().Data()
    );
  }
  cset_string.replace(cset_string.find("######"), 6, "PtResolution");
  cset = cset_jerc->at(cset_string);
  return safeEvaluate(cset, "GetJER", {eta, pt, rho});
}

float MyCorrection::GetJERSF(const float eta, const float pt,
                             const variation syst,
                             const TString &source) const {
  // 1. Get the base string for the current era
  string base_cset_string = JME_JER_GT.at("2024");
  if (std::string(GetEra().Data()) != "2024") {
    throw std::runtime_error(
        std::string("[MyCorrection::GetJER] Error: Expected era 2024, but got ") + 
        GetEra().Data()
    );
  }

  // ==========================================
  // Run 2: Legacy combined tag schema
  // ==========================================
  if (Run == 2) {
    string cset_string = base_cset_string;
    cset_string.replace(cset_string.find("######"), 6, "ScaleFactor");
    correction::Correction::Ref cset = cset_jerc->at(cset_string);
    // Old format takes {eta, syst}
    return safeEvaluate(cset, "GetJERSF_Run2", {eta, getSystString_JME(syst)});
  }

  // ==========================================
  // Run 3: New separated tags schema
  // ==========================================
  if (Run == 3) {
    // 2. Evaluate Nominal ScaleFactor
    string sf_cset_string = base_cset_string;
    sf_cset_string.replace(sf_cset_string.find("######"), 6, "ScaleFactor");
    correction::Correction::Ref cset_sf = cset_jerc->at(sf_cset_string);
    
    // New SF format drops the systematic string argument and adds pT
    float sf_nom = safeEvaluate(cset_sf, "GetJERSF_nom", {eta, pt});

    if (syst == variation::nom) {
      return sf_nom;
    }

    // 3. Evaluate SFUncertainty for variations
    string unc_cset_string = base_cset_string;
    unc_cset_string.replace(unc_cset_string.find("######"), 6, "SFUncertainty");
    correction::Correction::Ref cset_unc = cset_jerc->at(unc_cset_string);
    
    // SFUncertainty tag also expects {eta, pt}
    float sf_unc = safeEvaluate(cset_unc, "GetJERSF_unc", {eta, pt});

    // 4. Apply the updated formula
    if (syst == variation::up) {
      return sf_nom * (1.0f + sf_unc);
    } else if (syst == variation::down) {
      return sf_nom * (1.0f - sf_unc);
    }
  }

  return 1.0f;
}

// JESC
float MyCorrection::GetJESSF(const float area, const float eta, const float pt,
                             const float phi, const float rho,
                             const unsigned int runNumber) const {
  correction::CompoundCorrection::Ref cset = nullptr;
  string cset_string = JME_JES_GT.at(GetEra().Data());
  cset_string.replace(cset_string.find("######"), 6, "L1L2L3Res");
  cset = cset_jerc->compound().at(cset_string);
  vector<correction::Variable::Type> args;
  float JESSF = 1.;
  if (GetEra() == "2023BPix" || GetEra() == "2024") {
    args = {area, eta, pt, rho, phi};
    if (IsDATA)
      args = {area, eta, pt, rho, phi, static_cast<float>(runNumber)};
  } else if (GetEra() == "2023") {
    args = {area, eta, pt, rho};
    if (IsDATA)
      args = {area, eta, pt, rho, static_cast<float>(runNumber)};
  } else {
    args = {area, eta, pt, rho};
  }
  return safeEvaluate(cset, "GetJERSF", args);
}

// for UParTAK4 JESC
float MyCorrection::GetJESSF(const float area, const float eta, const float pt,
                             const float phi, const float rho,
                             const unsigned int runNumber, const TString &source) const {
  correction::Correction::Ref cset = nullptr;
  string cset_string = JME_JES_GT.at(GetEra().Data());
  cset_string.replace(cset_string.find("######"), 6, source);
  try {
    cset = cset_jerc->at(cset_string);
  } catch (const std::out_of_range& e) {
    // 에러 발생 시 cset_string의 정확한 값을 에러 메시지와 함께 출력합니다.
    std::cerr << "\n[Error] Correction set을 찾을 수 없습니다!" << std::endl;
    std::cerr << "찾으려는 cset_string: [" << cset_string << "]" << std::endl;
    std::cerr << "에러 내용: " << e.what() << "\n" << std::endl;
    
    // 디버깅 정보를 출력한 후, 프로그램을 안전하게 중단하기 위해 에러를 다시 던집니다.
    throw; 
  }
  vector<correction::Variable::Type> args;
  float JESSF = 1.;
  if (GetEra() == "2023BPix" || GetEra() == "2024") {
    args = {eta, pt};
    if (IsDATA)
      args = {static_cast<float>(runNumber), eta, pt};
  } else if (GetEra() == "2023") {
    args = {area, eta, pt, rho};
    if (IsDATA)
      args = {area, eta, pt, rho, static_cast<float>(runNumber)};
  } else {
    args = {eta, pt};
  }
  return safeEvaluate(cset, "GetJERSF", args);
}

float MyCorrection::GetJESUncertainty(const float eta, const float pt,
                                      const TString &source) const {
  correction::Correction::Ref cset = nullptr;
  string cset_string = JME_JES_GT.at(GetEra().Data());
  cset_string.replace(cset_string.find("######"), 6, source);
  cset = cset_jerc->at(cset_string);
  return safeEvaluate(cset, "GetJESUncertainty", {eta, pt});
}

float MyCorrection::GetJESUncertaintySF(const float eta, const float pt,
                                        const variation syst,
                                        const TString &source) const {
  int int_syst = 0;
  if (syst == variation::up)
    int_syst = 1;
  else if (syst == variation::down)
    int_syst = -1;
  else
    int_syst = 0;

  float this_factor = 1.;
  this_factor += int_syst * GetJESUncertainty(eta, pt, source);
  return this_factor;
}

bool MyCorrection::IsJetVetoZone(const float eta, const float phi,
                                 TString mapCategory) const {
  correction::Correction::Ref cset = nullptr;
  string cset_string = JME_vetomap_keys.at(GetEra().Data());
  cset = cset_jetvetomap->at(cset_string);
  if (safeEvaluate(cset, "IsJetVetoZone", {mapCategory.Data(), eta, phi}) > 0)
    return true;
  return false;
}

void MyCorrection::METXYCorrection(Particle &Met, const int RunNumber,
                                   const int npvs,
                                   const XYCorrection_MetType MetType) {
  if (Run == 3) {
    // No supprot in Run3
    return;
  }
  correction::Correction::Ref cset_pt = nullptr;
  correction::Correction::Ref cset_phi = nullptr;
  switch (MetType) {
  case XYCorrection_MetType::Type1PFMET:
    if (IsDATA) {
      cset_pt = cset_met->at("pt_metphicorr_pfmet_data");
      cset_phi = cset_met->at("phi_metphicorr_pfmet_data");
    } else {
      cset_pt = cset_met->at("pt_metphicorr_pfmet_mc");
      cset_phi = cset_met->at("phi_metphicorr_pfmet_mc");
    }
    break;
  case XYCorrection_MetType::Type1PuppiMET:
    if (IsDATA) {
      cset_pt = cset_met->at("pt_metphicorr_puppimet_data");
      cset_phi = cset_met->at("phi_metphicorr_puppimet_data");
    } else {
      cset_pt = cset_met->at("pt_metphicorr_puppimet_mc");
      cset_phi = cset_met->at("phi_metphicorr_puppimet_mc");
    }
    break;
  }
  float this_pt =
      cset_pt->evaluate({Met.Pt(), Met.Phi(), static_cast<float>(npvs),
                         static_cast<float>(RunNumber)});
  float this_phi =
      cset_phi->evaluate({Met.Pt(), Met.Phi(), static_cast<float>(npvs),
                          static_cast<float>(RunNumber)});
  float this_eta = Met.Eta();
  float this_m = Met.M();
  Met.SetPtEtaPhiM(this_pt, this_eta, this_phi, this_m);
}

float MyCorrection::GetTopPtReweight(
    const TLorentzVector &LastCopyTop,
    const TLorentzVector &LastCopyAntiTop) const {

  const float t_pt = LastCopyTop.Pt();
  const float t_y = LastCopyTop.Rapidity();
  const float t_phi = LastCopyTop.Phi();
  const float t_mass = LastCopyTop.M();

  const float at_pt = LastCopyAntiTop.Pt();
  const float at_y = LastCopyAntiTop.Rapidity();
  const float at_phi = LastCopyAntiTop.Phi();
  const float at_mass = LastCopyAntiTop.M();

  const TLorentzVector tt = LastCopyTop + LastCopyAntiTop;
  const float tt_pt = tt.Pt();
  const float tt_y = tt.Rapidity();
  const float tt_phi = tt.Phi();
  const float tt_mass = tt.M();

  float feats[12] = {tt_pt, tt_y,   tt_phi, tt_mass, t_pt,   t_y,
                     t_phi, t_mass, at_pt,  at_y,    at_phi, at_mass};

  struct NormSpec {
    double mean;
    double std;
    bool use_log;
  };

  static const NormSpec nrm_tt[4] = {
      {3.6520673599656903, 1.0123402362573612, true},
      {0.0001718810581680775, 1.0362455506718102, false},
      {2.8943571877384285e-05, 1.8139038706413384, false},
      {6.21729978047307, 0.2771419580231537, true},
  };
  static const NormSpec nrm_top[4] = {
      {4.595855742518925, 0.7101176940989488, true},
      {0.00022746366634849002, 1.213207643109532, false},
      {-0.00028213870737636996, 1.8136544140703632, false},
      {171.93706459943778, 6.9652037622153, false},
  };
  static const NormSpec nrm_atop[4] = {
      {4.5986175957604045, 0.7103218938891299, true},
      {0.00011712322394057398, 1.2076422016031159, false},
      {0.0003628069129526392, 1.8139415747773364, false},
      {171.93691192651536, 6.9500586980501575, false},
  };

  auto normalize_block = [](float *x, const NormSpec *specs) {
    for (int j = 0; j < 4; ++j) {
      double v = static_cast<double>(x[j]);
      if (specs[j].use_log) {
        v = std::log(std::max(v, 1e-6));
      }
      v -= specs[j].mean;
      if (specs[j].std >= 1e-2) {
        v /= specs[j].std;
      }
      x[j] = static_cast<float>(v);
    }
  };

  normalize_block(&feats[0], nrm_tt);
  normalize_block(&feats[4], nrm_top);
  normalize_block(&feats[8], nrm_atop);

  std::array<float, 15> input_minnlo;
  // tt
  input_minnlo[0] = feats[0];
  input_minnlo[1] = feats[1];
  input_minnlo[2] = feats[2];
  input_minnlo[3] = feats[3];
  input_minnlo[4] = 0.0f;
  // top
  input_minnlo[5] = feats[4];
  input_minnlo[6] = feats[5];
  input_minnlo[7] = feats[6];
  input_minnlo[8] = feats[7];
  input_minnlo[9] = 0.6f;
  // antitop
  input_minnlo[10] = feats[8];
  input_minnlo[11] = feats[9];
  input_minnlo[12] = feats[10];
  input_minnlo[13] = feats[11];
  input_minnlo[14] = -0.6f;

  std::unordered_map<std::string, VariousArray> inputDataMap;
  static const std::unordered_map<std::string, IntArray> inputShapeMap_toppt = {
      {"input", IntArray{1, 3, 5}}};

  inputDataMap["input"] = FloatArray(input_minnlo.data(),
                                     input_minnlo.data() + input_minnlo.size());
  auto outputDataMap =
      MLHelper_TopPtReweight->Run_ONNX_Model(inputDataMap, inputShapeMap_toppt);
  return outputDataMap["activation_6"][1] / outputDataMap["activation_6"][0];
}

float MyCorrection::GethDampReweight(const TLorentzVector &FirstCopyTop,
                                     const TLorentzVector &FirstCopyAntiTop,
                                     const variation &syst) const {
  if (syst == variation::nom)
    return 1.f;

  const float t_pt = FirstCopyTop.Pt();
  const float t_y = FirstCopyTop.Rapidity();
  const float t_phi = FirstCopyTop.Phi();
  const float t_mass = FirstCopyTop.M();

  const float at_pt = FirstCopyAntiTop.Pt();
  const float at_y = FirstCopyAntiTop.Rapidity();
  const float at_phi = FirstCopyAntiTop.Phi();
  const float at_mass = FirstCopyAntiTop.M();

  const TLorentzVector tt = FirstCopyTop + FirstCopyAntiTop;
  const float tt_pt = tt.Pt();
  if (tt_pt > 1000.f)
    return 1.f;

  FloatArray input_hdamp(12);
  input_hdamp = {
      std::log10(t_pt),  t_y,  t_phi,  t_mass / 243.95f,  0.1, 1.379,
      std::log10(at_pt), at_y, at_phi, at_mass / 243.95f, 0.2, 1.379};

  std::unordered_map<std::string, VariousArray> inputDataMap;
  static const std::unordered_map<std::string, IntArray> inputShapeMap_hdamp = {
      {"input", IntArray{1, 2, 6}}};

  inputDataMap["input"] = input_hdamp;

  MLHelper *this_model = (syst == variation::up) ? MLHelper_hDampUp.get()
                                                 : MLHelper_hDampDown.get();

  auto outputDataMap =
      this_model->Run_ONNX_Model(inputDataMap, inputShapeMap_hdamp);
  return outputDataMap["activation_6"][0] / outputDataMap["activation_6"][1];
}

float MyCorrection::GetBFragReweight(
    const TLorentzVector &LastCopyTop, const TLorentzVector &LastCopyAntiTop,
    const TLorentzVector &LastCopyWPlus, const TLorentzVector &LastCopyWMinus,
    const TLorentzVector &FirstCopyBHadronFromTop,
    const TLorentzVector &FirstCopyBHadronFromAntiTop, const variation &syst)
    const {
  if (syst == variation::down)
    return 1.f;
  const float x_e_top =
      2 * FirstCopyBHadronFromTop * LastCopyTop / LastCopyTop.M2();
  const float x_e_antitop =
      2 * FirstCopyBHadronFromAntiTop * LastCopyAntiTop / LastCopyAntiTop.M2();
  const float w_top = LastCopyWPlus.M2() / LastCopyTop.M2();
  const float w_antitop = LastCopyWMinus.M2() / LastCopyAntiTop.M2();
  const float clip_value = 1.2f;
  const float x_b_top = std::min(x_e_top / (1 - w_top), clip_value);
  const float x_b_antitop = std::min(x_e_antitop / (1 - w_antitop), clip_value);
  FloatArray input_bfrag(4);
  MLHelper *this_model =
      (syst == variation::up) ? MLHelper_rBUp.get() : MLHelper_rBnom.get();
  input_bfrag = {x_b_top, 0.855, x_b_antitop, 0.855};

  std::unordered_map<std::string, VariousArray> inputDataMap;
  static const std::unordered_map<std::string, IntArray> inputShapeMap_bfrag = {
      {"input", IntArray{1, 2, 2}}};

  inputDataMap["input"] = input_bfrag;
  auto outputDataMap =
      this_model->Run_ONNX_Model(inputDataMap, inputShapeMap_bfrag);
  return outputDataMap["activation_6"][0] / outputDataMap["activation_6"][1];
}

std::array<std::size_t, 6>
MyCorrection::GetGenIdxofTopDecayProducts(const GenViewCollection &gens) const {
  constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  std::array<std::size_t, 6> out{npos, npos, npos, npos, npos, npos};

  const auto &storage = gens.storage();
  if (!storage)
    return out;

  const std::size_t n = storage->size();
  if (n == 0)
    return out;

  const auto &pt = storage->pt;
  const auto &eta = storage->eta;
  const auto &phi = storage->phi;
  const auto &mass = storage->mass;
  const auto &pdgId = storage->pdgId;
  const auto &statusFlags = storage->statusFlags;
  const auto &motherIdx = storage->motherIdx;

  constexpr unsigned short FIRST_COPY = 1u << 12;
  constexpr unsigned short LAST_COPY = 1u << 13;

  auto isBottomHadron = [](int pdg) {
    int ap = std::abs(pdg);
    if (ap < 500)
      return false;
    int hundreds = (ap / 100) % 10;
    int thousands = (ap / 1000) % 10;
    return (hundreds == 5) || (thousands == 5);
  };

  std::size_t topIdx = npos;
  std::size_t antiTopIdx = npos;
  std::size_t WTopIdx = npos;
  std::size_t WAntiTopIdx = npos;
  std::size_t bFromTopIdx = npos;
  std::size_t bFromAntiTopIdx = npos;

  std::vector<std::size_t> Bcands;
  Bcands.reserve(8);

  for (std::size_t i = 0; i < n; ++i) {
    const int id = pdgId[i];
    const unsigned short flg = statusFlags[i];

    if (flg & LAST_COPY) {
      if (id == 6 && topIdx == npos)
        topIdx = i;
      else if (id == -6 && antiTopIdx == npos)
        antiTopIdx = i;
      else if (id == 24 && WTopIdx == npos)
        WTopIdx = i;
      else if (id == -24 && WAntiTopIdx == npos)
        WAntiTopIdx = i;
      else if (id == 5 && bFromTopIdx == npos)
        bFromTopIdx = i;
      else if (id == -5 && bFromAntiTopIdx == npos)
        bFromAntiTopIdx = i;
    }

    if ((flg & FIRST_COPY) && isBottomHadron(id)) {
      Bcands.push_back(i);
    }
  }

  if (topIdx == npos || antiTopIdx == npos || WTopIdx == npos ||
      WAntiTopIdx == npos) {
    throw runtime_error("[MyCorrection::GetGenIdxofTopDecayProducts] Unable to "
                        "find top/anti-top or W bosons in the event.");
  }

  if (Bcands.size() < 2) {
    throw runtime_error("[MyCorrection::GetGenIdxofTopDecayProducts] Unable to "
                        "find sufficient b-hadron candidates in the event.");
  }

  auto matchBHad = [&](std::size_t bIdx, bool isTopLeg) -> std::size_t {
    if (bIdx == npos || Bcands.empty())
      return npos;

    const float eta_b = eta[bIdx];
    const float phi_b = phi[bIdx];

    double bestDR2 = 1e9;
    std::size_t best = npos;
    const double maxDR2 = 0.2 * 0.2; // DR < 0.2

    for (auto cand : Bcands) {
      const int pid = pdgId[cand];

      // PDG sign 패턴: top vs anti-top
      const int pid100 = (pid / 100) % 10;
      const int pid1000 = (pid / 1000) % 10;

      bool signOK = false;
      if (isTopLeg) {
        // ((pid/100)%10 == -5 || (pid/1000)%10 == 5)
        signOK = (pid100 == -5) || (pid1000 == 5);
      } else {
        // ((pid/100)%10 == 5 || (pid/1000)%10 == -5)
        signOK = (pid100 == 5) || (pid1000 == -5);
      }
      if (!signOK)
        continue;

      // ΔR<0.2 조건
      const double dEta = static_cast<double>(eta[cand]) - eta_b;
      double dPhi = static_cast<double>(phi[cand]) - phi_b;
      if (dPhi > M_PI)
        dPhi -= 2 * M_PI;
      else if (dPhi < -M_PI)
        dPhi += 2 * M_PI;

      const double dr2 = dEta * dEta + dPhi * dPhi;
      if (dr2 < bestDR2 && dr2 < maxDR2) {
        bestDR2 = dr2;
        best = cand;
      }
    }
    return best;
  };

  std::size_t BHadTopIdx = matchBHad(bFromTopIdx, true);
  std::size_t BHadAntiIdx = matchBHad(bFromAntiTopIdx, false);

  out[0] = topIdx;
  out[1] = WTopIdx;
  out[2] = BHadTopIdx;
  out[3] = antiTopIdx;
  out[4] = WAntiTopIdx;
  out[5] = BHadAntiIdx;

  return out;
}
