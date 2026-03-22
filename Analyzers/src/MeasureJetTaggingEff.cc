#include "MeasureJetTaggingEff.h"

MeasureJetTaggingEff::MeasureJetTaggingEff() {}
MeasureJetTaggingEff::~MeasureJetTaggingEff() {}

void MeasureJetTaggingEff::initializeAnalyzer() {

  MuonIDs = { Muon::MuonID::POG_TIGHT };
  MuonIDSFKeys = { "NUM_TightID_DEN_TrackerMuons" };

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
  WP_Loose = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Loose);
  WP_Medium = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Medium);
  WP_Tight = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::Tight);
  WP_VeryTight = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::VeryTight);
  WP_SuperTight = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerWP::SuperTight);

  // init SystematicHelper
  string SKNANO_HOME = getenv("SKNANO_HOME");
  if (IsDATA) {
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
  } else {
    //systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", MCSample, DataEra);
  }
}

void MeasureJetTaggingEff::executeEvent() {

  AllMuonViews = GetAllMuonViews();
  AllElectronViews = GetAllElectronViews();
  AllJetViews = GetAllJetViews();

  ev = GetEvent();

  // Check this for Run3
  weight_Prefire = 1.; 

  // Systematic sources from YAML

  for (const auto &syst_dummy : *systHelper) {
    executeEventFromParameter();
  }
}

void MeasureJetTaggingEff::executeEventFromParameter() {

  const TString this_syst = systHelper->getCurrentSysName();

  Muon::MuonID this_muon_id = MuonIDs[0];
  TString this_muon_id_sf_key = MuonIDSFKeys[0];

  //==== MET Filter & Trigger
  if (!PassMetFilter(AllJetViews, ev)) return;

  //==== Lepton Selection
  std::vector<size_t> SelectedMuonIndices = SelectMuonIndices(AllMuonViews, this_muon_id, 15., 2.4);
  std::vector<size_t> SelectedElectronIndices = SelectElectronIndices(AllElectronViews, Electron::ElectronID::POG_LOOSE, 15., 2.5);

  if (SelectedMuonIndices.size() + SelectedElectronIndices.size() != 1) return;

  RVec<Muon> muons = MaterializeMuons(AllMuonViews, SelectedMuonIndices);
  RVec<Electron> electrons = MaterializeElectrons(AllElectronViews, SelectedElectronIndices);

  //==== Jet Selection (Updated pt 20., eta 2.5)
  std::vector<size_t> SelectedJetIndices = SelectJetIndices(AllJetViews, Jet::JetID::TIGHT, 20., 2.5);
  RVec<Jet> jets = MaterializeJets(AllJetViews, SelectedJetIndices);
  
  // Veto Lepton in Jets
  jets = JetsVetoLeptonInside(jets, electrons, muons, 0.3);

  //==== Sorting
  sort(muons.begin(), muons.end(), PtComparing);
  sort(jets.begin(), jets.end(), PtComparing);

  //==== B-Tagging Efficiency Setup
  // Define binning for abseta and pt using RVec
  const RVec<float> abseta_bins = {0., 0.8, 1.6, 2.5};
  const RVec<float> pt_bins = {20., 30., 50., 70., 100., 140., 200., 300., 600., 1000.};

  // Event weight
  float weight = 1.0; 

  // Store WP labels and their corresponding cut values in a vector for easy iteration
  std::vector<std::pair<TString, float>> wp_cuts = {
      {"L", WP_Loose},
      {"M", WP_Medium},
      {"T", WP_Tight},
      {"XT", WP_VeryTight},
      {"XXT", WP_SuperTight}
  };

  //==== Loop over all selected jets to fill efficiency histograms
  for (unsigned int ij = 0; ij < jets.size(); ij++) {
    float pt = jets.at(ij).Pt();
    float abseta = std::abs(jets.at(ij).Eta());
    int flavor = jets.at(ij).hadronFlavour();

    // Handle flavor: group everything other than b(5) and c(4) into udsg(0)
    if (flavor != 5 && flavor != 4) {
      flavor = 0;
    }

    // 1. Fill denominator (All jets regardless of WP. Filled ONLY ONCE per jet)
    TString den_key = TString::Format("tagging#b##era#%s##flavor#%d##systematic#central##den", 
                                      DataEra.Data(), flavor);
    FillHist(den_key.Data(), abseta, pt, weight, abseta_bins, pt_bins);

    // Get Tagger Score for the current jet
    double this_discr = jets.at(ij).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);

    // 2. Loop over the WPs to fill the appropriate numerator histograms
    for (const auto& wp : wp_cuts) {
      TString wp_label = wp.first;
      float wp_cut = wp.second;

      if (this_discr > wp_cut) {
        // Dynamically insert the WP label ("L", "M", etc.) into the key
        TString num_key = TString::Format("tagging#b##era#%s##tagger#UParTAK4##working_point#%s##flavor#%d##systematic#central##num", 
                                          DataEra.Data(), wp_label.Data(), flavor);
        FillHist(num_key.Data(), abseta, pt, weight, abseta_bins, pt_bins);
      }
    }
  }

}