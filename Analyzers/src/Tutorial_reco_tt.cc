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

  MuonIDs = { Muon::MuonID::POG_TIGHT };
  MuonIDISOSFKeys = { "NUM_TightID_DEN_TrackerMuons", "NUM_TightPFIso_DEN_TightID" };

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
    systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", MCSample, DataEra);
  }
}

void Tutorial_reco_tt::executeEvent() {

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

void Tutorial_reco_tt::executeEventFromParameter() {

  const TString this_syst = systHelper->getCurrentSysName();

  Muon::MuonID this_muon_id = MuonIDs[0];
  TString this_muon_id_sf_key = MuonIDISOSFKeys[0];
  TString this_muon_iso_sf_key = MuonIDISOSFKeys[1];

  FillHist(this_syst + "/cutflow_" + this_syst, 0.5, 1., 6, 0., 6.);

  //==== MET Filter & Trigger
  if (!PassMetFilter(AllJetViews, ev)) return;
  if (!(ev.PassTrigger(IsoMuTriggerName))) return;

  Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI); 

  //==== Lepton Selection
  std::vector<size_t> SelectedMuonIndices = SelectMuonIndices(AllMuonViews, this_muon_id, 15., 2.4);
  std::vector<size_t> SelectedElectronIndices = SelectElectronIndices(AllElectronViews, Electron::ElectronID::POG_LOOSE, 15., 2.5);

  if (SelectedMuonIndices.size() + SelectedElectronIndices.size() != 1) return;

  RVec<Muon> muons = MaterializeMuons(AllMuonViews, SelectedMuonIndices);
  RVec<Electron> electrons = MaterializeElectrons(AllElectronViews, SelectedElectronIndices);

  //==== Jet Selection
  std::vector<size_t> SelectedJetIndices = SelectJetIndices(AllJetViews, Jet::JetID::TIGHT, 30., 2.4);
  RVec<Jet> jets = MaterializeJets(AllJetViews, SelectedJetIndices);
  
  // Veto Lepton in Jets
  jets = JetsVetoLeptonInside(jets, electrons, muons, 0.3);

  //==== Sorting
  sort(muons.begin(), muons.end(), PtComparing);
  sort(jets.begin(), jets.end(), PtComparing);

  //==== Event selections
  if (muons.size() != 1) return;
  if (muons.at(0).Pt() <= TriggerSafePtCut) return;
  if (jets.size() < 4) return;
  if (METv.Pt() <= 20) return;
  FillHist(this_syst + "/cutflow_" + this_syst, 1.5, 1., 6, 0., 6.);

  if (!PassJetVetoMap(AllJetViews, AllMuonViews, "jetvetomap_fpix")) return;
  FillHist(this_syst + "/cutflow_" + this_syst, 2.5, 1., 6, 0., 6.);

  //==== B-Tagging (DeepJet Medium WP example)
  int NBJets = 0;
  float btag_wp_cut = myCorr->GetBTaggingWP(); 
  std::vector<bool> btag_vector;

  for (unsigned int ij = 0; ij < jets.size(); ij++) {
    double this_discr = jets.at(ij).GetTaggerResult(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B);
    if (this_discr > btag_wp_cut) {
      NBJets++;
      btag_vector.push_back(true);
    } else {
      btag_vector.push_back(false);
    }
  }

  if (NBJets != 2) return;
  FillHist(this_syst + "/cutflow_" + this_syst, 3.5, 1., 6, 0., 6.);

  //==== Take leading four jets in pT
  std::vector<unsigned int> top_b_jet_candidates;
  std::vector<unsigned int> had_W_candidates;

  for(unsigned int ij(0); ij < 4; ij++){
    if(btag_vector.at(ij)){
      if(top_b_jet_candidates.size() < 2) top_b_jet_candidates.push_back(ij);
    } else {
      if(had_W_candidates.size() < 2) had_W_candidates.push_back(ij);
    }
  }

  if(had_W_candidates.size() < 2 || top_b_jet_candidates.size() < 2) return;
  FillHist(this_syst + "/cutflow_" + this_syst, 4.5, 1., 6, 0., 6.);

  //==== Combinatorics
  Tutorial_reco_tt::ttCombinatoric tt_combinatoric_1, tt_combinatoric_2;

  tt_combinatoric_1.lepton            = &(muons.at(0));
  tt_combinatoric_1.jets              = &jets;
  tt_combinatoric_1.met               = &METv;
  tt_combinatoric_1.had_W_jet_idx_1   = had_W_candidates.at(0);
  tt_combinatoric_1.had_W_jet_idx_2   = had_W_candidates.at(1);
  tt_combinatoric_1.had_top_b_jet_idx = top_b_jet_candidates.at(0);
  tt_combinatoric_1.lep_top_b_jet_idx = top_b_jet_candidates.at(1);

  tt_combinatoric_2.lepton            = &(muons.at(0));
  tt_combinatoric_2.jets              = &jets;
  tt_combinatoric_2.met               = &METv;
  tt_combinatoric_2.had_W_jet_idx_1   = had_W_candidates.at(0);
  tt_combinatoric_2.had_W_jet_idx_2   = had_W_candidates.at(1);
  tt_combinatoric_2.had_top_b_jet_idx = top_b_jet_candidates.at(1);
  tt_combinatoric_2.lep_top_b_jet_idx = top_b_jet_candidates.at(0);

  this->EvalChi2(tt_combinatoric_1);
  this->EvalChi2(tt_combinatoric_2);

  Tutorial_reco_tt::ttCombinatoric* best_combinatoric;
  if(tt_combinatoric_1.best_chi2 < tt_combinatoric_2.best_chi2){
    best_combinatoric = &tt_combinatoric_1;
  } else {
    best_combinatoric = &tt_combinatoric_2;
  }

  //==== Event Weight
  float weight = 1.;
  if (!IsDATA) {
    weight *= MCweight();
    weight *= ev.GetTriggerLumi("Full");
    // muon official trigger SF is not available yet.
    float muon_id_sf = myCorr->GetMuonIDSF(this_muon_id_sf_key, muons, MyCorrection::variation::nom);
    weight *= muon_id_sf;
    float muon_iso_sf = myCorr->GetMuonIDSF(this_muon_iso_sf_key, muons, MyCorrection::variation::nom);
    weight *= muon_iso_sf;
    float pu_weight = myCorr->GetPUWeight(ev.nTrueInt(), MyCorrection::variation::nom);
    weight *= pu_weight;
    float btag_sf = myCorr->GetBTaggingSF(jets, 
            JetTagging::JetFlavTagger::ParT, 
            JetTagging::JetFlavTaggerWP::Medium,
            JetTagging::JetTaggingSFMethod::comb,
            MyCorrection::variation::nom
        );
    weight *= btag_sf;

    FillHist(this_syst + "/noChi2Cut/had_W_mass_" + this_syst, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
    FillHist(this_syst + "/noChi2Cut/had_top_mass_" + this_syst, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
    FillHist(this_syst + "/noChi2Cut/chi2_" + this_syst, best_combinatoric->best_chi2, weight, 1000, 0., 100000.);
  } else {
    FillHist(this_syst + "/noChi2Cut/had_W_mass_" + this_syst, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
    FillHist(this_syst + "/noChi2Cut/had_top_mass_" + this_syst, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
    FillHist(this_syst + "/noChi2Cut/chi2_" + this_syst, best_combinatoric->best_chi2, weight, 1000, 0., 100000.);
  }

  if(best_combinatoric->best_chi2 >= 2e3) return;
  FillHist(this_syst + "/cutflow_" + this_syst, 5.5, 1., 6, 0., 6.);
  
  if (!IsDATA) {
      FillHist(this_syst + "/Chi2Cut/had_W_mass_" + this_syst, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
      FillHist(this_syst + "/Chi2Cut/had_top_mass_" + this_syst, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
  } else {
      FillHist(this_syst + "/Chi2Cut/had_W_mass_" + this_syst, best_combinatoric->had_W_mass, weight, 40, 0., 200.);
      FillHist(this_syst + "/Chi2Cut/had_top_mass_" + this_syst, best_combinatoric->had_top_mass, weight, 80, 0., 400.);
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