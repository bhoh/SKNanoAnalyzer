#ifndef Tutorial_reco_tt_h
#define Tutorial_reco_tt_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class Tutorial_reco_tt : public AnalyzerCore {

  // particle mass in MC setting
  // will be initialized in constructor
  const double const_top_mass;
  const double const_top_width;
  const double const_w_mass;
  const double const_w_width;

public:

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter(const TString &override_syst = "",
                                 MyCorrection::variation met_variation =
                                     MyCorrection::variation::nom);

  TString IsoMuTriggerName;
  float TriggerSafePtCut;
  TString EleTriggerName;
  float EleTriggerSafePtCut;

  RVec<Muon::MuonID> MuonIDs;
  RVec<TString> MuonIDISOSFKeys;
  RVec<Electron::ElectronID> ElectronIDs;
  RVec<TString> ElectronIDSFKeys;
  RVec<TString> ElectronTriggerSFKeys;

  MuonViewCollection AllMuonViews;
  ElectronViewCollection AllElectronViews;
  JetViewCollection AllJetViews;
  GenViewCollection AllGenViews;
  GenJetViewCollection AllGenJetViews;

  Event ev;
  float weight_Prefire;

  unique_ptr<SystematicHelper> systHelper;

  array<std::size_t, 4> GetTopAndAntiTopIndices(const GenViewCollection &gens);

  Tutorial_reco_tt();
  ~Tutorial_reco_tt();

private:

  struct ttCombinatoric {
    Lepton* lepton;
    RVec<Jet>* jets; 
    Particle* met;

    unsigned int had_W_jet_idx_1;
    unsigned int had_W_jet_idx_2;
    unsigned int had_top_b_jet_idx;
    unsigned int lep_top_b_jet_idx;

    double had_W_mass;
    double had_top_mass;
    double had_W_b_tagger_score_1; 
    double had_W_b_tagger_score_2;

    std::vector<double> neu_pz;
    std::vector<double> lep_top_mass;
    std::vector<double> lep_W_mass;
    std::vector<double> chi2;

    double best_neu_pz;
    double best_lep_top_mass;
    double best_lep_W_mass;
    double best_chi2;

    void EvalHadronicPart();
    void EvalLeptonicPart();
  };

  void EvalChi2(ttCombinatoric& tt_combinatoric);
  double Chi2Function(double had_top_mass, double had_W_mass, double lep_top_mass, double lep_W_mass);
};

#endif
