#ifndef AtobbMLTree_h
#define AtobbMLTree_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class AtobbMLTree : public AnalyzerCore {
public:
  void initializeAnalyzer();
  void executeEvent();

  AtobbMLTree();
  ~AtobbMLTree();

private:
  TString IsoMuTriggerName;
  float TriggerSafePtCut;

  RVec<Muon::MuonID> MuonIDs;
  RVec<TString> MuonIDISOSFKeys;

  MuonViewCollection AllMuonViews;
  ElectronViewCollection AllElectronViews;
  JetViewCollection AllJetViews;
  GenViewCollection AllGenViews;

  Event ev;
  unique_ptr<SystematicHelper> systHelper;

  const double const_top_mass;
  const double const_top_width;
  const double const_w_mass;
  const double const_w_width;

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

  void FillJetBranches(const TString& prefix, const Jet& jet, float btag_wp_cut);
  void FillDummyJetBranches(const TString& prefix);

  void FillTreeBranches(
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
  );

  void EvalChi2(ttCombinatoric& tt_combinatoric);
  double Chi2Function(double had_top_mass, double had_W_mass,
                      double lep_top_mass, double lep_W_mass);
};

#endif
