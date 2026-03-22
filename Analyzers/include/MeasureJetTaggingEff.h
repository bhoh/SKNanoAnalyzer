#ifndef MeasureJetTaggingEff_h
#define MeasureJetTaggingEff_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class MeasureJetTaggingEff : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter();

  TString IsoMuTriggerName;
  float TriggerSafePtCut;

  RVec<Muon::MuonID> MuonIDs;
  RVec<TString> MuonIDSFKeys;

  MuonViewCollection AllMuonViews;
  ElectronViewCollection AllElectronViews;
  JetViewCollection AllJetViews;

  Event ev;
  float weight_Prefire;

  float WP_Loose;
  float WP_Medium;
  float WP_Tight;
  float WP_VeryTight;
  float WP_SuperTight;

  unique_ptr<SystematicHelper> systHelper;

  MeasureJetTaggingEff();
  ~MeasureJetTaggingEff();

};

#endif