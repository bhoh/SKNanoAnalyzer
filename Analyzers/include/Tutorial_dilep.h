#ifndef Tutorial_dilep_h
#define Tutorial_dilep_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class Tutorial_dilep : public AnalyzerCore {
public:
  void initializeAnalyzer();
  void executeEvent();
  void executeEventFromParameter();

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
  GenViewCollection AllGenViews;

  Event ev;
  unique_ptr<SystematicHelper> systHelper;

  Tutorial_dilep();
  ~Tutorial_dilep();
};

#endif
