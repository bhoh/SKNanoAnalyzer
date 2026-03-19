#ifndef ExampleRun_base_h
#define ExampleRun_base_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class ExampleRun_base : public AnalyzerCore {
public:
    ExampleRun_base();
    ~ExampleRun_base();

    void initializeAnalyzer();
    virtual void executeEvent();

    bool RunSyst;
    bool RunNewPDF;
    bool RunXSecSyst;

    TString IsoMuTriggerName;
    float TriggerSafePtCut;

    RVec<TString> MuonIDSFKeys;
    RVec<Muon::MuonID> MuonIDs;
    
    // Variables used in event execution
    RVec<Muon> AllMuons;
    RVec<Jet> AllJets;
    float weight_Prefire;

    unique_ptr<SystematicHelper> systHelper;
};

#endif