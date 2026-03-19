#include "ExampleRun_base.h"

ExampleRun_base::ExampleRun_base() {}
ExampleRun_base::~ExampleRun_base() {}

void ExampleRun_base::initializeAnalyzer() {
    //==== Example 1
    //==== Dimuon Z-peak events with two muon IDs, with systematics
    
    RunSyst = HasFlag("RunSyst");
    cout << "[ExampleRun_base::initializeAnalyzer] RunSyst = " << RunSyst << endl;

    MuonIDs = {Muon::MuonID::POG_MEDIUM, Muon::MuonID::POG_TIGHT};
    MuonIDSFKeys = {"NUM_MediumID_DEN_TrackerMuons", "NUM_TightID_DEN_TrackerMuons"};

    if (DataEra == "2016preVFP" || DataEra == "2016postVFP" || DataEra == "2018") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2017") {
        IsoMuTriggerName = "HLT_IsoMu27";
        TriggerSafePtCut = 29.;
    } else if (DataEra == "2022" || DataEra == "2022EE") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2023") {
        IsoMuTriggerName = "";
        TriggerSafePtCut = 26.;
    } else {
        cerr << "[ExampleRun_base::initializeAnalyzer] DataEra is not set properly" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "[ExampleRun_base::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[ExampleRun_base::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample ,IsDATA);

    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    }

    //==== Example 2
    //==== Using new PDF
    RunNewPDF = HasFlag("RunNewPDF");
    cout << "[ExampleRun_base::initializeAnalyzer] RunNewPDF = " << RunNewPDF << endl;
    if (RunNewPDF && !IsDATA) {
        LHAPDFHandler LHAPDFHandler_Prod;
        LHAPDFHandler_Prod.CentralPDFName = "NNPDF31_nnlo_hessian_pdfas";
        LHAPDFHandler_Prod.init();

        LHAPDFHandler LHAPDFHandler_New;
        LHAPDFHandler_New.CentralPDFName = "NNPDF31_nlo_hessian_pdfas";
        LHAPDFHandler_New.ErrorSetMember_Start = 1;
        LHAPDFHandler_New.ErrorSetMember_End = 100;
        LHAPDFHandler_New.AlphaSMember_Down = 101;
        LHAPDFHandler_New.AlphaSMember_Up = 102;
        LHAPDFHandler_New.init();
        
        pdfReweight->SetProdPDF( LHAPDFHandler_Prod.PDFCentral );
        pdfReweight->SetNewPDF( LHAPDFHandler_New.PDFCentral );
        pdfReweight->SetNewPDFErrorSet( LHAPDFHandler_New.PDFErrorSet );
        pdfReweight->SetNewPDFAlphaS( LHAPDFHandler_New.PDFAlphaS_Down, LHAPDFHandler_New.PDFAlphaS_Up );
    }

    //==== Example 3
    //==== How to estimate xsec errors (PDF & Scale)
    RunXSecSyst = false;
    cout << "[ExampleRun_base::initializeAnalyzer] RunXSecSyst = " << RunXSecSyst << endl;
}

void ExampleRun_base::executeEvent() {
    // Virtual function to be overridden in Python
    return;
}