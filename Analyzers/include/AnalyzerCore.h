#ifndef AnalyzerCore_h
#define AnalyzerCore_h

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/hash/hash.h"
#include <nlohmann/json.hpp>

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TObjString.h"
#include "TMath.h"
#include "Compression.h"  

#include "SKNanoLoader.h"
#include "Event.h"
#include "Particle.h"
#include "Lepton.h"
#include "Gen.h"
#include "GenView.h"
#include "LHE.h"
#include "Muon.h"
#include "MuonView.h"
#include "Electron.h"
#include "ElectronView.h"
#include "Jet.h"
#include "JetView.h"
#include "FatJet.h"
#include "Tau.h"
#include "Photon.h"
#include "GenJet.h"
#include "GenJetView.h"
#include "GenDressedLepton.h"
#include "GenIsolatedPhoton.h"
#include "GenVisTau.h"
#include "TrigObj.h"
#include "TrigObjView.h"

#include "LHAPDFHandler.h"
#include "PDFReweight.h"
#include "MyCorrection.h"
#include "JetTaggingParameter.h"
#include "PhysicalConstants.h"

struct TransparentStringHash {
    using is_transparent = void;
    std::size_t operator()(std::string_view s) const noexcept {
        return absl::Hash<std::string_view>{}(s);
    }
    std::size_t operator()(const std::string& s) const noexcept {
        return absl::Hash<std::string_view>{}(s);
    }
};

struct TransparentStringEq {
    using is_transparent = void;
    bool operator()(std::string_view a, std::string_view b) const noexcept {
        return a == b;
    }
    bool operator()(const std::string& a, const std::string& b) const noexcept {
        return a == b;
    }
    bool operator()(const std::string& a, std::string_view b) const noexcept {
        return a == b;
    }
    bool operator()(std::string_view a, const std::string& b) const noexcept {
        return a == b;
    }
}; 

class IDContainer {
public:
    IDContainer() {}
    IDContainer(const TString &tight, const TString &loose):
        j_tight(tight), j_loose(loose) {}

    TString GetID(const TString &wp) const {
        if (wp == "tight") return j_tight;
        else if (wp == "loose") return j_loose;
        else throw runtime_error("Invalid WP: " + wp);
    }

private:
    TString j_tight, j_loose;
};

class AnalyzerCore: public SKNanoLoader {
public:
    AnalyzerCore();
    ~AnalyzerCore();

    virtual void initializeAnalyzer() {};
    virtual void executeEvent() {};

    inline bool HasFlag(const TString &flag) { return std::find(Userflags.begin(), Userflags.end(), flag) != Userflags.end(); }

    inline static bool PtComparing(const Particle& p1, const Particle& p2) { return p1.Pt() > p2.Pt();}
    inline static bool PtComparingPtr(const Particle* p1, const Particle* p2) { return p1->Pt() > p2->Pt();}


    //MetFilter
    bool PassNoiseFilter(const RVec<Jet> &AllJets, const Event &ev, Event::MET_Type met_type = Event::MET_Type::PUPPI);
    bool PassNoiseFilter(const JetViewCollection &AllJets, const Event &ev, Event::MET_Type met_type = Event::MET_Type::PUPPI);
    bool PassMetFilter(const RVec<Jet> &AllJets, const Event &ev, Event::MET_Type met_type = Event::MET_Type::PUPPI);
    bool PassMetFilter(const JetViewCollection &AllJets, const Event &ev, Event::MET_Type met_type = Event::MET_Type::PUPPI);
    // PDF reweight
    PDFReweight *pdfReweight;
    float GetPDFWeight(LHAPDF::PDF *pdf_);
    float GetPDFReweight();
    float GetPDFReweight(int member);    
    // Correction
    MyCorrection *myCorr;
    //unique_ptr<CorrectionSet> csetMuon;
    //unique_ptr<CorrectionSet> csetElectron;;

    // MC weights
    float MCweight(bool usesign = true, bool norm_1invpb = true) const;

    struct ModellingPatch {
        std::vector<float> patch_ScaleVariation;
        std::vector<float> patch_PSVariation;
        float patch_hdamp_up = 0.f;
        float patch_hdamp_down = 0.f;
        float patch_minnlo = 0.f;
    };

    void load_modelling_json(const TString &filename);

    // Get objects
    Event GetEvent();
    MuonViewCollection GetAllMuonViews();
    GenViewCollection GetAllGenViews();
    JetViewCollection GetAllJetViews();
    void SmearJetViews(const JetViewCollection &jets, const float rho);
    GenJetViewCollection GetAllGenJetViews();
    RVec<Muon> GetAllMuons();
    std::vector<std::size_t> SelectMuonIndices(const MuonViewCollection &muons, const std::vector<std::size_t> &seed_indices, const Muon::MuonID ID, const float ptmin, const float fetamax) const;
    std::vector<std::size_t> SelectMuonIndices(const MuonViewCollection &muons, const Muon::MuonID ID, const float ptmin, const float fetamax) const;
    RVec<Muon> MaterializeMuons(const MuonViewCollection &muons) const;
    RVec<Muon> MaterializeMuons(const MuonViewCollection &muons, const std::vector<std::size_t> &indices) const;
    ElectronViewCollection GetAllElectronViews();
    RVec<Electron> GetAllElectrons();
    std::vector<std::size_t> SelectElectronIndices(const ElectronViewCollection &electrons, const std::vector<size_t> &seed_indices, const Electron::ElectronID ID, const float ptmin, const float fetamax, bool vetoHEM = false) const;
    std::vector<std::size_t> SelectElectronIndices(const ElectronViewCollection &electrons, const Electron::ElectronID ID, const float ptmin, const float fetamax, bool vetoHEM = false) const;
    RVec<Electron> MaterializeElectrons(const ElectronViewCollection &electrons) const;
    RVec<Electron> MaterializeElectrons(const ElectronViewCollection &electrons, const std::vector<std::size_t> &indices) const;
    Jet MaterializeJet(const JetViewCollection &jets, std::size_t index, const MyCorrection::variation &JESVariation = MyCorrection::variation::nom, const MyCorrection::variation &JERVariation = MyCorrection::variation::nom) const;
    JetViewCollection CloneJetViews(const JetViewCollection &jets) const;
    RVec<Jet> MaterializeJets(const JetViewCollection &jets, const MyCorrection::variation &JESVariation = MyCorrection::variation::nom, const MyCorrection::variation &JERVariation = MyCorrection::variation::nom) const;
    RVec<Jet> MaterializeJets(const JetViewCollection &jets, const std::vector<std::size_t> &indices, const MyCorrection::variation &JESVariation = MyCorrection::variation::nom, const MyCorrection::variation &JERVariation = MyCorrection::variation::nom) const;
    RVec<GenJet> MaterializeGenJets(const GenJetViewCollection &genjets) const;

    RVec<Muon> SelectMuons(const MuonViewCollection &muons, const Muon::MuonID ID, const float ptmin, const float absetamax) const {
        auto indices = SelectMuonIndices(muons, ID, ptmin, absetamax);
        return MaterializeMuons(muons, indices);
    }

    RVec<Electron> SelectElectrons(const ElectronViewCollection &electrons, const Electron::ElectronID ID, const float ptmin, const float absetamax, bool vetoHEM = false) const {
        auto indices = SelectElectronIndices(electrons, ID, ptmin, absetamax, vetoHEM);
        return MaterializeElectrons(electrons, indices);
    }
    RVec<Jet> GetAllJets();
    RVec<Gen> GetAllGens();
    RVec<LHE> GetAllLHEs();
    RVec<Jet> GetJets(const TString id, const float ptmin, const float fetamax);
    RVec<Tau> GetAllTaus();
    RVec<FatJet> GetAllFatJets();
    RVec<GenJet> GetAllGenJets();
    RVec<GenDressedLepton> GetAllGenDressedLeptons();
    RVec<GenIsolatedPhoton> GetAllGenIsolatedPhotons();
    RVec<GenVisTau> GetAllGenVisTaus();
    static void MuonEnsureThunk(void *ctx, Muon &muon, Muon::Property property);
    static void ElectronEnsureThunk(void *ctx, Electron &electron, Electron::Property property);
    static void JetEnsureThunk(void *ctx, Jet &jet, Jet::Property property);

    void EnsureMuonProperty(Muon &muon, Muon::Property property) const;
    void EnsureElectronProperty(Electron &electron, Electron::Property property) const;
    void EnsureJetProperty(Jet &jet, Jet::Property property) const;

    RVec<Photon> GetAllPhotons();
    RVec<Photon> GetPhotons(TString id, double ptmin, double fetamax);
    TrigObjViewCollection GetAllTrigObjViews();
    std::vector<std::size_t> SelectTrigObjIndices(const TrigObjViewCollection &trigobjs, const std::vector<std::size_t> &seed_indices, const int id, const float ptmin, const float fetamax) const;
    std::vector<std::size_t> SelectTrigObjIndices(const TrigObjViewCollection &trigobjs, const int id, const float ptmin, const float fetamax) const;
    RVec<TrigObj> MaterializeTrigObjs(const TrigObjViewCollection &trigobjs) const;
    RVec<TrigObj> MaterializeTrigObjs(const TrigObjViewCollection &trigobjs, const std::vector<std::size_t> &indices) const;
    RVec<TrigObj> GetAllTrigObjs();

    // Select objects
    RVec<Muon> SelectMuons(const RVec<Muon> &muons, TString ID, const float ptmin, const float absetamax) const;
    RVec<Muon> SelectMuons(const RVec<Muon> &muons, Muon::MuonID ID, const float ptmin, const float absetamax) const;
    RVec<Jet> SelectJets(const RVec<Jet> &jets, const Jet::JetID, const float ptmin, const float fetamax) const;
    std::vector<std::size_t> SelectJetIndices(const JetViewCollection &jets, const std::vector<size_t> &seed_indices, const Jet::JetID, const float ptmin, const float fetamax, const MyCorrection::variation &JESVariation = MyCorrection::variation::nom, const MyCorrection::variation &JERVariation = MyCorrection::variation::nom) const;
    std::vector<std::size_t> SelectJetIndices(const JetViewCollection &jets, const Jet::JetID, const float ptmin, const float fetamax, const MyCorrection::variation &JESVariation = MyCorrection::variation::nom, const MyCorrection::variation &JERVariation = MyCorrection::variation::nom) const;
    RVec<Jet> SelectJets(const JetViewCollection &jets, const Jet::JetID, const float ptmin, const float fetamax) const;
    RVec<FatJet> SelectFatJets(const RVec<FatJet> &fatjets, const FatJet::FatJetID ID, const float ptmin, const float fetamax) const;
    RVec<Jet> JetsVetoLeptonInside(const RVec<Jet> &jets, const ElectronViewCollection &electrons, const MuonViewCollection &muons, const float dR = 0.3) const;
    std::vector<std::size_t> JetsVetoLeptonInside(const JetViewCollection& jets,
                                   const std::vector<std::size_t>& jet_indices,
                                   const ElectronViewCollection& electrons,
                                   const std::vector<std::size_t>& electron_indices,
                                   const MuonViewCollection& muons,
                                   const std::vector<std::size_t>& muon_indices,
                                   const float dR = 0.3) const;
    RVec<Jet> JetsVetoLeptonInside(const RVec<Jet> &jets, const RVec<Electron> &electrons, const RVec<Muon> &muons, const float dR = 0.3) const;
    RVec<Electron> SelectElectrons(const RVec<Electron> &electrons, const TString id, const float ptmin, const float absetamax, bool vetoHEM = false) const;
    RVec<Electron> SelectElectrons(const RVec<Electron> &electrons, const Electron::ElectronID ID, const float ptmin, const float absetamax, bool vetoHEM = false) const;
    RVec<Tau> SelectTaus(const RVec<Tau> &taus, const TString ID, const float ptmin, const float absetamax) const;
    // Functions
    float GetScaleVariation(const MyCorrection::variation &muF_syst, const MyCorrection::variation &muR_syst);
    float GetPSWeight(const MyCorrection::variation &ISR_syst, const MyCorrection::variation &FSR_syst);
    inline float GetBTaggingWP(const JetTagging::JetFlavTagger &tagger, const JetTagging::JetFlavTaggerWP &wp) { return myCorr->GetBTaggingWP(tagger, wp); }
    inline pair<float, float> GetCTaggingWP(const JetTagging::JetFlavTagger &tagger, const JetTagging::JetFlavTaggerWP &wp) { return myCorr->GetCTaggingWP(tagger, wp); }
    inline float GetBTaggingWP(){ return myCorr->GetBTaggingWP(); }
    inline pair<float, float> GetCTaggingWP(){ return myCorr->GetCTaggingWP(); }
    float GetHT(const RVec<Jet> &jets);
    bool IsHEMElectron(const Electron& electron) const;
    bool IsHEMElectron(const ElectronView &electron) const;

    // Gen Matching
    void PrintGen(const RVec<Gen> &gens);
    static RVec<int> TrackGenSelfHistory(const Gen& me, const RVec<Gen>& gens);
    static Gen GetGenMatchedLepton(const Lepton& lep, const RVec<Gen>& gens);
    static Gen GetGenMatchedMuon(const Muon& muon, const RVec<Gen>& gens);
    static Gen GetGenMatchedPhoton(const Lepton& lep, const RVec<Gen>& gens);
    static bool IsFinalPhotonSt23_Public(const RVec<Gen>& gens);
    bool IsFromHadron(const Gen& me, const RVec<Gen>& gens);
    bool IsSignalPID(const int &pid);
    int GetLeptonType(const Lepton& lep, const RVec<Gen>& gens);
    int GetLeptonType(const Gen& gen, const RVec<Gen>& gens);
    int GetLeptonType_Public(const int& genIdx, const RVec<Gen>& gens);
    int GetGenPhotonType(const Gen& genph, const RVec<Gen>& gens);
    int GetPrElType_InSameSCRange_Public(int genIdx, const RVec<Gen>& gens);

    // Scale and smear
    void METType1Propagation(Particle &MET, RVec<Particle> &original_objects, RVec<Particle> &corrected_objects);
    struct Type1METInfo {
        Particle met;
        MyCorrection::variation jesVar = MyCorrection::variation::nom;
        MyCorrection::variation jerVar = MyCorrection::variation::nom;
        bool skip = false;
    };
    Type1METInfo PropagateType1MET(const Event &ev, JetViewCollection &jets, const std::string &systTarget, MyCorrection::variation variation, const TString &systSource, bool doBreakdown = false, Event::MET_Type met_type = Event::MET_Type::PUPPI) const;
    float GetL1PrefireWeight(MyCorrection::variation syst = MyCorrection::variation::nom);
    unordered_map<int, int> GenJetMatching(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const float &rho, const float dR = 0.2, const float pTJerCut = 3.);
    unordered_map<int, int> deltaRMatching(const RVec<Particle> &objs1, const RVec<Particle> &objs2, const float dR = 0.4);
    template <typename ViewCollection>
    std::vector<int> MatchViewsToTrigObjs(const ViewCollection &objects, const TrigObjViewCollection &trigobjs, float dR = 0.4f, int trigId = -1) const;
    RVec<Muon> ScaleMuons(const RVec<Muon> &muons, const MyCorrection::variation &syst=MyCorrection::variation::nom);
    RVec<Electron> ScaleElectrons(const Event &ev, const RVec<Electron> &electrons, const MyCorrection::variation &syst=MyCorrection::variation::nom);
    RVec<Electron> SmearElectrons(const RVec<Electron> &electrons, const MyCorrection::variation &syst=MyCorrection::variation::nom);

    RVec<Jet> SmearJets(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const MyCorrection::variation &syst=MyCorrection::variation::nom, const TString &source = "total");
    RVec<Jet> ScaleJets(const RVec<Jet> &jets, const MyCorrection::variation &syst=MyCorrection::variation::nom, const TString &source = "total");
    RVec<Jet> SmearJets(const JetViewCollection &jets, const std::vector<std::size_t> &indices, const RVec<GenJet> &genjets, const MyCorrection::variation &syst=MyCorrection::variation::nom, const TString &source = "total");
    RVec<Jet> ScaleJets(const JetViewCollection &jets, const std::vector<std::size_t> &indices, const MyCorrection::variation &syst=MyCorrection::variation::nom, const TString &source = "total");
    void ApplyJetScaleVariation(JetViewCollection &jets, const TString &source = "total") const;
    void ApplyJetSmearVariation(JetViewCollection &jets, const RVec<GenJet> &genjets, const TString &source = "total") const;
    
    // Histogram Handlers
    TFile* GetOutfile() { return outfile; }
    inline void SetOutfilePath(const TString &outpath) { outfile = new TFile(outpath, "RECREATE"); }
    TH1D* GetHist1D(const string &histname);
    bool PassJetVetoMap(const Jet &jet, const MuonViewCollection &AllMuons, const TString mapCategory="jetvetomap");
    bool PassJetVetoMap(const Jet &jet, const RVec<Muon> &AllMuons, const TString mapCategory="jetvetomap");
    bool PassJetVetoMap(const JetViewCollection &AllJets, const MuonViewCollection &AllMuons, const TString mapCategory="jetvetomap");
    bool PassJetVetoMap(const RVec<Jet> &AllJets, const MuonViewCollection &AllMuons, const TString mapCategory="jetvetomap");
    bool PassJetVetoMap(const RVec<Jet> &AllJets, const RVec<Muon> &AllMuons, const TString mapCategory="jetvetomap");
    inline void FillCutFlow(const int &val,const int &maxCutN=10){
        static int storedMaxCutN = maxCutN;
        FillHist("CutFlow", val, 1., storedMaxCutN, 0, storedMaxCutN);
    }
    void FillHist(std::string_view histname, float value, float weight, int n_bin, float x_min, float x_max);
    void FillHist(std::string_view histname, float value, float weight, int n_bin, float *xbins);
    void FillHist(std::string_view histname, float value_x, float value_y, float weight, 
                                          int n_binx, float x_min, float x_max, 
                                          int n_biny, float y_min, float y_max);
    void FillHist(std::string_view histname, float value_x, float value_y, float weight,
                                          int n_binx, float *xbins,
                                          int n_biny, float *ybins);
    void FillHist(std::string_view histname, float value_x, float value_y, float value_z, float weight,
                                          int n_binx, float x_min, float x_max,
                                          int n_biny, float y_min, float y_max,
                                          int n_binz, float z_min, float z_max);
    void FillHist(std::string_view histname, float value_x, float value_y, float value_z, float weight,
                                          int n_binx, float *xbins,
                                          int n_biny, float *ybins,
                                          int n_binz, float *zbins);
    inline void FillHist(std::string_view histname, float value, float weight, const RVec<float> &xbins) { FillHist(histname, value, weight, xbins.size()-1, const_cast<float*>(xbins.data())); } 
    inline void FillHist(std::string_view histname, float value_x, float value_y, float weight, const RVec<float> &xbins, const RVec<float> &ybins) {FillHist(histname, value_x, value_y, weight, xbins.size() - 1, const_cast<float *>(xbins.data()), ybins.size() - 1, const_cast<float *>(ybins.data())); }
    inline void FillHist(std::string_view histname, float value_x, float value_y, float value_z, float weight, const RVec<float> &xbins, const RVec<float> &ybins, const RVec<float> &zbins) {FillHist(histname, value_x, value_y, value_z, weight, xbins.size() - 1, const_cast<float *>(xbins.data()), ybins.size() - 1, const_cast<float *>(ybins.data()), zbins.size() - 1, const_cast<float *>(zbins.data())); }


    TTree* NewTree(const TString &treename, const RVec<TString> &keeps = {}, const RVec<TString> &drops = {});
    TTree* GetTree(const TString &treename);
    inline void SetBranch(const TString &treename, const TString &branchname, float val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = scalar_float_storage[key];
        if (!storage)
            storage = std::make_unique<float>();
        *storage = val;
        SetBranch(treename, branchname, (void *)(storage.get()), branchname + "/F");
    };
    inline void SetBranch(const TString &treename, const TString &branchname, double val) {
        SetBranch(treename, branchname, static_cast<float>(val));
    };
    inline void SetBranch(const TString &treename, const TString &branchname, int val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = scalar_int_storage[key];
        if (!storage)
            storage = std::make_unique<int>();
        *storage = val;
        SetBranch(treename, branchname, (void *)(storage.get()), branchname + "/I");
    };
    inline void SetBranch(const TString &treename, const TString &branchname, bool val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = scalar_bool_storage[key];
        if (!storage)
            storage = std::make_unique<bool>();
        *storage = val;
        SetBranch(treename, branchname, (void *)(storage.get()), branchname + "/O");
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const std::vector<int> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_int_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<int>>();
        *storage = val;
        SetBranch_Vector(treename, branchname, *storage);
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const std::vector<float> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_float_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<float>>();
        *storage = val;
        SetBranch_Vector(treename, branchname, *storage);
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const std::vector<double> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_double_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<double>>();
        *storage = val;
        SetBranch_Vector(treename, branchname, *storage);
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const std::vector<bool> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_bool_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<bool>>();
        *storage = val;
        SetBranch_Vector(treename, branchname, *storage);
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const RVec<int> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_int_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<int>>();
        storage->assign(val.begin(), val.end());
        SetBranch_Vector(treename, branchname, *storage);
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const RVec<float> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_float_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<float>>();
        storage->assign(val.begin(), val.end());
        SetBranch_Vector(treename, branchname, *storage);
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const RVec<double> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_double_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<double>>();
        storage->assign(val.begin(), val.end());
        SetBranch_Vector(treename, branchname, *storage);
    }
    inline void SetBranch(const TString &treename, const TString &branchname, const RVec<bool> &val) {
        const std::string key =
            std::string(treename.Data()) + "/" + std::string(branchname.Data());
        auto &storage = vector_bool_storage[key];
        if (!storage)
            storage = std::make_unique<std::vector<bool>>();
        storage->assign(val.begin(), val.end());
        SetBranch_Vector(treename, branchname, *storage);
    }

    void FillTrees(const TString &treename="");
    virtual void WriteHist();

protected:
    template <typename JetCollection>
    bool PassNoiseFilterCommon(const JetCollection &AllJets, const Event &ev, Event::MET_Type met_type) const;
    std::shared_ptr<JetSoA> CreateJetSoA() const;
    void InitialiseJetSystematics(JetSoA &storage) const;
    void PopulateJetStorageWithoutCorrections(JetSoA &storage) const;
    void ApplyJetEnergyCorrections(JetSoA &storage, float rho);
    std::vector<int> MatchJetsToGenJets(const JetViewCollection& jets,
    const GenJetViewCollection& genjets,
    float rho) const;
    void ApplyJetSmearingAndUncertainties(const JetViewCollection& jets,
    const std::vector<int>& matchedGenIdx,
    const GenJetViewCollection& genjets,
    bool isMC,
    float rho);
    const std::vector<TString> &JetEnergyScaleSources() const;
    bool useTH1F = false;
    absl::flat_hash_map<std::string, TH1*, TransparentStringHash, TransparentStringEq> histmap1d;
    absl::flat_hash_map<std::string, TH2*, TransparentStringHash, TransparentStringEq> histmap2d;
    absl::flat_hash_map<std::string, TH3*, TransparentStringHash, TransparentStringEq> histmap3d;
    unordered_map<string, TTree*> treemap;
    unordered_map<TTree*, unordered_map<string, TBranch*>> branchmaps; 
    std::unordered_map<std::string, std::unique_ptr<float>> scalar_float_storage;
    std::unordered_map<std::string, std::unique_ptr<int>> scalar_int_storage;
    std::unordered_map<std::string, std::unique_ptr<bool>> scalar_bool_storage;
    std::unordered_map<std::string, std::unique_ptr<std::vector<int>>> vector_int_storage;
    std::unordered_map<std::string, std::unique_ptr<std::vector<float>>> vector_float_storage;
    std::unordered_map<std::string, std::unique_ptr<std::vector<double>>> vector_double_storage;
    std::unordered_map<std::string, std::unique_ptr<std::vector<bool>>> vector_bool_storage;
    std::unordered_map<std::string, std::vector<int>*> vector_int_ptr_storage;
    std::unordered_map<std::string, std::vector<float>*> vector_float_ptr_storage;
    std::unordered_map<std::string, std::vector<double>*> vector_double_ptr_storage;
    std::unordered_map<std::string, std::vector<bool>*> vector_bool_ptr_storage;
    std::unordered_map<std::string, ModellingPatch> modelling_patches;
    nlohmann::json modelling_json;
    TFile *outfile;
    void SetBranch(const TString &treename, const TString &branchname, void *address, const TString &leaflist);
    template <typename T>
    void SetBranch_Vector(const TString &treename, const TString &branchname, std::vector<T> &address) {
        try {
            TTree *tree = GetTree(treename);
            const std::string key = std::string(treename.Data()) + "/" + std::string(branchname.Data());

            unordered_map<string, TBranch *> *this_branchmap = &branchmaps[tree];
            auto it = this_branchmap->find(string(branchname));

            std::vector<T>** ptr_to_ptr = nullptr;
            if constexpr (std::is_same_v<T, int>) {
                vector_int_ptr_storage[key] = &address;
                ptr_to_ptr = &vector_int_ptr_storage[key];
            } else if constexpr (std::is_same_v<T, float>) {
                vector_float_ptr_storage[key] = &address;
                ptr_to_ptr = &vector_float_ptr_storage[key];
            } else if constexpr (std::is_same_v<T, double>) {
                vector_double_ptr_storage[key] = &address;
                ptr_to_ptr = &vector_double_ptr_storage[key];
            } else if constexpr (std::is_same_v<T, bool>) {
                vector_bool_ptr_storage[key] = &address;
                ptr_to_ptr = &vector_bool_ptr_storage[key];
            } else {
                static_assert(std::is_same_v<T, int> || std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, bool>, "Unsupported vector type in SetBranch_Vector");
            }

            if (it == this_branchmap->end()) {
                auto br = tree->Branch(branchname, ptr_to_ptr);
                this_branchmap->insert({string(branchname), br});
            } else {
                it->second->SetAddress(ptr_to_ptr);
            }
        } catch (int e) {
            cout << "[AnalyzerCore::SetBranch] Error get tree: " << treename.Data() << endl;
            exit(e);
        }
    }
};

template <typename JetCollection>
bool AnalyzerCore::PassNoiseFilterCommon(const JetCollection &AllJets, const Event &ev, Event::MET_Type met_type) const
{
    if (Run == 2)
    {
        bool passNoiseFilter = Flag_goodVertices
                            && Flag_globalSuperTightHalo2016Filter
                            && Flag_HBHENoiseFilter
                            && Flag_HBHENoiseIsoFilter
                            && Flag_EcalDeadCellTriggerPrimitiveFilter
                            && Flag_BadPFMuonFilter
                            && Flag_BadPFMuonDzFilter
                            && Flag_hfNoisyHitsFilter
                            && Flag_eeBadScFilter;
        if (DataEra.Contains("2017") || DataEra.Contains("2018"))
            passNoiseFilter = passNoiseFilter && Flag_ecalBadCalibFilter;
        return passNoiseFilter;
    }

    bool passNoiseFilter = Flag_goodVertices
                        && Flag_globalSuperTightHalo2016Filter
                        && Flag_EcalDeadCellTriggerPrimitiveFilter
                        && Flag_BadPFMuonFilter
                        && Flag_BadPFMuonDzFilter
                        && Flag_hfNoisyHitsFilter
                        && Flag_eeBadScFilter;
    if (!passNoiseFilter)
        return false;

    if (!IsDATA)
        return passNoiseFilter;
    if (!(362433 <= RunNumber && RunNumber <= 367144))
        return passNoiseFilter;

    const Particle METv = ev.GetMETVector(met_type);
    if (METv.Pt() <= 100.)
        return passNoiseFilter;

    RVec<Jet> problematicJets = SelectJets(AllJets, Jet::JetID::NOCUT, 50., 0.5);
    for (const auto &jet : problematicJets)
    {
        bool badEcal = (jet.Pt() > 50.);
        badEcal = badEcal && (-0.5 < jet.Eta() && jet.Eta() < -0.1);
        badEcal = badEcal && (-2.1 < jet.Phi() && jet.Phi() < -1.8);
        badEcal = badEcal && (jet.neEmEF() > 0.9 || jet.chEmEF() > 0.9);
        badEcal = badEcal && jet.DeltaPhi(METv) > 2.9;
        if (badEcal)
            return false;
    }
    return true;
}

template <typename ViewCollection>
std::vector<int> AnalyzerCore::MatchViewsToTrigObjs(
    const ViewCollection &objects, const TrigObjViewCollection &trigobjs,
    float dR, int trigId) const {
  const std::size_t nObj = objects.size();
  std::vector<int> matchedIdx(nObj, -999);
  if (nObj == 0 || trigobjs.size() == 0)
    return matchedIdx;

  const auto &storage = trigobjs.storage();
  if (!storage)
    return matchedIdx;

  const std::size_t nTrig = storage->size();
  if (nTrig == 0)
    return matchedIdx;

  const float maxDR2 = dR * dR;
  constexpr float kPi = 3.14159265358979323846f;
  constexpr float kTwoPi = 6.28318530717958647692f;

  std::vector<std::tuple<std::size_t, std::size_t, float>> candidates;
  candidates.reserve(nObj * 2);

  for (std::size_t i = 0; i < nObj; ++i) {
    const auto &obj = objects[i];
    const float objEta = obj.Eta();
    const float objPhi = obj.Phi();
    for (std::size_t j = 0; j < nTrig; ++j) {
      if (trigId >= 0 && static_cast<int>(storage->id[j]) != trigId)
        continue;

      const float dEta = objEta - storage->eta[j];
      float dPhi = objPhi - storage->phi[j];
      if (dPhi > kPi)
        dPhi -= kTwoPi;
      else if (dPhi <= -kPi)
        dPhi += kTwoPi;

      const float dr2 = dEta * dEta + dPhi * dPhi;
      if (dr2 < maxDR2)
        candidates.emplace_back(i, j, dr2);
    }
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const auto &a, const auto &b) {
              return std::get<2>(a) < std::get<2>(b);
            });

  std::vector<bool> usedObj(nObj, false);
  std::vector<bool> usedTrig(nTrig, false);
  for (const auto &candidate : candidates) {
    const std::size_t objIdx = std::get<0>(candidate);
    const std::size_t trigIdx = std::get<1>(candidate);
    if (usedObj[objIdx] || usedTrig[trigIdx])
      continue;

    matchedIdx[objIdx] = static_cast<int>(trigIdx);
    usedObj[objIdx] = true;
    usedTrig[trigIdx] = true;
  }

  return matchedIdx;
}

#endif
