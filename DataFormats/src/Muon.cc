#include "Muon.h"
#include <cmath>
#include <limits>

ClassImp(Muon)

void Muon::initializeMembers() {
    storage_.reset();
    index_ = std::numeric_limits<std::size_t>::max();
    cachedProperties_.reset();
    lazy_.reset();

    this->SetLeptonFlavour(MUON);

    j_isTracker = false;
    j_isStandalone = false;
    j_isGlobal = false;

    j_looseId = false;
    j_mediumId = false;
    j_mediumPromptId = false;
    j_tightId = false;
    j_softId = false;
    j_softMvaId = false;
    j_triggerIdLoose = false;

    j_highPtId = 0;
    j_miniIsoId = 0;
    j_multiIsoId = 0;
    j_mvaMuId = 0;
    j_pfIsoId = 0;
    j_puppiIsoId = 0;
    j_tkIsoId = 0;

    j_softMva = -999.f;
    j_mvaLowPt = -999.f;
    j_mvaPrompt = -999.f;
    j_softMvaRun3 = -999.f;

    j_nTrackerLayers = -999;

    j_miniAODPt = -999.f;
    j_momentumScaleUp = -999.f;
    j_momentumScaleDown = -999.f;

    j_genPartIdx = -999;
    j_genPartFlav = 0;
    j_jetIdx = -1;
}

Muon::Muon() {
    initializeMembers();
}

Muon::~Muon() {}

Muon::Muon(std::shared_ptr<const MuonSoA> storage, std::size_t index) {
    initializeMembers();
    storage_ = std::move(storage);
    index_ = index;
    if (storage_) {
        const float correctedPt = storage_->correctedPt.empty() ? storage_->pt[index_] : storage_->correctedPt[index_];
        SetPtEtaPhiM(correctedPt, storage_->eta[index_], storage_->phi[index_], storage_->mass[index_]);
        SetCharge(storage_->charge[index_]);
        if (!storage_->miniAODPt.empty())
            SetMiniAODPt(storage_->miniAODPt[index_]);
        if (!storage_->momentumScaleUp.empty())
            SetMomentumScaleUpDown(storage_->momentumScaleUp[index_], storage_->momentumScaleDown[index_]);
    }
}


void Muon::AttachLazyPayload(void *context, EnsureCallback callback, int index) const {
    materialize();
    if (!lazy_) {
        lazy_ = std::make_shared<LazyPayload>();
    }
    lazy_->context = context;
    lazy_->callback = callback;
    lazy_->index = index;
    lazy_->loaded = cachedProperties_;
}

void Muon::DetachLazyPayload() const {
    lazy_.reset();
}

void Muon::ensure(Property property) const {
    const auto bit = static_cast<std::size_t>(property);
    if (cachedProperties_.test(bit))
        return;

    if (storage_) {
        loadFromStorage(property);
        return;
    }

    if (!lazy_ || !lazy_->callback)
        return;

    lazy_->callback(lazy_->context, const_cast<Muon &>(*this), property);
    markLoaded(property);
}

void Muon::markLoaded(Property property) const {
    const auto bit = static_cast<std::size_t>(property);
    cachedProperties_.set(bit);
    if (lazy_)
        lazy_->loaded.set(bit);
}

void Muon::materialize() const {
    if (!storage_)
        return;

    auto self = const_cast<Muon *>(this);
    const auto idx = index_;
    const float correctedPt = storage_->correctedPt.empty() ? storage_->pt[idx] : storage_->correctedPt[idx];
    self->SetPtEtaPhiM(correctedPt, storage_->eta[idx], storage_->phi[idx], storage_->mass[idx]);
    self->SetCharge(storage_->charge[idx]);
    if (!storage_->miniAODPt.empty())
        self->SetMiniAODPt(storage_->miniAODPt[idx]);
    if (!storage_->momentumScaleUp.empty())
        self->SetMomentumScaleUpDown(storage_->momentumScaleUp[idx], storage_->momentumScaleDown[idx]);

    for (std::size_t bit = 0; bit < static_cast<std::size_t>(Property::Count); ++bit) {
        if (!cachedProperties_.test(bit))
            loadFromStorage(static_cast<Property>(bit));
    }

    self->storage_.reset();
    self->index_ = std::numeric_limits<std::size_t>::max();
}

void Muon::loadFromStorage(Property property) const {
    if (!storage_)
        return;

    const auto idx = index_;
    auto self = const_cast<Muon *>(this);

    switch (property) {
    case Property::TkRelIso:
        self->SetTkRelIso(storage_->tkRelIso[idx]);
        break;
    case Property::PfRelIso03:
        self->SetPfRelIso03(storage_->pfRelIso03[idx]);
        break;
    case Property::PfRelIso04:
        self->SetPfRelIso04(storage_->pfRelIso04[idx]);
        break;
    case Property::MiniPFRelIso:
        self->SetMiniPFRelIso(storage_->miniPFRelIsoAll[idx]);
        break;
    case Property::Dxy:
        self->SetdXY(storage_->dxy[idx]);
        break;
    case Property::DxyErr:
        self->SetdXYErr(storage_->dxyErr[idx]);
        break;
    case Property::Dz:
        self->SetdZ(storage_->dz[idx]);
        break;
    case Property::DzErr:
        self->SetdZErr(storage_->dzErr[idx]);
        break;
    case Property::Ip3d:
        self->SetIP3D(storage_->ip3d[idx]);
        break;
    case Property::Sip3d:
        self->SetSIP3D(storage_->sip3d[idx]);
        break;
    case Property::TrackerLayers:
        self->SetNTrackerLayers(storage_->nTrackerLayers[idx]);
        break;
    case Property::GenPartFlav:
        self->SetGenPartFlav(storage_->genPartFlav[idx]);
        break;
    case Property::GenPartIdx:
        self->SetGenPartIdx(storage_->genPartIdx[idx]);
        break;
    case Property::JetIdx:
        self->SetJetIdx(storage_->jetIdx[idx]);
        break;
    case Property::LooseId:
        self->SetBIDBit(BooleanID::LOOSE, storage_->looseId[idx]);
        break;
    case Property::MediumId:
        self->SetBIDBit(BooleanID::MEDIUM, storage_->mediumId[idx]);
        break;
    case Property::MediumPromptId:
        self->SetBIDBit(BooleanID::MEDIUMPROMPT, storage_->mediumPromptId[idx]);
        break;
    case Property::TightId:
        self->SetBIDBit(BooleanID::TIGHT, storage_->tightId[idx]);
        break;
    case Property::SoftId:
        self->SetBIDBit(BooleanID::SOFT, storage_->softId[idx]);
        break;
    case Property::SoftMvaId:
        self->SetBIDBit(BooleanID::SOFTMVA, storage_->softMvaId[idx]);
        break;
    case Property::TriggerLooseId:
        self->SetBIDBit(BooleanID::TRIGGERLOOSE, storage_->triggerLooseId[idx]);
        break;
    case Property::HighPtId:
        self->SetWIDBit(WorkingPointID::HIGHPT, storage_->highPtId[idx]);
        break;
    case Property::MiniIsoId:
        self->SetWIDBit(WorkingPointID::MINIISO, storage_->miniIsoId[idx]);
        break;
    case Property::MultiIsoId:
        self->SetWIDBit(WorkingPointID::MULTIISO, storage_->multiIsoId[idx]);
        break;
    case Property::MvaMuId:
        self->SetWIDBit(WorkingPointID::MVAMU, storage_->mvaMuId[idx]);
        break;
    case Property::PfIsoId:
        self->SetWIDBit(WorkingPointID::PFISO, storage_->pfIsoId[idx]);
        break;
    case Property::PuppiIsoId:
        self->SetWIDBit(WorkingPointID::PUPPIISO, storage_->puppiIsoId[idx]);
        break;
    case Property::TkIsoId:
        self->SetWIDBit(WorkingPointID::TKISO, storage_->tkIsoId[idx]);
        break;
    case Property::SoftMva:
        self->SetMVAID(MVAID::SOFTMVA, storage_->softMva[idx]);
        break;
    case Property::MvaLowPt:
        self->SetMVAID(MVAID::MVALOWPT, storage_->mvaLowPt[idx]);
        break;
    case Property::MvaPrompt:
        self->SetMVAID(MVAID::MVAPROMPT, storage_->mvaPrompt[idx]);
        break;
    case Property::SoftMvaRun3:
        self->SetMVAID(MVAID::SOFTMVA, storage_->softMvaRun3[idx]);
        break;
    default:
        break;
    }

    markLoaded(property);
}

void Muon::SetBIDBit(BooleanID id, bool idbit) {
    switch (id) {
        case BooleanID::LOOSE:
            j_looseId = idbit;
            markLoaded(Property::LooseId);
            break;
        case BooleanID::MEDIUM:
            j_mediumId = idbit;
            markLoaded(Property::MediumId);
            break;
        case BooleanID::MEDIUMPROMPT:
            j_mediumPromptId = idbit;
            markLoaded(Property::MediumPromptId);
            break;
        case BooleanID::TIGHT:
            j_tightId = idbit;
            markLoaded(Property::TightId);
            break;
        case BooleanID::SOFT:
            j_softId = idbit;
            markLoaded(Property::SoftId);
            break;
        case BooleanID::SOFTMVA:
            j_softMvaId = idbit;
            markLoaded(Property::SoftMvaId);
            break;
        case BooleanID::TRIGGERLOOSE:
            j_triggerIdLoose = idbit;
            markLoaded(Property::TriggerLooseId);
            break;
        default:
            break;
    }
}

void Muon::SetWIDBit(WorkingPointID id, unsigned char wp) {
    switch (id) {
        case WorkingPointID::HIGHPT:
            j_highPtId = static_cast<unsigned char>(wp + 1);
            markLoaded(Property::HighPtId);
            break;
        case WorkingPointID::MINIISO:
            j_miniIsoId = static_cast<unsigned char>(wp + 1);
            markLoaded(Property::MiniIsoId);
            break;
        case WorkingPointID::MULTIISO:
            j_multiIsoId = static_cast<unsigned char>(wp + 1);
            markLoaded(Property::MultiIsoId);
            break;
        case WorkingPointID::MVAMU:
            j_mvaMuId = static_cast<unsigned char>(wp + 2);
            markLoaded(Property::MvaMuId);
            break;
        case WorkingPointID::PFISO:
            j_pfIsoId = static_cast<unsigned char>(wp);
            markLoaded(Property::PfIsoId);
            break;
        case WorkingPointID::PUPPIISO:
            j_puppiIsoId = static_cast<unsigned char>(wp + 1);
            markLoaded(Property::PuppiIsoId);
            break;
        case WorkingPointID::TKISO:
            j_tkIsoId = static_cast<unsigned char>(wp + 1);
            markLoaded(Property::TkIsoId);
            break;
        default:
            break;
    }
}

void Muon::SetMVAID(MVAID id, float score) {
    switch (id) {
        case MVAID::SOFTMVA:
            j_softMva = score;
            markLoaded(Property::SoftMva);
            break;
        case MVAID::MVALOWPT:
            j_mvaLowPt = score;
            markLoaded(Property::MvaLowPt);
            break;
        case MVAID::MVAPROMPT:
            j_mvaPrompt = score;
            markLoaded(Property::MvaPrompt);
            break;
        case MVAID::SOFTMVARUN3:
            j_softMvaRun3 = score;
            markLoaded(Property::SoftMvaRun3);
            break;
        default:
            break;
    }
}

bool Muon::PassID(const TString ID) const {
    if (ID == "")                 return true;
    if (ID == "NOCUT")            return true;
    if (ID == "POGTight")         return isPOGTightId();
    if (ID == "POGMedium")        return isPOGMediumId();
    if (ID == "POGMediumPrompt")  return isPOGMediumPromptId();
    if (ID == "POGLoose")         return isPOGLooseId();
    if (ID == "POGSoft")          return isPOGSoftId();
    if (ID == "POGSoftMVA")       return isPOGSoftMvaId();
    if (ID == "POGTriggerLoose")  return isPOGTriggerIdLoose();
    if (ID == "POGTrackerHighPt") return static_cast<int>(HighPtId()) == 1;
    if (ID == "POGGlobalHighPt")  return static_cast<int>(HighPtId()) == 2;
    if (ID == "POGMiniIsoLoose")  return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
    if (ID == "POGMiniIsoMedium") return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
    if (ID == "POGMiniIsoTight")  return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::TIGHT);
    if (ID == "POGMiniIsoVTight") return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::VTIGHT);
    if (ID == "POGMultiIsoLoose") return static_cast<int>(MultiIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
    if (ID == "POGMultiIsoMedium")return static_cast<int>(MultiIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
    if (ID == "POGMvaMuMedium")   return static_cast<int>(MvaMuId()) >= static_cast<int>(WorkingPoint::MEDIUM);
    if (ID == "POGMvaMuTight")    return static_cast<int>(MvaMuId()) >= static_cast<int>(WorkingPoint::TIGHT);
    if (ID == "POGPfIsoVLoose")   return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::VLOOSE);
    if (ID == "POGPfIsoLoose")    return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
    if (ID == "POGPfIsoMedium")   return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
    if (ID == "POGPfIsoTight")    return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::TIGHT);
    if (ID == "POGPfIsoVTight")   return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::VTIGHT);
    if (ID == "POGPfIsoVVTight")  return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::VVTIGHT);
    if (ID == "POGPuppiIsoLoose") return static_cast<int>(PuppiIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
    if (ID == "POGPuppiIsoMedium")return static_cast<int>(PuppiIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
    if (ID == "POGPuppiIsoTight") return static_cast<int>(PuppiIsoId()) >= static_cast<int>(WorkingPoint::TIGHT);
    if (ID == "POGTkIsoLoose")    return static_cast<int>(TkIsoId()) == 1;
    if (ID == "POGTkIsoTight")    return static_cast<int>(TkIsoId()) == 2;
    if (ID == "HcToWATight")      return Pass_HcToWATight();
    if (ID == "HcToWALoose")      return Pass_HcToWALoose();
    if (ID == "POGPromptMVA_WP0p64") return PromptMVA() > 0.64;
    if (ID == "POGSoftMVA_Run3_Tight")   return SoftMvaRun3() > 0.83;
    if (ID == "POGSoftMVA_Run3_Medium")  return SoftMvaRun3() > 0.74;
    if (ID == "POGSoftMVA_Run3_Loose")   return SoftMvaRun3() > 0.63;
    if (ID == "POGSoftMVA_Run3_VLoose")  return SoftMvaRun3() > 0.47;
    throw std::runtime_error("[Muon::PassID] " + std::string(ID.Data()) + " is not implemented.");

    return false;
}

bool Muon::PassID(const MuonID ID) const {
    switch(ID){
        case MuonID::NOCUT:
            return true;
        case MuonID::POG_TIGHT:
            return isPOGTightId();
        case MuonID::POG_MEDIUM:
            return isPOGMediumId();
        case MuonID::POG_MEDIUM_PROMPT:
            return isPOGMediumPromptId();
        case MuonID::POG_LOOSE:
            return isPOGLooseId();
        case MuonID::POG_SOFT:
            return isPOGSoftId();
        case MuonID::POG_SOFT_MVA:
            return isPOGSoftMvaId();
        case MuonID::POG_TRIGGER_LOOSE:
            return isPOGTriggerIdLoose();
        case MuonID::POG_TRACKER_HIGH_PT:
            return static_cast<int>(HighPtId()) == 1;
        case MuonID::POG_GLOBAL_HIGH_PT:
            return static_cast<int>(HighPtId()) == 2;
        case MuonID::POG_MINISO_LOOSE:
            return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
        case MuonID::POG_MINISO_MEDIUM:
            return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
        case MuonID::POG_MINISO_TIGHT:
            return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::TIGHT);
        case MuonID::POG_MINISO_VTIGHT:
            return static_cast<int>(MiniIsoId()) >= static_cast<int>(WorkingPoint::VTIGHT);
        case MuonID::POG_MULTISO_LOOSE:
            return static_cast<int>(MultiIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
        case MuonID::POG_MULTISO_MEDIUM:
            return static_cast<int>(MultiIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
        case MuonID::POG_MVA_MU_MEDIUM:
            return static_cast<int>(MvaMuId()) >= static_cast<int>(WorkingPoint::MEDIUM);
        case MuonID::POG_MVA_MU_TIGHT:
            return static_cast<int>(MvaMuId()) >= static_cast<int>(WorkingPoint::TIGHT);
        case MuonID::POG_PFISO_VLOOSE:
            return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::VLOOSE);
        case MuonID::POG_PFISO_LOOSE:
            return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
        case MuonID::POG_PFISO_MEDIUM:
            return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
        case MuonID::POG_PFISO_TIGHT:
            return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::TIGHT);
        case MuonID::POG_PFISO_VTIGHT:
            return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::VTIGHT);
        case MuonID::POG_PFISO_VVTIGHT:
            return static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::VVTIGHT);
        case MuonID::POG_PUPPIISO_LOOSE:
            return static_cast<int>(PuppiIsoId()) >= static_cast<int>(WorkingPoint::LOOSE);
        case MuonID::POG_PUPPIISO_MEDIUM:
            return static_cast<int>(PuppiIsoId()) >= static_cast<int>(WorkingPoint::MEDIUM);
        case MuonID::POG_PUPPIISO_TIGHT:
            return static_cast<int>(PuppiIsoId()) >= static_cast<int>(WorkingPoint::TIGHT);
        case MuonID::POG_TKISO_LOOSE:
            return static_cast<int>(TkIsoId()) == 1;
        case MuonID::POG_TKISO_TIGHT:
            return static_cast<int>(TkIsoId()) == 2;
        case MuonID::POG_PROMPTMVA_WP0p64:
            return PromptMVA() > 0.64;
        case MuonID::POG_SOFTMVA_RUN3_TIGHT:
            return SoftMvaRun3() > 0.83;
        case MuonID::POG_SOFTMVA_RUN3_MEDIUM:
            return SoftMvaRun3() > 0.74;
        case MuonID::POG_SOFTMVA_RUN3_LOOSE:
            return SoftMvaRun3() > 0.63;
        case MuonID::POG_SOFTMVA_RUN3_VLOOSE:
            return SoftMvaRun3() > 0.47;
        case MuonID::POG_TIGHT_PFISO_TIGHT:
            return isPOGTightId() && (static_cast<int>(PfIsoId()) >= static_cast<int>(WorkingPoint::TIGHT));
        default:
            throw std::runtime_error("[Muon::PassID] MuonID not implemented.");
    }
    return false;
}

bool Muon::Pass_HcToWATight() const {
    if (! isPOGMediumId()) return false;
    if (! (fabs(dZ()) < 0.1)) return false;
    if (! (SIP3D() < 3.)) return false;
    if (! (TkRelIso() < 0.4*Pt())) return false;
    if (! (MiniPFRelIso() < 0.1)) return false;
    return true;
}

bool Muon::Pass_HcToWALoose() const {
    if (! isPOGMediumId()) return false;
    if (! (fabs(dZ()) < 0.1)) return false;
    if (! (SIP3D() < 5.)) return false;
    if (! (TkRelIso() < 0.4*Pt())) return false;
    if (! (MiniPFRelIso() < 0.6)) return false;
    return true;
}
