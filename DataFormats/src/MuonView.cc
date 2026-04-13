#include "MuonView.h"

MuonViewCollection::MuonViewCollection(std::shared_ptr<MuonSoA> data)
    : payload(std::move(data)) {
    if (payload) {
        views.reserve(payload->size());
        for (std::size_t i = 0; i < payload->size(); ++i) {
            views.emplace_back(payload, i);
        }
    }
}

bool MuonView::PassID(MuonID id) const {
    switch (id) {
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
        return store->highPtId[idx] == 1;
    case MuonID::POG_GLOBAL_HIGH_PT:
        return store->highPtId[idx] == 2;
    case MuonID::POG_MINISO_LOOSE:
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::LOOSE);
    case MuonID::POG_MINISO_MEDIUM:
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::MEDIUM);
    case MuonID::POG_MINISO_TIGHT:
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::TIGHT);
    case MuonID::POG_MINISO_VTIGHT:
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::VTIGHT);
    case MuonID::POG_MULTISO_LOOSE:
        return workingPointAtLeast(store->multiIsoId, WorkingPoint::LOOSE);
    case MuonID::POG_MULTISO_MEDIUM:
        return workingPointAtLeast(store->multiIsoId, WorkingPoint::MEDIUM);
    case MuonID::POG_MVA_MU_MEDIUM:
        return workingPointAtLeast(store->mvaMuId, WorkingPoint::MEDIUM);
    case MuonID::POG_MVA_MU_TIGHT:
        return workingPointAtLeast(store->mvaMuId, WorkingPoint::TIGHT);
    case MuonID::POG_PFISO_VLOOSE:
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::VLOOSE);
    case MuonID::POG_PFISO_LOOSE:
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::LOOSE);
    case MuonID::POG_PFISO_MEDIUM:
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::MEDIUM);
    case MuonID::POG_PFISO_TIGHT:
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::TIGHT);
    case MuonID::POG_PFISO_VTIGHT:
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::VTIGHT);
    case MuonID::POG_PFISO_VVTIGHT:
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::VVTIGHT);
    case MuonID::POG_PUPPIISO_LOOSE:
        return workingPointAtLeast(store->puppiIsoId, WorkingPoint::LOOSE);
    case MuonID::POG_PUPPIISO_MEDIUM:
        return workingPointAtLeast(store->puppiIsoId, WorkingPoint::MEDIUM);
    case MuonID::POG_PUPPIISO_TIGHT:
        return workingPointAtLeast(store->puppiIsoId, WorkingPoint::TIGHT);
    case MuonID::POG_TKISO_LOOSE:
        return store->tkIsoId[idx] == 1;
    case MuonID::POG_TKISO_TIGHT:
        return store->tkIsoId[idx] == 2;
    case MuonID::POG_PROMPTMVA_WP0p64:
        return MvaPrompt() > 0.64;
    case MuonID::POG_SOFTMVA_RUN3_TIGHT:
        return SoftMvaRun3() > 0.83;
    case MuonID::POG_SOFTMVA_RUN3_MEDIUM:
        return SoftMvaRun3() > 0.74;
    case MuonID::POG_SOFTMVA_RUN3_LOOSE:
        return SoftMvaRun3() > 0.63;
    case MuonID::POG_SOFTMVA_RUN3_VLOOSE:
        return SoftMvaRun3() > 0.47;
    case MuonID::POG_TIGHT_PFISO_TIGHT:
        return isPOGTightId() && workingPointAtLeast(store->pfIsoId, WorkingPoint::TIGHT);
    default:
        break;
    }
    return false;
}

bool MuonView::PassID(const TString &id) const {
    if (id == "" || id == "NOCUT")
        return true;
    if (id == "POGTight")
        return isPOGTightId();
    if (id == "POGMedium")
        return isPOGMediumId();
    if (id == "POGMediumPrompt")
        return isPOGMediumPromptId();
    if (id == "POGLoose")
        return isPOGLooseId();
    if (id == "POGSoft")
        return isPOGSoftId();
    if (id == "POGSoftMVA")
        return isPOGSoftMvaId();
    if (id == "POGTriggerLoose")
        return isPOGTriggerIdLoose();
    if (id == "POGTrackerHighPt")
        return store->highPtId[idx] == 1;
    if (id == "POGGlobalHighPt")
        return store->highPtId[idx] == 2;
    if (id == "POGMiniIsoLoose")
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::LOOSE);
    if (id == "POGMiniIsoMedium")
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::MEDIUM);
    if (id == "POGMiniIsoTight")
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::TIGHT);
    if (id == "POGMiniIsoVTight")
        return workingPointAtLeast(store->miniIsoId, WorkingPoint::VTIGHT);
    if (id == "POGMultiIsoLoose")
        return workingPointAtLeast(store->multiIsoId, WorkingPoint::LOOSE);
    if (id == "POGMultiIsoMedium")
        return workingPointAtLeast(store->multiIsoId, WorkingPoint::MEDIUM);
    if (id == "POGMvaMuMedium")
        return workingPointAtLeast(store->mvaMuId, WorkingPoint::MEDIUM);
    if (id == "POGMvaMuTight")
        return workingPointAtLeast(store->mvaMuId, WorkingPoint::TIGHT);
    if (id == "POGPfIsoVLoose")
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::VLOOSE);
    if (id == "POGPfIsoLoose")
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::LOOSE);
    if (id == "POGPfIsoMedium")
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::MEDIUM);
    if (id == "POGPfIsoTight")
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::TIGHT);
    if (id == "POGPfIsoVTight")
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::VTIGHT);
    if (id == "POGPfIsoVVTight")
        return workingPointAtLeast(store->pfIsoId, WorkingPoint::VVTIGHT);
    if (id == "POGPuppiIsoLoose")
        return workingPointAtLeast(store->puppiIsoId, WorkingPoint::LOOSE);
    if (id == "POGPuppiIsoMedium")
        return workingPointAtLeast(store->puppiIsoId, WorkingPoint::MEDIUM);
    if (id == "POGPuppiIsoTight")
        return workingPointAtLeast(store->puppiIsoId, WorkingPoint::TIGHT);
    if (id == "POGTkIsoLoose")
        return store->tkIsoId[idx] == 1;
    if (id == "POGTkIsoTight")
        return store->tkIsoId[idx] == 2;
    if (id == "HcToWATight")
        return Pass_HcToWATight();
    if (id == "HcToWALoose")
        return Pass_HcToWALoose();
    if (id == "POGPromptMVA_WP0p64")
        return MvaPrompt() > 0.64;
    return false;
}

bool MuonView::Pass_HcToWATight() const {
    if (!isPOGMediumId())
        return false;
    if (!(std::fabs(dZ()) < 0.1))
        return false;
    if (!(SIP3D() < 3.))
        return false;
    if (!(TkRelIso() < 0.4 * Pt()))
        return false;
    if (!(MiniPFRelIso() < 0.1))
        return false;
    return true;
}

bool MuonView::Pass_HcToWALoose() const {
    if (!isPOGMediumId())
        return false;
    if (!(std::fabs(dZ()) < 0.1))
        return false;
    if (!(SIP3D() < 5.))
        return false;
    if (!(TkRelIso() < 0.4 * Pt()))
        return false;
    if (!(MiniPFRelIso() < 0.6))
        return false;
    return true;
}
