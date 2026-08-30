#include "AnalyzerCore.h"
#include "JetView.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <Compression.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>

AnalyzerCore::AnalyzerCore() {
  myCorr = nullptr;
  outfile = nullptr;
  if (HasFlag("useTH1F")) {
    cout << "[AnalyzerCore::AnalyzerCore] Using TH1F" << endl;
    useTH1F = true;
  } else {
    cout << "[AnalyzerCore::AnalyzerCore] Using TH1D" << endl;
    useTH1F = false;
  }
  histmap1d.reserve(512);
  histmap2d.reserve(256);
  histmap3d.reserve(128);
  // pdfReweight = new PDFReweight();
}

AnalyzerCore::~AnalyzerCore() {
  for (const auto &pair : histmap1d)
    delete pair.second;
  histmap1d.clear();
  for (const auto &pair : histmap2d)
    delete pair.second;
  histmap2d.clear();
  for (const auto &pair : histmap3d)
    delete pair.second;
  histmap3d.clear();
  if (outfile)
    delete outfile;
  if (myCorr)
    delete myCorr;
  // if (pdfReweight) delete pdfReweight;
}

void AnalyzerCore::load_modelling_json(const TString &filename) {
  using json = nlohmann::json;

  std::cout << "Loading modelling json file: " << filename.Data() << std::endl;
  std::ifstream file_stream(filename.Data());
  if (!file_stream.is_open()) {
    throw std::runtime_error("Could not open modelling json file: " +
                             std::string(filename.Data()));
  }

  json j;
  file_stream >> j;
  file_stream.close();

  modelling_json = j;

  if (!j.contains("PD")) {
    throw std::runtime_error("modelling_json does not contain 'PD' field");
  }
  std::string pd_key = j["PD"].get<std::string>();

  if (!j.contains("subprocesses")) {
    throw std::runtime_error(
        "modelling_json does not contain 'subprocesses' field");
  }
  const auto &subprocs = j["subprocesses"];

  modelling_patches.clear();

  for (auto it = subprocs.begin(); it != subprocs.end(); ++it) {
    const std::string sub_name = it.key();
    const auto &sub = it.value();

    ModellingPatch mp;

    if (sub.contains("patch_ScaleVariation") &&
        sub["patch_ScaleVariation"].is_array()) {
      for (const auto &v : sub["patch_ScaleVariation"]) {
        if (v.is_null())
          continue;
        mp.patch_ScaleVariation.push_back(static_cast<float>(v.get<double>()));
      }
    }

    if (sub.contains("patch_PSVariation") &&
        sub["patch_PSVariation"].is_array()) {
      for (const auto &v : sub["patch_PSVariation"]) {
        if (v.is_null())
          continue;
        mp.patch_PSVariation.push_back(static_cast<float>(v.get<double>()));
      }
    }

    if (sub.contains("patch_hdamp_up") && !sub["patch_hdamp_up"].is_null())
      mp.patch_hdamp_up = sub["patch_hdamp_up"].get<float>();
    else
      mp.patch_hdamp_up = 0.0f;

    if (sub.contains("patch_hdamp_down") && !sub["patch_hdamp_down"].is_null())
      mp.patch_hdamp_down = sub["patch_hdamp_down"].get<float>();
    else
      mp.patch_hdamp_down = 0.0f;

    if (sub.contains("patch_minnlo") && !sub["patch_minnlo"].is_null())
      mp.patch_minnlo = sub["patch_minnlo"].get<float>();
    else
      mp.patch_minnlo = 0.0f;

    modelling_patches.emplace(sub_name, std::move(mp));
  }

  std::cout << "[modelling] PD = " << pd_key << std::endl;
  std::cout << "[modelling] loaded subprocesses: ";
  for (const auto &kv : modelling_patches) {
    std::cout << kv.first << " ";
  }
  std::cout << std::endl;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
bool AnalyzerCore::PassNoiseFilter(const RVec<Jet> &Alljets, const Event &ev,
                                   Event::MET_Type met_type) {
  return PassNoiseFilterCommon(Alljets, ev, met_type);
}

bool AnalyzerCore::PassNoiseFilter(const JetViewCollection &Alljets,
                                   const Event &ev, Event::MET_Type met_type) {
  return PassNoiseFilterCommon(Alljets, ev, met_type);
}

bool AnalyzerCore::PassMetFilter(const RVec<Jet> &AllJets, const Event &ev,
                                 Event::MET_Type met_type) {
  return PassNoiseFilter(AllJets, ev, met_type);
}

bool AnalyzerCore::PassMetFilter(const JetViewCollection &AllJets,
                                 const Event &ev, Event::MET_Type met_type) {
  return PassNoiseFilter(AllJets, ev, met_type);
}

float AnalyzerCore::GetScaleVariation(const MyCorrection::variation &muF_syst,
                                      const MyCorrection::variation &muR_syst) {
  if (nLHEScaleWeight == 0)
    return 1.;
  if (muF_syst == MyCorrection::variation::down &&
      muR_syst == MyCorrection::variation::down)
    return LHEScaleWeight[0];
  else if (muF_syst == MyCorrection::variation::nom &&
           muR_syst == MyCorrection::variation::down)
    return LHEScaleWeight[1];
  else if (muF_syst == MyCorrection::variation::up &&
           muR_syst == MyCorrection::variation::down)
    return LHEScaleWeight[2];
  else if (muF_syst == MyCorrection::variation::down &&
           muR_syst == MyCorrection::variation::nom)
    return LHEScaleWeight[3];
  else if (muF_syst == MyCorrection::variation::up &&
           muR_syst == MyCorrection::variation::nom)
    return LHEScaleWeight[4];
  else if (muF_syst == MyCorrection::variation::down &&
           muR_syst == MyCorrection::variation::up)
    return LHEScaleWeight[5];
  else if (muF_syst == MyCorrection::variation::nom &&
           muR_syst == MyCorrection::variation::up)
    return LHEScaleWeight[6];
  else if (muF_syst == MyCorrection::variation::up &&
           muR_syst == MyCorrection::variation::up)
    return LHEScaleWeight[7];
  else if (muF_syst == MyCorrection::variation::nom &&
           muR_syst == MyCorrection::variation::nom)
    return 1.f;
  else {
    cout << "[AnalyzerCore::GetScaleVariation] muF_syst = " << int(muF_syst)
         << " and muR_syst = " << int(muR_syst) << " is not implemented"
         << endl;
    exit(ENODATA);
  }
}

float AnalyzerCore::GetPSWeight(const MyCorrection::variation &ISR_syst,
                                const MyCorrection::variation &FSR_syst) {
  //[0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1
  // FSR=0.5;
  if (nPSWeight == 1)
    return 1.;
  if (ISR_syst == MyCorrection::variation::up &&
      FSR_syst == MyCorrection::variation::nom)
    return PSWeight[0];
  else if (ISR_syst == MyCorrection::variation::nom &&
           FSR_syst == MyCorrection::variation::up)
    return PSWeight[1];
  else if (ISR_syst == MyCorrection::variation::down &&
           FSR_syst == MyCorrection::variation::nom)
    return PSWeight[2];
  else if (ISR_syst == MyCorrection::variation::nom &&
           FSR_syst == MyCorrection::variation::down)
    return PSWeight[3];
  else {
    cout << "[AnalyzerCore::GetPSWeight] ISR_syst = " << int(ISR_syst)
         << " and FSR_syst = " << int(FSR_syst) << " is not implemented"
         << endl;
    exit(ENODATA);
  }
}

float AnalyzerCore::GetHT(const RVec<Jet> &jets) {
  return std::transform_reduce(std::execution::unseq, jets.begin(), jets.end(),
                               0.f, std::plus<float>(),
                               [](const Jet &jet) { return jet.Pt(); });
}

void AnalyzerCore::METType1Propagation(Particle &MET,
                                       RVec<Particle> &original_objects,
                                       RVec<Particle> &corrected_objects) {
  const auto particle_reducer = [](Particle lhs, const Particle &rhs) {
    lhs += rhs;
    return lhs;
  };

  const Particle original_sum =
      std::reduce(std::execution::unseq, original_objects.begin(),
                  original_objects.end(), Particle(), particle_reducer);

  const Particle corrected_sum =
      std::reduce(std::execution::unseq, corrected_objects.begin(),
                  corrected_objects.end(), Particle(), particle_reducer);

  MET += (original_sum - corrected_sum);
}

float AnalyzerCore::GetL1PrefireWeight(MyCorrection::variation syst) {
  if (Run == 3)
    return 1.;
  if (IsDATA)
    return 1.;
  float weight = 1.;
  throw std::logic_error(
      "[TODO] I cannot find L1Prefiring weight in the v15, please update");
  // switch(syst){
  //     case MyCorrection::variation::nom:
  //         weight = L1PreFiringWeight_Nom;
  //         break;
  //     case MyCorrection::variation::up:
  //         weight = L1PreFiringWeight_Up;
  //         break;
  //     case MyCorrection::variation::down:
  //         weight = L1PreFiringWeight_Dn;
  //         break;
  //     default:
  //         exit(ENODATA);
  // }
  return weight;
}

unordered_map<int, int>
AnalyzerCore::GenJetMatching(const RVec<Jet> &jets, const RVec<GenJet> &genjets,
                             const float &rho, const float dR,
                             const float pTJerCut) {
  unordered_map<int, int> matched_genjet_idx;
  RVec<tuple<int, int, float, float>> possible_matches;

  // All possible matches that pass the dR and pTJerCut
  for (size_t i = 0; i < jets.size(); i++) {
    for (size_t j = 0; j < genjets.size(); j++) {
      float this_DeltaR = jets[i].DeltaR(genjets[j]);
      float this_pt_diff = fabs(jets[i].Pt() - genjets[j].Pt());
      float this_jer = myCorr->GetJER(jets[i].Eta(), jets[i].Pt(), rho);
      this_jer *= jets[i].Pt();
      if (jets[i].DeltaR(genjets[j]) < dR &&
          this_pt_diff < pTJerCut * this_jer) {
        possible_matches.emplace_back(i, j, this_DeltaR, this_pt_diff);
      }
    }
  }
  // sort them by dR and then pt_diff
  sort(possible_matches.begin(), possible_matches.end(),
       [](const tuple<int, int, float, float> &a,
          const tuple<int, int, float, float> &b) {
         if (get<2>(a) == get<2>(b))
           return get<3>(a) < get<3>(b);
         return get<2>(a) < get<2>(b);
       });
  // greeedy matching
  RVec<bool> used_jet(jets.size(), false);
  RVec<bool> used_genjet(genjets.size(), false);
  for (const auto &match : possible_matches) {
    int jet_idx = get<0>(match);
    int genjet_idx = get<1>(match);
    if (used_jet[jet_idx] || used_genjet[genjet_idx])
      continue;
    matched_genjet_idx[jet_idx] = genjet_idx;
    used_jet[jet_idx] = true;
    used_genjet[genjet_idx] = true;
  }
  // assign -999 to unmatched jet
  for (size_t i = 0; i < jets.size(); i++) {
    if (used_jet[i])
      continue;
    else
      matched_genjet_idx[i] = -999;
  }
  return matched_genjet_idx;
}

unordered_map<int, int>
AnalyzerCore::deltaRMatching(const RVec<Particle> &objs1,
                             const RVec<Particle> &objs2, const float dR) {
  RVec<tuple<size_t, size_t, float>> possible_matches;
  RVec<bool> used_obj1(objs1.size(), false);
  RVec<bool> used_obj2(objs2.size(), false);
  for (size_t i = 0; i < objs1.size(); i++) {
    for (size_t j = 0; j < objs2.size(); j++) {
      float this_dR = objs1[i].DeltaR(objs2[j]);
      if (this_dR < dR) {
        possible_matches.emplace_back(i, j, this_dR);
      }
    }
  }
  sort(possible_matches.begin(), possible_matches.end(),
       [](const tuple<size_t, size_t, float> &a,
          const tuple<size_t, size_t, float> &b) {
         return get<2>(a) < get<2>(b);
       });
  unordered_map<int, int> matched_idx;
  // print all
  for (size_t i = 0; i < possible_matches.size(); i++) {
    int obj1_idx = get<0>(possible_matches[i]);
    int obj2_idx = get<1>(possible_matches[i]);
    if (used_obj1[obj1_idx] || used_obj2[obj2_idx])
      continue;
    matched_idx[obj1_idx] = obj2_idx;
    used_obj1[obj1_idx] = true;
    used_obj2[obj2_idx] = true;
  }

  for (size_t i = 0; i < objs1.size(); i++) {
    if (used_obj1[i])
      continue;
    matched_idx[i] = -999;
  }

  return matched_idx;
}

RVec<Jet> AnalyzerCore::SmearJets(const RVec<Jet> &jets,
                                  const RVec<GenJet> &genjets,
                                  const MyCorrection::variation &syst,
                                  const TString &source) {
  unordered_map<int, int> matched_idx =
      GenJetMatching(jets, genjets, Rho_fixedGridRhoFastjetAll);
  RVec<Jet> smeared_jets;
  gRandom->SetSeed(int(PuppiMET_pt * 1e6));
  float MIN_JET_ENERGY = 1e-2;
  for (size_t i = 0; i < jets.size(); i++) {
    Jet this_jet = jets.at(i);

    // backward smearing for systematic variation
    if (this_jet.GetUnsmearedP4().E() < 0) {
      // never smeared yet
      this_jet.SetUnsmearedP4(this_jet);
    } else {
      // already smeared
      const float backward_factor =
          this_jet.GetUnsmearedP4().Pt() / this_jet.Pt();
      this_jet *= backward_factor;
    }

    // Following the procedures in
    // https://github.com/choij1589/SKFlatMaker/blob/master/SKFlatMaker/src/SKFlatMaker.cc#L3378-L3419
    float this_corr = 1.;
    const float this_jer = myCorr->GetJER(this_jet.Eta(), this_jet.Pt(),
                                          Rho_fixedGridRhoFastjetAll);
    const float this_sf =
        myCorr->GetJERSF(this_jet.Eta(), this_jet.Pt(), syst, source);
    if (matched_idx[i] > 0 && matched_idx[i] < genjets.size()) {
      // found matched jet
      const float matched_genjet_pt = genjets[matched_idx[i]].Pt();
      this_corr += (this_sf - 1.) * (1. - matched_genjet_pt / this_jet.Pt());
    } else {
      this_corr +=
          (gRandom->Gaus(0., this_jer)) * sqrt(max(this_sf * this_sf - 1., 0.));
    }
    // To avoid flipping direction (this_corr < 0)
    const float min_corr = MIN_JET_ENERGY / this_jet.E();
    this_corr = max(this_corr, min_corr);

    /*
    if(matched_idx[i] < 0) {
        // if the jet is not matched to any genjet, do stochastic smearing
        if(this_sf < 1.) {
            this_corr = MIN_JET_ENERGY / this_jet.E();
        } else{
            this_corr += (gRandom->Gaus(0, this_jer)) * sqrt(max(this_sf *
    this_sf - 1., 0.0)); float new_corr = MIN_JET_ENERGY / this_jet.E();
            this_corr = max(this_corr, new_corr);
        }
    } else {
        float this_genjet_pt = genjets[matched_idx[i]].Pt();
        this_corr += (this_sf - 1.) * (1. - this_genjet_pt / this_jet.Pt());
        float new_corr = MIN_JET_ENERGY / this_jet.E();
        this_corr = max(this_corr, new_corr);
    }
    */
    this_jet *= this_corr;
    smeared_jets.push_back(this_jet);
  }
  return smeared_jets;
}

RVec<Jet> AnalyzerCore::ScaleJets(const RVec<Jet> &jets,
                                  const MyCorrection::variation &syst,
                                  const TString &source) {
  if (syst == MyCorrection::variation::nom)
    return jets;

  RVec<Jet> scaled_jets;
  RVec<TString> syst_sources = {
      "AbsoluteMPFBias", "AbsoluteScale",   "AbsoluteStat",   "FlavorQCD",
      "Fragmentation",   "PileUpDataMC",    "PileUpPtBB",     "PileUpPtEC1",
      "PileUpPtEC2",     "PileUpPtHF",      "PileUpPtRef",    "RelativeJEREC1",
      "RelativeJEREC2",  "RelativeJERHF",   "RelativePtBB",   "RelativePtEC1",
      "RelativePtEC2",   "RelativePtHF",    "RelativeBal",    "RelativeSample",
      "RelativeFSR",     "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF",
      "SinglePionECAL",  "SinglePionHCAL",  "TimePtEta"};

  if (syst == MyCorrection::variation::nom)
    return jets;

  for (const auto &jet : jets) {
    Jet this_jet = jet;

    TLorentzVector this_unsmearedP4 = jet.GetUnsmearedP4();
    if (this_unsmearedP4.E() < 0) {
      cerr << "[AnalyzerCore::ScaleJets] It seems the attribute UnsmearedP4 "
              "has never been set for this Jet."
           << endl;
      cerr << "[AnalyzerCore::ScaleJets] Since all jets in MC are smeared when "
              "calling GetAllJets(), UnsmearedP4 should have been set by "
              "AnalyzerCore::SmearJets at the same time."
           << endl;
      cerr << "[AnalyzerCore::ScaleJets] This may indicate spurious behavior "
              "in your code. For now, Jet.pt() is used for scaling, but I "
              "recommend inspecting your code."
           << endl;

      this_unsmearedP4 = jet;
    }

    if (source == "total") {
      for (const auto &it : syst_sources) {
        this_jet *= myCorr->GetJESUncertaintySF(
            this_unsmearedP4.Eta(), this_unsmearedP4.Pt(), syst, it);
      }
      scaled_jets.push_back(this_jet);
    } else {
      this_jet *= myCorr->GetJESUncertaintySF(
          this_unsmearedP4.Eta(), this_unsmearedP4.Pt(), syst, source);
      scaled_jets.push_back(this_jet);
    }
  }
  return scaled_jets;
}

// pdfs
float AnalyzerCore::GetPDFWeight(LHAPDF::PDF *pdf) {
  float pdf1 = pdf->xfxQ2(Generator_id1, Generator_x1, Generator_scalePDF);
  float pdf2 = pdf->xfxQ2(Generator_id2, Generator_x2, Generator_scalePDF);
  return pdf1 * pdf2;
}

float AnalyzerCore::GetPDFReweight() {
  return GetPDFWeight(pdfReweight->NewPDF) / GetPDFWeight(pdfReweight->ProdPDF);
}

float AnalyzerCore::GetPDFReweight(int member) {
  return GetPDFWeight(pdfReweight->PDFErrorSet.at(member)) /
         GetPDFWeight(pdfReweight->ProdPDF);
}

// MC weights
// TODO: Treat MiNNLO and Sherpa weights
float AnalyzerCore::MCweight(bool usesign, bool norm_1invpb) const {
  if (IsDATA)
    return 1.;
  float weight = genWeight;
  if (usesign) {
    if (weight > 0)
      weight = 1.0;
    else if (weight < 0)
      weight = -1.0;
    else
      weight = 0.0;
  }

  if (norm_1invpb) {
    if (usesign)
      weight *= xsec / sumSign;
    else
      weight *= xsec / sumW;
  }
  return weight;
}

// Objects
Event AnalyzerCore::GetEvent() {
  Event ev;
  ev.SetRunLumiEvent(RunNumber, luminosityBlock, event);
  ev.SetnPileUp(Pileup_nPU);
  ev.SetnTrueInt(Pileup_nTrueInt);
  ev.SetnPVsGood(PV_npvsGood);

  ev.SetTrigger(TriggerMap);
  ev.SetEra(GetEra());
  RVec<float> MET_pts = {PuppiMET_pt, PuppiMET_ptUnclusteredUp,
                         PuppiMET_ptUnclusteredDown};
  RVec<float> MET_phis = {PuppiMET_phi, PuppiMET_phiUnclusteredUp,
                          PuppiMET_phiUnclusteredDown};
  ev.SetMET(MET_pts, MET_phis);
  ev.setRho(Rho_fixedGridRhoFastjetAll);
  return ev;
}

GenViewCollection AnalyzerCore::GetAllGenViews() {
  if (IsDATA)
    return {};

  auto storage = std::make_shared<GenSoA>();
  storage->pt.bind(&GenPart_pt);
  storage->eta.bind(&GenPart_eta);
  storage->phi.bind(&GenPart_phi);
  storage->mass.bind(&GenPart_mass);
  storage->pdgId.bind(&GenPart_pdgId);
  storage->status.bind(&GenPart_status);
  storage->motherIdx.bind(&GenPart_genPartIdxMother);
  storage->statusFlags.bind(&GenPart_statusFlags);
  return GenViewCollection(std::move(storage));
}

GenJetViewCollection AnalyzerCore::GetAllGenJetViews() {
  if (IsDATA)
    return {};

  auto storage = std::make_shared<GenJetSoA>();
  storage->pt.bind(&GenJet_pt);
  storage->eta.bind(&GenJet_eta);
  storage->phi.bind(&GenJet_phi);
  storage->mass.bind(&GenJet_mass);
  storage->partonFlavour.bind(&GenJet_partonFlavour);
  storage->hadronFlavour.bind(&GenJet_hadronFlavour);
  return GenJetViewCollection(std::move(storage));
}

MuonViewCollection AnalyzerCore::GetAllMuonViews() {
  auto storage = std::make_shared<MuonSoA>();
  storage->pt.bind(&Muon_pt);
  storage->eta.bind(&Muon_eta);
  storage->phi.bind(&Muon_phi);
  storage->mass.bind(&Muon_mass);
  storage->charge.bind(&Muon_charge);
  storage->tkRelIso.bind(&Muon_tkRelIso);
  storage->pfRelIso03.bind(&Muon_pfRelIso03_all);
  storage->pfRelIso04.bind(&Muon_pfRelIso04_all);
  storage->miniPFRelIsoAll.bind(&Muon_miniPFRelIso_all);
  storage->dxy.bind(&Muon_dxy);
  storage->dxyErr.bind(&Muon_dxyErr);
  storage->dz.bind(&Muon_dz);
  storage->dzErr.bind(&Muon_dzErr);
  storage->ip3d.bind(&Muon_ip3d);
  storage->sip3d.bind(&Muon_sip3d);
  storage->highPtId.bind(&Muon_highPtId);
  storage->looseId.bind(&Muon_looseId);
  storage->mediumId.bind(&Muon_mediumId);
  storage->mediumPromptId.bind(&Muon_mediumPromptId);
  storage->tightId.bind(&Muon_tightId);
  storage->softId.bind(&Muon_softId);
  storage->softMvaId.bind(&Muon_softMvaId);
  storage->triggerLooseId.bind(&Muon_triggerIdLoose);
  storage->miniIsoId.bind(&Muon_miniIsoId);
  storage->multiIsoId.bind(&Muon_multiIsoId);
  storage->mvaMuId.bind(&Muon_mvaMuID_WP);
  storage->pfIsoId.bind(&Muon_pfIsoId);
  storage->puppiIsoId.bind(&Muon_puppiIsoId);
  storage->tkIsoId.bind(&Muon_tkIsoId);
  storage->nTrackerLayers.bind(&Muon_nTrackerLayers);
  storage->softMva.bind(&Muon_softMva);
  storage->softMvaRun3.bind(&Muon_softMvaRun3);
  storage->mvaLowPt.bind(&Muon_mvaLowPt);
  storage->mvaPrompt.bind(&Muon_promptMVA);
  storage->genPartFlav.bind(&Muon_genPartFlav);
  storage->genPartIdx.bind(&Muon_genPartIdx);
  storage->jetIdx.bind(&Muon_jetIdx);

  const std::size_t n = storage->pt.size();
  storage->correctedPt.resize(n, 0.f);
  storage->miniAODPt.resize(n, 0.f);
  storage->momentumScaleUp.resize(n, 0.f);
  storage->momentumScaleDown.resize(n, 0.f);

  if (n == 0)
    return MuonViewCollection(std::move(storage));

  RVec<Gen> truth;
  if (!IsDATA)
    truth = GetAllGens();

  for (std::size_t i = 0; i < n; ++i) {
    Muon probe;
    probe.SetPtEtaPhiM(storage->pt[i], storage->eta[i], storage->phi[i],
                       storage->mass[i]);
    probe.SetCharge(storage->charge[i]);

    float roccor = 1.f;
    float roccor_err = 0.f;
    if (IsDATA) {
      roccor = myCorr->GetMuonScaleSF(probe, MyCorrection::variation::nom);
      roccor_err =
          myCorr->GetMuonScaleSF(probe, MyCorrection::variation::up) - roccor;
    } else {
      Gen matched_gen = GetGenMatchedMuon(probe, truth);
      float matched_pt = matched_gen.Pt();
      roccor = myCorr->GetMuonScaleSF(probe, MyCorrection::variation::nom,
                                      matched_pt);
      roccor_err = myCorr->GetMuonScaleSF(probe, MyCorrection::variation::up,
                                          matched_pt) -
                   roccor;
    }

    const float miniAODPt = probe.Pt();
    storage->miniAODPt[i] = miniAODPt;
    storage->momentumScaleUp[i] = miniAODPt * (roccor + roccor_err);
    storage->momentumScaleDown[i] = miniAODPt * (roccor - roccor_err);
    storage->correctedPt[i] = miniAODPt * roccor;
  }

  return MuonViewCollection(std::move(storage));
}

RVec<Muon> AnalyzerCore::GetAllMuons() {
  MuonViewCollection view = GetAllMuonViews();
  return MaterializeMuons(view);
}

RVec<Muon>
AnalyzerCore::MaterializeMuons(const MuonViewCollection &muons) const {
  RVec<Muon> materialised;
  if (muons.size() == 0)
    return materialised;

  auto storage = muons.storage();
  if (!storage)
    return materialised;

  materialised.reserve(muons.size());
  for (std::size_t idx = 0; idx < muons.size(); ++idx) {
    materialised.emplace_back(storage, idx);
  }
  return materialised;
}

RVec<Muon>
AnalyzerCore::MaterializeMuons(const MuonViewCollection &muons,
                               const std::vector<std::size_t> &indices) const {
  RVec<Muon> materialised;
  if (indices.empty() || muons.size() == 0)
    return materialised;

  auto storage = muons.storage();
  if (!storage)
    return materialised;

  materialised.reserve(indices.size());
  const std::size_t n = muons.size();
  for (const auto idx : indices) {
    if (idx >= n)
      continue;
    materialised.emplace_back(storage, idx);
  }
  return materialised;
}

std::vector<std::size_t>
AnalyzerCore::SelectMuonIndices(const MuonViewCollection &muons,
                                const std::vector<std::size_t> &seed_indices,
                                const Muon::MuonID ID, const float ptmin,
                                const float fetamax) const {
  const auto viewID = static_cast<MuonView::MuonID>(ID);
  std::vector<std::size_t> selected;
  selected.reserve(muons.size());
  for (auto i : seed_indices) {
    const auto &mu = muons[i];
    if (!(mu.Pt() > ptmin))
      continue;
    if (!(std::abs(mu.Eta()) < fetamax))
      continue;
    if (!mu.PassID(viewID))
      continue;
    selected.push_back(i);
  }
  return selected;
}

std::vector<std::size_t>
AnalyzerCore::SelectMuonIndices(const MuonViewCollection &muons,
                                const Muon::MuonID ID, const float ptmin,
                                const float fetamax) const {
  size_t n = muons.size();
  std::vector<std::size_t> seed_indices(n);
  ;
  std::iota(seed_indices.begin(), seed_indices.end(), 0);
  return SelectMuonIndices(muons, seed_indices, ID, ptmin, fetamax);
}

RVec<Muon> AnalyzerCore::ScaleMuons(const RVec<Muon> &muons,
                                    const MyCorrection::variation &syst) {
  RVec<Muon> scaled_muons;
  for (const auto &muon : muons) {
    Muon scaled_muon = muon;
    if (syst == MyCorrection::variation::up) {
      scaled_muon.SetPtEtaPhiM(muon.MomentumScaleUp(), muon.Eta(), muon.Phi(),
                               muon.M());
    } else if (syst == MyCorrection::variation::down) {
      scaled_muon.SetPtEtaPhiM(muon.MomentumScaleDown(), muon.Eta(), muon.Phi(),
                               muon.M());
    } else {
      throw runtime_error("[AnalyzerCore::ScaleMuons] Invalid variation");
    }
    scaled_muons.emplace_back(scaled_muon);
  }
  return scaled_muons;
}

RVec<Muon> AnalyzerCore::SelectMuons(const RVec<Muon> &muons, const TString ID,
                                     const float ptmin,
                                     const float fetamax) const {
  RVec<Muon> selected_muons;

  for (const auto &muon : muons) {
    if (!(muon.Pt() > ptmin))
      continue;
    if (!(fabs(muon.Eta()) < fetamax))
      continue;
    if (ID.Contains("&&")) {
      TObjArray *IDs = ID.Tokenize("&&");
      bool pass = true;
      for (int i = 0; i < IDs->GetEntries(); i++) {
        TString thisID = ((TObjString *)IDs->At(i))->GetString();
        if (!muon.PassID(thisID)) {
          pass = false;
          break;
        }
      }
      delete IDs;
      if (!pass)
        continue;
    } else {
      if (!muon.PassID(ID))
        continue;
    }
    selected_muons.push_back(muon);
  }
  return selected_muons;
}

RVec<Muon> AnalyzerCore::SelectMuons(const RVec<Muon> &muons,
                                     const Muon::MuonID ID, const float ptmin,
                                     const float fetamax) const {
  RVec<Muon> selected_muons;

  for (const auto &muon : muons) {
    if (!(muon.Pt() > ptmin))
      continue;
    if (!(fabs(muon.Eta()) < fetamax))
      continue;
    if (!muon.PassID(ID))
      continue;
    selected_muons.push_back(muon);
  }
  return selected_muons;
}

ElectronViewCollection AnalyzerCore::GetAllElectronViews() {
  auto storage = std::make_shared<ElectronSoA>();
  storage->pt.bind(&Electron_pt);
  storage->eta.bind(&Electron_eta);
  storage->phi.bind(&Electron_phi);
  storage->mass.bind(&Electron_mass);
  storage->charge.bind(&Electron_charge);
  storage->pfRelIso03.bind(&Electron_pfRelIso03_all);
  storage->miniPFRelIso.bind(&Electron_miniPFRelIso_all);
  storage->dxy.bind(&Electron_dxy);
  storage->dxyErr.bind(&Electron_dxyErr);
  storage->dz.bind(&Electron_dz);
  storage->dzErr.bind(&Electron_dzErr);
  storage->ip3d.bind(&Electron_ip3d);
  storage->sip3d.bind(&Electron_sip3d);
  storage->convVeto.bind(&Electron_convVeto);
  storage->lostHits.bind(&Electron_lostHits);
  storage->seedGain.bind(&Electron_seedGain);
  storage->tightCharge.bind(&Electron_tightCharge);
  storage->sieie.bind(&Electron_sieie);
  storage->hoe.bind(&Electron_hoe);
  storage->eInvMinusPInv.bind(&Electron_eInvMinusPInv);
  storage->dr03EcalRecHitSumEt.bind(&Electron_dr03EcalRecHitSumEt);
  storage->dr03HcalDepth1TowerSumEt.bind(&Electron_dr03HcalDepth1TowerSumEt);
  storage->dr03TkSumPt.bind(&Electron_dr03TkSumPt);
  storage->dr03TkSumPtHEEP.bind(&Electron_dr03TkSumPtHEEP);
  storage->r9.bind(&Electron_r9);
  storage->energyErr.bind(&Electron_energyErr);
  storage->cutBasedHEEP.bind(&Electron_cutBased_HEEP);
  storage->promptMVA.bind(&Electron_promptMVA);
  storage->mvaIsoWP80.bind(&Electron_mvaIso_WP80);
  storage->mvaIsoWP90.bind(&Electron_mvaIso_WP90);
  storage->mvaNoIsoWP80.bind(&Electron_mvaNoIso_WP80);
  storage->mvaNoIsoWP90.bind(&Electron_mvaNoIso_WP90);
  storage->mvaIso.bind(&Electron_mvaIso);
  storage->mvaNoIso.bind(&Electron_mvaNoIso);
  storage->cutBased.bind(&Electron_cutBased);
  storage->genPartFlav.bind(&Electron_genPartFlav);
  storage->genPartIdx.bind(&Electron_genPartIdx);
  storage->jetIdx.bind(&Electron_jetIdx);
  storage->scEta.bind(&Electron_superclusterEta);
  storage->deltaEtaSC.bind(&Electron_deltaEtaSC);

  const std::size_t n = storage->size();
  storage->rho.assign(n, Rho_fixedGridRhoFastjetAll);
  storage->dEsigmaUp.assign(n, -999.f);
  storage->dEsigmaDown.assign(n, -999.f);
  storage->ecalPFClusterIso.assign(n, -999.f);
  storage->hcalPFClusterIso.assign(n, -999.f);
  storage->deltaEtaSeed.assign(n, -999.f);
  storage->deltaPhiSC.assign(n, -999.f);
  storage->deltaPhiSeed.assign(n, -999.f);

  return ElectronViewCollection(std::move(storage));
}

RVec<Electron> AnalyzerCore::GetAllElectrons() {
  ElectronViewCollection view = GetAllElectronViews();
  return MaterializeElectrons(view);
}

RVec<Electron> AnalyzerCore::SelectElectrons(const RVec<Electron> &electrons,
                                             const TString ID,
                                             const float ptmin,
                                             const float fetamax,
                                             bool vetoHEM) const {
  RVec<Electron> selected_electrons;
  for (const auto &electron : electrons) {
    if (!(electron.Pt() > ptmin))
      continue;
    if (!(fabs(electron.Eta()) < fetamax))
      continue;
    if (!electron.PassID(ID))
      continue;
    if (vetoHEM && IsHEMElectron(electron))
      continue;
    selected_electrons.push_back(electron);
  }
  return selected_electrons;
}

RVec<Electron> AnalyzerCore::SelectElectrons(const RVec<Electron> &electrons,
                                             const Electron::ElectronID ID,
                                             const float ptmin,
                                             const float fetamax,
                                             bool vetoHEM) const {
  RVec<Electron> selected_electrons;
  for (const auto &electron : electrons) {
    if (!(electron.Pt() > ptmin))
      continue;
    if (!(fabs(electron.Eta()) < fetamax))
      continue;
    if (!electron.PassID(ID))
      continue;
    if (vetoHEM && IsHEMElectron(electron))
      continue;
    selected_electrons.push_back(electron);
  }
  return selected_electrons;
}

std::vector<std::size_t> AnalyzerCore::SelectElectronIndices(
    const ElectronViewCollection &electrons,
    const std::vector<size_t> &seed_indices, const Electron::ElectronID ID,
    const float ptmin, const float fetamax, bool vetoHEM) const {
  std::vector<std::size_t> selected;
  selected.reserve(electrons.size());
  const auto viewID = static_cast<ElectronView::ElectronID>(ID);
  for (auto i : seed_indices) {
    const auto &electron = electrons[i];
    if (!(electron.Pt() > ptmin))
      continue;
    if (!(std::abs(electron.Eta()) < fetamax))
      continue;
    bool passID = true;
    passID = electron.PassID(static_cast<ElectronView::ElectronID>(ID));
    if (!passID)
      continue;
    if (vetoHEM && IsHEMElectron(electron))
      continue;
    selected.push_back(i);
  }
  return selected;
}

std::vector<std::size_t> AnalyzerCore::SelectElectronIndices(
    const ElectronViewCollection &electrons, const Electron::ElectronID ID,
    const float ptmin, const float fetamax, bool vetoHEM) const {
  size_t n = electrons.size();
  std::vector<std::size_t> seed_indices(n);
  ;
  std::iota(seed_indices.begin(), seed_indices.end(), 0);
  return SelectElectronIndices(electrons, seed_indices, ID, ptmin, fetamax,
                               vetoHEM);
}

RVec<Electron> AnalyzerCore::MaterializeElectrons(
    const ElectronViewCollection &electrons) const {
  RVec<Electron> materialised;
  if (electrons.empty())
    return materialised;

  materialised.reserve(electrons.size());
  auto *ctx = const_cast<AnalyzerCore *>(this);
  for (std::size_t idx = 0; idx < electrons.size(); ++idx) {
    const auto &view = electrons[idx];
    Electron electron;
    const auto sourceIndex = view.rawIndex();
    electron.SetPtEtaPhiM(view.Pt(), view.Eta(), view.Phi(), view.M());
    electron.SetCharge(view.Charge());
    electron.SetScEta(view.ScEta());
    electron.SetDeltaEtaInSC(view.deltaEtaInSC());
    electron.SetDeltaPhiInSC(view.deltaPhiInSC());
    electron.SetDeltaEtaInSeed(view.deltaEtaInSeed());
    electron.SetDeltaPhiInSeed(view.deltaPhiInSeed());
    electron.SetPFClusterIso(view.ecalPFClusterIso(), view.hcalPFClusterIso());
    electron.SetRho(view.Rho());
    electron.AttachLazyPayload(ctx, &AnalyzerCore::ElectronEnsureThunk,
                               static_cast<int>(sourceIndex));

    float nominalPt = view.Pt();
    float correctedPt = nominalPt;
    float energyErr_corr = view.energyErr();

    if (Run == 3) {
      electron.SetEnergyResUnc(nominalPt, view.energyErr());

      float scEta = view.ScEta();
      if (scEta < -990.f) scEta = view.Eta() + view.deltaEtaInSC();
      electron.SetScEta(scEta);

      if (IsDATA) {
        float scale = myCorr->GetElectronScale(scEta, view.r9(), nominalPt, view.SeedGain(), RunNumber);
        correctedPt = nominalPt * scale;

        float smear_width = myCorr->GetElectronSmearWidth(electron, MyCorrection::variation::nom);
        float raw_energy = nominalPt * std::cosh(view.Eta());
        float energyErr_sq = view.energyErr() * view.energyErr();
        float smear_term = raw_energy * smear_width;
        energyErr_corr = std::sqrt(energyErr_sq + smear_term * smear_term) * scale;
      } else {
        float smearing = myCorr->GetElectronSmearUnc(electron, MyCorrection::variation::nom, int(Rho_fixedGridRhoFastjetAll));
        correctedPt = nominalPt * smearing;

        float smear_width = myCorr->GetElectronSmearWidth(electron, MyCorrection::variation::nom);
        float raw_energy = nominalPt * std::cosh(view.Eta());
        float energyErr_sq = view.energyErr() * view.energyErr();
        float smear_term = raw_energy * smear_width;
        energyErr_corr = std::sqrt(energyErr_sq + smear_term * smear_term) * smearing;
      }
      
      electron.SetPtEtaPhiM(correctedPt, view.Eta(), view.Phi(), view.M());
      electron.SetEnergyErr(energyErr_corr);
    } else {
      electron.SetEnergyErr(view.energyErr());
    }

    materialised.emplace_back(std::move(electron));
  }
  return materialised;
}

RVec<Electron> AnalyzerCore::MaterializeElectrons(
    const ElectronViewCollection &electrons,
    const std::vector<std::size_t> &indices) const {
  RVec<Electron> materialised;
  if (indices.empty() || electrons.empty())
    return materialised;

  materialised.reserve(indices.size());
  auto *ctx = const_cast<AnalyzerCore *>(this);
  const std::size_t n = electrons.size();
  for (const auto idx : indices) {
    if (idx >= n)
      continue;
    const auto &view = electrons[idx];
    Electron electron;
    const auto sourceIndex = view.rawIndex();
    electron.SetPtEtaPhiM(view.Pt(), view.Eta(), view.Phi(), view.M());
    electron.SetCharge(view.Charge());
    electron.SetScEta(view.ScEta());
    electron.SetDeltaEtaInSC(view.deltaEtaInSC());
    electron.SetDeltaPhiInSC(view.deltaPhiInSC());
    electron.SetDeltaEtaInSeed(view.deltaEtaInSeed());
    electron.SetDeltaPhiInSeed(view.deltaPhiInSeed());
    electron.SetPFClusterIso(view.ecalPFClusterIso(), view.hcalPFClusterIso());
    electron.SetRho(view.Rho());
    electron.AttachLazyPayload(ctx, &AnalyzerCore::ElectronEnsureThunk,
                               static_cast<int>(sourceIndex));

    float nominalPt = view.Pt();
    float correctedPt = nominalPt;
    float energyErr_corr = view.energyErr();

    if (Run == 3) {
      electron.SetEnergyResUnc(nominalPt, view.energyErr());

      float scEta = view.ScEta();
      if (scEta < -990.f) scEta = view.Eta() + view.deltaEtaInSC();
      electron.SetScEta(scEta);

      if (IsDATA) {
        float scale = myCorr->GetElectronScale(scEta, view.r9(), nominalPt, view.SeedGain(), RunNumber);
        correctedPt = nominalPt * scale;

        float smear_width = myCorr->GetElectronSmearWidth(electron, MyCorrection::variation::nom);
        float raw_energy = nominalPt * std::cosh(view.Eta());
        float energyErr_sq = view.energyErr() * view.energyErr();
        float smear_term = raw_energy * smear_width;
        energyErr_corr = std::sqrt(energyErr_sq + smear_term * smear_term) * scale;
      } else {
        float smearing = myCorr->GetElectronSmearUnc(electron, MyCorrection::variation::nom, int(Rho_fixedGridRhoFastjetAll));
        correctedPt = nominalPt * smearing;

        float smear_width = myCorr->GetElectronSmearWidth(electron, MyCorrection::variation::nom);
        float raw_energy = nominalPt * std::cosh(view.Eta());
        float energyErr_sq = view.energyErr() * view.energyErr();
        float smear_term = raw_energy * smear_width;
        energyErr_corr = std::sqrt(energyErr_sq + smear_term * smear_term) * smearing;
      }
      
      electron.SetPtEtaPhiM(correctedPt, view.Eta(), view.Phi(), view.M());
      electron.SetEnergyErr(energyErr_corr);
    } else {
      electron.SetEnergyErr(view.energyErr());
    }

    materialised.emplace_back(std::move(electron));
  }
  return materialised;
}

RVec<Electron>
AnalyzerCore::ScaleElectrons(const Event &ev, const RVec<Electron> &electrons,
                             const MyCorrection::variation &syst) {
  if (IsDATA || syst == MyCorrection::variation::nom)
    return electrons;

  RVec<Electron> scaled_electrons;
  for (const auto &electron : electrons) {
    float scale_variation = 1.;
    float raw_pt = (Run == 3 && electron.dEsigmaUp() > 0.0) ? electron.dEsigmaUp() : electron.Pt();
    float scEta = electron.scEta();
    if (scEta < -990.f) scEta = electron.Eta() + electron.deltaEtaInSC();

    if (syst == MyCorrection::variation::up) {
      scale_variation = myCorr->GetElectronScaleUnc(
          scEta, electron.SeedGain(), ev.run(), electron.r9(),
          raw_pt, MyCorrection::variation::up);
    } else if (syst == MyCorrection::variation::down) {
      scale_variation = myCorr->GetElectronScaleUnc(
          scEta, electron.SeedGain(), ev.run(), electron.r9(),
          raw_pt, MyCorrection::variation::down);
    } else {
      throw runtime_error("[AnalyzerCore::ScaleElectrons] Invalid variation");
    }
    Electron scaled_electron = electron;
    scaled_electron.SetPtEtaPhiM(electron.Pt() * scale_variation,
                                 electron.Eta(), electron.Phi(), electron.M());
    if (Run == 3) {
      scaled_electron.SetEnergyErr(electron.energyErr() * scale_variation);
    }
    scaled_electrons.emplace_back(scaled_electron);
  }
  return scaled_electrons;
}

RVec<Electron>
AnalyzerCore::SmearElectrons(const RVec<Electron> &electrons,
                             const MyCorrection::variation &syst) {
  RVec<Electron> smeared_electrons;
  switch (Run) {
  case 2:
    for (const auto &electron : electrons) {
      float smeared_pt = electron.Pt();
      if (syst == MyCorrection::variation::up) {
        smeared_pt *= (electron.E() - electron.dEsigmaUp()) / electron.E();
      } else if (syst == MyCorrection::variation::down) {
        smeared_pt *= (electron.E() - electron.dEsigmaDown()) / electron.E();
      } else {
        throw runtime_error("[AnalyzerCore::SmearElectrons] Invalid variation");
      }
      Electron smeared_electron = electron;
      smeared_electron.SetPtEtaPhiM(smeared_pt, electron.Eta(), electron.Phi(),
                                    electron.M());
      smeared_electrons.emplace_back(smeared_electron);
    }
    break;
  case 3:
    for (const auto &electron : electrons) {
      float raw_pt = (electron.dEsigmaUp() > 0.0) ? electron.dEsigmaUp() : electron.Pt();
      float raw_energyErr = (electron.dEsigmaDown() > 0.0) ? electron.dEsigmaDown() : electron.energyErr();
      float raw_energy = raw_pt * std::cosh(electron.Eta());
      
      float smear_width = myCorr->GetElectronSmearWidth(electron, syst);
      float smearing = myCorr->GetElectronSmearUnc(electron, syst, int(Rho_fixedGridRhoFastjetAll));
      
      float smeared_pt = raw_pt * smearing;
      float energyErr_sq = raw_energyErr * raw_energyErr;
      float smear_term = raw_energy * smear_width;
      float smeared_energyErr = std::sqrt(energyErr_sq + smear_term * smear_term) * smearing;

      Electron smeared_electron = electron;
      smeared_electron.SetPtEtaPhiM(smeared_pt, electron.Eta(), electron.Phi(),
                                    electron.M());
      smeared_electron.SetEnergyErr(smeared_energyErr);
      smeared_electrons.emplace_back(smeared_electron);
    }
    break;
  default:
    throw runtime_error("[AnalyzerCore::SmearElectrons] Run " + to_string(Run) +
                        " is not supported");
  }
  return smeared_electrons;
}

RVec<GenJet>
AnalyzerCore::MaterializeGenJets(const GenJetViewCollection &genjets) const {
  RVec<GenJet> materialised;
  if (genjets.size() == 0)
    return materialised;

  materialised.reserve(genjets.size());
  for (std::size_t idx = 0; idx < genjets.size(); ++idx) {
    const auto &view = genjets[idx];
    GenJet jet;
    jet.SetPtEtaPhiM(view.Pt(), view.Eta(), view.Phi(), view.Mass());
    jet.SetGenFlavours(view.PartonFlavour(), view.HadronFlavour());
    materialised.emplace_back(jet);
  }
  return materialised;
}

RVec<Gen> AnalyzerCore::GetAllGens() {
  RVec<Gen> gens;
  if (IsDATA)
    return gens;

  GenViewCollection view = GetAllGenViews();
  gens.reserve(view.size());
  auto storage = view.storage();
  for (std::size_t i = 0; i < view.size(); ++i) {
    gens.emplace_back(storage, i);
  }

  return gens;
}

RVec<LHE> AnalyzerCore::GetAllLHEs() {
  RVec<LHE> lhes;
  if (IsDATA)
    return lhes;

  // Safety check for array sizes
  if (nLHEPart <= 0)
    return lhes;

  for (int i = 0; i < nLHEPart; i++) {
    LHE lhe;
    lhe.SetIsEmpty(false);
    lhe.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i],
                     LHEPart_mass[i]);
    lhe.SetStatus(LHEPart_status[i]);
    lhe.SetSpin(LHEPart_spin[i]);
    lhe.SetIncomingPz(LHEPart_incomingpz[i]);
    lhe.SetPdgId(LHEPart_pdgId[i]);
    lhes.push_back(lhe);
  }
  return lhes;
}

RVec<Tau> AnalyzerCore::GetAllTaus() {

  RVec<Tau> taus;

  for (int i = 0; i < nTau; i++) {

    Tau tau;
    tau.SetPtEtaPhiM(Tau_pt[i], Tau_eta[i], Tau_phi[i], Tau_mass[i]);
    tau.SetdXY(Tau_dxy[i]);
    tau.SetdZ(Tau_dz[i]);
    tau.SetGenPartFlav(Tau_genPartFlav[i]);
    tau.SetCharge(Tau_charge[i]);
    tau.SetDecayMode(Tau_decayMode[i]);
    tau.SetIdDeepTau2018v2p5VSjet(Tau_idDeepTau2018v2p5VSjet[i]);
    tau.SetIdDeepTau2018v2p5VSmu(Tau_idDeepTau2018v2p5VSmu[i]);
    tau.SetIdDeepTau2018v2p5VSe(Tau_idDeepTau2018v2p5VSe[i]);
    if (!IsDATA)
      tau.SetGenPartIdx(Tau_genPartIdx[i]);

    tau.SetIdDecayModeNewDMs(Tau_idDecayModeNewDMs[i]);
    taus.push_back(tau);
  }

  return taus;
}

RVec<Tau> AnalyzerCore::SelectTaus(const RVec<Tau> &taus, const TString ID,
                                   const float ptmin,
                                   const float absetamax) const {

  RVec<Tau> selected_taus;

  for (const auto &tau : taus) {

    if (!(tau.Pt() > ptmin))
      continue;
    if (!(fabs(tau.Eta()) < absetamax))
      continue;
    if (!tau.PassID(ID))
      continue;
    selected_taus.push_back(tau);
  }

  return selected_taus;
}

std::shared_ptr<JetSoA> AnalyzerCore::CreateJetSoA() const {
  auto storage = std::make_shared<JetSoA>();
  storage->pt.bind(&Jet_pt);
  storage->eta.bind(&Jet_eta);
  storage->phi.bind(&Jet_phi);
  storage->mass.bind(&Jet_mass);
  storage->rawFactor.bind(&Jet_rawFactor);
  storage->area.bind(&Jet_area);
  storage->chHEF.bind(&Jet_chHEF);
  storage->neHEF.bind(&Jet_neHEF);
  storage->neEmEF.bind(&Jet_neEmEF);
  storage->chEmEF.bind(&Jet_chEmEF);
  storage->muEF.bind(&Jet_muEF);
  storage->partonFlavour.bind(&Jet_partonFlavour);
  storage->hadronFlavour.bind(&Jet_hadronFlavour);
  storage->chMultiplicity.bind(&Jet_chMultiplicity);
  storage->neMultiplicity.bind(&Jet_neMultiplicity);
  storage->nConstituents.bind(&Jet_nConstituents);
  storage->nElectrons.bind(&Jet_nElectrons);
  storage->nMuons.bind(&Jet_nMuons);
  storage->nSVs.bind(&Jet_nSVs);
  storage->electronIdx1.bind(&Jet_electronIdx1);
  storage->electronIdx2.bind(&Jet_electronIdx2);
  storage->muonIdx1.bind(&Jet_muonIdx1);
  storage->muonIdx2.bind(&Jet_muonIdx2);
  storage->svIdx1.bind(&Jet_svIdx1);
  storage->svIdx2.bind(&Jet_svIdx2);
  storage->genJetIdx.bind(&Jet_genJetIdx);

  storage->deepFlavB.bind(&Jet_btagDeepFlavB);
  storage->deepFlavCvB.bind(&Jet_btagDeepFlavCvB);
  storage->deepFlavCvL.bind(&Jet_btagDeepFlavCvL);
  storage->deepFlavQG.bind(&Jet_btagDeepFlavQG);
  storage->pnetB.bind(&Jet_btagPNetB);
  storage->pnetCvB.bind(&Jet_btagPNetCvB);
  storage->pnetCvL.bind(&Jet_btagPNetCvL);
  storage->pnetCvNotB.bind(&Jet_btagPNetCvNotB);
  storage->pnetQvG.bind(&Jet_btagPNetQvG);
  storage->pnetTauVJet.bind(&Jet_btagPNetTauVJet);
  storage->uparTAK4B.bind(&Jet_btagUParTAK4B);
  storage->uparTAK4CvB.bind(&Jet_btagUParTAK4CvB);
  storage->uparTAK4CvL.bind(&Jet_btagUParTAK4CvL);
  storage->uparTAK4CvNotB.bind(&Jet_btagUParTAK4CvNotB);
  storage->uparTAK4Ele.bind(&Jet_btagUParTAK4Ele);
  storage->uparTAK4Mu.bind(&Jet_btagUParTAK4Mu);
  storage->uparTAK4QvG.bind(&Jet_btagUParTAK4QvG);
  storage->uparTAK4SvCB.bind(&Jet_btagUParTAK4SvCB);
  storage->uparTAK4SvUDG.bind(&Jet_btagUParTAK4SvUDG);
  storage->uparTAK4TauVJet.bind(&Jet_btagUParTAK4TauVJet);
  storage->uparTAK4UDG.bind(&Jet_btagUParTAK4UDG);
  storage->uparTAK4ProbB.bind(&Jet_btagUParTAK4probb);
  storage->uparTAK4ProbBB.bind(&Jet_btagUParTAK4probbb);
  storage->pnetRegPtRawCorr.bind(&Jet_PNetRegPtRawCorr);
  storage->pnetRegPtRawCorrNeutrino.bind(&Jet_PNetRegPtRawCorrNeutrino);
  storage->pnetRegPtRawRes.bind(&Jet_PNetRegPtRawRes);
  storage->uparTAK4RegPtRawCorr.bind(&Jet_UParTAK4RegPtRawCorr);
  storage->uparTAK4RegPtRawCorrNeutrino.bind(&Jet_UParTAK4RegPtRawCorrNeutrino);
  storage->uparTAK4RegPtRawRes.bind(&Jet_UParTAK4RegPtRawRes);
  storage->uparTAK4V1RegPtRawCorr.bind(&Jet_UParTAK4V1RegPtRawCorr);
  storage->uparTAK4V1RegPtRawCorrNeutrino.bind(
      &Jet_UParTAK4V1RegPtRawCorrNeutrino);
  storage->uparTAK4V1RegPtRawRes.bind(&Jet_UParTAK4V1RegPtRawRes);
  storage->puIdDisc.bind(&Jet_puIdDisc);
  return storage;
}

void AnalyzerCore::InitialiseJetSystematics(JetSoA &storage) const {
  const std::size_t n = storage.size();
  storage.jecFactor.assign(n, 1.f);
  storage.correctedPt.assign(n, 0.f);
  storage.correctedMass.assign(n, 0.f);
  storage.smearedPtNominal.assign(n, 0.f);
  storage.smearedPtUp.assign(n, 0.f);
  storage.smearedPtDown.assign(n, 0.f);
  storage.smearedMassNominal.assign(n, 0.f);
  storage.smearedMassUp.assign(n, 0.f);
  storage.smearedMassDown.assign(n, 0.f);
  storage.jesPtUp.assign(n, 0.f);
  storage.jesPtDown.assign(n, 0.f);
  storage.jesMassUp.assign(n, 0.f);
  storage.jesMassDown.assign(n, 0.f);
}

void AnalyzerCore::PopulateJetStorageWithoutCorrections(JetSoA &storage) const {
  const std::size_t n = storage.size();
  for (std::size_t i = 0; i < n; ++i) {
    const float rawPt = storage.pt[i] * (1.f - storage.rawFactor[i]);
    const float rawMass = storage.mass[i] * (1.f - storage.rawFactor[i]);
    storage.jecFactor[i] = 1.f;
    storage.correctedPt[i] = -999.f;
    storage.correctedMass[i] = -999.f;
    storage.smearedPtNominal[i] = -999.f;
    storage.smearedPtUp[i] = -999.f;
    storage.smearedPtDown[i] = -999.f;
    storage.smearedMassNominal[i] = -999.f;
    storage.smearedMassUp[i] = -999.f;
    storage.smearedMassDown[i] = -999.f;
    storage.jesPtUp[i] = -999.f;
    storage.jesPtDown[i] = -999.f;
    storage.jesMassUp[i] = -999.f;
    storage.jesMassDown[i] = -999.f;
  }
}

void AnalyzerCore::ApplyJetEnergyCorrections(JetSoA &storage, float rho) {
  const std::size_t n = storage.size();
  for (std::size_t i = 0; i < n; ++i) {
    const float rawPt = storage.pt[i] * (1.f - storage.rawFactor[i]);
    const float rawMass = storage.mass[i] * (1.f - storage.rawFactor[i]);
    const float jecSF =
        myCorr ? myCorr->GetJESSF(storage.area[i], storage.eta[i], rawPt,
                                  storage.phi[i], rho, RunNumber)
               : 1.f;
    storage.jecFactor[i] = jecSF;
    storage.correctedPt[i] = rawPt * jecSF;
    storage.correctedMass[i] = rawMass * jecSF;
  }
}

std::vector<int>
AnalyzerCore::MatchJetsToGenJets(const JetViewCollection &jets,
                                 const GenJetViewCollection &genjets,
                                 float rho) const {
  const auto &jetSoAPtr = jets.storage();
  const auto &genSoAPtr = genjets.storage();
  const std::size_t njet = jetSoAPtr ? jetSoAPtr->size() : 0;
  const std::size_t ngen = genSoAPtr ? genSoAPtr->size() : 0;

  std::vector<int> matchedGenIdx(njet, -999);
  if (njet == 0 || ngen == 0)
    return matchedGenIdx;

  const JetSoA &jsoa = *jetSoAPtr;
  const GenJetSoA &gsoa = *genSoAPtr;

  std::vector<std::tuple<std::size_t, std::size_t, float, float>> candidates;
  candidates.reserve(njet * 2);

  const float maxDR = 0.2f;
  const float maxDR2 = maxDR * maxDR;
  const float ptJerCut = 3.f;
  constexpr float pi = 3.14159265358979323846f;
  constexpr float twoPi = 6.28318530717958647692f;

  std::vector<float> jetJerSigma(njet, 0.f);
  for (std::size_t i = 0; i < njet; ++i) {
    const float pt =
        jsoa.correctedPt.empty() ? jsoa.pt[i] : jsoa.correctedPt[i];
    if (pt > 0.f)
      jetJerSigma[i] = myCorr->GetJER(jsoa.eta[i], pt, rho) * pt;
  }

  // 후보 수집
  for (std::size_t i = 0; i < njet; ++i) {
    const float jetPt =
        jsoa.correctedPt.empty() ? jsoa.pt[i] : jsoa.correctedPt[i];
    const float jetEta = jsoa.eta[i];
    const float jetPhi = jsoa.phi[i];
    const float sigma = jetJerSigma[i];
    if (jetPt <= 0.f || sigma <= 0.f)
      continue;

    const float maxPtDiff = ptJerCut * sigma;

    for (std::size_t j = 0; j < ngen; ++j) {
      const float genEta = gsoa.eta[j];
      const float genPhi = gsoa.phi[j];
      const float genPt = gsoa.pt[j];

      const float dEta = jetEta - genEta;
      float dPhi = jetPhi - genPhi;
      if (dPhi > pi)
        dPhi -= twoPi;
      else if (dPhi <= -pi)
        dPhi += twoPi;

      const float dr2 = dEta * dEta + dPhi * dPhi;
      if (dr2 >= maxDR2)
        continue;

      const float ptDiff = std::fabs(jetPt - genPt);
      if (ptDiff < maxPtDiff)
        candidates.emplace_back(i, j, dr2, ptDiff);
    }
  }

  // 가까운 순서로 정렬
  std::sort(candidates.begin(), candidates.end(),
            [](const auto &a, const auto &b) {
              if (std::get<2>(a) == std::get<2>(b))
                return std::get<3>(a) < std::get<3>(b);
              return std::get<2>(a) < std::get<2>(b);
            });

  // 1:1 할당
  std::vector<bool> usedJet(njet, false);
  std::vector<bool> usedGen(ngen, false);
  for (const auto &c : candidates) {
    const auto i = std::get<0>(c);
    const auto j = std::get<1>(c);
    if (usedJet[i] || usedGen[j])
      continue;
    matchedGenIdx[i] = static_cast<int>(j);
    usedJet[i] = true;
    usedGen[j] = true;
  }

  return matchedGenIdx;
}

void AnalyzerCore::ApplyJetSmearingAndUncertainties(
    const JetViewCollection &jets, const std::vector<int> &matchedGenIdx,
    const GenJetViewCollection &genjets, bool isMC, float rho) {
  if (IsDATA) {
    auto storage = jets.storage();
    if (storage) {
      for (size_t i = 0; i < storage->size(); ++i) {
        storage->smearedPtNominal[i] = storage->correctedPt[i];
        storage->smearedMassNominal[i] = storage->correctedMass[i];
      }
    }
    return;
  }
  const auto &jesSources = JetEnergyScaleSources();
  constexpr float MIN_JET_ENERGY = 1e-2f;

  auto &storage = *jets.storage();
  const std::size_t n = storage.size();
  for (std::size_t i = 0; i < n; ++i) {
    const float correctedPt = storage.correctedPt[i];
    const float correctedMass = storage.correctedMass[i];
    const float eta = storage.eta[i];
    const float phi = storage.phi[i];

    TLorentzVector unsmearedP4;
    unsmearedP4.SetPtEtaPhiM(correctedPt, eta, phi, correctedMass);
    const float minCorr =
        MIN_JET_ENERGY /
        std::max(static_cast<float>(unsmearedP4.E()), MIN_JET_ENERGY);

    float smearNom = 1.f;
    float smearUp = 1.f;
    float smearDown = 1.f;

    if (isMC && correctedPt > 0.f) {
      const int matched = (!matchedGenIdx.empty() && i < matchedGenIdx.size())
                              ? matchedGenIdx[i]
                              : -1;
      const float sfNom =
          myCorr->GetJERSF(eta, correctedPt, MyCorrection::variation::nom);
      const float sfUp =
          myCorr->GetJERSF(eta, correctedPt, MyCorrection::variation::up);
      const float sfDown =
          myCorr->GetJERSF(eta, correctedPt, MyCorrection::variation::down);

      const float jerSigma =
          (matched < 0) ? myCorr->GetJER(eta, correctedPt, rho) : 0.f;
      const float gausSample =
          (matched < 0) ? gRandom->Gaus(0., jerSigma) : 0.f;

      const auto computeCorr = [&](float sf) {
        float corr = 1.f;
        if (matched >= 0 && matched < genjets.size()) {
          const float genPt = genjets[matched].Pt();
          const float scale = 1.f - genPt / std::max(correctedPt, 1e-6f);
          corr += (sf - 1.f) * scale;
        } else {
          const float variance = std::max(sf * sf - 1.f, 0.f);
          corr += gausSample * std::sqrt(variance);
        }
        return std::max(corr, minCorr);
      };

      smearNom = computeCorr(sfNom);
      smearUp = computeCorr(sfUp);
      smearDown = computeCorr(sfDown);
    }

    storage.smearedPtNominal[i] = correctedPt * smearNom;
    storage.smearedPtUp[i] = correctedPt * smearUp;
    storage.smearedPtDown[i] = correctedPt * smearDown;
    storage.smearedMassNominal[i] = correctedMass * smearNom;
    storage.smearedMassUp[i] = correctedMass * smearUp;
    storage.smearedMassDown[i] = correctedMass * smearDown;
  }
}

void AnalyzerCore::ApplyJetScaleVariation(JetViewCollection &jets,
                                          const TString &source) const {
  if (!myCorr || IsDATA)
    return;

  auto storagePtr = jets.storage();
  if (!storagePtr)
    return;

  JetSoA &storage = *storagePtr;
  const std::size_t n = storage.size();
  if (n == 0)
    return;

  const auto &jesSources = JetEnergyScaleSources();
  const bool useTotalSource =
      source.IsNull() || source.EqualTo("total", TString::kIgnoreCase);

  for (std::size_t i = 0; i < n; ++i) {
    const float baseSmearedPt = storage.smearedPtNominal[i];
    float unc = 0.f;

    if (useTotalSource) {
      for (const auto &src : jesSources) {
        const float srcUnc = myCorr->GetJESUncertainty(
            storage.eta[i], storage.correctedPt[i], src);
        unc += (srcUnc) * (srcUnc);
      }
      unc = std::sqrt(unc);
    } else {
      unc = myCorr->GetJESUncertainty(storage.eta[i], storage.correctedPt[i],
                                      source);
    }

    // storage.jecFactor[i]; untouched
    // storage.correctedPt[i]; untouched
    /// storage.correctedMass[i]; untouched
    // storage.smearedPtNominal[i]  = -1.f; // to indicate that nominal has
    // changed storage.smearedPtUp[i] = -1.f; storage.smearedPtDown[i] = -1.f;
    // storage.smearedMassNominal[i] = -1.f;
    // storage.smearedMassUp[i] = -1.f;
    // storage.smearedMassDown[i] = -1.f;

    storage.jesPtUp[i] = storage.smearedPtNominal[i] * (1.f + unc);
    storage.jesPtDown[i] = storage.smearedPtNominal[i] * (1.f - unc);
    storage.jesMassUp[i] = storage.smearedMassNominal[i] * (1.f + unc);
    storage.jesMassDown[i] = storage.smearedMassNominal[i] * (1.f - unc);
  }
}

void AnalyzerCore::ApplyJetSmearVariation(JetViewCollection &jets,
                                          const RVec<GenJet> &genjets,
                                          const TString &source) const {
  throw runtime_error(
      "AnalyzerCore::ApplyJetSmearVariation not implemented"); // not need
                                                               // untill someone
                                                               // needs JER
                                                               // breakdown
}

const std::vector<TString> &AnalyzerCore::JetEnergyScaleSources() const {
  static const std::vector<TString> sources = {
      "AbsoluteMPFBias", "AbsoluteScale",   "AbsoluteStat",   "FlavorQCD",
      "Fragmentation",   "PileUpDataMC",    "PileUpPtBB",     "PileUpPtEC1",
      "PileUpPtEC2",     "PileUpPtHF",      "PileUpPtRef",    "RelativeJEREC1",
      "RelativeJEREC2",  "RelativeJERHF",   "RelativePtBB",   "RelativePtEC1",
      "RelativePtEC2",   "RelativePtHF",    "RelativeBal",    "RelativeSample",
      "RelativeFSR",     "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF",
      "SinglePionECAL",  "SinglePionHCAL",  "TimePtEta"};
  ;
  return sources;
}

JetViewCollection AnalyzerCore::GetAllJetViews() {
  auto storage = CreateJetSoA();
  InitialiseJetSystematics(*storage);

  const std::size_t n = storage->size();
  if (n == 0)
    return JetViewCollection(std::move(storage));

  if (!myCorr) {
    PopulateJetStorageWithoutCorrections(*storage);
    return JetViewCollection(std::move(storage));
  }

  const bool isMC = !IsDATA;
  const float rho = Rho_fixedGridRhoFastjetAll;
  ApplyJetEnergyCorrections(*storage, rho);

  JetViewCollection jetViews = JetViewCollection(std::move(storage));
  SmearJetViews(jetViews, rho);

  return jetViews;
}

void AnalyzerCore::SmearJetViews(const JetViewCollection &jets,
                                 const float rho) {

  GenJetViewCollection genjetviews = GetAllGenJetViews();
  std::vector<int> matchedGenIdx;
  matchedGenIdx = MatchJetsToGenJets(jets, genjetviews, rho);
  gRandom->SetSeed(static_cast<int>(PuppiMET_pt * 1e6));

  ApplyJetSmearingAndUncertainties(jets, matchedGenIdx, genjetviews, true, rho);
}

JetViewCollection
AnalyzerCore::CloneJetViews(const JetViewCollection &jets) const {
  const auto &storage = jets.storage();
  if (!storage)
    return JetViewCollection();
  auto clone = std::make_shared<JetSoA>(*storage);
  return JetViewCollection(std::move(clone));
}

Jet AnalyzerCore::MaterializeJet(const JetViewCollection &jets,
                                 std::size_t index,
                                 const MyCorrection::variation &jesVar,
                                 const MyCorrection::variation &jerVar) const {
  if (jesVar != MyCorrection::variation::nom &&
      jerVar != MyCorrection::variation::nom) {
    throw runtime_error(
        "[AnalyzerCore::MaterializeJet] Both JES and JER variations are "
        "non-nominal, not supported");
  }

  Jet jet;
  auto storage = jets.storage();
  const std::size_t nStorage = storage->size();
  const std::size_t nViews = jets.size();
  if (index >= nStorage || index >= nViews)
    throw runtime_error(
        "[AnalyzerCore::MaterializeJet] Invalid storage or index");

  const auto &jetView = jets[index];
  if (!jetView.valid())
    return jet;

  jet = Jet(storage, index);

  const float rawPt = storage->pt[index] * (1.f - storage->rawFactor[index]);
  const float rawMass =
      storage->mass[index] * (1.f - storage->rawFactor[index]);
  const float correctedPt =
      !storage->correctedPt.empty() ? storage->correctedPt[index] : rawPt;
  const float correctedMass =
      !storage->correctedMass.empty() ? storage->correctedMass[index] : rawMass;
  const float smearedPt = storage->smearedPtNominal[index];
  const float smearedMass = storage->smearedMassNominal[index];

  float selectedPt = smearedPt;
  float selectedMass = smearedMass;
  if (!IsDATA) {
    if (jesVar == MyCorrection::variation::up) {
      selectedPt = storage->jesPtUp[index];
      selectedMass = storage->jesMassUp[index];
    } else if (jesVar == MyCorrection::variation::down) {
      selectedPt = storage->jesPtDown[index];
      selectedMass = storage->jesMassDown[index];
    } else if (jerVar == MyCorrection::variation::up) {
      selectedPt = storage->smearedPtUp[index];
      selectedMass = storage->smearedMassUp[index];
    } else if (jerVar == MyCorrection::variation::down) {
      selectedPt = storage->smearedPtDown[index];
      selectedMass = storage->smearedMassDown[index];
    }
  }

  jet.SetPtEtaPhiM(IsDATA ? correctedPt : selectedPt, jetView.Eta(),
                   jetView.Phi(), IsDATA ? correctedMass : selectedMass);
  jet.SetRawPt(rawPt);
  jet.SetOriginalPt(jetView.Pt());
  jet.SetArea(jetView.Area());
  jet.SetOriginalIndex(static_cast<int>(index));
  jet.SetEnergyFractions(jetView.ChHEF(), jetView.NeHEF(), jetView.NeEmEF(),
                         jetView.ChEmEF(), jetView.MuEF());
  if (!IsDATA) {
    jet.SetJetFlavours(jetView.PartonFlavour(), jetView.HadronFlavour());
  }

  jet.SetHadronMultiplicities(jetView.ChMultiplicity(),
                              jetView.NeMultiplicity());
  jet.SetMultiplicities(jetView.NConstituents(), jetView.NElectrons(),
                        jetView.NMuons(), jetView.NSVs());
  if (!IsDATA) {
    jet.SetMatchingIndices(jetView.ElectronIdx1(), jetView.ElectronIdx2(),
                           jetView.MuonIdx1(), jetView.MuonIdx2(),
                           jetView.SvIdx1(), jetView.SvIdx2(),
                           jetView.GenJetIdx());
  } else {
    jet.SetMatchingIndices(jetView.ElectronIdx1(), jetView.ElectronIdx2(),
                           jetView.MuonIdx1(), jetView.MuonIdx2(),
                           jetView.SvIdx1(), jetView.SvIdx2());
  }

  RVec<float> corrs = {jetView.PNetRegPtRawCorr(),
                       jetView.PNetRegPtRawCorrNeutrino(),
                       jetView.PNetRegPtRawRes(),
                       jetView.UParTAK4RegPtRawCorr(),
                       jetView.UParTAK4RegPtRawCorrNeutrino(),
                       jetView.UParTAK4RegPtRawRes(),
                       jetView.UParTAK4V1RegPtRawCorr(),
                       jetView.UParTAK4V1RegPtRawCorrNeutrino(),
                       jetView.UParTAK4V1RegPtRawRes(),
                       jetView.RawFactor()};

  jet.SetJetPuIDScore(jetView.PuIdDisc());
  jet.SetCorrections(corrs);

  // const int nPFCand = static_cast<int>(nJetPFCand);
  // for (int idx = 0; idx < nPFCand; ++idx) {
  //   if (JetPFCand_jetIdx[idx] == static_cast<int>(index)) {
  //     JetConstituent constituent;
  //     const auto pfCandTableidx =
  //         static_cast<std::size_t>(JetPFCand_pfCandIdx[idx]);
  //     if (pfCandTableidx >= PFCand_pt.size() ||
  //         pfCandTableidx >= PFCand_eta.size() ||
  //         pfCandTableidx >= PFCand_phi.size() ||
  //         pfCandTableidx >= PFCand_mass.size()) {
  //       continue; // skip malformed constituent entries
  //     }
  //     constituent.SetPtEtaPhiM(
  //         PFCand_pt[pfCandTableidx], PFCand_eta[pfCandTableidx],
  //         PFCand_phi[pfCandTableidx], PFCand_mass[pfCandTableidx]);
  //     constituent.SetPdgId(PFCand_pdgId[pfCandTableidx]);
  //     if (PFCand_puppiWeight.valid() &&
  //         pfCandTableidx < PFCand_puppiWeight.size())
  //       constituent.SetPUPPIWeight(PFCand_puppiWeight[pfCandTableidx]);
  //   }
  // }

  Jet unsmeared;
  unsmeared.SetPtEtaPhiM(correctedPt, jetView.Eta(), jetView.Phi(),
                         correctedMass);
  unsmeared.SetArea(jetView.Area());
  unsmeared.SetRawPt(rawPt);
  unsmeared.SetOriginalPt(jetView.Pt());
  jet.SetUnsmearedP4(unsmeared);

  return jet;
}

namespace {
inline void ApplyJetVariationScale(const JetView &view,
                                   const MyCorrection::variation &jesVar,
                                   const MyCorrection::variation &jerVar,
                                   Jet &jet) {
  if (jesVar != MyCorrection::variation::nom &&
      jerVar != MyCorrection::variation::nom) {
    throw runtime_error(
        "[AnalyzerCore::MaterializeJets] Both JES and JER variations are "
        "non-nominal, not supported");
  }

  const float nominalPt = view.SmearedPtNominal();
  if (nominalPt <= 0.f)
    return;

  if (jesVar == MyCorrection::variation::up)
    jet *= view.JesPtUp() / nominalPt;
  else if (jesVar == MyCorrection::variation::down)
    jet *= view.JesPtDown() / nominalPt;
  else if (jerVar == MyCorrection::variation::up)
    jet *= view.SmearedPtUp() / nominalPt;
  else if (jerVar == MyCorrection::variation::down)
    jet *= view.SmearedPtDown() / nominalPt;
}
} // namespace

RVec<Jet>
AnalyzerCore::MaterializeJets(const JetViewCollection &jets,
                              const MyCorrection::variation &jesVar,
                              const MyCorrection::variation &jerVar) const {
  RVec<Jet> materialised;
  if (jets.empty())
    return materialised;

  materialised.reserve(jets.size());
  for (std::size_t idx = 0; idx < jets.size(); ++idx) {
    const auto &view = jets[idx];
    if (!view.valid())
      continue;
    materialised.emplace_back(MaterializeJet(jets, idx, jesVar, jerVar));
  }
  return materialised;
}

RVec<Jet>
AnalyzerCore::MaterializeJets(const JetViewCollection &jets,
                              const std::vector<std::size_t> &indices,
                              const MyCorrection::variation &jesVar,
                              const MyCorrection::variation &jerVar) const {
  RVec<Jet> materialised;
  if (indices.empty() || jets.empty())
    return materialised;

  materialised.reserve(indices.size());
  const std::size_t n = jets.size();
  for (const auto idx : indices) {
    if (idx >= n)
      throw runtime_error(
          "[AnalyzerCore::MaterializeJets] Index out of range: " +
          to_string(idx) + " >= " + to_string(n));
    const auto &view = jets[idx];
    if (!view.valid())
      continue;
    materialised.emplace_back(MaterializeJet(jets, idx, jesVar, jerVar));
  }
  return materialised;
}

RVec<Jet> AnalyzerCore::GetAllJets() {
  JetViewCollection view = GetAllJetViews();
  return MaterializeJets(view);
}

void AnalyzerCore::MuonEnsureThunk(void *ctx, Muon &muon,
                                   Muon::Property property) {
  if (!ctx)
    return;
  static_cast<AnalyzerCore *>(ctx)->EnsureMuonProperty(muon, property);
}

void AnalyzerCore::ElectronEnsureThunk(void *ctx, Electron &electron,
                                       Electron::Property property) {
  if (!ctx)
    return;
  static_cast<AnalyzerCore *>(ctx)->EnsureElectronProperty(electron, property);
}

void AnalyzerCore::JetEnsureThunk(void *ctx, Jet &jet, Jet::Property property) {
  if (!ctx)
    return;
  static_cast<AnalyzerCore *>(ctx)->EnsureJetProperty(jet, property);
}

void AnalyzerCore::EnsureMuonProperty(Muon &muon,
                                      Muon::Property property) const {
  if (!muon.lazy_)
    return;
  const int idx = muon.lazy_->index;
  switch (property) {
  case Muon::Property::TkRelIso:
    muon.SetTkRelIso(Muon_tkRelIso[idx]);
    break;
  case Muon::Property::PfRelIso03:
    muon.SetPfRelIso03(Muon_pfRelIso03_all[idx]);
    break;
  case Muon::Property::PfRelIso04:
    muon.SetPfRelIso04(Muon_pfRelIso04_all[idx]);
    break;
  case Muon::Property::MiniPFRelIso:
    muon.SetMiniPFRelIso(Muon_miniPFRelIso_all[idx]);
    break;
  case Muon::Property::Dxy:
    muon.SetdXY(Muon_dxy[idx]);
    break;
  case Muon::Property::DxyErr:
    muon.SetdXYErr(Muon_dxyErr[idx]);
    break;
  case Muon::Property::Dz:
    muon.SetdZ(Muon_dz[idx]);
    break;
  case Muon::Property::DzErr:
    muon.SetdZErr(Muon_dzErr[idx]);
    break;
  case Muon::Property::Ip3d:
    muon.SetIP3D(Muon_ip3d[idx]);
    break;
  case Muon::Property::Sip3d:
    muon.SetSIP3D(Muon_sip3d[idx]);
    break;
  case Muon::Property::TrackerLayers:
    muon.SetNTrackerLayers(Muon_nTrackerLayers[idx]);
    break;
  case Muon::Property::GenPartFlav:
    muon.SetGenPartFlav(Muon_genPartFlav[idx]);
    break;
  case Muon::Property::GenPartIdx:
    muon.SetGenPartIdx(Muon_genPartIdx[idx]);
    break;
  case Muon::Property::JetIdx:
    muon.SetJetIdx(Muon_jetIdx[idx]);
    break;
  case Muon::Property::LooseId:
    muon.SetBIDBit(Muon::BooleanID::LOOSE, Muon_looseId[idx]);
    break;
  case Muon::Property::MediumId:
    muon.SetBIDBit(Muon::BooleanID::MEDIUM, Muon_mediumId[idx]);
    break;
  case Muon::Property::MediumPromptId:
    muon.SetBIDBit(Muon::BooleanID::MEDIUMPROMPT, Muon_mediumPromptId[idx]);
    break;
  case Muon::Property::TightId:
    muon.SetBIDBit(Muon::BooleanID::TIGHT, Muon_tightId[idx]);
    break;
  case Muon::Property::SoftId:
    muon.SetBIDBit(Muon::BooleanID::SOFT, Muon_softId[idx]);
    break;
  case Muon::Property::SoftMvaId:
    muon.SetBIDBit(Muon::BooleanID::SOFTMVA, Muon_softMvaId[idx]);
    break;
  case Muon::Property::TriggerLooseId:
    muon.SetBIDBit(Muon::BooleanID::TRIGGERLOOSE, Muon_triggerIdLoose[idx]);
    break;
  case Muon::Property::HighPtId:
    muon.SetWIDBit(Muon::WorkingPointID::HIGHPT, Muon_highPtId[idx]);
    break;
  case Muon::Property::MiniIsoId:
    muon.SetWIDBit(Muon::WorkingPointID::MINIISO, Muon_miniIsoId[idx]);
    break;
  case Muon::Property::MultiIsoId:
    muon.SetWIDBit(Muon::WorkingPointID::MULTIISO, Muon_multiIsoId[idx]);
    break;
  case Muon::Property::MvaMuId:
    muon.SetWIDBit(Muon::WorkingPointID::MVAMU, Muon_mvaMuID_WP[idx]);
    break;
  case Muon::Property::PfIsoId:
    muon.SetWIDBit(Muon::WorkingPointID::PFISO, Muon_pfIsoId[idx]);
    break;
  case Muon::Property::PuppiIsoId:
    muon.SetWIDBit(Muon::WorkingPointID::PUPPIISO, Muon_puppiIsoId[idx]);
    break;
  case Muon::Property::TkIsoId:
    muon.SetWIDBit(Muon::WorkingPointID::TKISO, Muon_tkIsoId[idx]);
    break;
  case Muon::Property::SoftMva:
    muon.SetMVAID(Muon::MVAID::SOFTMVA, Muon_softMva[idx]);
    break;
  case Muon::Property::MvaLowPt:
    muon.SetMVAID(Muon::MVAID::MVALOWPT, Muon_mvaLowPt[idx]);
    break;
  case Muon::Property::MvaPrompt:
    muon.SetMVAID(Muon::MVAID::MVAPROMPT, Muon_promptMVA[idx]);
    break;
  case Muon::Property::SoftMvaRun3:
    muon.SetMVAID(Muon::MVAID::SOFTMVARUN3, Muon_softMvaRun3[idx]);
    break;
  default:
    throw std::runtime_error(
        "[AnalyzerCore::EnsureMuonProperty] Unsupported property");
  }
}

void AnalyzerCore::EnsureElectronProperty(Electron &electron,
                                          Electron::Property property) const {
  if (!electron.lazy_)
    return;
  const int idx = electron.lazy_->index;
  switch (property) {
  case Electron::Property::PfRelIso03:
    electron.SetPfRelIso03(Electron_pfRelIso03_all[idx]);
    break;
  case Electron::Property::MiniPFRelIso:
    electron.SetMiniPFRelIso(Electron_miniPFRelIso_all[idx]);
    break;
  case Electron::Property::Dxy:
    electron.SetdXY(Electron_dxy[idx]);
    break;
  case Electron::Property::DxyErr:
    electron.SetdXYErr(Electron_dxyErr[idx]);
    break;
  case Electron::Property::Dz:
    electron.SetdZ(Electron_dz[idx]);
    break;
  case Electron::Property::DzErr:
    electron.SetdZErr(Electron_dzErr[idx]);
    break;
  case Electron::Property::Ip3d:
    electron.SetIP3D(Electron_ip3d[idx]);
    break;
  case Electron::Property::Sip3d:
    electron.SetSIP3D(Electron_sip3d[idx]);
    break;
  case Electron::Property::ConvVeto:
    electron.SetConvVeto(Electron_convVeto[idx]);
    break;
  case Electron::Property::LostHits:
    electron.SetLostHits(Electron_lostHits[idx]);
    break;
  case Electron::Property::SeedGain:
    electron.SetSeedGain(Electron_seedGain[idx]);
    break;
  case Electron::Property::TightCharge:
    electron.SetTightCharge(Electron_tightCharge[idx]);
    break;
  case Electron::Property::Sieie:
    electron.SetSieie(Electron_sieie[idx]);
    break;
  case Electron::Property::Hoe:
    electron.SetHoe(Electron_hoe[idx]);
    break;
  case Electron::Property::EInvMinusPInv:
    electron.SetEInvMinusPInv(Electron_eInvMinusPInv[idx]);
    break;
  case Electron::Property::Dr03EcalRecHitSumEt:
    electron.SetDr03EcalRecHitSumEt(Electron_dr03EcalRecHitSumEt[idx]);
    break;
  case Electron::Property::Dr03HcalDepth1TowerSumEt:
    electron.SetDr03HcalDepth1TowerSumEt(
        Electron_dr03HcalDepth1TowerSumEt[idx]);
    break;
  case Electron::Property::Dr03TkSumPt:
    electron.SetDr03TkSumPt(Electron_dr03TkSumPt[idx]);
    break;
  case Electron::Property::Dr03TkSumPtHEEP:
    electron.SetDr03TkSumPtHEEP(Electron_dr03TkSumPtHEEP[idx]);
    break;
  case Electron::Property::R9:
    electron.SetR9(Electron_r9[idx]);
    break;
  case Electron::Property::EnergyErr:
    electron.SetEnergyErr(Electron_energyErr[idx]);
    break;
  case Electron::Property::CutBasedHEEP:
    electron.SetBIDBit(Electron::BooleanID::CUTBASEDHEEP,
                       Electron_cutBased_HEEP[idx]);
    break;
  case Electron::Property::PromptMVA:
    electron.SetMVA(Electron::MVATYPE::MVAPROMPT, Electron_promptMVA[idx]);
    break;
  case Electron::Property::MvaIsoWP80:
    electron.SetBIDBit(Electron::BooleanID::MVAISOWP80,
                       Electron_mvaIso_WP80[idx]);
    break;
  case Electron::Property::MvaIsoWP90:
    electron.SetBIDBit(Electron::BooleanID::MVAISOWP90,
                       Electron_mvaIso_WP90[idx]);
    break;
  case Electron::Property::MvaNoIsoWP80:
    electron.SetBIDBit(Electron::BooleanID::MVANOISOWP80,
                       Electron_mvaNoIso_WP80[idx]);
    break;
  case Electron::Property::MvaNoIsoWP90:
    electron.SetBIDBit(Electron::BooleanID::MVANOISOWP90,
                       Electron_mvaNoIso_WP90[idx]);
    break;
  case Electron::Property::MvaIso:
    electron.SetMVA(Electron::MVATYPE::MVAISO, Electron_mvaIso[idx]);
    break;
  case Electron::Property::MvaNoIso:
    electron.SetMVA(Electron::MVATYPE::MVANOISO, Electron_mvaNoIso[idx]);
    break;
  case Electron::Property::CutBased:
    electron.SetCBIDBit(Electron::CutBasedID::CUTBASED, Electron_cutBased[idx]);
    break;
  case Electron::Property::GenPartFlav:
    electron.SetGenPartFlav(Electron_genPartFlav[idx]);
    break;
  case Electron::Property::GenPartIdx:
    electron.SetGenPartIdx(Electron_genPartIdx[idx]);
    break;
  case Electron::Property::JetIdx:
    electron.SetJetIdx(Electron_jetIdx[idx]);
    break;
  default:
    throw std::runtime_error(
        "[AnalyzerCore::EnsureElectronProperty] Unsupported property");
  }
}

void AnalyzerCore::EnsureJetProperty(Jet &jet, Jet::Property property) const {
  if (!jet.lazy_)
    return;
  const int idx = jet.lazy_->index;
  using Tagger = JetTagging::JetFlavTagger;
  using Score = JetTagging::JetFlavTaggerScoreType;
  switch (property) {
  case Jet::Property::DeepFlavB:
    jet.SetTaggerScore(Tagger::DeepJet, Score::B, Jet_btagDeepFlavB[idx]);
    break;
  case Jet::Property::DeepFlavCvB:
    jet.SetTaggerScore(Tagger::DeepJet, Score::CvB, Jet_btagDeepFlavCvB[idx]);
    break;
  case Jet::Property::DeepFlavCvL:
    jet.SetTaggerScore(Tagger::DeepJet, Score::CvL, Jet_btagDeepFlavCvL[idx]);
    break;
  case Jet::Property::DeepFlavQG:
    jet.SetTaggerScore(Tagger::DeepJet, Score::QvG, Jet_btagDeepFlavQG[idx]);
    break;
  case Jet::Property::PNetB:
    jet.SetTaggerScore(Tagger::ParticleNet, Score::B, Jet_btagPNetB[idx]);
    break;
  case Jet::Property::PNetCvB:
    jet.SetTaggerScore(Tagger::ParticleNet, Score::CvB, Jet_btagPNetCvB[idx]);
    break;
  case Jet::Property::PNetCvL:
    jet.SetTaggerScore(Tagger::ParticleNet, Score::CvL, Jet_btagPNetCvL[idx]);
    break;
  case Jet::Property::PNetCvNotB:
    jet.SetTaggerScore(Tagger::ParticleNet, Score::CvNotB,
                       Jet_btagPNetCvNotB[idx]);
    break;
  case Jet::Property::PNetQvG:
    jet.SetTaggerScore(Tagger::ParticleNet, Score::QvG, Jet_btagPNetQvG[idx]);
    break;
  case Jet::Property::PNetTauVJet:
    jet.SetTaggerScore(Tagger::ParticleNet, Score::TauVJet,
                       Jet_btagPNetTauVJet[idx]);
    break;
  case Jet::Property::ParTB:
    jet.SetTaggerScore(Tagger::ParT, Score::B, Jet_btagUParTAK4B[idx]);
    break;
  case Jet::Property::ParTCvB:
    jet.SetTaggerScore(Tagger::ParT, Score::CvB, Jet_btagUParTAK4CvB[idx]);
    break;
  case Jet::Property::ParTCvL:
    jet.SetTaggerScore(Tagger::ParT, Score::CvL, Jet_btagUParTAK4CvL[idx]);
    break;
  case Jet::Property::ParTCvNotB:
    jet.SetTaggerScore(Tagger::ParT, Score::CvNotB,
                       Jet_btagUParTAK4CvNotB[idx]);
    break;
  case Jet::Property::ParTEle:
    jet.SetTaggerScore(Tagger::ParT, Score::Ele, Jet_btagUParTAK4Ele[idx]);
    break;
  case Jet::Property::ParTMu:
    jet.SetTaggerScore(Tagger::ParT, Score::Mu, Jet_btagUParTAK4Mu[idx]);
    break;
  case Jet::Property::ParTQvG:
    jet.SetTaggerScore(Tagger::ParT, Score::QvG, Jet_btagUParTAK4QvG[idx]);
    break;
  case Jet::Property::ParTSvCB:
    jet.SetTaggerScore(Tagger::ParT, Score::SvCB, Jet_btagUParTAK4SvCB[idx]);
    break;
  case Jet::Property::ParTSvUDG:
    jet.SetTaggerScore(Tagger::ParT, Score::SvUDG, Jet_btagUParTAK4SvUDG[idx]);
    break;
  case Jet::Property::ParTTauVJet:
    jet.SetTaggerScore(Tagger::ParT, Score::TauVJet,
                       Jet_btagUParTAK4TauVJet[idx]);
    break;
  case Jet::Property::ParTUDG:
    jet.SetTaggerScore(Tagger::ParT, Score::probUDG, Jet_btagUParTAK4UDG[idx]);
    break;
  case Jet::Property::ParTProbB:
    jet.SetTaggerScore(Tagger::ParT, Score::probB, Jet_btagUParTAK4probb[idx]);
    break;
  case Jet::Property::ParTProbBB:
    jet.SetTaggerScore(Tagger::ParT, Score::probBB,
                       Jet_btagUParTAK4probbb[idx]);
    break;
  default:
    throw std::runtime_error(
        "[AnalyzerCore::EnsureJetProperty] Unsupported property");
  }
}

RVec<Photon> AnalyzerCore::GetAllPhotons() {
  RVec<Photon> Photons;
  for (int i = 0; i < nPhoton; i++) {
    Photon photon;
    photon.SetPtEtaPhiE(Photon_pt[i], Photon_eta[i], Photon_phi[i], 1.0);
    float photon_theta = photon.Theta();
    float photon_E = Photon_pt[i] / sin(photon_theta);
    photon.SetPtEtaPhiE(Photon_pt[i], Photon_eta[i], Photon_phi[i], photon_E);
    photon.SetEnergy(photon_E);
    photon.SetSigmaIetaIeta(Photon_sieie[i]);
    photon.SetHoe(Photon_hoe[i]);
    if (Run == 3) {
      photon.SetEnergyRaw(Photon_energyRaw[i]);
    }
    photon.SetPixelSeed(Photon_pixelSeed[i]);
    photon.SetisScEtaEB(Photon_isScEtaEB[i]);
    photon.SetisScEtaEE(Photon_isScEtaEE[i]);
    //        Need to Add PV Info
    //        Photon.SetSCEta(Photon_eta[i], Photon_phi[i], PV_x, PV_y, PV_z,
    //        Photon_isScEtaEB[i], Photon_isScEtaEE[i]);
    photon.SetSCEta(Photon_eta[i], Photon_phi[i], 0, 0, 0, false, false);
    photon.SetBIDBit(Photon::BooleanID::MVAIDWP80, Photon_mvaID_WP80[i]);
    photon.SetBIDBit(Photon::BooleanID::MVAIDWP90, Photon_mvaID_WP90[i]);
    photon.SetCBIDBit(Photon::CutBasedID::CUTBASED, Photon_cutBased[i]);

    Photons.push_back(photon);
  }
  return Photons;
}

RVec<Photon> AnalyzerCore::GetPhotons(TString id, double ptmin,
                                      double fetamax) {
  RVec<Photon> photons = GetAllPhotons();
  RVec<Photon> out;
  for (unsigned int i = 0; i < photons.size(); i++) {
    if (!(photons.at(i).Pt() > ptmin)) {
      continue;
    }
    if (!(fabs(photons.at(i).scEta()) < fetamax)) {
      continue;
    }
    if (!(photons.at(i).PassID(id))) {
      continue;
    }
    out.push_back(photons.at(i));
  }
  return out;
}

RVec<Jet> AnalyzerCore::SelectJets(const JetViewCollection &jets,
                                   const Jet::JetID ID, const float ptmin,
                                   const float fetamax) const {
  RVec<Jet> selected_jets;
  if (jets.empty())
    return selected_jets;

  selected_jets.reserve(jets.size());
  for (std::size_t idx = 0; idx < jets.size(); ++idx) {
    const auto &jetView = jets[idx];
    if (!jetView.valid())
      continue;
    const float jetPt = jetView.SmearedPtNominal();
    if (jetPt < ptmin)
      continue;
    if (std::fabs(jetView.Eta()) > fetamax)
      continue;
    Jet jet = MaterializeJet(jets, idx);
    if (!myCorr->PassJetID(jet, ID))
      continue;
    selected_jets.push_back(std::move(jet));
  }
  return selected_jets;
}

std::vector<std::size_t> AnalyzerCore::SelectJetIndices(
    const JetViewCollection &jets, const std::vector<std::size_t> &seed_indices,
    const Jet::JetID id, const float ptmin, const float fetamax,
    const MyCorrection::variation &jesVar,
    const MyCorrection::variation &jerVar) const {
  std::vector<std::size_t> selected;
  if (jets.empty() || seed_indices.empty())
    return selected;

  selected.reserve(seed_indices.size());
  const float etaMax2 = fetamax * fetamax;

  const bool isData = IsDATA;
  const bool jesIsNom = (jesVar == MyCorrection::variation::nom);
  const bool jerIsNom = (jerVar == MyCorrection::variation::nom);

  // 둘 다 non-nominal은 아직 지원 안 한다
  if (!jesIsNom && !jerIsNom) {
    throw std::runtime_error("[AnalyzerCore::SelectJetIndices] Both JES and "
                             "JER are non-nominal, not supported");
  }

// 벡터라이저에게 힌트
#pragma GCC ivdep
  for (std::size_t k = 0; k < seed_indices.size(); ++k) {
    const auto idx = seed_indices[k];
    const auto &jetView = jets[idx];
    if (!jetView.valid())
      continue;

    float jetPt = 0.f;
    if (isData) {
      jetPt = jetView.CorrectedPt();
    } else {
      if (jerIsNom) {
        if (jesIsNom) {
          jetPt = jetView.SmearedPtNominal();
        } else if (jesVar == MyCorrection::variation::up) {
          jetPt = jetView.JesPtUp();
        } else if (jesVar == MyCorrection::variation::down) {
          jetPt = jetView.JesPtDown();
        } else {
          continue;
        }
      } else {
        if (jerVar == MyCorrection::variation::up) {
          jetPt = jetView.SmearedPtUp();
        } else if (jerVar == MyCorrection::variation::down) {
          jetPt = jetView.SmearedPtDown();
        } else {
          continue;
        }
      }
    }

    if (jetPt < ptmin)
      continue;

    const float eta = jetView.Eta();
    if (eta * eta > etaMax2)
      continue;

    if (!myCorr->PassJetID(jetView, id))
      continue;

    selected.push_back(idx);
  }

  return selected;
}

std::vector<std::size_t> AnalyzerCore::SelectJetIndices(
    const JetViewCollection &jets, const Jet::JetID ID, const float ptmin,
    const float fetamax, const MyCorrection::variation &JESVariation,
    const MyCorrection::variation &JERVariation) const {
  std::vector<std::size_t> seed_indices = std::vector<std::size_t>(jets.size());
  std::iota(seed_indices.begin(), seed_indices.end(), 0);
  return SelectJetIndices(jets, seed_indices, ID, ptmin, fetamax, JESVariation,
                          JERVariation);
}

RVec<Jet> AnalyzerCore::ScaleJets(const JetViewCollection &jets,
                                  const std::vector<std::size_t> &indices,
                                  const MyCorrection::variation &syst,
                                  const TString &source) {
  RVec<Jet> materialised = MaterializeJets(jets, indices);
  return ScaleJets(materialised, syst, source);
}

RVec<Jet> AnalyzerCore::SmearJets(const JetViewCollection &jets,
                                  const std::vector<std::size_t> &indices,
                                  const RVec<GenJet> &genjets,
                                  const MyCorrection::variation &syst,
                                  const TString &source) {
  RVec<Jet> materialised = MaterializeJets(jets, indices);
  return SmearJets(materialised, genjets, syst, source);
}

RVec<Jet> AnalyzerCore::SelectJets(const RVec<Jet> &jets, const Jet::JetID ID,
                                   const float ptmin,
                                   const float fetamax) const {
  RVec<Jet> selected_jets;
  for (const auto &jet : jets) {
    if (jet.Pt() < ptmin)
      continue;
    if (fabs(jet.Eta()) > fetamax)
      continue;
    if (!myCorr->PassJetID(jet, ID))
      continue;
    selected_jets.push_back(jet);
  }
  return selected_jets;
}

RVec<FatJet> AnalyzerCore::SelectFatJets(const RVec<FatJet> &fatjets,
                                         const FatJet::FatJetID ID,
                                         const float ptmin,
                                         const float fetamax) const {
  RVec<FatJet> selected_jets;
  for (const auto &fatjet : fatjets) {
    if (fatjet.Pt() < ptmin)
      continue;
    if (fabs(fatjet.Eta()) > fetamax)
      continue;
    if (!myCorr->PassFatJetID(fatjet, ID))
      continue;
    selected_jets.push_back(fatjet);
  }
  return selected_jets;
}

RVec<Jet> AnalyzerCore::JetsVetoLeptonInside(
    const RVec<Jet> &jets, const ElectronViewCollection &electrons,
    const MuonViewCollection &muons, const float dR) const {
  RVec<Jet> selected_jets;
  selected_jets.reserve(jets.size());
  for (const auto &jet : jets) {
    bool isLeptonInside = false;
    for (const auto &electron : electrons) {
      if (jet.DeltaR(electron.P4()) < dR) {
        isLeptonInside = true;
        break;
      }
    }
    if (isLeptonInside)
      continue;
    for (const auto &muon : muons) {
      if (jet.DeltaR(muon.P4()) < dR) {
        isLeptonInside = true;
        break;
      }
    }
    if (!isLeptonInside)
      selected_jets.push_back(jet);
  }
  return selected_jets;
}

std::vector<std::size_t> AnalyzerCore::JetsVetoLeptonInside(
    const JetViewCollection &jets, const std::vector<std::size_t> &jet_indices,
    const ElectronViewCollection &electrons,
    const std::vector<std::size_t> &electron_indices,
    const MuonViewCollection &muons,
    const std::vector<std::size_t> &muon_indices, const float dR) const {
  std::vector<std::size_t> filtered;
  if (jets.empty() || jet_indices.empty())
    return filtered;

  filtered.reserve(jet_indices.size());
  const float dR2 = dR * dR;

  constexpr float PI = 3.14159265358979323846f;
  constexpr float TWOPI = 2.0f * PI;

  auto deltaPhiWrap = [PI, TWOPI](float dphi) noexcept -> float {
    if (dphi > PI)
      dphi -= TWOPI;
    else if (dphi <= -PI)
      dphi += TWOPI;
    return dphi;
  };

  const std::size_t nEle = electron_indices.size();
  const std::size_t nMu = muon_indices.size();

  std::vector<float> ele_eta, ele_phi, mu_eta, mu_phi;
  ele_eta.reserve(nEle);
  ele_phi.reserve(nEle);
  mu_eta.reserve(nMu);
  mu_phi.reserve(nMu);

  for (std::size_t i = 0; i < nEle; ++i) {
    const auto &e = electrons[electron_indices[i]];
    ele_eta.push_back(e.Eta());
    ele_phi.push_back(e.Phi());
  }
  for (std::size_t i = 0; i < nMu; ++i) {
    const auto &m = muons[muon_indices[i]];
    mu_eta.push_back(m.Eta());
    mu_phi.push_back(m.Phi());
  }

  // 2) jet 루프
  for (const auto idx : jet_indices) {
    if (idx >= jets.size())
      continue;

    const auto &j = jets[idx];
    const float j_eta = j.Eta();
    const float j_phi = j.Phi();

    bool veto = false;

    // electrons
    for (std::size_t i = 0; i < nEle; ++i) {
      const float dEta = j_eta - ele_eta[i];
      float dPhi = j_phi - ele_phi[i];
      dPhi = deltaPhiWrap(dPhi);

      const float dr2 = dEta * dEta + dPhi * dPhi;
      if (dr2 < dR2) {
        veto = true;
        break;
      }
    }
    if (veto)
      continue;

    // muons
    for (std::size_t i = 0; i < nMu; ++i) {
      const float dEta = j_eta - mu_eta[i];
      float dPhi = j_phi - mu_phi[i];
      dPhi = deltaPhiWrap(dPhi);

      const float dr2 = dEta * dEta + dPhi * dPhi;
      if (dr2 < dR2) {
        veto = true;
        break;
      }
    }

    if (!veto)
      filtered.push_back(idx);
  }

  return filtered;
}

RVec<Jet> AnalyzerCore::JetsVetoLeptonInside(const RVec<Jet> &jets,
                                             const RVec<Electron> &electrons,
                                             const RVec<Muon> &muons,
                                             const float dR) const {
  RVec<Jet> selected_jets;
  for (const auto &jet : jets) {
    bool isLeptonInside = false;
    for (const auto &electron : electrons) {
      if (jet.DeltaR(electron) < dR) {
        isLeptonInside = true;
        break;
      }
    }
    if (isLeptonInside)
      continue;
    for (const auto &muon : muons) {
      if (jet.DeltaR(muon) < dR) {
        isLeptonInside = true;
        break;
      }
    }
    if (isLeptonInside)
      continue;
    selected_jets.push_back(jet);
  }
  return selected_jets;
}

// For Run2, reject jets within the veto map
bool AnalyzerCore::PassJetVetoMap(const Jet &jet,
                                  const MuonViewCollection &AllMuons,
                                  const TString mapCategory) {
  if (!(Run == 2))
    return true;

  if (jet.chEmEF() + jet.neEmEF() > 0.9)
    return true;

  bool pass_loose_selection = jet.Pt() > 15.;
  pass_loose_selection =
      pass_loose_selection && myCorr->PassJetID(jet, Jet::JetID::TIGHT);
  pass_loose_selection =
      pass_loose_selection &&
      (jet.Pt() > 50. || myCorr->PassJetID(jet, Jet::JetID::PUID_LOOSE));
  for (const auto &muon : AllMuons) {
    if (jet.DeltaR(muon.P4()) <= 0.2) {
      pass_loose_selection = false;
      break;
    }
  }
  bool pass_veto_map =
      pass_loose_selection &&
      (!myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory));
  return pass_veto_map;
}

bool AnalyzerCore::PassJetVetoMap(const JetViewCollection &AllJets,
                                  const MuonViewCollection &AllMuons,
                                  const TString mapCategory) {
  if (!(Run == 3))
    return true;
  std::vector<size_t> this_jet_indices =
      SelectJetIndices(AllJets, Jet::JetID::TIGHT, 15., 5.0);
  ElectronViewCollection empty_electron_views;
  std::vector<size_t> empty_electron_indices;
  std::vector<size_t> all_muon_indices(AllMuons.size());
  std::iota(all_muon_indices.begin(), all_muon_indices.end(), 0);

  this_jet_indices = JetsVetoLeptonInside(
      AllJets, this_jet_indices, empty_electron_views, empty_electron_indices,
      AllMuons, all_muon_indices, 0.2);
  for (const auto idx : this_jet_indices) {
    if (idx >= AllJets.size())
      throw std::runtime_error(
          "[AnalyzerCore::PassJetVetoMap] Jet index out of range");

    const auto &jet = AllJets[idx]; // JetView
    if (myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory)) {
      return false;
    }
  }
  return true;
}

bool AnalyzerCore::PassJetVetoMap(const RVec<Jet> &AllJets,
                                  const MuonViewCollection &AllMuons,
                                  const TString mapCategory) {
  if (!(Run == 3))
    return true;
  RVec<Jet> this_jet = SelectJets(AllJets, Jet::JetID::TIGHT, 15., 5.0);
  this_jet = SelectJets(this_jet, Jet::JetID::PUID_LOOSE, 15., 5.0);
  ElectronViewCollection empty_electrons;
  this_jet = JetsVetoLeptonInside(this_jet, empty_electrons, AllMuons, 0.2);
  for (const auto &jet : this_jet) {
    if (myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory))
      return false;
  }
  return true;
}

bool AnalyzerCore::PassJetVetoMap(const Jet &jet, const RVec<Muon> &AllMuons,
                                  const TString mapCategory) {
  if (!(Run == 2))
    return true;
  // Only apply to the jets with em energy fraction less than 0.9
  if (jet.chEmEF() + jet.neEmEF() > 0.9)
    return true;

  // Selections should be looser than Analysis jet selections
  bool pass_loose_selection = jet.Pt() > 15.;
  pass_loose_selection =
      pass_loose_selection && myCorr->PassJetID(jet, Jet::JetID::TIGHT);
  pass_loose_selection =
      pass_loose_selection &&
      (jet.Pt() > 50. || myCorr->PassJetID(jet, Jet::JetID::PUID_LOOSE));
  for (const auto &muon : AllMuons) {
    pass_loose_selection = pass_loose_selection && (jet.DeltaR(muon) > 0.2);
  }
  bool pass_veto_map =
      pass_loose_selection &&
      (!myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory));
  return pass_veto_map;
}

// For Run3, reject events if any jet is within the veto map
bool AnalyzerCore::PassJetVetoMap(const RVec<Jet> &AllJets,
                                  const RVec<Muon> &AllMuons,
                                  const TString mapCategory) {
  if (!(Run == 3))
    return true;
  RVec<Jet> selected_jets;
  RVec<Electron> empty_electrons;

  RVec<Jet> this_jet = SelectJets(AllJets, Jet::JetID::TIGHT, 15., 5.0);
  this_jet = JetsVetoLeptonInside(this_jet, empty_electrons, AllMuons, 0.2);
  for (const auto &jet : this_jet) {
    if (jet.chEmEF() + jet.neEmEF() < 0.9)
      selected_jets.push_back(jet);
  }
  for (const auto &jet : selected_jets) {
    if (myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory))
      return false;
  }
  return true;
}

RVec<FatJet> AnalyzerCore::GetAllFatJets() {

  RVec<FatJet> FatJets;

  for (int i = 0; i < nFatJet; i++) {

    FatJet fatjet;
    fatjet.SetOriginalIndex(i);
    fatjet.SetPtEtaPhiM(FatJet_pt[i], FatJet_eta[i], FatJet_phi[i],
                        FatJet_mass[i]);
    fatjet.SetArea(FatJet_area[i]);
    fatjet.SetRawFactor(FatJet_rawFactor[i]);
    fatjet.SetEnergyFraction(FatJet_chEmEF[i], FatJet_chHEF[i],
                             FatJet_hfEmEF[i], FatJet_hfHEF[i], FatJet_muEF[i],
                             FatJet_neEmEF[i], FatJet_neHEF[i]);
    fatjet.SetConstituents(FatJet_chMultiplicity[i], FatJet_neMultiplicity[i],
                           FatJet_nConstituents[i]);
    fatjet.SetGenMatchIDs(FatJet_genJetAK8Idx[i], FatJet_subJetIdx1[i],
                          FatJet_subJetIdx2[i]);
    fatjet.SetSDMass(FatJet_msoftdrop[i]);
    fatjet.SetSubjettiness(FatJet_tau1[i], FatJet_tau2[i], FatJet_tau3[i],
                           FatJet_n2b1[i], FatJet_n3b1[i]);
    fatjet.SetLSF3(FatJet_lsf3[i]);

    const int nFatPFCand = static_cast<int>(nFatJetPFCand);
    for (int idx = 0; idx < nFatPFCand; ++idx) {
      if (FatJetPFCand_jetIdx[idx] == i) {
        JetConstituent constituent;
        const auto pfCandTableidx =
            static_cast<std::size_t>(FatJetPFCand_pfCandIdx[idx]);
        constituent.SetPtEtaPhiM(
            PFCand_pt[pfCandTableidx], PFCand_eta[pfCandTableidx],
            PFCand_phi[pfCandTableidx], PFCand_mass[pfCandTableidx]);
        constituent.SetPdgId(PFCand_pdgId[pfCandTableidx]);
        if (PFCand_puppiWeight.valid() &&
            pfCandTableidx < PFCand_puppiWeight.size())
          constituent.SetPUPPIWeight(PFCand_puppiWeight[pfCandTableidx]);
      }
    }
    RVec<float> pnet;
    pnet = {FatJet_particleNetWithMass_H4qvsQCD[i],
            FatJet_particleNetWithMass_HbbvsQCD[i],
            FatJet_particleNetWithMass_HccvsQCD[i],
            FatJet_particleNetWithMass_QCD[i],
            FatJet_particleNetWithMass_TvsQCD[i],
            FatJet_particleNetWithMass_WvsQCD[i],
            FatJet_particleNetWithMass_ZvsQCD[i]};

    RVec<float> pnet_nomass;
    pnet_nomass = {
        FatJet_particleNet_QCD[i],      FatJet_particleNet_QCD0HF[i],
        FatJet_particleNet_QCD1HF[i],   FatJet_particleNet_QCD2HF[i],
        FatJet_particleNet_WVsQCD[i],   FatJet_particleNet_XbbVsQCD[i],
        FatJet_particleNet_XccVsQCD[i], FatJet_particleNet_XggVsQCD[i],
        FatJet_particleNet_XqqVsQCD[i], FatJet_particleNet_XteVsQCD[i],
        FatJet_particleNet_XtmVsQCD[i], FatJet_particleNet_XttVsQCD[i],
        FatJet_particleNet_massCorr[i]};

    RVec<float> gloparT;
    gloparT = {FatJet_globalParT3_QCD[i],
               FatJet_globalParT3_TopbWev[i],
               FatJet_globalParT3_TopbWmv[i],
               FatJet_globalParT3_TopbWq[i],
               FatJet_globalParT3_TopbWqq[i],
               FatJet_globalParT3_TopbWtauhv[i],
               FatJet_globalParT3_WvsQCD[i],
               FatJet_globalParT3_XWW3q[i],
               FatJet_globalParT3_XWW4q[i],
               FatJet_globalParT3_XWWqqev[i],
               FatJet_globalParT3_XWWqqmv[i],
               FatJet_globalParT3_Xbb[i],
               FatJet_globalParT3_Xcc[i],
               FatJet_globalParT3_Xcs[i],
               FatJet_globalParT3_Xqq[i],
               FatJet_globalParT3_Xtauhtaue[i],
               FatJet_globalParT3_Xtauhtauh[i],
               FatJet_globalParT3_Xtauhtaum[i],
               FatJet_globalParT3_withMassTopvsQCD[i],
               FatJet_globalParT3_withMassWvsQCD[i],
               FatJet_globalParT3_withMassZvsQCD[i]};

    fatjet.SetParticleNetResult(pnet);
    fatjet.SetParticleNetNoMassResult(pnet_nomass);
    fatjet.SetGloParTResult(gloparT);
    FatJets.push_back(fatjet);
  }

  return FatJets;
}

RVec<GenJet> AnalyzerCore::GetAllGenJets() {
  GenJetViewCollection view = GetAllGenJetViews();
  return MaterializeGenJets(view);
}

RVec<GenDressedLepton> AnalyzerCore::GetAllGenDressedLeptons() {
  RVec<GenDressedLepton> GenDressedLeptons;
  if (IsDATA)
    return GenDressedLeptons;

  for (int i = 0; i < nGenDressedLepton; i++) {
    GenDressedLepton lep;
    lep.SetPtEtaPhiM(GenDressedLepton_pt[i], GenDressedLepton_eta[i],
                     GenDressedLepton_phi[i], GenDressedLepton_mass[i]);
    lep.SetPdgId(GenDressedLepton_pdgId[i]);
    lep.SetHasTauAnc(GenDressedLepton_hasTauAnc[i]);
    GenDressedLeptons.push_back(lep);
  }
  return GenDressedLeptons;
}

RVec<GenIsolatedPhoton> AnalyzerCore::GetAllGenIsolatedPhotons() {
  RVec<GenIsolatedPhoton> GenIsolatedPhotons;
  if (IsDATA)
    return GenIsolatedPhotons;

  for (int i = 0; i < nGenIsolatedPhoton; i++) {
    GenIsolatedPhoton photon;
    photon.SetPtEtaPhiM(GenIsolatedPhoton_pt[i], GenIsolatedPhoton_eta[i],
                        GenIsolatedPhoton_phi[i], GenIsolatedPhoton_mass[i]);
    GenIsolatedPhotons.push_back(photon);
  }
  return GenIsolatedPhotons;
}

RVec<GenVisTau> AnalyzerCore::GetAllGenVisTaus() {
  RVec<GenVisTau> GenVisTaus;
  if (IsDATA)
    return GenVisTaus;

  for (int i = 0; i < nGenVisTau; i++) {
    GenVisTau tau;
    tau.SetPtEtaPhiM(GenVisTau_pt[i], GenVisTau_eta[i], GenVisTau_phi[i],
                     GenVisTau_mass[i]);
    tau.SetCharge(GenVisTau_charge[i]);
    tau.SetGenPartIdxMother(GenVisTau_genPartIdxMother[i]);
    tau.SetStatus(GenVisTau_status[i]);
    GenVisTaus.push_back(tau);
  }
  return GenVisTaus;
}

TrigObjViewCollection AnalyzerCore::GetAllTrigObjViews() {
  auto storage = std::make_shared<TrigObjSoA>();
  storage->run = Run;
  storage->pt.bind(&TrigObj_pt);
  storage->eta.bind(&TrigObj_eta);
  storage->phi.bind(&TrigObj_phi);
  storage->id.bind(&TrigObj_id);
  storage->filterBits.bind(&TrigObj_filterBits);
  storage->l1charge.bind(&TrigObj_l1charge);
  storage->l1iso.bind(&TrigObj_l1iso);
  storage->l1pt.bind(&TrigObj_l1pt);
  storage->l1pt2.bind(&TrigObj_l1pt_2);
  storage->l2pt.bind(&TrigObj_l2pt);
  return TrigObjViewCollection(std::move(storage));
}

RVec<TrigObj>
AnalyzerCore::MaterializeTrigObjs(const TrigObjViewCollection &trigobjs) const {
  RVec<TrigObj> materialised;
  if (trigobjs.size() == 0)
    return materialised;

  auto storage = trigobjs.storage();
  if (!storage)
    return materialised;

  materialised.reserve(trigobjs.size());
  for (std::size_t idx = 0; idx < trigobjs.size(); ++idx) {
    materialised.emplace_back(storage, idx);
  }
  return materialised;
}

RVec<TrigObj> AnalyzerCore::MaterializeTrigObjs(
    const TrigObjViewCollection &trigobjs,
    const std::vector<std::size_t> &indices) const {
  RVec<TrigObj> materialised;
  if (indices.empty() || trigobjs.size() == 0)
    return materialised;

  auto storage = trigobjs.storage();
  if (!storage)
    return materialised;

  materialised.reserve(indices.size());
  const std::size_t n = trigobjs.size();
  for (const auto idx : indices) {
    if (idx >= n)
      continue;
    materialised.emplace_back(storage, idx);
  }
  return materialised;
}

std::vector<std::size_t> AnalyzerCore::SelectTrigObjIndices(
    const TrigObjViewCollection &trigobjs,
    const std::vector<std::size_t> &seed_indices, const int id,
    const float ptmin, const float fetamax) const {
  std::vector<std::size_t> selected;
  selected.reserve(seed_indices.size());
  for (const auto idx : seed_indices) {
    const auto &obj = trigobjs[idx];
    if (obj.Id() != id)
      continue;
    if (!(obj.Pt() > ptmin))
      continue;
    if (!(std::abs(obj.Eta()) < fetamax))
      continue;
    selected.push_back(idx);
  }
  return selected;
}

std::vector<std::size_t> AnalyzerCore::SelectTrigObjIndices(
    const TrigObjViewCollection &trigobjs, const int id, const float ptmin,
    const float fetamax) const {
  size_t n = trigobjs.size();
  std::vector<std::size_t> seed_indices(n);
  std::iota(seed_indices.begin(), seed_indices.end(), 0);
  return SelectTrigObjIndices(trigobjs, seed_indices, id, ptmin, fetamax);
}

RVec<TrigObj> AnalyzerCore::GetAllTrigObjs() {
  TrigObjViewCollection views = GetAllTrigObjViews();
  return MaterializeTrigObjs(views);
}

bool AnalyzerCore::IsHEMElectron(const Electron &electron) const {
  if (DataYear != 2018)
    return false;

  if (electron.Eta() < -1.25) {
    if ((electron.Phi() < -0.82) && (electron.Phi() > -1.62))
      return true;
  }
  return false;
}

bool AnalyzerCore::IsHEMElectron(const ElectronView &electron) const {
  if (DataYear != 2018)
    return false;

  if (electron.Eta() < -1.25) {
    if ((electron.Phi() < -0.82) && (electron.Phi() > -1.62))
      return true;
  }
  return false;
}

// Gen Matching
void AnalyzerCore::PrintGen(const RVec<Gen> &gens) {
  cout << "===========================================================" << endl;
  // cout << "RunNumber:EventNumber = " << Run << ":" << event << endl;
  cout << "index\tPID\tStatus\tMIdx\tMPID\tStart\tPt\tEta\tPhi\tM" << endl;
  for (unsigned int i = 2; i < gens.size(); i++) {
    Gen gen = gens.at(i);
    RVec<int> history = TrackGenSelfHistory(gen, gens);
    cout << i << "\t" << gen.PID() << "\t" << gen.Status() << "\t"
         << gen.MotherIndex() << "\t" << gens.at(gen.MotherIndex()).PID()
         << "\t" << history[0] << "\t" << fixed << setprecision(2) << gen.Pt()
         << "\t" << gen.Eta() << "\t" << gen.Phi() << "\t" << gen.M() << endl;
  }
}

RVec<int> AnalyzerCore::TrackGenSelfHistory(const Gen &me,
                                            const RVec<Gen> &gens) {
  // returns {index of the first history of the gen,
  //          index of the last history of the gen's mother}
  int myindex = me.Index();
  if (myindex < 2)
    return {myindex, -1};

  int mypid = gens.at(myindex).PID();
  int currentidx = myindex;
  int motherindex = me.MotherIndex();
  if (motherindex < 0)
    return {myindex, -1};

  while (gens.at(motherindex).PID() == mypid) {
    // Go one generation up
    currentidx = motherindex;
    motherindex = gens.at(motherindex).MotherIndex();
    if (motherindex < 0)
      break;
  }
  RVec<int> out = {currentidx, motherindex};
  return out;
}

Gen AnalyzerCore::GetGenMatchedLepton(const Lepton &lep,
                                      const RVec<Gen> &gens) {
  //==== find status 1 lepton
  int reco_PID = -999;
  if (lep.LeptonFlavour() == Lepton::ELECTRON) {
    reco_PID = 11;
  } else if (lep.LeptonFlavour() == Lepton::MUON) {
    reco_PID = 13;
  } else {
    cerr << "[AnalyzerCore::GetGenMatchedLepton] input lepton flavour not set"
         << endl;
    exit(EXIT_FAILURE);
  }

  float min_dR = 0.1;
  Gen gen_closest;
  bool found_match = false;
  for (unsigned int i = 2; i < gens.size(); i++) {
    const Gen &gen = gens.at(i);
    if (gen.Status() != 1)
      continue;
    if (abs(gen.PID()) != reco_PID)
      continue;
    if (gen.MotherIndex() < 0)
      continue;                     // reject ISR
    if (gen.DeltaR(lep) < min_dR) { // dR matching
      min_dR = gen.DeltaR(lep);
      gen_closest = gen;
      found_match = true;
    }
  }
  if (!found_match) {
    // Return a Gen with invalid index if no match found
    Gen dummy;
    dummy.SetIndexPIDStatus(-1, 0, -1);
    return dummy;
  }
  return gen_closest;
}

Gen AnalyzerCore::GetGenMatchedMuon(const Muon &muon, const RVec<Gen> &gens) {
  // Find status 1 muon
  int reco_PID = 13;

  float max_dR = 0.2;
  float min_dR = 0.1;
  Gen gen_closest;
  bool found_match = false;
  float distance = 100000;
  for (const auto &gen : gens) {
    if (gen.Status() != 1)
      continue;
    if (abs(gen.PID()) != reco_PID)
      continue;
    if (gen.MotherIndex() < 0)
      continue; // reject ISR
    if (gen.DeltaR(muon) > max_dR)
      continue;
    float this_distance = pow(gen.DeltaR(muon) / 0.005, 2) +
                          pow((muon.Pt() / gen.Pt() - 1) / 0.02, 2);
    if (this_distance < distance) {
      distance = this_distance;
      gen_closest = gen;
      found_match = true;
    }
  }
  if (!found_match) {
    Gen dummy;
    return dummy;
  }
  return gen_closest;
}

Gen AnalyzerCore::GetGenMatchedPhoton(const Lepton &lep,
                                      const RVec<Gen> &gens) {
  // Find if there is a photon candidate for a source of external conversion
  // (similar pt, direction of lepton) 1) Validity of cuts are checked for
  // electron with PT>10. Note that the previous dPtRel<0.2 cut in the CatNtuple
  // analysis was optimized for PT(e)>25.
  //    External conversion rate proportionate to M^{-2}, thus muon external
  //    conversion rate is negligible.
  // 2) As the algorithm targets stable photon near lepton, only status-1 photon
  // should be investigated theoretically.
  //    However, in some events, last status in photon history is 23. This is
  //    believed to be due to skimming of GEN history between pythia & MiniAOD.
  //    This case is also covered here. Note that IsFinalPhotonSt23 has to
  //    analyze whole history, i.e. computationally heavier part in this
  //    function. Therefore it is put at last part of criteria for efficient
  //    functioning.
  // For more details about the optimization, ask J.B..
  float min_dR = 0.2; // 1)
  Gen gen_closest;
  bool found_match = false;
  float pt_min = 10., dPtRelmax = 0.5, dRmax = 0.2; // 1)
  for (unsigned int i = 2; i < gens.size(); i++) {
    const Gen &gen = gens.at(i);
    if (gen.MotherIndex() < 0)
      continue;
    if (!(abs(gen.PID()) == 22 && (gen.Status() == 1 || gen.Status() == 23)))
      continue;
    if (gen.Pt() < pt_min)
      continue;
    if (!(lep.Pt() / gen.Pt() > (1. - dPtRelmax) &&
          lep.Pt() / gen.Pt() < (1. + dPtRelmax)))
      continue;
    if (gen.DeltaR(lep) > dRmax)
      continue;
    if (gen.Status() == 23 && !IsFinalPhotonSt23_Public(gens))
      continue;

    if (gen.DeltaR(lep) < min_dR) {
      min_dR = gen.DeltaR(lep);
      gen_closest = gen;
      found_match = true;
    }
  }
  if (!found_match) {
    // Return a Gen with invalid index if no match found
    Gen dummy;
    dummy.SetIndexPIDStatus(-1, 0, -1);
    return dummy;
  }
  return gen_closest;
}

bool AnalyzerCore::IsFinalPhotonSt23_Public(const RVec<Gen> &gens) {
  // In Some XG proc events, it is observed that some of photons' last status
  // is 23. Presumably due to skimming of generator history between pythia and
  // MiniAOD. The function returns if this is the case. And this is designed
  // only for 1 hard photon case as W+G or Z+G or TT+G
  bool IsFinalGammaStatus23 = false;
  bool HasStatus23Photon = false;
  for (unsigned int i = 2; i < gens.size(); i++) {
    const Gen &gen = gens.at(i);
    int fpid = fabs(gen.PID());
    int GenSt = gen.Status();
    int MPID_direct = gens.at(gen.MotherIndex()).PID();
    if (!((fpid != 22 && MPID_direct == 22) ||
          (fpid == 22 && (GenSt == 23 || GenSt == 1))))
      continue;

    RVec<int> my_history, mom_history;
    my_history = TrackGenSelfHistory(gen, gens);
    int LastSelfIdx = my_history[0];
    int MotherIdx = my_history[1];
    int LastSelfSt = gens.at(LastSelfIdx).Status();
    int LastSelfMIdx = -1, MStatus_orig = -1;
    if (MotherIdx != -1) {
      mom_history = TrackGenSelfHistory(gens.at(MotherIdx), gens);
      LastSelfMIdx = mom_history[0];
      MStatus_orig = gens.at(LastSelfMIdx).Status();
    }

    if (fpid == 22) {
      if (GenSt == 23) {
        HasStatus23Photon = true;
        IsFinalGammaStatus23 = true;
      } else if (GenSt == 1 && LastSelfSt == 23) {
        IsFinalGammaStatus23 = false;
        break; // a
      }
    } else if (MPID_direct == 22 && MStatus_orig == 23) {
      IsFinalGammaStatus23 = false;
      break; // b
    }
  }
  if (!HasStatus23Photon)
    return false;
  return IsFinalGammaStatus23;
  //**footnotes
  // a) Status-23 photon's last is 1. Thus status-23 photon is not the last
  // history. b) Daughter particle of status-23 photon is found. Thus status-23
  // photon is not the last history.
}

bool AnalyzerCore::IsFromHadron(const Gen &me, const RVec<Gen> &gens) {
  bool out = false;
  int myindex = me.Index();
  if (myindex < 2)
    return true;
  RVec<int> my_history = TrackGenSelfHistory(me, gens);
  if (my_history[1] < 0)
    return true;
  Gen Start = gens.at(my_history[0]);
  Gen MotherOf_Start = gens.at(my_history[1]);

  //==== Status 21~29 are from hardprocess
  //==== Means it's lepton from hardprocess
  //==== e.g., leptons from Z start their lives with status 23
  if (20 < Start.Status() && Start.Status() < 30)
    return false;

  Gen current_me = Start;     // me will always be Start
  Gen current_mother = Start; // initializing
  while (current_mother.Index() >= 2) {
    RVec<int> current_history = TrackGenSelfHistory(current_me, gens);
    if (current_history[1] == -1) { // current_me is the initial partons
      out = true;
      break;
    }
    //==== Go one generation up
    current_me = gens.at(current_history[1]);

    //==== Now look at mother of previous "me"
    current_mother = gens.at(current_history[1]);

    RVec<int> current_mother_history =
        TrackGenSelfHistory(current_mother, gens);
    Gen StartOf_current_mother = gens.at(current_mother_history[0]);
    int current_mother_PID = current_mother.PID();

    //==== From Z,W,H,t
    if (current_mother_PID == 23 || current_mother_PID == 24 ||
        current_mother_PID == 25 || current_mother_PID == 6) {
      out = false;
      break;
    }

    //==== From Signal
    if (IsSignalPID(current_mother_PID)) {
      out = false;
      break;
    }
    if ((current_mother_PID == 11 || current_mother_PID == 13 ||
         current_mother_PID == 15 || current_mother_PID == 22) &&
        (StartOf_current_mother.Status() > 20 &&
         StartOf_current_mother.Status() < 30)) {
      out = false;
      break;
    }
    if (current_mother_PID > 50) {
      out = true;
      break;
    }
    if ((current_mother_PID >= 1 && current_mother_PID <= 5) ||
        current_mother_PID == 21) {
      out = true;
      break;
    }
  }

  return out;
}

bool AnalyzerCore::IsSignalPID(const int &pid) {
  const unsigned int fpid = abs(pid);

  //==== HeavyNeutrino
  if (fpid >= 9900000)
    return true;
  //==== CP-odd Higgs
  if (fpid == 36)
    return true;
  return false;
}

//==== [Type]
//====  1 : EWPrompt
//====  2 : Signal Daughter
//====  3 : EWtau daughter
//====  4 : Internal Conversion daughter from t/EWV/EWlep(Implicit,Explicit)
//====  5 : Internal Conversion daughter from HardScatterPhoton
//====  6 : From Offshell W (mother == 37)
//==== -1 : Unmatched & not EW Conversion candidate
//==== -2 : Hadron daughter
//==== -3 : Daughter of tau from hadron or parton
//==== -4 : Internal conversion daughter(implicit,explicit) having hadronic
// origin
//==== -5 : External conversion candidate(Hard scattered photon)
//==== -6 : External conversion from t/EWV/EWlep
//==== (-4: Daughter of Non-hard scattered photon & has parton or hadron
// ancestor OR implicit Conv from quark)
//==== -------------------------------------------------------------
//====  0 : Error
//==== >0 : Non-fake; Non-hadronic origin
//==== <0 : Fakes   ; Hadronic origin or external conversion

int AnalyzerCore::GetLeptonType(const Gen &gen, const RVec<Gen> &gens) {
  int LeptonType = GetLeptonType_Public(gen.Index(), gens); // 4)
  if (LeptonType >= 4 && LeptonType <= 5) {                 // 5)
    int NearbyPrElType = GetPrElType_InSameSCRange_Public(gen.Index(), gens);
    if (NearbyPrElType > 0)
      LeptonType = NearbyPrElType;
  }
  return LeptonType;
}

int AnalyzerCore::GetLeptonType(const Lepton &lep, const RVec<Gen> &gens) {
  int LeptonType = 0, MatchedTruthIdx = -1;
  Gen gen_closest = GetGenMatchedLepton(lep, gens);
  // cout << "[AnalyzerCore::GetLeptonType] [Reco] pt = " << lep.Pt() << "\t,
  // eta = " << lep.Eta() << endl; cout << "[AnalyzerCore::GetLeptonType] [Gen]
  // Index = " << gen_closest.Index() << endl;

  //==== No matched gen lepton
  if (gen_closest.IsEmpty()) {
    Gen gen_photon_closest = GetGenMatchedPhoton(lep, gens);
    int NearPhotonType = GetGenPhotonType(gen_photon_closest, gens);
    if (NearPhotonType == 1)
      LeptonType = -5;
    else if (NearPhotonType == 2)
      LeptonType = -6;
    else
      LeptonType = -1; // NearPhotonType <= 0
  } else {             //==== Has macthed gen lepton
    MatchedTruthIdx = gen_closest.Index();
    LeptonType = GetLeptonType_Public(MatchedTruthIdx, gens);
    if (LeptonType >= 4 && LeptonType <= 5) {
      int NearbyPrElType =
          GetPrElType_InSameSCRange_Public(MatchedTruthIdx, gens);
      if (NearbyPrElType > 0)
        LeptonType = NearbyPrElType;
    }
  }
  return LeptonType;
  //**footnote
  // 1) matched to no gen-lepton nor gen-photon -> mis-reco.
  //   matched to no gen-lepton, but to photon with hadronic origin -> mis-reco.
  //   (e.g. pions->e) or external conversion from photon in jets
  // 2) matched to no gen-lepton, but to photon with non-hadronic origin (hard
  // process) -> external conversion from photon with non-hadronic source
  // (ME-level) 3) matched to no gen-lepton, but to photon with non-hadronic
  // origin (soft QED radiation) -> external conversion from photon with
  // non-hadronic source (PS-level) 4) matched to gen-lepton, categorize based
  // on the truth categorization algo.: AnalyzerCore::GetLeptonType_Public(int
  // TruthIdx, std::vector<Gen>& TruthColl) 5) collimated e/gm objects are
  // merged in SC, hence if there is prompt electron within SC-merging range,
  // reco-electron's properties actually represent pre-QED-FSR prompt-electron,
  //   rather than the closest internal conversion electron. Therefore shift the
  //   type to the prompt lepton's type.
  //- Note: distinction between type 4 vs. 5 and type -5 vs. -6 is unphysical.
  // it is only for debugging.
}

int AnalyzerCore::GetLeptonType_Public(const int &genIdx,
                                       const RVec<Gen> &gens) {
  // Type : 1:EW-Prompt / 2: BSM-Prompt / 3:EW/BSM-Prompt-Tau Daughter / 6: From
  // Offshell W (mother == 37)
  //        4:Internal Conversion from Soft QED Radiation (PS-level) /
  //        5:Internal Conversion from Hard Process Photon (ME-level)
  //       -1:Unmatched & not EW Conversion candidate (mis-reco. or external
  //       conversion within jets) -2:Hadron Daughter / -3:Daughter of Tau from
  //       Hadron or Parton / -4:Internal Conversion Daughter having hadronic
  //       origin -5:External Conversion from Hard process photon / -6:External
  //       conversion from t/EWV/EWlep
  //       (-4:Daughter of Non-hard scattered photon & has parton or hadron
  //       ancestor OR implicit Conv from quark)
  //        0:Error / >0: Non-fake: Non-hadronic origin / <0 : Fakes: Hadronic
  //        origin or external conversion
  // Note: There is no physical distinction between type 4 vs. 5 and type -5 vs.
  // -6. The distinction is only intended for straightforward debugging.

  // Only consider Status 1 lepton
  if (genIdx < 2)
    return 0;
  if (gens.at(genIdx).Status() != 1)
    return 0;
  if (!(fabs(gens.at(genIdx).PID()) == 11 || fabs(gens.at(genIdx).PID()) == 13))
    return 0;

  int LeptonType = 0;
  int MPID = 0, GrMPID = 0;
  RVec<int> my_history, mom_history, grmom_history;
  int LastSelfIdx = -1, MotherIdx = -1, LastSelfMIdx = -1, GrMotherIdx = -1,
      LastSelfGrMIdx = -1;
  int Status_orig = 0, MStatus_orig = 0, MStatus_last = 0, GrMStatus_orig = 0,
      GrMStatus_last = 0;
  bool HadronicOrigin = false;

  my_history = TrackGenSelfHistory(gens.at(genIdx), gens);
  LastSelfIdx = my_history[0];
  MotherIdx = my_history[1];
  Status_orig = gens.at(LastSelfIdx).Status();
  HadronicOrigin = IsFromHadron(gens.at(genIdx), gens);

  if (MotherIdx != -1) {
    mom_history = TrackGenSelfHistory(gens.at(MotherIdx), gens);
    LastSelfMIdx = mom_history[0];
    GrMotherIdx = mom_history[1];
    MPID = gens.at(MotherIdx).PID();
    MStatus_orig = gens.at(LastSelfMIdx).Status();
    MStatus_last = gens.at(MotherIdx).Status();
  }
  if (GrMotherIdx != -1) {
    grmom_history = TrackGenSelfHistory(gens.at(GrMotherIdx), gens);
    LastSelfGrMIdx = grmom_history[0];
    GrMPID = gens.at(GrMotherIdx).PID();
    GrMStatus_orig = gens.at(LastSelfGrMIdx).Status();
    GrMStatus_last = gens.at(GrMotherIdx).Status();
  }

  if (genIdx == -1)
    LeptonType = 0;
  else if (fabs(MPID) == 23 || fabs(MPID) == 24 || fabs(MPID) == 25)
    LeptonType = 1;
  else if (IsSignalPID(MPID))
    LeptonType = 2;
  else if (abs(MPID) == 37)
    LeptonType = 6;
  else if (Status_orig > 20 && Status_orig < 30)
    LeptonType = 1; // 1)
  else if (fabs(MPID) > 50)
    LeptonType = -2;
  else if (fabs(MPID) == 15 && MStatus_last == 2) {
    if (fabs(GrMPID) == 23 || fabs(GrMPID) == 24 || fabs(GrMPID) == 25)
      LeptonType = 3;
    else if (IsSignalPID(GrMPID))
      LeptonType = 3;
    else if (MStatus_orig > 20 && MStatus_orig < 30)
      LeptonType = 3; // 1)
    else if (HadronicOrigin)
      LeptonType = -3; // 2-a)
    else if (fabs(GrMPID) == 22 && GrMStatus_orig > 20 && GrMStatus_orig < 30)
      LeptonType = 5; // 2-b)
    else if (fabs(GrMPID) == 22)
      LeptonType = 4; // 2-c)
    else if ((fabs(GrMPID) == 11 || fabs(GrMPID) == 13 || fabs(GrMPID) == 15) &&
             GrMStatus_last != 2)
      LeptonType = 4; // 2-d)
    else
      LeptonType = 0;
  } else if (fabs(MPID) == 22) {
    if (MStatus_orig > 20 && MStatus_orig < 30)
      LeptonType = 5; // 3-a)
    else if (HadronicOrigin)
      LeptonType = -4; // 3-b)
    else if (fabs(GrMPID) == 24 || fabs(GrMPID) == 23 || fabs(GrMPID) == 6)
      LeptonType = 4; // 3-c)
    else if (fabs(GrMPID) == 11 || fabs(GrMPID) == 13 || fabs(GrMPID) == 15)
      LeptonType = 4; // 3-d)
    else
      LeptonType = 0;
  } else if ((fabs(MPID) == 11 || fabs(MPID) == 13 || fabs(MPID) == 15) &&
             MStatus_last != 2 && !HadronicOrigin)
    LeptonType = 4; // 4-a)
  else if (((fabs(MPID) >= 1 && fabs(MPID) <= 5) || fabs(MPID) == 21) &&
           MStatus_last != 2)
    LeptonType = -4; // 4-b)
  else if (fabs(MPID) == 6)
    LeptonType = 4; // 4-c)
  else
    LeptonType = 0;

  return LeptonType;
}

//==== [Type]
//====  0 : Invalid input or Error or HardScatter is input when hardscatter is
// not final state
//====  1 : HardScatter / 2: Else prompt daughter(l,V,t)
//==== -1 : Reserved for unmatched(Not used now) / -2: Hadronic origin
int AnalyzerCore::GetGenPhotonType(const Gen &genph, const RVec<Gen> &gens) {
  if (genph.IsEmpty())
    return 0;
  int genph_index = genph.Index();
  if (genph_index < 2)
    return 0;
  if (!(gens.at(genph_index).PID() == 22 &&
        (gens.at(genph_index).Status() == 1 ||
         gens.at(genph_index).Status() == 23)))
    return 0;

  if (gens.at(genph_index).Status() == 23) {
    if (IsFinalPhotonSt23_Public(gens))
      return 1;
    else
      return 0;
  } // From this pt, only St1 Photon is treated.

  RVec<int> phhist = TrackGenSelfHistory(genph, gens);
  const Gen &Start = gens.at(phhist[0]);
  const Gen &MotherOf_Start = gens.at(phhist[1]);
  const int &MotherOf_Start_PID = abs(MotherOf_Start.PID()); // |PID|
  const bool fromhadron = IsFromHadron(genph, gens);

  if (20 < Start.Status() && Start.Status() < 30)
    return 1; // 1)
  else if (MotherOf_Start_PID == 23 || MotherOf_Start_PID == 25)
    return 1; // 2)
  else if (MotherOf_Start_PID == 24 || MotherOf_Start_PID == 6 ||
           IsSignalPID(MotherOf_Start_PID))
    return 2; // 3)
  else if (fromhadron)
    return -2; // 4)
  else if (MotherOf_Start_PID == 11 || MotherOf_Start_PID == 13 ||
           MotherOf_Start_PID == 15)
    return 2; // 5)
  else
    return 0;

  //**footnote
  // 1) Hard process photon's original status is 20's, unless the history is
  // skimmed in MiniAOD. 2) If part of intermediate history is skimmed at
  // MinoAOD, then last photon status is 1 without any preceding history lines.
  //   In that case, intermediate mother is written in history, and this line
  //   catches the case.
  // 3) top and charged bosons can radiate photons.
  // - Note that distinction between 1 and 2 is not physical. You should use
  // both 1 & 2 for prompt photons. 4) this category does not include top, but
  // photons from hadrons and quarks. Predominantly, in most of the cases they
  // are daughter of pi0.
  //   Rarer, but other mesons as eta, B, or even some quarks can also radiate
  //   energetic photons.
  // 5) Photons radiated from lepton FSR. Sometimes they are quite energetic.
}

int AnalyzerCore::GetPrElType_InSameSCRange_Public(int genIdx,
                                                   const RVec<Gen> &gens) {
  // Abbreviation: Get Prompt Electron Type In Same Supercluster Range (Public
  // (shared) version of J.B.'s original gen-matching code) Prompt e>eee (int.
  // conv.) case, collimated electrons can be merged in one SC & track is
  // selected among them, and reconstructed as single electron. In this case,
  // still there will be a nearby prompt electron of LepType 1-3 within
  // supercluster merging range. Supercluster merging range: |dphi|<0.3/0.3
  // (EB/EE), |deta|<0/0.07 (EB/EE) (additionally put 0.03 for marginal
  // difference from size of crystals) Ref:JINST 10 (2015) P06005,
  // arXiv:1502.02701 Return Value: 1/2/3:LeptonType-1/2/3 electron found in
  // same SC range
  //              -1: No LeptonType-1/2/3 electron found in same SC range
  // Note: Hadronic origins are not included in this algo. as I am not sure of
  // effect of many nearby particles in jets on the ele-reco performance.
  //       + it is not of my interest yet.

  if (genIdx < 2)
    return -1;
  if (abs(gens.at(genIdx).PID()) != 11)
    return -1;
  if (gens.at(genIdx).Status() != 1)
    return -1;

  float dPhiMax = 0.3, dEtaMax = 0.1;
  int NearbyElType = -1;

  for (unsigned int it_gen = 2; it_gen < gens.size(); it_gen++) {
    if (gens.at(it_gen).Status() != 1)
      continue;
    if (abs(gens.at(it_gen).PID()) != 11)
      continue;
    if (fabs(gens.at(genIdx).Eta() - gens.at(it_gen).Eta()) > dEtaMax)
      continue;
    if (gens.at(genIdx).DeltaPhi(gens.at(it_gen)) > dPhiMax)
      continue;

    int LepType = GetLeptonType_Public(it_gen, gens);
    if (!(LepType >= 1 && LepType <= 3))
      continue;
    NearbyElType = LepType;
  }
  return NearbyElType;
}

// Histogram Handlers
namespace {
template <typename HistT> struct HistCache {
  static constexpr size_t kSize = 8;
  std::array<std::string, kSize> keys{};
  std::array<HistT *, kSize> hists{};
  size_t count = 0;
  size_t next = 0;

  HistT *Find(std::string_view key) {
    for (size_t i = 0; i < count; ++i) {
      if (keys[i] == key)
        return hists[i];
    }
    return nullptr;
  }

  void Insert(std::string_view key, HistT *hist) {
    if (count < kSize) {
      keys[count] = key;
      hists[count] = hist;
      ++count;
      return;
    }
    keys[next] = key;
    hists[next] = hist;
    next = (next + 1) % kSize;
  }
};

static thread_local HistCache<TH1> cache1d;
static thread_local HistCache<TH2> cache2d;
static thread_local HistCache<TH3> cache3d;
} // namespace

void AnalyzerCore::FillHist(std::string_view histname, float value,
                            float weight, int n_bin, float x_min, float x_max) {
  if (TH1 *cached = cache1d.Find(histname)) {
    cached->Fill(value, weight);
    return;
  }

  auto it = histmap1d.find(histname);
  if (it == histmap1d.end()) {
    std::string key(histname);
    TH1 *this_hist = nullptr;
    if (useTH1F)
      this_hist = new TH1F(key.c_str(), "", n_bin, x_min, x_max);
    else
      this_hist = new TH1D(key.c_str(), "", n_bin, x_min, x_max);

    this_hist->SetDirectory(nullptr);
    it = histmap1d.emplace(std::move(key), this_hist).first;
  }

  cache1d.Insert(histname, it->second);
  it->second->Fill(value, weight);
}

void AnalyzerCore::FillHist(std::string_view histname, float value,
                            float weight, int n_bin, float *xbins) {
  if (TH1 *cached = cache1d.Find(histname)) {
    cached->Fill(value, weight);
    return;
  }

  auto it = histmap1d.find(histname);
  if (it == histmap1d.end()) {
    std::string key(histname);
    TH1 *this_hist = nullptr;
    if (useTH1F)
      this_hist = new TH1F(key.c_str(), "", n_bin, xbins);
    else
      this_hist = new TH1D(key.c_str(), "", n_bin, xbins);

    this_hist->SetDirectory(nullptr);
    it = histmap1d.emplace(std::move(key), this_hist).first;
  }

  cache1d.Insert(histname, it->second);
  it->second->Fill(value, weight);
}

void AnalyzerCore::FillHist(std::string_view histname, float value_x,
                            float value_y, float weight, int n_binx,
                            float x_min, float x_max, int n_biny, float y_min,
                            float y_max) {
  if (TH2 *cached = cache2d.Find(histname)) {
    cached->Fill(value_x, value_y, weight);
    return;
  }

  auto it = histmap2d.find(histname);
  if (it == histmap2d.end()) {
    std::string key(histname);
    TH2 *this_hist = nullptr;
    if (useTH1F)
      this_hist =
          new TH2F(key.c_str(), "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    else
      this_hist =
          new TH2D(key.c_str(), "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    this_hist->SetDirectory(nullptr);
    it = histmap2d.emplace(std::move(key), this_hist).first;
  }

  cache2d.Insert(histname, it->second);
  it->second->Fill(value_x, value_y, weight);
}

void AnalyzerCore::FillHist(std::string_view histname, float value_x,
                            float value_y, float weight, int n_binx,
                            float *xbins, int n_biny, float *ybins) {
  if (TH2 *cached = cache2d.Find(histname)) {
    cached->Fill(value_x, value_y, weight);
    return;
  }

  auto it = histmap2d.find(histname);
  if (it == histmap2d.end()) {
    std::string key(histname);
    TH2 *this_hist = nullptr;
    if (useTH1F)
      this_hist = new TH2F(key.c_str(), "", n_binx, xbins, n_biny, ybins);
    else
      this_hist = new TH2D(key.c_str(), "", n_binx, xbins, n_biny, ybins);
    this_hist->SetDirectory(nullptr);
    it = histmap2d.emplace(std::move(key), this_hist).first;
  }

  cache2d.Insert(histname, it->second);
  it->second->Fill(value_x, value_y, weight);
}

void AnalyzerCore::FillHist(std::string_view histname, float value_x,
                            float value_y, float value_z, float weight,
                            int n_binx, float x_min, float x_max, int n_biny,
                            float y_min, float y_max, int n_binz, float z_min,
                            float z_max) {
  if (TH3 *cached = cache3d.Find(histname)) {
    cached->Fill(value_x, value_y, value_z, weight);
    return;
  }

  auto it = histmap3d.find(histname);
  if (it == histmap3d.end()) {
    std::string key(histname);
    TH3 *this_hist = nullptr;
    if (useTH1F)
      this_hist = new TH3F(key.c_str(), "", n_binx, x_min, x_max, n_biny, y_min,
                           y_max, n_binz, z_min, z_max);
    else
      this_hist = new TH3D(key.c_str(), "", n_binx, x_min, x_max, n_biny, y_min,
                           y_max, n_binz, z_min, z_max);
    this_hist->SetDirectory(nullptr);
    it = histmap3d.emplace(std::move(key), this_hist).first;
  }

  cache3d.Insert(histname, it->second);
  it->second->Fill(value_x, value_y, value_z, weight);
}

void AnalyzerCore::FillHist(std::string_view histname, float value_x,
                            float value_y, float value_z, float weight,
                            int n_binx, float *xbins, int n_biny, float *ybins,
                            int n_binz, float *zbins) {
  if (TH3 *cached = cache3d.Find(histname)) {
    cached->Fill(value_x, value_y, value_z, weight);
    return;
  }

  auto it = histmap3d.find(histname);
  if (it == histmap3d.end()) {
    std::string key(histname);
    TH3 *this_hist = nullptr;
    if (useTH1F)
      this_hist = new TH3F(key.c_str(), "", n_binx, xbins, n_biny, ybins,
                           n_binz, zbins);
    else
      this_hist = new TH3D(key.c_str(), "", n_binx, xbins, n_biny, ybins,
                           n_binz, zbins);
    this_hist->SetDirectory(nullptr);
    it = histmap3d.emplace(std::move(key), this_hist).first;
  }

  cache3d.Insert(histname, it->second);
  it->second->Fill(value_x, value_y, value_z, weight);
}

TTree *AnalyzerCore::NewTree(const TString &treename,
                             const RVec<TString> &keeps,
                             const RVec<TString> &drops) {
  auto treekey = string(treename);
  auto it = treemap.find(treekey);
  if (it == treemap.end()) {
    const bool makeEmptyTree =
        keeps.empty() &&
        (drops.empty() || (drops.size() == 1 && drops[0] == "*"));
    if (makeEmptyTree) {
      TTree *newtree = new TTree(treekey.c_str(), "");
      treemap[treekey] = newtree;
      return newtree;
    } else {
      // check tree is empty.
      if (fChain->GetEntries() == 0) {
        cout << "[AnalyzerCore::NewTree] fChain is empty." << endl;
        exit(0);
      }
      TTree *newtree = fChain->CloneTree(0);
      newtree->SetName(treekey.c_str());
      for (const auto &drop : drops) {
        newtree->SetBranchStatus(drop, 0);
      }
      for (const auto &keep : keeps) {
        newtree->SetBranchStatus(keep, 1);
      }
      treemap[treekey] = newtree;
      unordered_map<string, TBranch *> this_branchmap;
      branchmaps[newtree] = this_branchmap;
      return newtree;
    }
  } else {
    return it->second;
  }
}

TTree *AnalyzerCore::GetTree(const TString &treename) {
  auto treekey = string(treename);
  auto it = treemap.find(treekey);
  if (it == treemap.end()) {
    cout << "[AnalyzerCore::GetTree] Tree " << treename << " not found" << endl;
    exit(ENODATA);
  }
  return it->second;
}

void AnalyzerCore::SetBranch(const TString &treename, const TString &branchname,
                             void *address, const TString &leaflist) {
  try {
    void *this_address = address;
    TTree *tree = GetTree(treename);

    unordered_map<string, TBranch *> *this_branchmap = &branchmaps[tree];
    auto it = this_branchmap->find(string(branchname));

    if (it == this_branchmap->end()) {
      auto br = tree->Branch(branchname, this_address, leaflist);
      this_branchmap->insert({string(branchname), br});
    } else {
      it->second->SetAddress(this_address);
    }
  } catch (int e) {
    cout << "[AnalyzerCore::SetBranch] Error get tree: " << treename.Data()
         << endl;
    exit(e);
  }
}

template void AnalyzerCore::SetBranch_Vector<int>(const TString &,
                                                  const TString &,
                                                  std::vector<int> &);
template void AnalyzerCore::SetBranch_Vector<float>(const TString &,
                                                    const TString &,
                                                    std::vector<float> &);
template void AnalyzerCore::SetBranch_Vector<double>(const TString &,
                                                     const TString &,
                                                     std::vector<double> &);
template void AnalyzerCore::SetBranch_Vector<bool>(const TString &,
                                                   const TString &,
                                                   std::vector<bool> &);

void AnalyzerCore::FillTrees(const TString &treename) {
  if (treename == "") {
    for (const auto &pair : treemap) {
      const string &treename = pair.first;
      TTree *tree = pair.second;
      tree->Fill();
    }
  } else {
    // Convert treeName to std::string for comparison
    std::string treeNameStr(treename.Data());

    auto it = treemap.find(treeNameStr);
    if (it != treemap.end()) {
      // Tree with the given name exists, fill it
      TTree *tree = it->second;
      tree->Fill();
    } else {
      // Handle the case where the treeName is not found in the map
      throw std::runtime_error("[AnalyzerCore::FillTrees] Tree with name '" +
                               treeNameStr + "' not found in treemap.");
    }
  }
}

void AnalyzerCore::WriteHist() {
  const int compression_level = 4;
  const int compression_algorithm = ROOT::RCompressionSetting::EAlgorithm::kLZ4;
  cout << "[AnalyzerCore::WriteHist] Writing histograms to "
       << outfile->GetName() << endl;
  cout << "[AnalyzerCore::WriteHist] Set compression algorithm to LZ4 with "
          "level "
       << compression_level << endl;
  outfile->SetCompressionAlgorithm(compression_algorithm);
  outfile->SetCompressionLevel(compression_level);
  std::vector<std::pair<std::string, TH1 *>> sorted_histograms1d(
      histmap1d.begin(), histmap1d.end());
  std::vector<std::pair<std::string, TH2 *>> sorted_histograms2d(
      histmap2d.begin(), histmap2d.end());
  std::vector<std::pair<std::string, TH3 *>> sorted_histograms3d(
      histmap3d.begin(), histmap3d.end());
  std::sort(
      sorted_histograms1d.begin(), sorted_histograms1d.end(),
      [](const std::pair<std::string, TH1 *> &a,
         const std::pair<std::string, TH1 *> &b) { return a.first < b.first; });
  std::sort(
      sorted_histograms2d.begin(), sorted_histograms2d.end(),
      [](const std::pair<std::string, TH2 *> &a,
         const std::pair<std::string, TH2 *> &b) { return a.first < b.first; });
  std::sort(
      sorted_histograms3d.begin(), sorted_histograms3d.end(),
      [](const std::pair<std::string, TH3 *> &a,
         const std::pair<std::string, TH3 *> &b) { return a.first < b.first; });
  for (const auto &pair : sorted_histograms1d) {
    const string &histname = pair.first;
    TH1 *hist = pair.second;
    cout << "[AnalyzerCore::WriteHist] Writing 1D histogram: " << histname
         << endl;
    // Split the directory and name
    // e.g. "dir1/dir2/histname" -> "dir1/dir2", "histname"
    // e.g. "histname" -> "", "histname"
    size_t last_slash = histname.find_last_of('/');
    string this_prefix, this_name;
    last_slash == string::npos ? this_prefix = ""
                               : this_prefix = histname.substr(0, last_slash);
    last_slash == string::npos ? this_name = histname
                               : this_name = histname.substr(last_slash + 1);

    TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
    if (!this_dir)
      outfile->mkdir(this_prefix.c_str());
    outfile->cd(this_prefix.c_str());
    hist->Write(this_name.c_str());
  }
  for (const auto &pair : sorted_histograms2d) {
    const string &histname = pair.first;
    cout << "[AnalyzerCore::WriteHist] Writing 2D histogram: " << histname
         << endl;
    TH2 *hist = pair.second;
    // Split the directory and name
    // e.g. "dir1/dir2/histname" -> "dir1/dir2", "histname"
    // e.g. "histname" -> "", "histname"
    size_t last_slash = histname.find_last_of('/');
    string this_prefix, this_name;
    last_slash == string::npos ? this_prefix = ""
                               : this_prefix = histname.substr(0, last_slash);
    last_slash == string::npos ? this_name = histname
                               : this_name = histname.substr(last_slash + 1);

    TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
    if (!this_dir)
      outfile->mkdir(this_prefix.c_str());
    outfile->cd(this_prefix.c_str());
    hist->Write(this_name.c_str());
  }
  for (const auto &pair : sorted_histograms3d) {
    const string &histname = pair.first;
    cout << "[AnalyzerCore::WriteHist] Writing 3D histogram: " << histname
         << endl;
    TH3 *hist = pair.second;
    // Split the directory and name
    // e.g. "dir1/dir2/histname" -> "dir1/dir2", "histname"
    // e.g. "histname" -> "", "histname"
    size_t last_slash = histname.find_last_of('/');
    string this_prefix, this_name;
    last_slash == string::npos ? this_prefix = ""
                               : this_prefix = histname.substr(0, last_slash);
    last_slash == string::npos ? this_name = histname
                               : this_name = histname.substr(last_slash + 1);

    TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
    if (!this_dir)
      outfile->mkdir(this_prefix.c_str());
    outfile->cd(this_prefix.c_str());
    hist->Write(this_name.c_str());
  }
  for (const auto &pair : treemap) {
    const string &treename = pair.first;
    cout << "[AnalyzerCore::WriteHist] Writing tree: " << treename << endl;
    TTree *tree = pair.second;

    size_t last_slash = treename.find_last_of('/');
    string this_prefix, this_name;
    last_slash == string::npos ? this_prefix = ""
                               : this_prefix = treename.substr(0, last_slash);
    last_slash == string::npos ? this_name = treename
                               : this_name = treename.substr(last_slash + 1);

    TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
    if (!this_dir)
      outfile->mkdir(this_prefix.c_str());
    outfile->cd(this_prefix.c_str());
    tree->Write(this_name.c_str());
    delete tree;
  }
  cout << "[AnalyzerCore::WriteHist] Writing histograms done" << endl;
  outfile->Close();
}
