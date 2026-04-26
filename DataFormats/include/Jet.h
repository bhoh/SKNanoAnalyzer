#ifndef Jet_h
#define Jet_h

#include <array>
#include <bitset>
#include <memory>
#include <limits>

#include "Particle.h"
#include "JetTaggingParameter.h"
#include "JetConstituent.h"
#include "JetView.h"

class AnalyzerCore;

class Jet : public Particle
{
public:
  Jet();
  Jet(std::shared_ptr<const JetSoA> storage, std::size_t index);
  ~Jet();

  enum class Property {
    DeepFlavB,
    DeepFlavCvB,
    DeepFlavCvL,
    DeepFlavQG,
    PNetB,
    PNetCvB,
    PNetCvL,
    PNetCvNotB,
    PNetQvG,
    PNetTauVJet,
    ParTB,
    ParTCvB,
    ParTCvL,
    ParTCvNotB,
    ParTEle,
    ParTMu,
    ParTQvG,
    ParTSvCB,
    ParTSvUDG,
    ParTTauVJet,
    ParTUDG,
    ParTProbB,
    ParTProbBB,
    Count
  };

  using EnsureCallback = void (*)(void*, Jet&, Property);
  void AttachLazyPayload(void *context, EnsureCallback callback, int index) const;
  void DetachLazyPayload() const;
  bool HasLazyPayload() const { return static_cast<bool>(lazy_); }

  enum class JetID
  {
    NOCUT,
    LOOSE,
    TIGHT,
    TIGHTLEPVETO,
    PUID_LOOSE,
    PUID_MEDIUM,
    PUID_TIGHT,
  };
  inline void SetRawPt(float pt) { jet_rawPt = pt; };
  inline void SetOriginalPt(float pt) { jet_originalPt = pt; };
  float GetRawPt() const { return jet_rawPt; };
  float GetOriginalPt() const { return jet_originalPt; };

  inline void SetArea(double area) { j_area = area; };
  inline void SetJetFlavours(short pf, unsigned char hf)
  {
    j_partonFlavour = static_cast<short>(pf);
    j_hadronFlavour = static_cast<short>(hf);
  };

  void SetTaggerScore(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerScoreType scoreType, float value);
  inline void SetTaggerResults(RVec<float> deepjet, RVec<float> pnet, RVec<float> parT)
  {
    SetTaggerScore(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerScoreType::B, deepjet[0]);
    SetTaggerScore(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerScoreType::CvB, deepjet[1]);
    SetTaggerScore(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerScoreType::CvL, deepjet[2]);
    SetTaggerScore(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerScoreType::QvG, deepjet[3]);

    SetTaggerScore(JetTagging::JetFlavTagger::ParticleNet, JetTagging::JetFlavTaggerScoreType::B, pnet[0]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParticleNet, JetTagging::JetFlavTaggerScoreType::CvB, pnet[1]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParticleNet, JetTagging::JetFlavTaggerScoreType::CvL, pnet[2]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParticleNet, JetTagging::JetFlavTaggerScoreType::CvNotB, pnet[3]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParticleNet, JetTagging::JetFlavTaggerScoreType::QvG, pnet[4]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParticleNet, JetTagging::JetFlavTaggerScoreType::TauVJet, pnet[5]);

    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::B, parT[0]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::CvB, parT[1]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::CvL, parT[2]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::CvNotB, parT[3]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::Ele, parT[4]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::Mu, parT[5]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::QvG, parT[6]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::SvCB, parT[7]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::SvUDG, parT[8]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::TauVJet, parT[9]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::probUDG, parT[10]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::probB, parT[11]);
    SetTaggerScore(JetTagging::JetFlavTagger::ParT, JetTagging::JetFlavTaggerScoreType::probBB, parT[12]);
  };
  inline void SetEnergyFractions(float cH, float nH, float nEM, float cEM, float muE)
  {
    j_chHEF = cH;
    j_neHEF = nH;
    j_neEmEF = nEM;
    j_chEmEF = cEM;
    j_muEF = muE;
  };
  inline float chHEF() const { return j_chHEF; }
  inline float neHEF() const { return j_neHEF; }
  inline float neEmEF() const { return j_neEmEF; }
  inline float chEmEF() const { return j_chEmEF; }
  inline float muEF() const { return j_muEF; }

  inline void SetMultiplicities(unsigned char nC, unsigned char nEl, unsigned char nM, unsigned char nSV)
  {
    j_nConstituents = short(nC);
    j_nElectrons = short(nEl);
    j_nMuons = short(nM);
    j_nSVs = short(nSV);
  };
  inline short nConstituents() const { return j_nConstituents; }
  inline short nElectrons() const { return j_nElectrons; }
  inline short nMuons() const { return j_nMuons; }
  inline short nSVs() const { return j_nSVs; }

  inline void SetHadronMultiplicities(unsigned char chMult, unsigned char neMult)
  {
    j_chMultiplicity = chMult;
    j_neMultiplicity = neMult;
  };
  inline unsigned char chMultiplicity() const { return j_chMultiplicity; }
  inline unsigned char neMultiplicity() const { return j_neMultiplicity; }

  inline void SetMatchingIndices(short e1, short e2, short m1, short m2, short sv1, short sv2, short gj)
  {
    j_electronIdx1 = e1;
    j_electronIdx2 = e2;
    j_muonIdx1 = m1;
    j_muonIdx2 = m2;
    j_svIdx1 = sv1;
    j_svIdx2 = sv2;
    j_genJetIdx = gj;
  };

  inline void SetMatchingIndices(short e1, short e2, short m1, short m2, short sv1, short sv2)
  {
    j_electronIdx1 = e1;
    j_electronIdx2 = e2;
    j_muonIdx1 = m1;
    j_muonIdx2 = m2;
    j_svIdx1 = sv1;
    j_svIdx2 = sv2;
  };

  inline void SetOriginalIndex(int idx) { j_originalIndex = idx; };
  inline int OriginalIndex() const { return j_originalIndex; };

  inline void SetJetPuIDScore(float puIDScore)
  {
    j_puIDScore = puIDScore;
  };

  inline void SetCorrections(RVec<float> corrs)
  {
    j_PNetRegPtRawCorr = corrs[0];
    j_PNetRegPtRawCorrNeutrino = corrs[1];
    j_PNetRegPtRawRes = corrs[2];
    j_UParTAK4RegPtRawCorr = corrs[3];
    j_UParTAK4RegPtRawCorrNeutrino = corrs[4];
    j_UParTAK4RegPtRawRes = corrs[5];
    j_UParTAK4V1RegPtRawCorr = corrs[6];
    j_UParTAK4V1RegPtRawCorrNeutrino = corrs[7];
    j_UParTAK4V1RegPtRawRes = corrs[8];
    j_rawFactor = corrs[9];
  };

  inline void SetM(double jet_m)
  {
    j_m = jet_m;
  };

  inline void SetUnsmearedP4(Jet jet)
  {
    j_unsmearedP4 = jet;
  };
  inline double GetM() { return j_m; }
  inline int partonFlavour() const { return j_partonFlavour; };
  inline int hadronFlavour() const { return j_hadronFlavour; };
  inline int genJetIdx() const { return j_genJetIdx; };
  inline float chargedEMFraction() const { return j_chEmEF; }
  inline float chargedHadronFraction() const { return j_chHEF; }
  inline float neutralEMFraction() const { return j_neEmEF; }
  inline float neutralHadronFraction() const { return j_neHEF; }
  inline float EMFraction() const { return j_chEmEF + j_neEmEF; }
  void SetPFConstituents(RVec<JetConstituent> constituents)
  {
    j_constituents = std::move(constituents);
  }

  const RVec<JetConstituent> &PFConstituents() const noexcept
  {
    return j_constituents;
  }
  float GetTaggerResult(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerScoreType) const;
  TLorentzVector GetUnsmearedP4() const;

  float UParTAK4RegPtRawCorr() const { return j_UParTAK4RegPtRawCorr; }
  float UParTAK4RegPtRawCorrNeutrino() const { return j_UParTAK4RegPtRawCorrNeutrino; }
  float PuIdDisc() const { return j_puIDScore; }

private:
  void ensure(Property property) const;
  void markLoaded(Property property) const;
  void ensureScore(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerScoreType scoreType) const;
  static Property propertyFor(JetTagging::JetFlavTagger tagger, JetTagging::JetFlavTaggerScoreType scoreType);

  void initializeMembers();
  void materialize() const;
  void loadFromStorage(Property property) const;

  struct LazyPayload {
    void *context = nullptr;
    EnsureCallback callback = nullptr;
    int index = -1;
    mutable std::bitset<static_cast<std::size_t>(Property::Count)> loaded;
  };

  static constexpr float kInvalidTagScore = -1.f;

  std::shared_ptr<const JetSoA> storage_;
  std::size_t index_ = std::numeric_limits<std::size_t>::max();

  mutable std::bitset<static_cast<std::size_t>(Property::Count)> cachedProperties_;
  mutable std::shared_ptr<LazyPayload> lazy_;

  float jet_rawPt;
  float jet_originalPt;
  int j_originalIndex;
  float j_puIDScore;
  float j_PNetRegPtRawCorr;
  float j_PNetRegPtRawCorrNeutrino;
  float j_PNetRegPtRawRes;
  float j_UParTAK4RegPtRawCorr;
  float j_UParTAK4RegPtRawCorrNeutrino;
  float j_UParTAK4RegPtRawRes;
  float j_UParTAK4V1RegPtRawCorr;
  float j_UParTAK4V1RegPtRawCorrNeutrino;
  float j_UParTAK4V1RegPtRawRes;
  float j_rawFactor;
  std::array<float, JetTagging::JetFlavTaggerScoreCount> j_btagDeepFlav;
  std::array<float, JetTagging::JetFlavTaggerScoreCount> j_btagPNet;
  std::array<float, JetTagging::JetFlavTaggerScoreCount> j_btagUParTAK4;
  float j_chEmEF;
  float j_chHEF;
  float j_neEmEF;
  float j_neHEF;
  float j_muEF;
  short j_nConstituents;
  short j_nElectrons;
  short j_nMuons;
  short j_nSVs;
  short j_chMultiplicity;
  short j_neMultiplicity;
  short j_electronIdx1;
  short j_electronIdx2;
  short j_muonIdx1;
  short j_muonIdx2;
  short j_svIdx1;
  short j_svIdx2;
  short j_genJetIdx;
  short j_hadronFlavour;
  short j_partonFlavour;
  float j_area;
  float j_m;
  TLorentzVector j_unsmearedP4;
  RVec<JetConstituent> j_constituents;
  friend class AnalyzerCore;
  ClassDef(Jet, 1)
};

#endif
