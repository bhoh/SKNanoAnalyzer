from ROOT import TString
from ROOT.VecOps import RVec, Sort, Reverse
from ROOT import ExampleRun_base
from ROOT import MyCorrection
from ROOT import Event, Muon, Particle

class ExampleRun(ExampleRun_base):
    def __init__(self):
        super(ExampleRun, self).__init__()

    def initializePyAnalyzer(self):
        self.initializeAnalyzer()

    def executeEvent(self):
        # Fetch all NANOAOD muons once per event to save CPU time
        self.AllMuons = self.GetAllMuons()
        
        # No prefire weight for Run3?
        self.weight_Prefire = 1.0  
        
        # Loop over systematic sources that require separate event loop
        for syst_dummy in self.systHelper:
            self.executeEventFromParameter()

        # Using new PDF
        if self.RunNewPDF and not self.IsDATA:
            self.FillHist("NewPDF_PDFReweight", self.GetPDFReweight(), 1., 2000, 0.90, 1.10)
            for i in range(self.pdfReweight.NErrorSet):
                hist_name = "NewPDF_PDFErrorSet/PDFReweight_Member_{}".format(i)
                self.FillHist(hist_name, self.GetPDFReweight(i), 1., 2000, 0.90, 1.10)

    def executeEventFromParameter(self):
        this_syst = str(self.systHelper.getCurrentSysName())
        
        if this_syst == "Central":
            this_muon_id = self.MuonIDs.at(0)
            this_muon_id_sf_key = self.MuonIDSFKeys.at(0)
        elif this_syst == "Muon_ID_Tight_Up":
            this_muon_id = self.MuonIDs.at(1)
            this_muon_id_sf_key = self.MuonIDSFKeys.at(1)
        else:
            return

        # No cut
        self.FillHist(this_syst + "/NoCut", 0., 1., 1, 0., 1.)

        ev = self.GetEvent()
        if not ev.PassTrigger(self.IsoMuTriggerName):
            return

        # Apply ID selections
        muons = self.SelectMuons(self.AllMuons, this_muon_id, 20., 2.4)
        
        # Sort in pt-order
        muons = Reverse(Sort(muons))

        # Event selection: dimuon
        if muons.size() != 2:
            return
        # Leading muon trigger-safe cut
        if muons.at(0).Pt() <= self.TriggerSafePtCut:
            return
        # On-Z mass window
        ZCand = muons.at(0) + muons.at(1)
        if not (abs(ZCand.M() - 91.2) < 15.):
            return

        # Event weight
        weight = 1.0
        if not self.IsDATA:
            weight *= self.MCweight()
            weight *= ev.GetTriggerLumi("Full")
            weight *= self.weight_Prefire
            
            default_weight = weight
            weight_map = self.systHelper.calculateWeight()
            
            for w_pair in weight_map:
                hist_path = "{}/ZCand_Mass_{}".format(this_syst, str(w_pair.first))
                self.FillHist(hist_path, ZCand.M(), default_weight * w_pair.second, 50, 70., 110.)
        else:
            hist_path = "{}/ZCand_Mass_{}".format(this_syst, this_syst)
            self.FillHist(hist_path, ZCand.M(), weight, 50, 70., 110.)


if __name__ == "__main__":
    module = ExampleRun()
    module.SetTreeName("Events")
    module.LogEvery = 5000
    module.IsDATA = False
    module.MCSample = "DYJets"
    module.xsec = 6077.22 
    module.sumW = 1.0 
    module.sumSign = 1.0
    module.SetEra("2022")
    module.Userflags = RVec(TString)(["RunSyst"])
    module.AddFile("NANOAOD_EXAMPLE.root")
    module.MaxEvent = max(1, int(module.fChain.GetEntries()))
    module.SetOutfilePath("hist.root")
    module.Init()
    module.initializePyAnalyzer()
    module.Loop()
    module.WriteHist()