import ROOT
import cmsstyle

# 1. Configuration for MC and Data files
# Map each MC file to a specific ROOT color and legend label
# TTbar uses Red shades, Single Top uses Azure/Teal shades

base_dir = "/data9/Users/bhoh/SKNanoOutput/Tutorial_reco_tt/2024"


DATA_FILES = [
    "Muon0_C", "Muon1_C",
    "Muon0_D", "Muon1_D",
    "Muon0_E", "Muon1_E",
    "Muon0_F", "Muon1_F",
    "Muon0_G", "Muon1_G",
    "Muon0_H", "Muon1_H",
    "Muon0_I", "Muon1_I",
]

plotConfigs = {
    '0': {
        'name': DATA_FILES, 
        'legend': 'Data', 
        'color': ROOT.kBlack
    },
    '1': {
        'name': ['TTLJ_powheg', 'TTLL_powheg'], 
        'legend': 't#bar{t}', 
        'color': ROOT.kRed - 4 # or cmsstyle.p10.kRed
    },
    '2': {
        'name': [
            'ST_tW_antitop_Lep', 'ST_tW_top_Lep', 'ST_tW_top_Semilep', 
            'ST_tch_antitop_lep', 'ST_tW_antitop_Semilep', 'ST_tch_top_lep'
        ], 
        'legend': 'Single t', 
        'color': ROOT.kAzure + 1 # or cmsstyle.p10.kBlue
    }
}

HIST_PATHS = [
    "Central/Chi2Cut/had_W_mass_Central",
    "Central/Chi2Cut/had_top_mass_Central",
    "Central/cutflow_Central"
]

reduction = 10 # scale MC hist by this factor

def get_plots_from_files(histo_path):
    """Extract and combine histograms based on plotConfigs."""
    histos = []
    colors = []
    labels = []
    
    htotal = None

    for key, config in sorted(plotConfigs.items()):
        print(f"  > Processing sample group: {config['legend']}")
        names = config['name']
        
        h_combined = None
        
        # Loop through files in the current sample group
        for name in names:
            file_name = f"{base_dir}/{name}.root"
            infile = ROOT.TFile.Open(file_name)
            
            if not infile or infile.IsZombie():
                print(f"Warning: Cannot open {file_name}")
                continue
                
            h_temp = infile.Get(histo_path)
            if not h_temp:
                print(f"Warning: Histogram {histo_path} not found in {file_name}")
                infile.Close()
                continue
            
            # Clone and detach from directory
            if h_combined is None:
                h_combined = h_temp.Clone(f"combined_{key}_{histo_path.replace('/', '_')}")
                h_combined.SetDirectory(0)
            else:
                h_temp_clone = h_temp.Clone()
                h_combined.Add(h_temp_clone)
                
            infile.Close()

        if h_combined:
            histos.append(h_combined)
            colors.append(config['color'])
            labels.append(config['legend'])
            
            # Calculate total MC for uncertainty and ratio (Skip DATA which is key '0')
            if key != '0':
                if htotal is None:
                    htotal = h_combined.Clone(f"htotal_{histo_path.replace('/', '_')}")
                    htotal.SetDirectory(0)
                    htotal.SetNameTitle('Uncertainty', 'Uncertainty')
                else:
                    htotal.Add(h_combined)

    # Note: For strict stat. uncertainty handling on htotal,
    # the Add() function already propagates sum-of-weights-squared if Sumw2() was called during creation.
    return histos, colors, labels, htotal


def draw_stack_with_ratio(histos, colors, labels, htotal_prediction, label, out_name, is_log=False):
    """Plotting function using cmsstyle guidelines."""
    if not histos or not htotal_prediction:
        print("Error: No histograms to draw.")
        return

    cmsstyle.setCMSStyle()
    cmsstyle.getCMSStyle().SetNdivisions(5, "X")

    cmsstyle.SetLumi(110)
    cmsstyle.SetEnergy(13)
    cmsstyle.SetExtraText('Preliminary')

    # Data is at index 0, MC starts from index 1
    ref_hist = histos[0]
    x_axis = ref_hist.GetXaxis()
    x_min = x_axis.GetXmin()
    x_max = x_axis.GetXmax()
    
    # Title from histogram path
    x_title = out_name.replace(".png", "").split("_")[-1]

    # Build THStack directly to prevent AttributeError
    hs = ROOT.THStack("hs", "")
    for h, color in zip(histos[1:], colors[1:]):
        h.SetFillColor(color)
        h.SetLineColor(1) # Black border
        h.SetLineWidth(1)
        h.SetFillStyle(1001)
        h.Scale(reduction) # Scale MC histograms for visual comparison with data    
        hs.Add(h)
    htotal_prediction.Scale(reduction) # Scale total MC histogram for uncertainty and ratio 

    # Y-axis scaling
    y_max = cmsstyle.cmsReturnMaxY(ref_hist)
    y_max = max(y_max, cmsstyle.cmsReturnMaxY(htotal_prediction))
    
    y_min = 0.08 if is_log else 0.0
    y_scale = 100 if is_log else 1.5
    
    # Use cmsDiCanvas for clean Pad split
    c = cmsstyle.cmsDiCanvas(f"c_{out_name}", x_min, x_max, y_min, y_max * y_scale, 
                             0.5, 1.5, x_title, "Events / bin", "Data / MC")

    # Legend Setup (Using AddEntry to prevent AttributeError)
    plotlegend = cmsstyle.cmsLeg(0.50, 0.65, 0.92, 0.88, textSize=0.045, columns=2)  
    
    # Add Data
    plotlegend.AddEntry(histos[0], labels[0], 'pe')
    # Add MC (Reverse order so visual stack matches legend top-to-bottom)
    for h, l in zip(reversed(histos[1:]), reversed(labels[1:])):
        plotlegend.AddEntry(h, l, 'f')
    plotlegend.AddEntry(htotal_prediction, 'Stat. Unc.', 'f')

    # --- Pad 1: Main Plot ---
    c.cd(1)
    if is_log:
        ROOT.gPad.SetLogy(1) 
        
    cmsstyle.cmsObjectDraw(hs, "HIST")
    cmsstyle.GetcmsCanvasHist(ROOT.gPad).GetYaxis().SetMaxDigits(3)
    cmsstyle.cmsObjectDraw(htotal_prediction, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)
    cmsstyle.cmsObjectDraw(histos[0], "E", MarkerStyle=ROOT.kFullCircle)

    # Extra Label (e.g., specific region or cut)
    extraLabel = ROOT.TLatex(0.18, 0.75, label)
    extraLabel.SetNDC()
    cmsstyle.cmsObjectDraw(extraLabel, TextFont=cmsstyle.additionalInfoFont)
    
    # --- Pad 2: Ratio Plot ---
    c.cd(2)
    data_ratio = histos[0].Clone("data_ratio")
    prediction_ratio = htotal_prediction.Clone("prediction_ratio")
    
    for i in range(1, data_ratio.GetNbinsX() + 1):
        mc_val = htotal_prediction.GetBinContent(i)
        mc_err = htotal_prediction.GetBinError(i)
        
        data_val = histos[0].GetBinContent(i)
        data_err = histos[0].GetBinError(i)
        
        if mc_val > 0:
            prediction_ratio.SetBinContent(i, 1.0)
            prediction_ratio.SetBinError(i, mc_err / mc_val)
            
            data_ratio.SetBinContent(i, data_val / mc_val)
            data_ratio.SetBinError(i, data_err / mc_val)
        else:
            prediction_ratio.SetBinContent(i, 1.0)
            prediction_ratio.SetBinError(i, 0.0)
            data_ratio.SetBinContent(i, 0.0)
            data_ratio.SetBinError(i, 0.0)

    cmsstyle.cmsObjectDraw(prediction_ratio, "E2", FillStyle=3345, LineWidth=0, FillColor=12, MarkerSize=0)
    cmsstyle.cmsObjectDraw(data_ratio, "E", MarkerStyle=ROOT.kFullCircle)

    cmsstyle.UpdatePad(c)
    c.SaveAs(out_name)


if __name__ == "__main__":
    print("Starting Histogram Processing...")
    
    for path in HIST_PATHS:
        print(f"\n--- Processing {path} ---")
        histos, colors, labels, htotal = get_plots_from_files(path)
        
        hist_name = path.split("/")[-1]
        out_name = f"{hist_name}_plot.png"
        
        # cutflow의 경우 y축 값의 차이가 크므로 Log scale을 적용하는 것이 좋습니다.
        is_log_scale = True if "cutflow" in hist_name else False
        
        draw_stack_with_ratio(
            histos, colors, labels, htotal, 
            label="Single Muon", 
            out_name=out_name,
            is_log=is_log_scale
        )
        print(f"Saved successfully: {out_name}")