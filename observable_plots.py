import ROOT
ROOT.EnableImplicitMT() # Tell ROOT you want to go parallel
RDF = ROOT.ROOT.RDataFrame


########################## Importing datasets ##########################
rootFile_1 = ROOT.TFile("../JetToyHIResultConstituentSubtraction.root")
rootTree_1 = rootFile_1.Get("jetTree")
RDF_1 = RDF(rootTree_1)



########################## Making histograms ##########################
H1 = RDF_1.Histo1D(( "sigJetPt", "Signal pT",  100, 0,500) ,"sigJetPt") # These are the histograms of variable bin size
H1.SetStats(0)
H1.SetLineColor(2)
H1.SetLineWidth(2)
#H1.Scale(1/H1.GetEntries())
H1.GetYaxis().SetTitle("counts")
H1.GetXaxis().SetTitle("Jet transverse momentum in GeV")
H3 = RDF_1.Histo1D(( "rawJetPt", "raw pT",  100, 0,500) ,"rawJetPt") # These are the histograms of variable bin size
H3.SetLineWidth(2)

H2 = RDF_1.Histo1D(( "sigJetM", "Signal M",  100, 0,80) ,"sigJetM") # These are the histograms of variable bin size
H2.SetStats(0)
H2.SetLineColor(2)
H2.SetLineWidth(2)
H2.GetYaxis().SetTitle("counts")
H2.GetXaxis().SetTitle("Jet mass in GeV")
H4 = RDF_1.Histo1D(( "rawJetM", "raw pT",  100, 0,80) ,"rawJetM") # These are the histograms of variable bin size
H4.SetLineWidth(2)

c1 = ROOT.TCanvas()
H1.Draw()
H3.Draw("same")
legend1 = ROOT.TLegend (0.55 ,0.6 ,0.85 ,0.75)
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetFillStyle(0)
legend1.SetTextFont(42)
legend1.SetTextSize(0.035)
legend1.SetLineWidth (0)
legend1.AddEntry("sigJetPt","Signal, mean = {}".format("%0.2f" % H1.GetMean()),"l")
legend1.AddEntry("rawJetPt","Raw, mean = {}".format("%0.2f" % H3.GetMean()),"l")
legend1.Draw("same")
c1.SaveAs("~/cernbox/plots/sigJetPt.png")

c2 = ROOT.TCanvas()
H2.Draw()
H4.Draw("same")
legend2 = ROOT.TLegend (0.55 ,0.6 ,0.85 ,0.75)
legend2.SetBorderSize(0)
legend2.SetFillColor(0)
legend2.SetFillStyle(0)
legend2.SetTextFont(42)
legend2.SetTextSize(0.035)
legend2.SetLineWidth (0)
legend2.AddEntry("sigJetPt","Signal, mean = {}".format("%0.2f" % H2.GetMean()),"l")
legend2.AddEntry("rawJetPt","Raw, mean = {}".format("%0.2f" % H4.GetMean()),"l")
legend2.Draw("same")
c2.SaveAs("~/cernbox/plots/sigJetM.png")



########################## Making legend ##########################


