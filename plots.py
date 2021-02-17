import ROOT
ROOT.EnableImplicitMT() # Tell ROOT you want to go parallel
RDF = ROOT.ROOT.RDataFrame


########################## Importing datasets ##########################
rootFile_1 = ROOT.TFile("JetToyHIResultConstituentSubtraction.root")
rootTree_1 = rootFile_1.Get("jetTree")
RDF_1 = RDF(rootTree_1)



########################## Making histograms ##########################
H1 = RDF_1.Histo1D(( "Transverse momentum", "Signal pT",  100, 25,250) ,"sigJetPt") # These are the histograms of variable bin size
H1.SetStats(0)
H1.SetLineColor(1)
H1.SetLineWidth(1)
#H1.Scale(1/H1.GetEntries())
H1.GetYaxis().SetTitle("counts")
H1.GetXaxis().SetTitle("Jet transverse momentum in GeV")
H3 = RDF_1.Histo1D(( "Transverse momentum", "raw pT",  100, 25,250) ,"rawJetPt") # These are the histograms of variable bin size


H2 = RDF_1.Histo1D(( "Mass", "Signal M",  100, 0,30) ,"sigJetM") # These are the histograms of variable bin size
#H2.SetStats(0)
H2.SetLineColor(2)
H2.SetLineWidth(1)
H2.GetYaxis().SetTitle("counts")
H2.GetXaxis().SetTitle("Jet mass in GeV")
H4 = RDF_1.Histo1D(( "Transverse momentum", "raw pT",  100, 25,250) ,"rawJetPt") # These are the histograms of variable bin size

c1 = ROOT.TCanvas()
H1.Draw()
H3.Draw("same")
c1.SaveAs("~/cernbox/plots/sigJetPt.png")

c2 = ROOT.TCanvas()
H2.Draw()
H4.Draw("same")
c2.SaveAs("~/cernbox/plots/sigJetM.png")



########################## Making legend ##########################
"""
legend = ROOT.TLegend (0.7 ,0.6 ,0.85 ,0.75)

legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetLineWidth (0)

legend.AddEntry(H1.DrawCopy(),"a0r025")
legend.AddEntry(H2.DrawCopy(),"a1r025")

legend.Draw("same")
"""