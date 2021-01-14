import ROOT
ROOT.EnableImplicitMT() # Tell ROOT you want to go parallel
RDF = ROOT.ROOT.RDataFrame


########################## Importing datasets ##########################
rootFile_1 = ROOT.TFile("~/cernbox/a0r025ghA005.root")
rootTree_1 = rootFile_1.Get("jetTree")
RDF_1 = RDF(rootTree_1)

rootFile_2 = ROOT.TFile("~/cernbox/a1r025ghA005.root")
rootTree_2 = rootFile_2.Get("jetTree")
RDF_2 = RDF(rootTree_2)


########################## Making histograms ##########################
H1 = RDF_1.Histo1D(( "H1", "(Reconstructed - Signal) / Signal",  100, -1,1) ,"ptPull") # These are the histograms of variable bin size
H1.SetStats(0)
H1.SetLineColor(1)
H1.SetLineWidth(1)

H2 = RDF_2.Histo1D(( "H1", "(Reconstructed - Signal) / Signal",  100, -1,1) ,"ptPull") # These are the histograms of variable bin size
H2.SetStats(0)
H2.SetLineColor(2)
H2.SetLineWidth(1)


c1 = ROOT.TCanvas()

H1.Draw("same")
H2.Draw("same")

#c1.SaveAs("df002_trN.png")

########################## Making legend ##########################
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