import ROOT
import numpy as np
import uproot3
import array
"""
# Pythia
# weight = sig / tries
ptWeightList=[
[9,  1.107e+01,    73952],
[12, 4.265e+00,    72896],
[16, 1.589e+00,    72717],
[21, 5.982e-01,    73571],
[28, 2.097e-01,    75433],
[36, 8.096e-02,    82621],
[45, 3.451e-02,    81655],
[57, 1.344e-02,    78917],
[70, 5.817e-03,    78952],
[85, 2.593e-03,    78632],
[99, 1.365e-03,    78467],
[115, 7.206e-04,    77883],
[132, 3.912e-04,    78854]]
"""
# Herwig
ptWeightList=[
[1,  0.6845e+06,  1],
[2, 0.3442e+06,   1],
[3, 0.13859e+06,  1],
[4, 60.75e+03,    1],
[5, 21.21e+03,    1],
[6, 8.296e+03,    1],
[7, 3.933e+03,    1],
[8, 1.531e+03,    1],
[9, 0.6995e+03,   1],
[10, 0.2842e+03,  1],
[11, 0.1571e+03,  1],
[12, 82.96e+00,   1],
[13, 45.93e+00,   1]]

def makePartialResponse(pt_min,weight,job):
    #jets = uproot3.open("results_Pythia/Analysis_R_0_"+R+"_PT_"+str(pt_min)+"_job_"+str(job)+".root")["jetTree"]
    jets = uproot3.open("results_Herwig/AnalysisHerwig_R_0_"+R+"_PT_"+str(pt_min)+"_job_"+str(job)+".root")["jetTree"]

    cutJets = ((jets.array("sigJet_Detector_smallerPt").flatten()!=0) & (jets.array("sigJet_Detector_biggerPt").flatten()!=0) & (jets.array("sigJet_Truth_biggerPt").flatten()!=0))
    cutEtaLower = ((jets.array("sigJet_Truth_smallerEta").flatten()>-.5) & (jets.array("sigJet_Truth_biggerEta").flatten()>-.5) & (jets.array("sigJet_Detector_smallerEta").flatten()>-.5) & (jets.array("sigJet_Detector_biggerEta").flatten()>-.5))
    cutEtaHigher = ((jets.array("sigJet_Truth_smallerEta").flatten()<.5) & (jets.array("sigJet_Truth_biggerEta").flatten()<.5) & (jets.array("sigJet_Detector_smallerEta").flatten()<.5) & (jets.array("sigJet_Detector_biggerEta").flatten()<.5))
    
    #cutOutlier = ((jets.array("sigJet_Truth_smallerPt").flatten()<3*pt_min) & (jets.array("sigJet_Truth_biggerPt").flatten()<3*pt_min))
    pt_cut = [9,12,16,21,28,36,45,57,70,85,99,115,132]
    cutOutlier = ((jets.array("sigJet_Truth_smallerPt").flatten()<3*pt_cut[pt_min]) & (jets.array("sigJet_Truth_biggerPt").flatten()<3*pt_cut[pt_min]))

    cutJets = cutJets & cutEtaLower & cutEtaHigher & cutOutlier

    sigJet_Truth_smallerPt  = jets.array("sigJet_Truth_smallerPt").flatten()[cutJets]
    sigJet_Detector_smallerPt = jets.array("sigJet_Detector_smallerPt").flatten()[cutJets]

    sigJet_Truth_biggerPt   = jets.array("sigJet_Truth_biggerPt").flatten()[cutJets]       
    sigJet_Detector_biggerPt = jets.array("sigJet_Detector_biggerPt").flatten()[cutJets]

    dPt_Truth       = sigJet_Truth_biggerPt-sigJet_Truth_smallerPt
    dPt_Detector    = sigJet_Detector_biggerPt - sigJet_Detector_smallerPt 

    sigJet_Truth_smallerEta  = jets.array("sigJet_Truth_smallerEta").flatten()[cutJets]
    sigJet_Truth_biggerEta  = jets.array("sigJet_Truth_biggerEta").flatten()[cutJets]
    sigJet_Detector_smallerEta  = jets.array("sigJet_Detector_smallerEta").flatten()[cutJets]
    sigJet_Detector_biggerEta  = jets.array("sigJet_Detector_biggerEta").flatten()[cutJets]

    for i in range(0,len(dPt_Truth)-1):
        h_pt.Fill(sigJet_Detector_smallerPt[i],sigJet_Truth_smallerPt[i],weight)
        h_pt_Bigger.Fill(sigJet_Detector_biggerPt[i],sigJet_Truth_biggerPt[i],weight)
        h_4D.Fill(array.array('d',[sigJet_Truth_smallerPt[i],sigJet_Detector_smallerPt[i],dPt_Truth[i],dPt_Detector[i]]),weight)
        if (sigJet_Truth_smallerPt[i] > 40 and sigJet_Truth_smallerPt[i] < 60):
            h_Dpt.Fill(dPt_Detector[i],dPt_Truth[i],weight)

        if (dPt_Truth[i] != 0):
            xDh_pt.Fill((dPt_Truth[i]-dPt_Detector[i])/dPt_Truth[i],sigJet_Truth_smallerPt[i],weight)
        
        xh_pt.Fill((sigJet_Truth_smallerPt[i]-sigJet_Detector_smallerPt[i])/sigJet_Truth_smallerPt[i],sigJet_Truth_smallerPt[i],weight)
        xh_pt_Bigger.Fill((sigJet_Truth_biggerPt[i]-sigJet_Detector_biggerPt[i])/sigJet_Truth_biggerPt[i],sigJet_Truth_biggerPt[i],weight)

        xh_eta.Fill(sigJet_Truth_smallerEta[i])
        xh_eta.Fill(sigJet_Truth_biggerEta[i])
        xh_eta.Fill(sigJet_Detector_smallerEta[i])
        xh_eta.Fill(sigJet_Detector_biggerEta[i])

#for R in (["05","10","15","20","25","30","35"]):
for R in (["25"]):
    #f = ROOT.TFile(R+"_responsePythia.root", "RECREATE")
    f = ROOT.TFile(R+"_responseHerwig.root", "RECREATE")

    h_pt = ROOT.TH2D("h_pt_R_"+R, "2D jet pT response for R = 0."+R, 40, 20, 60, 50, 20, 70)
    h_pt.GetXaxis().SetTitle("Detector pT")
    h_pt.GetYaxis().SetTitle("Truth pT")
    h_pt.Sumw2()

    h_pt_Bigger = ROOT.TH2D("Bigger_h_pt_R_"+R, "2D jet pT response for R bigger", 40, 20, 60, 50, 20, 70)
    h_pt_Bigger.GetXaxis().SetTitle("Detector pT")
    h_pt_Bigger.GetYaxis().SetTitle("Truth pT")    
    h_pt_Bigger.Sumw2()

    h_Dpt = ROOT.TH2D("h_Dpt_R_"+R, "pT response for R = 0."+R, 30, 0, 30, 35, 0, 35)
    h_Dpt.GetXaxis().SetTitle("Detector DpT")
    h_Dpt.GetYaxis().SetTitle("Truth DpT") 
    h_Dpt.Sumw2()

    xh_pt = ROOT.TH2D("xh_pt_R_"+R, " DpT response for R = 0."+R, 100, -1, 1,40,0,200)
    xh_pt.GetXaxis().SetTitle("Truth pT - Detector pT / Truth pT")       
    xh_pt.Sumw2()

    xh_pt_Bigger = ROOT.TH2D("xh_pt_BiggerR_"+R, " DpT response for Bigger R = 0."+R, 100, -1, 1,40,0,200)
    xh_pt_Bigger.GetXaxis().SetTitle("Truth pT - Detector pT / Truth pT")       
    xh_pt_Bigger.Sumw2()

    xDh_pt = ROOT.TH2D("xh_Dpt_R_"+R, " DpT response for R = 0."+R, 100, -1, 1,40,0,200)
    xDh_pt.GetXaxis().SetTitle("Truth DpT - Detector DpT / Truth DpT")       
    xDh_pt.Sumw2()

    h_4D = ROOT.THnD("h_4D_R_"+R, "4D jet response for R = 0."+R+"_", 4, array.array('i', [40, 40, 100, 100]), array.array('d', [0, 0, -20, -20]), array.array('d', [200, 200, 80, 80]))
    h_4D.GetAxis(0).SetTitle("Truth pT")
    h_4D.GetAxis(1).SetTitle("Detector pT")
    h_4D.GetAxis(2).SetTitle("Truth DpT")
    h_4D.GetAxis(3).SetTitle("Detector DpT")
    h_4D.Sumw2()  

    xh_eta = ROOT.TH1D("xh_eta_R_"+R, "eta R = 0."+R, 100, -1, 1)
    xh_eta.Sumw2()  

    for i in range(0,len(ptWeightList)-1):
        for job in range(1,301):
            makePartialResponse(ptWeightList[i][0],ptWeightList[i][1]/ptWeightList[i][2],job)

    h_pt.Write()
    h_pt_Bigger.Write()
    h_Dpt.Write()
    h_4D.Write()

    xh_pt.Write()
    xh_pt_Bigger.Write()
    xDh_pt.Write()

    xh_eta.Write()

f.Close()
