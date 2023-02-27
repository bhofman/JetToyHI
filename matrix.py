import ROOT
import numpy as np
import uproot3
import array

#weight = sig / tries
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
[115, 1.365e-03,    0],
[132, 1.365e-03,    0],
[200,0,0]]

#for R in (["05","10","15","20","25","30","35"]):
for R in (["05"]):
    h_pt = ROOT.TH2D("h_pt_R_"+R, "2D jet pT response for R = 0."+R, 40, 20, 60, 50, 20, 70)
    h_pt.GetXaxis().SetTitle("Detector pT")
    h_pt.GetYaxis().SetTitle("Truth pT")

    h_pt_Bigger = ROOT.TH2D("Bigger_h_pt_R_"+R, "2D jet pT response for R bigger", 40, 20, 60, 50, 20, 70)
    h_pt_Bigger.GetXaxis().SetTitle("Detector pT")
    h_pt_Bigger.GetYaxis().SetTitle("Truth pT")    

    h_Dpt = ROOT.TH2D("h_Dpt_R_"+R, "Jet pT [40,60]: DpT response for R = 0."+R, 30, 0, 30, 35, 0, 35)
    h_Dpt.GetXaxis().SetTitle("Detector DpT")
    h_Dpt.GetYaxis().SetTitle("Truth DpT")    
        
    h_4D = ROOT.THnD("h_4D_R_"+R, "4D jet response for R = 0."+R, 4, array.array('i', [40, 40, 100, 100]), array.array('d', [0, 0, -20, -20]), array.array('d', [200, 200, 80, 80]))
    h_4D.GetAxis(0).SetTitle("Truth pT")
    h_4D.GetAxis(1).SetTitle("Detector pT")
    h_4D.GetAxis(2).SetTitle("Truth DpT")
    h_4D.GetAxis(3).SetTitle("Detector DpT")

    def makePartialResponse(pt_min,pt_max,weight,job):
        jets = uproot3.open("results_new/Analysis_R_0_"+R+"_PT_"+str(pt_min)+"_job_"+str(job)+".root")["jetTree"]
        #jets = uproot3.open("test.root")["jetTree"]

        cutJets = ((jets.array("sigJet_Detector_smallerPt").flatten()!=0) & (jets.array("sigJet_Detector_biggerPt").flatten()!=0) & (jets.array("sigJet_Truth_biggerPt").flatten()!=0))

        sigJet_Truth_smallerPt  = jets.array("sigJet_Truth_smallerPt").flatten()[cutJets]
        sigJet_Detector_smallerPt = jets.array("sigJet_Detector_smallerPt").flatten()[cutJets]

        sigJet_Truth_biggerPt   = jets.array("sigJet_Truth_biggerPt").flatten()[cutJets]       
        sigJet_Detector_biggerPt = jets.array("sigJet_Detector_biggerPt").flatten()[cutJets]

        dPt_Truth       = sigJet_Truth_biggerPt-sigJet_Truth_smallerPt
        dPt_Detector    = sigJet_Detector_biggerPt - sigJet_Detector_smallerPt 

        for i in range(0,len(dPt_Truth)-1):
            h_pt.Fill(sigJet_Detector_smallerPt[i],sigJet_Truth_smallerPt[i],weight)
            h_pt_Bigger.Fill(sigJet_Detector_biggerPt[i],sigJet_Truth_biggerPt[i],weight)
            h_4D.Fill(array.array('d',[sigJet_Truth_smallerPt[i],sigJet_Detector_smallerPt[i],dPt_Truth[i],dPt_Detector[i]]),weight)
            if (sigJet_Truth_smallerPt[i] > 40 and sigJet_Truth_smallerPt[i] < 60):
                h_Dpt.Fill(dPt_Detector[i],dPt_Truth[i],weight)
        
    
    for i in range(0,len(ptWeightList)-1):
        for job in range(1,301):
            makePartialResponse(ptWeightList[i][0],ptWeightList[i+1][0],ptWeightList[i][1]/ptWeightList[i][2],job)
    
    #makePartialResponse(ptWeightList[0][0],ptWeightList[1][0],ptWeightList[0][1]/ptWeightList[0][2],1)

    h_pt.SetStats(0)
    h_pt_Bigger.SetStats(0)
    h_Dpt.SetStats(0)

    # Save the histogram to a ROOT file
    f = ROOT.TFile(R+"_response.root", "RECREATE")
    h_pt.Write()
    h_pt_Bigger.Write()
    h_Dpt.Write()
    h_4D.Write()
    f.Close()
