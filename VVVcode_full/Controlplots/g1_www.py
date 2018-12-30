#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser


from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack, TGraph, TGraphErrors,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, TVectorD, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite


############################################
#              Job steering                #
############################################

parser = OptionParser()

#parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="mu")


parser.add_option('--control', action='store_true', dest='control', default=False, help='control plot')
parser.add_option('--numak8jets', action='store',type="string", dest='numak8jets', default="1", help='preselection numak8jets; "1":numak8jets==1 ;"2":numak8jets==2;"1or2":numak8jets==1or2')

parser.add_option('--stcut', action="store",type="int",dest="stcut",default=1000)
parser.add_option('--scale1', action="store",type="int",dest="scale1",default=100)
parser.add_option('--scale2', action="store",type="int",dest="scale2",default=10)

#parser.add_option('--category', action="store",type="string",dest="category",default="ALLP")


(options, args) = parser.parse_args()




###############################
## doFit Class Implemetation ##
###############################

def UnderOverFlow1D(h):
    Bins=h.GetNbinsX();
    h.SetBinContent( 1,  h.GetBinContent(1)+h.GetBinContent(0) );
    h.SetBinError(   1,  math.sqrt( h.GetBinError(1)*h.GetBinError(1) + h.GetBinError(0)*h.GetBinError(0)) );
    h.SetBinContent( Bins,  h.GetBinContent(Bins)+h.GetBinContent(Bins+1) );
    h.SetBinError(   Bins,  math.sqrt( h.GetBinError(Bins)*h.GetBinError(Bins) + h.GetBinError(Bins+1)*h.GetBinError(Bins+1)) );
    return h;

class doFit_wj_and_wlvj:
    


    def __init__(self, in_channel):

        self.setTDRStyle();
        self.channel=in_channel;

        self.file_Directory="/eos/cms/store/user/xulyu/data2016/2016/";
            


        #if options.closuretest == 0:
        #    self.file_data = (self.channel+"_PKUTree_data_xww.root");#keep blind!!!!
        #else:
        #    self.file_data = (self.channel+"_PKUTree_data_xww.root");#keep blind!!!!
        
        if options.control==1:
            self.file_data = (self.channel+"_PKUTree_2017vvv.root");
                #self.file_data = (self.channel+"_out_case2_1500.root");


        self.file_signal     = (self.channel+"_out_case1_off_3000.root");
        self.file_signal2     = (self.channel+"_out_case1_off_1500.root");
        self.file_signal3     = (self.channel+"_out_case2_off_3000.root");
        self.file_signal4     = (self.channel+"_out_case2_off_1500.root");

        self.file_WJets0_mc  = (self.channel+"_PKUTree_WJetsPt180_xww.root");
        self.file_VV_mc      = (self.channel+"_PKUTree_VV_xww.root");# WW+WZ
        self.file_TTbar_mc   = (self.channel+"_PKUTree_TTBARpowheg_xww.root");
        self.file_STop_mc    = (self.channel+"_PKUTree_SingleTop_xww.root");


        ## color palette 
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : kRed,
            'VV' : 6,
            'STop' : 7,
            'TTbar' : kGray,
            'Signal': 1,
            'Uncertainty' : kBlack,
        
        }

        ## for basic selection         
        self.vpt_cut   = 200;
        self.MET_cut = 40;
        self.lpt_cut   = 55;

        if self.channel=="el":
            self.MET_cut= 80; 
            self.lpt_cut = 45;
        self.deltaPhi_METj_cut =2.0;
        if options.numak8jets == "1":
            self.nak8jet1 =1;
            self.nak8jet2 =1;
        if options.numak8jets == "2":
            self.nak8jet1 =2;
            self.nak8jet2 =2;
        if options.numak8jets == "1or2":
            self.nak8jet1 =1;
            self.nak8jet2 =2;



    def setTDRStyle(self):
        self.tdrStyle =TStyle("tdrStyle","Style for P-TDR");
        #For the canvas:
        self.tdrStyle.SetCanvasBorderMode(0);
        self.tdrStyle.SetCanvasColor(kWhite);
        self.tdrStyle.SetCanvasDefH(600); #Height of canvas
        self.tdrStyle.SetCanvasDefW(600); #Width of canvas
        self.tdrStyle.SetCanvasDefX(0); #POsition on screen
        self.tdrStyle.SetCanvasDefY(0);
      
        #For the Pad:
        self.tdrStyle.SetPadBorderMode(0);
        self.tdrStyle.SetPadColor(kWhite);
        self.tdrStyle.SetPadGridX(False);
        self.tdrStyle.SetPadGridY(False);
        self.tdrStyle.SetGridColor(0);
        self.tdrStyle.SetGridStyle(3);
        self.tdrStyle.SetGridWidth(1);
      
        #For the frame:
        self.tdrStyle.SetFrameBorderMode(0);
        self.tdrStyle.SetFrameBorderSize(1);
        self.tdrStyle.SetFrameFillColor(0);
        self.tdrStyle.SetFrameFillStyle(0);
        self.tdrStyle.SetFrameLineColor(1);
        self.tdrStyle.SetFrameLineStyle(1);
        self.tdrStyle.SetFrameLineWidth(1);
      
        #For the histo:
        self.tdrStyle.SetHistLineColor(1);
        self.tdrStyle.SetHistLineStyle(0);
        self.tdrStyle.SetHistLineWidth(1);
        self.tdrStyle.SetEndErrorSize(2);
        self.tdrStyle.SetErrorX(0.);
        self.tdrStyle.SetMarkerStyle(20);
      
        #For the fit/function:
        self.tdrStyle.SetOptFit(1);
        self.tdrStyle.SetFitFormat("5.4g");
        self.tdrStyle.SetFuncColor(2);
        self.tdrStyle.SetFuncStyle(1);
        self.tdrStyle.SetFuncWidth(1);
      
        #For the date:
        self.tdrStyle.SetOptDate(0);
      
        #For the statistics box:
        self.tdrStyle.SetOptFile(0);
        self.tdrStyle.SetOptStat(0); #To display the mean and RMS:
        self.tdrStyle.SetStatColor(kWhite);
        self.tdrStyle.SetStatFont(42);
        self.tdrStyle.SetStatFontSize(0.025);
        self.tdrStyle.SetStatTextColor(1);
        self.tdrStyle.SetStatFormat("6.4g");
        self.tdrStyle.SetStatBorderSize(1);
        self.tdrStyle.SetStatH(0.1);
        self.tdrStyle.SetStatW(0.15);
      
        #Margins:
        self.tdrStyle.SetPadTopMargin(0.05);
        self.tdrStyle.SetPadBottomMargin(0.13);
        self.tdrStyle.SetPadLeftMargin(0.18);
        self.tdrStyle.SetPadRightMargin(0.06);
      
        #For the Global title:
        self.tdrStyle.SetOptTitle(0);
        self.tdrStyle.SetTitleFont(42);
        self.tdrStyle.SetTitleColor(1);
        self.tdrStyle.SetTitleTextColor(1);
        self.tdrStyle.SetTitleFillColor(10);
        self.tdrStyle.SetTitleFontSize(0.05);
      
        #For the axis titles:
        self.tdrStyle.SetTitleColor(1, "XYZ");
        self.tdrStyle.SetTitleFont(42, "XYZ");
        self.tdrStyle.SetTitleSize(0.03, "XYZ");
        self.tdrStyle.SetTitleXOffset(0.9);
        self.tdrStyle.SetTitleYOffset(1.5);
      
        #For the axis labels:
        self.tdrStyle.SetLabelColor(1, "XYZ");
        self.tdrStyle.SetLabelFont(42, "XYZ");
        self.tdrStyle.SetLabelOffset(0.007, "XYZ");
        self.tdrStyle.SetLabelSize(0.03, "XYZ");
      
        #For the axis:
        self.tdrStyle.SetAxisColor(1, "XYZ");
        self.tdrStyle.SetStripDecimals(kTRUE);
        self.tdrStyle.SetTickLength(0.03, "XYZ");
        self.tdrStyle.SetNdivisions(510, "XYZ");
        self.tdrStyle.SetPadTickX(1); #To get tick marks on the opposite side of the frame
        self.tdrStyle.SetPadTickY(1);
      
        #Change for log plots:
        self.tdrStyle.SetOptLogx(0);
        self.tdrStyle.SetOptLogy(0);
        self.tdrStyle.SetOptLogz(0);
      
        #Postscript options:
        self.tdrStyle.SetPaperSize(20.,20.);
        self.tdrStyle.cd();

    ######## ++++++++++++++
    def ControlPlots(self):
        
        ######definition of self.nak8jet1 and self.nak8jet2 can be seen at Line 109

### common selection (for Nj=1 and Nj=2)
###### 𝜏_41 max<0.3,for M=40-100
###### Mjmax > 40
###### M_j1 > 40
###### ST > 1000 GeV
###### 𝜏_42 max<0.55 if Mjmax > 100
###### Nb-jet = 0
###### Nak4-jet < 5

### if Nj=2:
###### 𝜏_41 min<0.3 if Mjmin > 40
###### Mjmin < 90

        selection_nj1="jetAK8puppi_sdcorr>40&&ST>%s&&Mj_max>40&&(Mj_max<100?tau41_Mj_max<0.3:1)&&nbtag==0&&num_ak4jets<5&&(((%s==2)?(num_ak8jets==%s&&jet_pt_puppi_2>0&& abs(jetAK8puppi_eta_2)<2.4&&Mj_min<90&&((Mj_max>100)?(tau42_Mj_max<0.55):1)&&((Mj_min>40)?(tau41_Mj_min<0.3):1)):(num_ak8jets==%s))||((%s==2)?(num_ak8jets==%s&&jet_pt_puppi_2>0&& abs(jetAK8puppi_eta_2)<2.4&Mj_min<90&&((Mj_max>100)?(tau42_Mj_max<0.55):1)&&((Mj_min>40)?(tau41_Mj_min<0.3):1)):(num_ak8jets==%s)))&&(((%s==2)?(m_lvj>1000):(m_jlv>1000))||((%s==2)?(m_lvj>1000):(m_jlv>1000))) &&abs(l_eta)<2.4 &&  W_pt > 200 && abs(jetAK8puppi_eta)<2.4  &&jet_pt_puppi>200&& passFilter_HBHE>0 && passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0   && l_pt> %s && MET_et>%s"%(options.stcut,self.nak8jet1,self.nak8jet1,self.nak8jet1,self.nak8jet2, self.nak8jet2, self.nak8jet2,self.nak8jet1,self.nak8jet2, self.lpt_cut, self.MET_cut)
        selection_nj1_case2a="status_1!=4&&status_2==4&&status_3==4&&jetAK8puppi_sdcorr>40&&ST>%s&&Mj_max>40&&(Mj_max<100?tau41_Mj_max<0.3:1)&&nbtag==0&&num_ak4jets<5&&(((%s==2)?(num_ak8jets==%s&&jet_pt_puppi_2>0&& abs(jetAK8puppi_eta_2)<2.4&Mj_min<90&&((Mj_max>100)?(tau42_Mj_max<0.55):1)&&((Mj_min>40)?(tau41_Mj_min<0.3):1)):(num_ak8jets==%s))||((%s==2)?(num_ak8jets==%s&&jet_pt_puppi_2>0&& abs(jetAK8puppi_eta_2)<2.4&Mj_min<90&&((Mj_max>100)?(tau42_Mj_max<0.55):1)&&((Mj_min>40)?(tau41_Mj_min<0.3):1)):(num_ak8jets==%s)))&&(((%s==2)?(m_lvj>1000):(m_jlv>1000))||((%s==2)?(m_lvj>1000):(m_jlv>1000))) &&abs(l_eta)<2.4 &&  W_pt > 200 && abs(jetAK8puppi_eta)<2.4  &&jet_pt_puppi>200&& passFilter_HBHE>0 && passFilter_GlobalHalo>0 && passFilter_HBHEIso>0 && passFilter_ECALDeadCell>0 && passFilter_GoodVtx>0 && passFilter_EEBadSc>0   && l_pt> %s && MET_et>%s"%(options.stcut,self.nak8jet1,self.nak8jet1,self.nak8jet1,self.nak8jet2, self.nak8jet2, self.nak8jet2,self.nak8jet1,self.nak8jet2, self.lpt_cut, self.MET_cut)


        self.Make_Controlplots(selection_nj1,selection_nj1_case2a,"preselection",options.numak8jets);

    ######## ++++++++++++++
    def Make_Controlplots(self,selection_nj1,selection_nj1_case2a,tag,numak8jets,TTBarControl=0):
       #       if numak8jets==2:
        self.make_controlplot(numak8jets,"tau42_Mj_max",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau42_Mj_max","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau43_Mj_max",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau43_Mj_max","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau32_Mj_max",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau32_Mj_max","Events/(0.05)", 0, TTBarControl);
        
        self.make_controlplot(numak8jets,"tau42_Mj_min",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau42_Mj_min","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau43_Mj_min",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau43_Mj_min","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau32_Mj_min",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau32_Mj_min","Events/(0.05)", 0, TTBarControl);
        

    	self.make_controlplot(numak8jets,"m_lvj",selection_nj1,selection_nj1_case2a,tag,20,0,4000,"M_{jjlv} (GeV)","Events/(200 GeV)",0, TTBarControl );
        self.make_controlplot(numak8jets,"num_ak8jets",selection_nj1,selection_nj1_case2a,tag,2,0.5, 2.5,"N_{AK8jets}","Events", 0, TTBarControl);
        self.make_controlplot(numak8jets,"m_jlv",selection_nj1,selection_nj1_case2a,tag,35,0,3500,"M_{jlv} (GeV)","Events/(100 GeV)",0, TTBarControl );
        if(numak8jets=="2" or numak8jets=="1or2"):
            self.make_controlplot(numak8jets,"MassVV[0]",selection_nj1,selection_nj1_case2a,tag,35,0,3500,"M_{W_{1}W_{2}} (GeV)","Events/(100 GeV)",0, TTBarControl );

        self.make_controlplot(numak8jets,"jet_mass_puppi_corr",selection_nj1,selection_nj1_case2a,tag,28,0,280,"M_{j1(Ak8)} (GeV)","Events/(10 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"Mj_min",selection_nj1,selection_nj1_case2a,tag,28,0,280,"Mj_min","Events/(10 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"Mj_max",selection_nj1,selection_nj1_case2a,tag,28,0,280,"Mj_max","Events/(10 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau31_Mj_max",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau31_Mj_max","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau21_Mj_max",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau21_Mj_max","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau41_Mj_max",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau41_Mj_max","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau31_Mj_min",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau31_Mj_min","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau21_Mj_min",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau21_Mj_min","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"tau41_Mj_min",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau41_Mj_min","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"PTw2divPTw1",selection_nj1,selection_nj1_case2a,tag,20,0,1,"PTw2/PTw1","Events", 0, TTBarControl);
        self.make_controlplot(numak8jets,"PTw3divPTw1",selection_nj1,selection_nj1_case2a,tag,20,0,1,"PTw3/PTw1","Events", 0, TTBarControl);
        self.make_controlplot(numak8jets,"PTw3divPTw2",selection_nj1,selection_nj1_case2a,tag,20,0,1,"PTw3/PTw2","Events", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jet_mass_puppi_corr",selection_nj1,selection_nj1_case2a,tag,28,0,280,"M_{j1(Ak8)} (GeV)","Events/(10 GeV)", 0, TTBarControl);
        #        if numak8jets==1:
        self.make_controlplot(numak8jets,"MassVV[1]",selection_nj1,selection_nj1_case2a,tag,35,0,3500,"M_{W_{1}W_{3}} (GeV)","Events/(100 GeV)",0, TTBarControl );
        self.make_controlplot(numak8jets,"MassVV[2]",selection_nj1,selection_nj1_case2a,tag,35,0,3500,"M_{W_{2}W_{3}} (GeV)","Events/(100 GeV)",0, TTBarControl );
        
        #       if numak8jets==2:
        self.make_controlplot(numak8jets,"W_pt",selection_nj1,selection_nj1_case2a,tag,30,0, 1500,"PT(l,MET) (GeV)","Events/(50 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jet_pt_puppi",selection_nj1,selection_nj1_case2a,tag,30,0, 3000,"PT_{j1(AK8)} (GeV)","Events/(100 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jet_pt_puppi_2",selection_nj1,selection_nj1_case2a,tag,20,0, 1000,"PT_{j2(AK8)} (GeV)","Events/(50 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"ST",selection_nj1,selection_nj1_case2a,tag,35,0, 3500,"S_{T} (GeV)","Events/(100 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"num_ak8jets",selection_nj1,selection_nj1_case2a,tag,2,0.5, 2.5,"N_{AK8jets}","Events", 0, TTBarControl);
        self.make_controlplot(numak8jets,"num_ak4jets",selection_nj1,selection_nj1_case2a,tag,8,0.5, 8.5,"N_{AK4jets}","Events", 0, TTBarControl);
        self.make_controlplot(numak8jets,"mtVlepnew",selection_nj1,selection_nj1_case2a,tag,25,0, 500,"MT(l,MET) (GeV)","Events/(20 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"l_pt",selection_nj1,selection_nj1_case2a,tag,24,0, 1200,"PT_{l} (GeV)","Events/(50 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"trackIso",selection_nj1,selection_nj1_case2a,tag,20,0,2.5,"trackIso(l) (GeV)","Events/0.125", 0, TTBarControl);
        self.make_controlplot(numak8jets,"muisolation",selection_nj1,selection_nj1_case2a,tag,20,0,0.1,"muisolation","Events/0.05", 0, TTBarControl);
        self.make_controlplot(numak8jets,"l_eta",selection_nj1,selection_nj1_case2a,tag,20,0,2.5,"#eta_{l}","Events/(0.125)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"deltaRlepjet",selection_nj1,selection_nj1_case2a,tag,20,0,5,"#DeltaR(lep,jet1)","Events/(0.25)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"deltaRlepjet_2",selection_nj1,selection_nj1_case2a,tag,20,0,5,"#DeltaR(lep,jet2)","Events/(0.25)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"delPhijetmet",selection_nj1,selection_nj1_case2a,tag,20,0,3.14159,"#Delta#phi(jet1,MET)","Events/(0.262)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"delPhijetmet_2",selection_nj1,selection_nj1_case2a,tag,12,0,3.14159,"#Delta#phi(jet2,MET)","Events/(0.262)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"delPhijetlep",selection_nj1,selection_nj1_case2a,tag,12,0,3.14159,"#Delta#phi(jet1,lep)","Events/(0.262)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"delPhijetlep_2",selection_nj1,selection_nj1_case2a,tag,12,0,3.14159,"#Delta#phi(jet2,lep)","Events/(0.262)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"deltaRjet1jet2",selection_nj1,selection_nj1_case2a,tag,20,0,6,"#DeltaR(j1,j2)","Events/(0.3)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"MET_et",selection_nj1,selection_nj1_case2a,tag,24,0,1200,"MET (GeV)","Events/(50 GeV)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"nPV",selection_nj1,selection_nj1_case2a,tag,20,0,40,"nPV","Events/(2)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"nbtag",selection_nj1,selection_nj1_case2a,tag,5,-0.5,4.5,"number of b-jets","Events", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jetAK8puppi_sdcorr_2",selection_nj1,selection_nj1_case2a,tag,22,0,220,"M_{j2(AK8)} (GeV)","Events/(10 GeV)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jet_tau2tau1_puppi",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{21_j1}","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jet_tau2tau1_puppi_2",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{21_j2}","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jet_tau4tau2_puppi",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{42_j1}","Events/(0.05)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"jet_tau4tau2_puppi_2",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{42_j2}","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jetAK8puppi_tau31",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{31_j1}","Events/(0.05)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"jetAK8puppi_tau32",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{32_j1}","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jetAK8puppi_tau41",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{41_j1}","Events/(0.05)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"jetAK8puppi_tau43",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{43_j1}","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jetAK8puppi_tau31_2",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{31_j2}","Events/(0.05)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"jetAK8puppi_tau32_2",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{32_j2}","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"jetAK8puppi_tau41_2",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{41_j2}","Events/(0.05)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"jetAK8puppi_tau43_2",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#tau_{43_j2}","Events/(0.05)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"deltaetajet1jet2",selection_nj1,selection_nj1_case2a,tag,20,0,4,"#Delta#eta_{j1,j2}","Events(0.2)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"deltaetajet1lep",selection_nj1,selection_nj1_case2a,tag,20,0,4,"#Delta#eta_{j1,lep}","Events(0.2)", 0, TTBarControl);
        #self.make_controlplot(numak8jets,"deltaetajet2lep",selection_nj1,selection_nj1_case2a,tag,20,0,4,"#Delta#eta_{j2,lep}","Events(0.2)", 0, TTBarControl);
        self.make_controlplot(numak8jets,"delPhilepmet",selection_nj1,selection_nj1_case2a,tag,12,0,3.14159,"#Delta#phi(lep,met)","Events/(0.262)", 0, TTBarControl);
        self.make_controlplot_productof3tau(numak8jets,"jet_tau2tau1_puppi","jetAK8puppi_tau31","jetAK8puppi_tau41",selection_nj1,selection_nj1_case2a,tag,20,0,1,"#(tau_{21_j1}*#tau_{31_j1}*#tau_{41_j1})^{1/3}","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot_averageof3tau(numak8jets,"jet_tau2tau1_puppi","jetAK8puppi_tau31","jetAK8puppi_tau41",selection_nj1,selection_nj1_case2a,tag,20,0,1,"(#tau_{21_j1}+#tau_{31_j1}+#tau_{41_j1})/3","Events/(0.05)", 0, TTBarControl);
        self.make_controlplot_productof3tau(numak8jets,"W_pt","jet_pt_puppi","jet_pt_puppi_2",selection_nj1,selection_nj1_case2a,tag,30,0,1500,"(PT_{W_{l}}*PT_{j1}*PT_{j2})^{1/3} [GeV]","Events/(50 GeV)", 0, TTBarControl);
        self.make_controlplot_productof2Wpt(numak8jets,"W_pt","jet_pt_puppi",selection_nj1,selection_nj1_case2a,tag,30,0,1800,"(PT_{W_{l}}*PT_{j1})^{1/2} [GeV]","Events/(60 GeV)", 0, TTBarControl);
        self.make_controlplot_averageof3tau(numak8jets,"MassVV[0]","MassVV[1]","MassVV[2]",selection_nj1,selection_nj1_case2a,tag,35,0,3500,"(M_{W_{1}W_{2}}+M_{W_{1}W_{3}}+M_{W_{2}W_{3}})/3","Events/(100 GeV)", 0, TTBarControl);
        self.make_controlplot_averageof2MassVV(numak8jets,"MassVV[1]","MassVV[2]",selection_nj1,selection_nj1_case2a,tag,35,0,3500,"(M_{W_{1}W_{3}}+M_{W_{2}W_{3}})/2","Events/(100 GeV)", 0, TTBarControl);


    ######## ++++++++++++++
    def make_controlplot(self,numak8jets,variable,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="", logy=0 ,TTBarControl=0,):
        tmp_lumi=self.GetLumi()
        tmp_signal_scale="%s"%(options.scale1)
        tmp_signal_scale1="%s"%(options.scale2)
        weight_mc_forSignal="weight*%s*%s"%(tmp_lumi, tmp_signal_scale);
        weight_mc_forSignal1="weight*%s*%s"%(tmp_lumi, tmp_signal_scale1);


        weight_mc_forV="weight*%s"%(tmp_lumi);
        weight_mc_forT="weight*%s"%(tmp_lumi);
        ##weight_mc_forV="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forV.getVal());#little error rrv_wtagger_eff_reweight_forV
        ##weight_mc_forT="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forT.getVal());#little error rrv_wtagger_eff_reweight_forT
        weight_mc_forG="weight*%s"%(tmp_lumi); #General

        tmp_WJets_scale=0.8
        tmp_TTBar_scale=1.0

        ##if TTBarControl==0: 
        ##    if self.channel=="mu": tmp_WJets_scale=1.2.18
        ##    if self.channel=="el": tmp_WJets_scale=1.2.01 
        ##else: 
        ##    if self.channel=="mu": tmp_TTBar_scale=0.85
        ##    if self.channel=="el": tmp_TTBar_scale=0.70

        if self.channel=="mu": tmp_WJets_scale=1.2
        if self.channel=="el": tmp_WJets_scale=0.97
        if self.channel=="mu": tmp_TTBar_scale=1
        if self.channel=="el": tmp_TTBar_scale=0.81

        weight_mc_forWJets="weight*%s*%s"%(tmp_lumi, tmp_WJets_scale); #General
        weight_mc_forTTBar="weight*%s*%s"%(tmp_lumi, tmp_TTBar_scale); #General

        weightcut_mc_forSignal="(%s)*(%s)"%(weight_mc_forSignal,cut);
        weightcut_mc_forSignala="(%s)*(%s)"%(weight_mc_forSignal,cut1);
        weightcut_mc_forSignal1="(%s)*(%s)"%(weight_mc_forSignal1,cut);
        weightcut_mc_forSignala1="(%s)*(%s)"%(weight_mc_forSignal1,cut1);

        weightcut_mc_forV="(%s)*(%s)"%(weight_mc_forV,cut);
        weightcut_mc_forT="(%s)*(%s)"%(weight_mc_forT,cut);
        weightcut_mc_forG="(%s)*(%s)"%(weight_mc_forG,cut);
        weightcut_mc_forWJets="(%s)*(%s)"%(weight_mc_forWJets,cut);
        weightcut_mc_forTTBar="(%s)*(%s)"%(weight_mc_forTTBar,cut);
        weightcut_data="%s"%(cut);
        print "weightcut_mc_forV="+weightcut_mc_forV;
        print "weightcut_mc_forT="+weightcut_mc_forT;
        print "weightcut_mc_forG="+weightcut_mc_forG;
        print "weightcut_mc_forWJets="+weightcut_mc_forWJets;
        print "weightcut_mc_forTTBar="+weightcut_mc_forTTBar;
        hist_data =TH1D("hist_data","hist_data"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal =TH1D("hist_Signal","hist_Signal"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal.Sumw2();
        hist_Signal2 =TH1D("hist_Signal2","hist_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal2.Sumw2();
        hist_Signal3 =TH1D("hist_Signal3","hist_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3.Sumw2();
        hist_Signal4 =TH1D("hist_Signal4","hist_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4.Sumw2();
        hist_Signal3a =TH1D("hist_Signal3a","hist_Signal3a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3a.Sumw2();
        hist_Signal4a =TH1D("hist_Signal4a","hist_Signal4a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4a.Sumw2();
        hist_WJets=TH1D("hist_WJets","hist_WJets"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_WJets.Sumw2();
        hist_TTbar=TH1D("hist_TTbar","hist_TTbar"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TTbar.Sumw2();
        hist_STop =TH1D("hist_STop","hist_STop"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_STop.Sumw2();
        hist_VV   =TH1D("hist_VV","hist_VV"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_VV.Sumw2();
        hist_TotalMC =TH1D("hist_TotalMC","hist_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TotalMC.Sumw2();


        hstack_TotalMC = THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle))
        if TTBarControl==0:
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);
            
            hstack_TotalMC.Add(hist_WJets); 
        else:
            hstack_TotalMC.Add(hist_WJets); 
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);


        hist_data.SetLineColor(self.color_palet["data"]); hist_data.SetFillColor(self.color_palet["data"]);
        hist_Signal.SetLineColor(self.color_palet["Signal"]); hist_Signal.SetFillColor(self.color_palet["Signal"]); hist_Signal.SetFillStyle(0);hist_Signal.SetLineWidth(2);
        hist_Signal2.SetLineColor(kBlue); hist_Signal2.SetFillColor(self.color_palet["Signal"]); hist_Signal2.SetFillStyle(0);hist_Signal2.SetLineWidth(2);
        hist_Signal3.SetLineColor(kGreen); hist_Signal3.SetFillColor(self.color_palet["Signal"]); hist_Signal3.SetFillStyle(0);hist_Signal3.SetLineWidth(2);
        hist_Signal4.SetLineColor(kOrange); hist_Signal4.SetFillColor(self.color_palet["Signal"]); hist_Signal4.SetFillStyle(0);hist_Signal4.SetLineWidth(2);
        hist_Signal3a.SetLineColor(kGreen); hist_Signal3a.SetFillColor(self.color_palet["Signal"]); hist_Signal3a.SetFillStyle(0);hist_Signal3a.SetLineWidth(2);hist_Signal3a.SetLineStyle(2);
        hist_Signal4a.SetLineColor(kOrange); hist_Signal4a.SetFillColor(self.color_palet["Signal"]); hist_Signal4a.SetFillStyle(0);hist_Signal4a.SetLineWidth(2);hist_Signal4a.SetLineStyle(2);


        hist_WJets.SetLineColor(kBlack); hist_WJets.SetFillColor(self.color_palet["WJets"]);
        hist_TTbar.SetLineColor(kBlack); hist_TTbar.SetFillColor(self.color_palet["TTbar"]);
        hist_STop.SetLineColor(kBlack); hist_STop.SetFillColor(self.color_palet["STop"]);
        hist_VV.SetLineColor(kBlack); hist_VV.SetFillColor(self.color_palet["VV"]);



        tree_data   =TChain("PKUTree");  tree_data.Add(self.file_Directory+self.file_data);
        tree_Signal =TChain("PKUTree");  tree_Signal.Add(self.file_Directory+self.file_signal);
        tree_Signal2 =TChain("PKUTree");  tree_Signal2.Add(self.file_Directory+self.file_signal2);
        tree_Signal3 =TChain("PKUTree");  tree_Signal3.Add(self.file_Directory+self.file_signal3);
        tree_Signal4 =TChain("PKUTree");  tree_Signal4.Add(self.file_Directory+self.file_signal4);

        tree_WJets  =TChain("PKUTree");tree_WJets.Add(self.file_Directory+self.file_WJets0_mc);
        tree_TTbar  =TChain("PKUTree");tree_TTbar.Add(self.file_Directory+self.file_TTbar_mc);
        tree_STop   =TChain("PKUTree");  tree_STop.Add(self.file_Directory+self.file_STop_mc);
        tree_VV     =TChain("PKUTree");      tree_VV.Add(self.file_Directory+self.file_VV_mc);

        
        tree_data.Draw("abs(%s) >> hist_data"%(variable), weightcut_data);
        tree_Signal.Draw("abs(%s) >> hist_Signal"%(variable), weightcut_mc_forSignal);
        tree_Signal2.Draw("abs(%s) >> hist_Signal2"%(variable), weightcut_mc_forSignal1);
        tree_Signal3.Draw("abs(%s) >> hist_Signal3"%(variable), weightcut_mc_forSignal);
        tree_Signal4.Draw("abs(%s) >> hist_Signal4"%(variable), weightcut_mc_forSignal1);
        tree_Signal3.Draw("abs(%s) >> hist_Signal3a"%(variable), weightcut_mc_forSignala);
        tree_Signal4.Draw("abs(%s) >> hist_Signal4a"%(variable), weightcut_mc_forSignala1);

        #tree_WJets.Draw("%s >> hist_WJets"%(variable), weightcut_mc_forG);
        tree_WJets.Draw("abs(%s) >> hist_WJets"%(variable), weightcut_mc_forWJets);
        #tree_TTbar.Draw("%s >> hist_TTbar"%(variable), weightcut_mc_forT);
        tree_TTbar.Draw("abs(%s) >> hist_TTbar"%(variable), weightcut_mc_forTTBar);
        tree_STop.Draw("abs(%s) >> hist_STop"%(variable), weightcut_mc_forT);
        tree_VV.Draw("abs(%s) >> hist_VV"%(variable), weightcut_mc_forV);

 
       
        hist_WJets=UnderOverFlow1D(hist_WJets);
        hist_TTbar=UnderOverFlow1D(hist_TTbar);
        hist_STop=UnderOverFlow1D(hist_STop);
        hist_VV=UnderOverFlow1D(hist_VV);
        hist_Signal=UnderOverFlow1D(hist_Signal);
        hist_Signal2=UnderOverFlow1D(hist_Signal2);
        hist_Signal3=UnderOverFlow1D(hist_Signal3);
        hist_Signal4=UnderOverFlow1D(hist_Signal4);
        hist_Signal3a=UnderOverFlow1D(hist_Signal3a);
        hist_Signal4a=UnderOverFlow1D(hist_Signal4a);

        hist_TotalMC.Add(hist_WJets); hist_TotalMC.Add(hist_TTbar); hist_TotalMC.Add(hist_STop); hist_TotalMC.Add(hist_VV);


        canvas_controlplot = TCanvas("canvas_controlplot"+variable,"canvas_controlplot"+variable, 500,500);

        fPads1 = TPad("pad1", "Run2", 0.0, 0.0, 1.00, 1.00);
        #fPads2 = TPad("pad2", "", 0.00, 0.00, 1.00, 0.28);
        #canvas_controlplot.cd();
        fPads1.SetBottomMargin(0.07);
        fPads1.SetLeftMargin(0.10);
        fPads1.SetRightMargin(0.03);
##fPads2.SetTopMargin(0);
##      fPads2.SetBottomMargin(0.25);

        fPads1.Draw();
##       fPads2.Draw();
        fPads1.cd();

            
        
        hist_data.SetBinErrorOption(TH1D.kPoisson);
        #hist_data.Draw("e");
        #hstack_TotalMC.GetXaxis().SetTitle("%s"%(xtitle));
        #hist_data.Draw("same e");
        hist_Signal.Sumw2();
        hist_Signal2.Sumw2();
        hist_Signal3.Sumw2();
        hist_Signal4.Sumw2();
        hist_Signal3a.Sumw2();
        hist_Signal4a.Sumw2();

        hist_WJets.Sumw2();
        


        #hist_Signal.Scale(200);
        #hist_Signal2.Scale(200);
        #hist_Signal3a.Scale(200);
        #hist_Signal4a.Scale(200);
        #hist_Signal3.Scale(200);
        #hist_Signal4.Scale(200);


            #hist_Signal3a.Scale(hist_Signal3.GetMaximum()/hist_Signal3a.GetMaximum());
            #hist_Signal4a.Scale(hist_Signal4.GetMaximum()/hist_Signal4a.GetMaximum());

        histsigmax=TMath.Max(hist_Signal.GetMaximum(),hist_Signal2.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal3.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal4.GetMaximum());

        hist_Signal.GetYaxis().SetRangeUser(1,TMath.Max(histsigmax,hist_TotalMC.GetMaximum())*1.2);

        hist_Signal.Draw("HIST");
        hstack_TotalMC.Draw("same HIST");
        hist_Signal.Draw("same HIST");

        hist_Signal2.Draw("same HIST");
        hist_Signal3.Draw("same HIST");
        hist_Signal4.Draw("same HIST");
        hist_Signal3a.Draw("same HIST");
        hist_Signal4a.Draw("same HIST");



        hist_data.GetXaxis().SetTitleOffset(1.2);
        hist_data.GetYaxis().SetTitleOffset(1.3);
        hist_data.GetYaxis().SetTitleSize(0.07);
        hist_data.GetXaxis().SetTitleSize(0.08);
        hist_data.GetXaxis().SetLabelSize(0.06);
        hist_data.GetYaxis().SetLabelSize(0.06);


        banner = TLatex(0.95, 0.96, "35.9 fb^{-1} (13 TeV)");
        banner.SetNDC(); banner.SetTextSize(0.038); banner.SetTextFont(42); banner.SetTextAlign(31); banner.SetLineWidth(2); banner.Draw();
        CMStext = TLatex(0.10,0.96,"CMS");
        CMStext.SetNDC(); CMStext.SetTextSize(0.041); CMStext.SetTextFont(61); CMStext.SetTextAlign(11); CMStext.SetLineWidth(2); CMStext.Draw();
#    if self.channel=="el":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow e#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
#    elif self.channel=="mu":
        Extratext = TLatex(0.20, 0.96, "Preliminary");        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if(numak8jets=="1"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="1or2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1or2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
#    elif self.channel=="em":

#    elif self.channel=="em":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow l#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();

        theLeg = TLegend(0.51, 0.75, 0.76, 0.91, "", "NDC");
        theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);
        theLeg.SetTextSize(.025);
        # theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.025);
        theLeg.SetTextFont(42);

#theLeg.AddEntry(hist_data, "Data W#rightarrow #mu#nu","ep");
        theLeg.AddEntry(hist_WJets, "W+jets","F");
        theLeg.AddEntry(hist_VV, "WW/WZ","F");
        theLeg.AddEntry(hist_TTbar, "t#bar{t}","F");
        theLeg.AddEntry(hist_STop, "Single Top","F");
        #theLeg.AddEntry(gr_MCStat, "Sys.","F");
        theLeg.AddEntry(hist_Signal, "Signal(3000,1500) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal2, "Signal(1500,750) (#times %s)"%(tmp_signal_scale1),"L");
        theLeg.AddEntry(hist_Signal3, "Signal(3000,180) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal3a, "(3000,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal4, "Signal(1500,180) (#times %s)"%(tmp_signal_scale1),"L");

        theLeg.AddEntry(hist_Signal4a, "(1500,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale1),"L");


        theLeg.SetY1NDC(0.9 - 0.05*6 - 0.005);
        theLeg.SetY1(theLeg.GetY1NDC());
        theLeg.Draw();
        



        lineAtZero = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.0, hist_TotalMC.GetXaxis().GetXmax(), 1.0);
        lineAtZero.SetLineColor(2);
#lineAtZero.Draw();
        lineAtPlusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.5, hist_TotalMC.GetXaxis().GetXmax(), 1.5);
        lineAtPlusTwo.SetLineColor(2);
        lineAtPlusTwo.SetLineStyle(2);
#       lineAtPlusTwo.Draw();
        lineAtMinusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 0.5, hist_TotalMC.GetXaxis().GetXmax(), 0.5);
        lineAtMinusTwo.SetLineColor(2);
        lineAtMinusTwo.SetLineStyle(2);
#      lineAtMinusTwo.Draw();


        Directory=TString("ST_more_than_%sGeV/plots_nj%s/"%(options.stcut,numak8jets));


        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()): os.system("mkdir -p  "+Directory.Data());

        #only draw png
        rlt_file=TString(Directory.Data()+"controlplot_"+variable+"_"+tag+".pdf");
        rlt_file1=TString(Directory.Data()+"controlplot_"+variable+"_"+tag+".png");
        #fPads1.SetLogy();
        fPads1.Update();
        canvas_controlplot.SaveAs(rlt_file.Data());
        canvas_controlplot.SaveAs(rlt_file1.Data());

        if logy:
            #canvas_controlplot.SetLogy() ;
            #fPads1.SetLogy();
            fPads1.Update();
            #fPads2.Update();
            canvas_controlplot.Update();
            rlt_file.ReplaceAll(".pdf","_log.pdf"); 
            canvas_controlplot.SaveAs(rlt_file.Data());




    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
                          
    def make_controlplot_averageof3tau(self,numak8jets,variable,variable1,variable2,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="", logy=0 ,TTBarControl=0,):
        tmp_lumi=self.GetLumi()
        tmp_signal_scale="%s"%(options.scale1)
        tmp_signal_scale1="%s"%(options.scale2)
        
        weight_mc_forSignal="weight*%s*%s"%(tmp_lumi, tmp_signal_scale);
        weight_mc_forSignal1="weight*%s*%s"%(tmp_lumi, tmp_signal_scale1);


        weight_mc_forV="weight*%s"%(tmp_lumi);
        weight_mc_forT="weight*%s"%(tmp_lumi);
        ##weight_mc_forV="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forV.getVal());#little error rrv_wtagger_eff_reweight_forV
        ##weight_mc_forT="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forT.getVal());#little error rrv_wtagger_eff_reweight_forT
        weight_mc_forG="weight*%s"%(tmp_lumi); #General

        tmp_WJets_scale=0.8
        tmp_TTBar_scale=1.0

        ##if TTBarControl==0: 
        ##    if self.channel=="mu": tmp_WJets_scale=1.2.18
        ##    if self.channel=="el": tmp_WJets_scale=1.2.01 
        ##else: 
        ##    if self.channel=="mu": tmp_TTBar_scale=0.85
        ##    if self.channel=="el": tmp_TTBar_scale=0.70

        if self.channel=="mu": tmp_WJets_scale=1.2
        if self.channel=="el": tmp_WJets_scale=0.97
        if self.channel=="mu": tmp_TTBar_scale=1
        if self.channel=="el": tmp_TTBar_scale=0.81

        weight_mc_forWJets="weight*%s*%s"%(tmp_lumi, tmp_WJets_scale); #General
        weight_mc_forTTBar="weight*%s*%s"%(tmp_lumi, tmp_TTBar_scale); #General

        weightcut_mc_forSignal="(%s)*(%s)"%(weight_mc_forSignal,cut);
        weightcut_mc_forSignala="(%s)*(%s)"%(weight_mc_forSignal,cut1);
        weightcut_mc_forSignal1="(%s)*(%s)"%(weight_mc_forSignal1,cut);
        weightcut_mc_forSignala1="(%s)*(%s)"%(weight_mc_forSignal1,cut1);

        weightcut_mc_forV="(%s)*(%s)"%(weight_mc_forV,cut);
        weightcut_mc_forT="(%s)*(%s)"%(weight_mc_forT,cut);
        weightcut_mc_forG="(%s)*(%s)"%(weight_mc_forG,cut);
        weightcut_mc_forWJets="(%s)*(%s)"%(weight_mc_forWJets,cut);
        weightcut_mc_forTTBar="(%s)*(%s)"%(weight_mc_forTTBar,cut);
        weightcut_data="%s"%(cut);
        print "weightcut_mc_forV="+weightcut_mc_forV;
        print "weightcut_mc_forT="+weightcut_mc_forT;
        print "weightcut_mc_forG="+weightcut_mc_forG;
        print "weightcut_mc_forWJets="+weightcut_mc_forWJets;
        print "weightcut_mc_forTTBar="+weightcut_mc_forTTBar;
        hist_data =TH1D("hist_data","hist_data"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal =TH1D("hist_Signal","hist_Signal"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal.Sumw2();
        hist_Signal2 =TH1D("hist_Signal2","hist_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal2.Sumw2();
        hist_Signal3 =TH1D("hist_Signal3","hist_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3.Sumw2();
        hist_Signal4 =TH1D("hist_Signal4","hist_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4.Sumw2();
        hist_Signal3a =TH1D("hist_Signal3a","hist_Signal3a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3a.Sumw2();
        hist_Signal4a =TH1D("hist_Signal4a","hist_Signal4a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4a.Sumw2();
        hist_WJets=TH1D("hist_WJets","hist_WJets"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_WJets.Sumw2();
        hist_TTbar=TH1D("hist_TTbar","hist_TTbar"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TTbar.Sumw2();
        hist_STop =TH1D("hist_STop","hist_STop"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_STop.Sumw2();
        hist_VV   =TH1D("hist_VV","hist_VV"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_VV.Sumw2();
        hist_TotalMC =TH1D("hist_TotalMC","hist_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TotalMC.Sumw2();


        hstack_TotalMC = THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle))
        if TTBarControl==0:
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);
            
            hstack_TotalMC.Add(hist_WJets); 
        else:
            hstack_TotalMC.Add(hist_WJets); 
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);


        hist_data.SetLineColor(self.color_palet["data"]); hist_data.SetFillColor(self.color_palet["data"]);
        hist_Signal.SetLineColor(self.color_palet["Signal"]); hist_Signal.SetFillColor(self.color_palet["Signal"]); hist_Signal.SetFillStyle(0);hist_Signal.SetLineWidth(2);
        hist_Signal2.SetLineColor(kBlue); hist_Signal2.SetFillColor(self.color_palet["Signal"]); hist_Signal2.SetFillStyle(0);hist_Signal2.SetLineWidth(2);
        hist_Signal3.SetLineColor(kGreen); hist_Signal3.SetFillColor(self.color_palet["Signal"]); hist_Signal3.SetFillStyle(0);hist_Signal3.SetLineWidth(2);
        hist_Signal4.SetLineColor(kOrange); hist_Signal4.SetFillColor(self.color_palet["Signal"]); hist_Signal4.SetFillStyle(0);hist_Signal4.SetLineWidth(2);
        hist_Signal3a.SetLineColor(kGreen); hist_Signal3a.SetFillColor(self.color_palet["Signal"]); hist_Signal3a.SetFillStyle(0);hist_Signal3a.SetLineWidth(2);hist_Signal3a.SetLineStyle(2);
        hist_Signal4a.SetLineColor(kOrange); hist_Signal4a.SetFillColor(self.color_palet["Signal"]); hist_Signal4a.SetFillStyle(0);hist_Signal4a.SetLineWidth(2);hist_Signal4a.SetLineStyle(2);


        hist_WJets.SetLineColor(kBlack); hist_WJets.SetFillColor(self.color_palet["WJets"]);
        hist_TTbar.SetLineColor(kBlack); hist_TTbar.SetFillColor(self.color_palet["TTbar"]);
        hist_STop.SetLineColor(kBlack); hist_STop.SetFillColor(self.color_palet["STop"]);
        hist_VV.SetLineColor(kBlack); hist_VV.SetFillColor(self.color_palet["VV"]);



        tree_data   =TChain("PKUTree");  tree_data.Add(self.file_Directory+self.file_data);
        tree_Signal =TChain("PKUTree");  tree_Signal.Add(self.file_Directory+self.file_signal);
        tree_Signal2 =TChain("PKUTree");  tree_Signal2.Add(self.file_Directory+self.file_signal2);
        tree_Signal3 =TChain("PKUTree");  tree_Signal3.Add(self.file_Directory+self.file_signal3);
        tree_Signal4 =TChain("PKUTree");  tree_Signal4.Add(self.file_Directory+self.file_signal4);

        tree_WJets  =TChain("PKUTree");tree_WJets.Add(self.file_Directory+self.file_WJets0_mc);
        tree_TTbar  =TChain("PKUTree");tree_TTbar.Add(self.file_Directory+self.file_TTbar_mc);
        tree_STop   =TChain("PKUTree");  tree_STop.Add(self.file_Directory+self.file_STop_mc);
        tree_VV     =TChain("PKUTree");      tree_VV.Add(self.file_Directory+self.file_VV_mc);

        
        tree_data.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_data"%(variable,variable1,variable2), weightcut_data);
        tree_Signal.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_Signal"%(variable,variable1,variable2), weightcut_mc_forSignal);
        tree_Signal2.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_Signal2"%(variable,variable1,variable2), weightcut_mc_forSignal1);
        tree_Signal3.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_Signal3"%(variable,variable1,variable2), weightcut_mc_forSignal);
        tree_Signal4.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_Signal4"%(variable,variable1,variable2), weightcut_mc_forSignal1);
        tree_Signal3.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_Signal3a"%(variable,variable1,variable2), weightcut_mc_forSignala);
        tree_Signal4.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_Signal4a"%(variable,variable1,variable2), weightcut_mc_forSignala1);
        
        #tree_WJets.Draw("%s >> hist_WJets"%(variable), weightcut_mc_forG);
        tree_WJets.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_WJets"%(variable,variable1,variable2), weightcut_mc_forWJets);
        #tree_TTbar.Draw("%s >> hist_TTbar"%(variable), weightcut_mc_forT);
        tree_TTbar.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_TTbar"%(variable,variable1,variable2), weightcut_mc_forTTBar);
        tree_STop.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_STop"%(variable,variable1,variable2), weightcut_mc_forT);
        tree_VV.Draw("(abs(%s)+abs(%s)+abs(%s))/3 >> hist_VV"%(variable,variable1,variable2), weightcut_mc_forV);

 
       
        hist_WJets=UnderOverFlow1D(hist_WJets);
        hist_TTbar=UnderOverFlow1D(hist_TTbar);
        hist_STop=UnderOverFlow1D(hist_STop);
        hist_VV=UnderOverFlow1D(hist_VV);
        hist_Signal=UnderOverFlow1D(hist_Signal);
        hist_Signal2=UnderOverFlow1D(hist_Signal2);
        hist_Signal3=UnderOverFlow1D(hist_Signal3);
        hist_Signal4=UnderOverFlow1D(hist_Signal4);
        hist_Signal3a=UnderOverFlow1D(hist_Signal3a);
        hist_Signal4a=UnderOverFlow1D(hist_Signal4a);

        hist_TotalMC.Add(hist_WJets); hist_TotalMC.Add(hist_TTbar); hist_TotalMC.Add(hist_STop); hist_TotalMC.Add(hist_VV);


        canvas_controlplot = TCanvas("canvas_controlplot"+variable,"canvas_controlplot"+variable, 500,500);

        fPads1 = TPad("pad1", "Run2", 0.0, 0.0, 1.00, 1.00);
        #fPads2 = TPad("pad2", "", 0.00, 0.00, 1.00, 0.28);
        #canvas_controlplot.cd();
        fPads1.SetBottomMargin(0.07);
        fPads1.SetLeftMargin(0.10);
        fPads1.SetRightMargin(0.03);
##fPads2.SetTopMargin(0);
##      fPads2.SetBottomMargin(0.25);

        fPads1.Draw();
##       fPads2.Draw();
        fPads1.cd();

            
        
        hist_data.SetBinErrorOption(TH1D.kPoisson);
        #hist_data.Draw("e");
        #hstack_TotalMC.GetXaxis().SetTitle("%s"%(xtitle));
        #hist_data.Draw("same e");
        hist_Signal.Sumw2();
        hist_Signal2.Sumw2();
        hist_Signal3.Sumw2();
        hist_Signal4.Sumw2();
        hist_Signal3a.Sumw2();
        hist_Signal4a.Sumw2();

        hist_WJets.Sumw2();
        


        #hist_Signal.Scale(200);
        #hist_Signal2.Scale(200);
        #hist_Signal3a.Scale(200);
        #hist_Signal4a.Scale(200);
        #hist_Signal3.Scale(200);
        #hist_Signal4.Scale(200);


            #hist_Signal3a.Scale(hist_Signal3.GetMaximum()/hist_Signal3a.GetMaximum());
            #hist_Signal4a.Scale(hist_Signal4.GetMaximum()/hist_Signal4a.GetMaximum());

        histsigmax=TMath.Max(hist_Signal.GetMaximum(),hist_Signal2.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal3.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal4.GetMaximum());

        hist_Signal.GetYaxis().SetRangeUser(1,TMath.Max(histsigmax,hist_TotalMC.GetMaximum())*1.2);

        hist_Signal.Draw("HIST");
        hstack_TotalMC.Draw("same HIST");
        hist_Signal.Draw("same HIST");

        hist_Signal2.Draw("same HIST");
        hist_Signal3.Draw("same HIST");
        hist_Signal4.Draw("same HIST");
        hist_Signal3a.Draw("same HIST");
        hist_Signal4a.Draw("same HIST");



        hist_data.GetXaxis().SetTitleOffset(1.2);
        hist_data.GetYaxis().SetTitleOffset(1.3);
        hist_data.GetYaxis().SetTitleSize(0.07);
        hist_data.GetXaxis().SetTitleSize(0.08);
        hist_data.GetXaxis().SetLabelSize(0.06);
        hist_data.GetYaxis().SetLabelSize(0.06);


        banner = TLatex(0.95, 0.96, "35.9 fb^{-1} (13 TeV)");
        banner.SetNDC(); banner.SetTextSize(0.038); banner.SetTextFont(42); banner.SetTextAlign(31); banner.SetLineWidth(2); banner.Draw();
        CMStext = TLatex(0.10,0.96,"CMS");
        CMStext.SetNDC(); CMStext.SetTextSize(0.041); CMStext.SetTextFont(61); CMStext.SetTextAlign(11); CMStext.SetLineWidth(2); CMStext.Draw();
#    if self.channel=="el":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow e#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
#    elif self.channel=="mu":
        Extratext = TLatex(0.20, 0.96, "Preliminary");        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if(numak8jets=="1"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="1or2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1or2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
#    elif self.channel=="em":

#    elif self.channel=="em":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow l#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();

        theLeg = TLegend(0.51, 0.75, 0.76, 0.91, "", "NDC");
        theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);
        theLeg.SetTextSize(.025);
        # theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.025);
        theLeg.SetTextFont(42);

#theLeg.AddEntry(hist_data, "Data W#rightarrow #mu#nu","ep");
        theLeg.AddEntry(hist_WJets, "W+jets","F");
        theLeg.AddEntry(hist_VV, "WW/WZ","F");
        theLeg.AddEntry(hist_TTbar, "t#bar{t}","F");
        theLeg.AddEntry(hist_STop, "Single Top","F");
        #theLeg.AddEntry(gr_MCStat, "Sys.","F");
        theLeg.AddEntry(hist_Signal, "Signal(3000,1500) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal2, "Signal(1500,750) (#times %s)"%(tmp_signal_scale1),"L");
        theLeg.AddEntry(hist_Signal3, "Signal(3000,180) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal3a, "(3000,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal4, "Signal(1500,180) (#times %s)"%(tmp_signal_scale1),"L");

        theLeg.AddEntry(hist_Signal4a, "(1500,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale1),"L");


        theLeg.SetY1NDC(0.9 - 0.05*6 - 0.005);
        theLeg.SetY1(theLeg.GetY1NDC());
        theLeg.Draw();
        



        lineAtZero = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.0, hist_TotalMC.GetXaxis().GetXmax(), 1.0);
        lineAtZero.SetLineColor(2);
#lineAtZero.Draw();
        lineAtPlusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.5, hist_TotalMC.GetXaxis().GetXmax(), 1.5);
        lineAtPlusTwo.SetLineColor(2);
        lineAtPlusTwo.SetLineStyle(2);
#       lineAtPlusTwo.Draw();
        lineAtMinusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 0.5, hist_TotalMC.GetXaxis().GetXmax(), 0.5);
        lineAtMinusTwo.SetLineColor(2);
        lineAtMinusTwo.SetLineStyle(2);
#      lineAtMinusTwo.Draw();


        Directory=TString("ST_more_than_%sGeV/plots_nj%s/"%(options.stcut,numak8jets));


        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()): os.system("mkdir -p  "+Directory.Data());
        #fPads1.SetLogy();
        fPads1.Update();
        #only draw png
        if variable=="jet_tau2tau1_puppi":
            rlt_file=TString(Directory.Data()+"controlplot_"+"averageof3tau"+"_"+tag+".pdf");
            rlt_file1=TString(Directory.Data()+"controlplot_"+"averageof3tau"+"_"+tag+".png");
        if variable=="MassVV[0]":
            rlt_file=TString(Directory.Data()+"controlplot_"+"averageof3MassVV"+"_"+tag+".pdf");
            rlt_file1=TString(Directory.Data()+"controlplot_"+"averageof3MassVV"+"_"+tag+".png");

        canvas_controlplot.SaveAs(rlt_file.Data());
        canvas_controlplot.SaveAs(rlt_file1.Data());
  
        if logy:
            #canvas_controlplot.SetLogy() ;
            #fPads1.SetLogy();
            fPads1.Update();
            #fPads2.Update();
            canvas_controlplot.Update();
            rlt_file.ReplaceAll(".pdf","_log.pdf"); 
            canvas_controlplot.SaveAs(rlt_file.Data());




    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def make_controlplot_averageof2MassVV(self,numak8jets,variable,variable1,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="", logy=0 ,TTBarControl=0,):
        tmp_lumi=self.GetLumi()
        tmp_signal_scale="%s"%(options.scale1)
        tmp_signal_scale1="%s"%(options.scale2)
        
        weight_mc_forSignal="weight*%s*%s"%(tmp_lumi, tmp_signal_scale);
        weight_mc_forSignal1="weight*%s*%s"%(tmp_lumi, tmp_signal_scale1);


        weight_mc_forV="weight*%s"%(tmp_lumi);
        weight_mc_forT="weight*%s"%(tmp_lumi);
        ##weight_mc_forV="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forV.getVal());#little error rrv_wtagger_eff_reweight_forV
        ##weight_mc_forT="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forT.getVal());#little error rrv_wtagger_eff_reweight_forT
        weight_mc_forG="weight*%s"%(tmp_lumi); #General

        tmp_WJets_scale=0.8
        tmp_TTBar_scale=1.0

        ##if TTBarControl==0: 
        ##    if self.channel=="mu": tmp_WJets_scale=1.2.18
        ##    if self.channel=="el": tmp_WJets_scale=1.2.01 
        ##else: 
        ##    if self.channel=="mu": tmp_TTBar_scale=0.85
        ##    if self.channel=="el": tmp_TTBar_scale=0.70

        if self.channel=="mu": tmp_WJets_scale=1.2
        if self.channel=="el": tmp_WJets_scale=0.97
        if self.channel=="mu": tmp_TTBar_scale=1
        if self.channel=="el": tmp_TTBar_scale=0.81

        weight_mc_forWJets="weight*%s*%s"%(tmp_lumi, tmp_WJets_scale); #General
        weight_mc_forTTBar="weight*%s*%s"%(tmp_lumi, tmp_TTBar_scale); #General

        weightcut_mc_forSignal="(%s)*(%s)"%(weight_mc_forSignal,cut);
        weightcut_mc_forSignala="(%s)*(%s)"%(weight_mc_forSignal,cut1);
        weightcut_mc_forSignal1="(%s)*(%s)"%(weight_mc_forSignal1,cut);
        weightcut_mc_forSignala1="(%s)*(%s)"%(weight_mc_forSignal1,cut1);

        weightcut_mc_forV="(%s)*(%s)"%(weight_mc_forV,cut);
        weightcut_mc_forT="(%s)*(%s)"%(weight_mc_forT,cut);
        weightcut_mc_forG="(%s)*(%s)"%(weight_mc_forG,cut);
        weightcut_mc_forWJets="(%s)*(%s)"%(weight_mc_forWJets,cut);
        weightcut_mc_forTTBar="(%s)*(%s)"%(weight_mc_forTTBar,cut);
        weightcut_data="%s"%(cut);
        print "weightcut_mc_forV="+weightcut_mc_forV;
        print "weightcut_mc_forT="+weightcut_mc_forT;
        print "weightcut_mc_forG="+weightcut_mc_forG;
        print "weightcut_mc_forWJets="+weightcut_mc_forWJets;
        print "weightcut_mc_forTTBar="+weightcut_mc_forTTBar;
        hist_data =TH1D("hist_data","hist_data"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal =TH1D("hist_Signal","hist_Signal"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal.Sumw2();
        hist_Signal2 =TH1D("hist_Signal2","hist_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal2.Sumw2();
        hist_Signal3 =TH1D("hist_Signal3","hist_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3.Sumw2();
        hist_Signal4 =TH1D("hist_Signal4","hist_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4.Sumw2();
        hist_Signal3a =TH1D("hist_Signal3a","hist_Signal3a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3a.Sumw2();
        hist_Signal4a =TH1D("hist_Signal4a","hist_Signal4a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4a.Sumw2();
        hist_WJets=TH1D("hist_WJets","hist_WJets"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_WJets.Sumw2();
        hist_TTbar=TH1D("hist_TTbar","hist_TTbar"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TTbar.Sumw2();
        hist_STop =TH1D("hist_STop","hist_STop"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_STop.Sumw2();
        hist_VV   =TH1D("hist_VV","hist_VV"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_VV.Sumw2();
        hist_TotalMC =TH1D("hist_TotalMC","hist_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TotalMC.Sumw2();


        hstack_TotalMC = THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle))
        if TTBarControl==0:
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);
            
            hstack_TotalMC.Add(hist_WJets); 
        else:
            hstack_TotalMC.Add(hist_WJets); 
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);


        hist_data.SetLineColor(self.color_palet["data"]); hist_data.SetFillColor(self.color_palet["data"]);
        hist_Signal.SetLineColor(self.color_palet["Signal"]); hist_Signal.SetFillColor(self.color_palet["Signal"]); hist_Signal.SetFillStyle(0);hist_Signal.SetLineWidth(2);
        hist_Signal2.SetLineColor(kBlue); hist_Signal2.SetFillColor(self.color_palet["Signal"]); hist_Signal2.SetFillStyle(0);hist_Signal2.SetLineWidth(2);
        hist_Signal3.SetLineColor(kGreen); hist_Signal3.SetFillColor(self.color_palet["Signal"]); hist_Signal3.SetFillStyle(0);hist_Signal3.SetLineWidth(2);
        hist_Signal4.SetLineColor(kOrange); hist_Signal4.SetFillColor(self.color_palet["Signal"]); hist_Signal4.SetFillStyle(0);hist_Signal4.SetLineWidth(2);
        hist_Signal3a.SetLineColor(kGreen); hist_Signal3a.SetFillColor(self.color_palet["Signal"]); hist_Signal3a.SetFillStyle(0);hist_Signal3a.SetLineWidth(2);hist_Signal3a.SetLineStyle(2);
        hist_Signal4a.SetLineColor(kOrange); hist_Signal4a.SetFillColor(self.color_palet["Signal"]); hist_Signal4a.SetFillStyle(0);hist_Signal4a.SetLineWidth(2);hist_Signal4a.SetLineStyle(2);


        hist_WJets.SetLineColor(kBlack); hist_WJets.SetFillColor(self.color_palet["WJets"]);
        hist_TTbar.SetLineColor(kBlack); hist_TTbar.SetFillColor(self.color_palet["TTbar"]);
        hist_STop.SetLineColor(kBlack); hist_STop.SetFillColor(self.color_palet["STop"]);
        hist_VV.SetLineColor(kBlack); hist_VV.SetFillColor(self.color_palet["VV"]);



        tree_data   =TChain("PKUTree");  tree_data.Add(self.file_Directory+self.file_data);
        tree_Signal =TChain("PKUTree");  tree_Signal.Add(self.file_Directory+self.file_signal);
        tree_Signal2 =TChain("PKUTree");  tree_Signal2.Add(self.file_Directory+self.file_signal2);
        tree_Signal3 =TChain("PKUTree");  tree_Signal3.Add(self.file_Directory+self.file_signal3);
        tree_Signal4 =TChain("PKUTree");  tree_Signal4.Add(self.file_Directory+self.file_signal4);

        tree_WJets  =TChain("PKUTree");tree_WJets.Add(self.file_Directory+self.file_WJets0_mc);
        tree_TTbar  =TChain("PKUTree");tree_TTbar.Add(self.file_Directory+self.file_TTbar_mc);
        tree_STop   =TChain("PKUTree");  tree_STop.Add(self.file_Directory+self.file_STop_mc);
        tree_VV     =TChain("PKUTree");      tree_VV.Add(self.file_Directory+self.file_VV_mc);

        
        tree_data.Draw("(abs(%s)+abs(%s))/2 >> hist_data"%(variable,variable1), weightcut_data);
        tree_Signal.Draw("(abs(%s)+abs(%s))/2 >> hist_Signal"%(variable,variable1), weightcut_mc_forSignal);
        tree_Signal2.Draw("(abs(%s)+abs(%s))/2 >> hist_Signal2"%(variable,variable1), weightcut_mc_forSignal1);
        tree_Signal3.Draw("(abs(%s)+abs(%s))/2 >> hist_Signal3"%(variable,variable1), weightcut_mc_forSignal);
        tree_Signal4.Draw("(abs(%s)+abs(%s))/2 >> hist_Signal4"%(variable,variable1), weightcut_mc_forSignal1);
        tree_Signal3.Draw("(abs(%s)+abs(%s))/2 >> hist_Signal3a"%(variable,variable1), weightcut_mc_forSignala);
        tree_Signal4.Draw("(abs(%s)+abs(%s))/2 >> hist_Signal4a"%(variable,variable1), weightcut_mc_forSignala1);
        
        #tree_WJets.Draw("%s >> hist_WJets"%(variable), weightcut_mc_forG);
        tree_WJets.Draw("(abs(%s)+abs(%s))/2 >> hist_WJets"%(variable,variable1), weightcut_mc_forWJets);
        #tree_TTbar.Draw("%s >> hist_TTbar"%(variable), weightcut_mc_forT);
        tree_TTbar.Draw("(abs(%s)+abs(%s))/2 >> hist_TTbar"%(variable,variable1), weightcut_mc_forTTBar);
        tree_STop.Draw("(abs(%s)+abs(%s))/2 >> hist_STop"%(variable,variable1), weightcut_mc_forT);
        tree_VV.Draw("(abs(%s)+abs(%s))/2 >> hist_VV"%(variable,variable1), weightcut_mc_forV);

 
       
        hist_WJets=UnderOverFlow1D(hist_WJets);
        hist_TTbar=UnderOverFlow1D(hist_TTbar);
        hist_STop=UnderOverFlow1D(hist_STop);
        hist_VV=UnderOverFlow1D(hist_VV);
        hist_Signal=UnderOverFlow1D(hist_Signal);
        hist_Signal2=UnderOverFlow1D(hist_Signal2);
        hist_Signal3=UnderOverFlow1D(hist_Signal3);
        hist_Signal4=UnderOverFlow1D(hist_Signal4);
        hist_Signal3a=UnderOverFlow1D(hist_Signal3a);
        hist_Signal4a=UnderOverFlow1D(hist_Signal4a);

        hist_TotalMC.Add(hist_WJets); hist_TotalMC.Add(hist_TTbar); hist_TotalMC.Add(hist_STop); hist_TotalMC.Add(hist_VV);


        canvas_controlplot = TCanvas("canvas_controlplot"+variable,"canvas_controlplot"+variable, 500,500);

        fPads1 = TPad("pad1", "Run2", 0.0, 0.0, 1.00, 1.00);
        #fPads2 = TPad("pad2", "", 0.00, 0.00, 1.00, 0.28);
        #canvas_controlplot.cd();
        fPads1.SetBottomMargin(0.07);
        fPads1.SetLeftMargin(0.10);
        fPads1.SetRightMargin(0.03);
##fPads2.SetTopMargin(0);
##      fPads2.SetBottomMargin(0.25);

        fPads1.Draw();
##       fPads2.Draw();
        fPads1.cd();

            
        
        hist_data.SetBinErrorOption(TH1D.kPoisson);
        #hist_data.Draw("e");
        #hstack_TotalMC.GetXaxis().SetTitle("%s"%(xtitle));
        #hist_data.Draw("same e");
        hist_Signal.Sumw2();
        hist_Signal2.Sumw2();
        hist_Signal3.Sumw2();
        hist_Signal4.Sumw2();
        hist_Signal3a.Sumw2();
        hist_Signal4a.Sumw2();

        hist_WJets.Sumw2();
        


        #hist_Signal.Scale(200);
        #hist_Signal2.Scale(200);
        #hist_Signal3a.Scale(200);
        #hist_Signal4a.Scale(200);
        #hist_Signal3.Scale(200);
        #hist_Signal4.Scale(200);


            #hist_Signal3a.Scale(hist_Signal3.GetMaximum()/hist_Signal3a.GetMaximum());
            #hist_Signal4a.Scale(hist_Signal4.GetMaximum()/hist_Signal4a.GetMaximum());

        histsigmax=TMath.Max(hist_Signal.GetMaximum(),hist_Signal2.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal3.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal4.GetMaximum());

        hist_Signal.GetYaxis().SetRangeUser(1,TMath.Max(histsigmax,hist_TotalMC.GetMaximum())*1.2);

        hist_Signal.Draw("HIST");
        hstack_TotalMC.Draw("same HIST");
        hist_Signal.Draw("same HIST");

        hist_Signal2.Draw("same HIST");
        hist_Signal3.Draw("same HIST");
        hist_Signal4.Draw("same HIST");
        hist_Signal3a.Draw("same HIST");
        hist_Signal4a.Draw("same HIST");



        hist_data.GetXaxis().SetTitleOffset(1.2);
        hist_data.GetYaxis().SetTitleOffset(1.3);
        hist_data.GetYaxis().SetTitleSize(0.07);
        hist_data.GetXaxis().SetTitleSize(0.08);
        hist_data.GetXaxis().SetLabelSize(0.06);
        hist_data.GetYaxis().SetLabelSize(0.06);


        banner = TLatex(0.95, 0.96, "35.9 fb^{-1} (13 TeV)");
        banner.SetNDC(); banner.SetTextSize(0.038); banner.SetTextFont(42); banner.SetTextAlign(31); banner.SetLineWidth(2); banner.Draw();
        CMStext = TLatex(0.10,0.96,"CMS");
        CMStext.SetNDC(); CMStext.SetTextSize(0.041); CMStext.SetTextFont(61); CMStext.SetTextAlign(11); CMStext.SetLineWidth(2); CMStext.Draw();
#    if self.channel=="el":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow e#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
#    elif self.channel=="mu":
        Extratext = TLatex(0.20, 0.96, "Preliminary");        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if(numak8jets=="1"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="1or2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1or2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
#    elif self.channel=="em":

#    elif self.channel=="em":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow l#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();

        theLeg = TLegend(0.51, 0.75, 0.76, 0.91, "", "NDC");
        theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);
        theLeg.SetTextSize(.025);
        # theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.025);
        theLeg.SetTextFont(42);

#theLeg.AddEntry(hist_data, "Data W#rightarrow #mu#nu","ep");
        theLeg.AddEntry(hist_WJets, "W+jets","F");
        theLeg.AddEntry(hist_VV, "WW/WZ","F");
        theLeg.AddEntry(hist_TTbar, "t#bar{t}","F");
        theLeg.AddEntry(hist_STop, "Single Top","F");
        #theLeg.AddEntry(gr_MCStat, "Sys.","F");
        theLeg.AddEntry(hist_Signal, "Signal(3000,1500) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal2, "Signal(1500,750) (#times %s)"%(tmp_signal_scale1),"L");
        theLeg.AddEntry(hist_Signal3, "Signal(3000,180) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal3a, "(3000,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal4, "Signal(1500,180) (#times %s)"%(tmp_signal_scale1),"L");

        theLeg.AddEntry(hist_Signal4a, "(1500,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale1),"L");


        theLeg.SetY1NDC(0.9 - 0.05*6 - 0.005);
        theLeg.SetY1(theLeg.GetY1NDC());
        theLeg.Draw();
        



        lineAtZero = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.0, hist_TotalMC.GetXaxis().GetXmax(), 1.0);
        lineAtZero.SetLineColor(2);
#lineAtZero.Draw();
        lineAtPlusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.5, hist_TotalMC.GetXaxis().GetXmax(), 1.5);
        lineAtPlusTwo.SetLineColor(2);
        lineAtPlusTwo.SetLineStyle(2);
#       lineAtPlusTwo.Draw();
        lineAtMinusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 0.5, hist_TotalMC.GetXaxis().GetXmax(), 0.5);
        lineAtMinusTwo.SetLineColor(2);
        lineAtMinusTwo.SetLineStyle(2);
#      lineAtMinusTwo.Draw();


        Directory=TString("ST_more_than_%sGeV/plots_nj%s/"%(options.stcut,numak8jets));


        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()): os.system("mkdir -p  "+Directory.Data());
        #fPads1.SetLogy();
        fPads1.Update();
        #only draw png
        if variable=="MassVV[1]":
            rlt_file=TString(Directory.Data()+"controlplot_"+"averageof2MassVV"+"_"+tag+".pdf");
            rlt_file1=TString(Directory.Data()+"controlplot_"+"averageof2MassVV"+"_"+tag+".png");

        canvas_controlplot.SaveAs(rlt_file.Data());
        canvas_controlplot.SaveAs(rlt_file1.Data());
  
        if logy:
            #canvas_controlplot.SetLogy() ;
            #fPads1.SetLogy();
            fPads1.Update();
            #fPads2.Update();
            canvas_controlplot.Update();
            rlt_file.ReplaceAll(".pdf","_log.pdf"); 
            canvas_controlplot.SaveAs(rlt_file.Data());




    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def make_controlplot_productof3tau(self,numak8jets,variable,variable1,variable2,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="", logy=0 ,TTBarControl=0,):
        tmp_lumi=self.GetLumi()
        tmp_signal_scale="%s"%(options.scale1)
        tmp_signal_scale1="%s"%(options.scale2)
        
        weight_mc_forSignal="weight*%s*%s"%(tmp_lumi, tmp_signal_scale);
        weight_mc_forSignal1="weight*%s*%s"%(tmp_lumi, tmp_signal_scale1);


        weight_mc_forV="weight*%s"%(tmp_lumi);
        weight_mc_forT="weight*%s"%(tmp_lumi);
        ##weight_mc_forV="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forV.getVal());#little error rrv_wtagger_eff_reweight_forV
        ##weight_mc_forT="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forT.getVal());#little error rrv_wtagger_eff_reweight_forT
        weight_mc_forG="weight*%s"%(tmp_lumi); #General

        tmp_WJets_scale=0.8
        tmp_TTBar_scale=1.0

        ##if TTBarControl==0: 
        ##    if self.channel=="mu": tmp_WJets_scale=1.2.18
        ##    if self.channel=="el": tmp_WJets_scale=1.2.01 
        ##else: 
        ##    if self.channel=="mu": tmp_TTBar_scale=0.85
        ##    if self.channel=="el": tmp_TTBar_scale=0.70

        if self.channel=="mu": tmp_WJets_scale=1.2
        if self.channel=="el": tmp_WJets_scale=0.97
        if self.channel=="mu": tmp_TTBar_scale=1
        if self.channel=="el": tmp_TTBar_scale=0.81

        weight_mc_forWJets="weight*%s*%s"%(tmp_lumi, tmp_WJets_scale); #General
        weight_mc_forTTBar="weight*%s*%s"%(tmp_lumi, tmp_TTBar_scale); #General

        weightcut_mc_forSignal="(%s)*(%s)"%(weight_mc_forSignal,cut);
        weightcut_mc_forSignala="(%s)*(%s)"%(weight_mc_forSignal,cut1);
        weightcut_mc_forSignal1="(%s)*(%s)"%(weight_mc_forSignal1,cut);
        weightcut_mc_forSignala1="(%s)*(%s)"%(weight_mc_forSignal1,cut1);
        
        weightcut_mc_forV="(%s)*(%s)"%(weight_mc_forV,cut);
        weightcut_mc_forT="(%s)*(%s)"%(weight_mc_forT,cut);
        weightcut_mc_forG="(%s)*(%s)"%(weight_mc_forG,cut);
        weightcut_mc_forWJets="(%s)*(%s)"%(weight_mc_forWJets,cut);
        weightcut_mc_forTTBar="(%s)*(%s)"%(weight_mc_forTTBar,cut);
        weightcut_data="%s"%(cut);
        print "weightcut_mc_forV="+weightcut_mc_forV;
        print "weightcut_mc_forT="+weightcut_mc_forT;
        print "weightcut_mc_forG="+weightcut_mc_forG;
        print "weightcut_mc_forWJets="+weightcut_mc_forWJets;
        print "weightcut_mc_forTTBar="+weightcut_mc_forTTBar;
        hist_data =TH1D("hist_data","hist_data"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal =TH1D("hist_Signal","hist_Signal"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal.Sumw2();
        hist_Signal2 =TH1D("hist_Signal2","hist_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal2.Sumw2();
        hist_Signal3 =TH1D("hist_Signal3","hist_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3.Sumw2();
        hist_Signal4 =TH1D("hist_Signal4","hist_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4.Sumw2();
        hist_Signal3a =TH1D("hist_Signal3a","hist_Signal3a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3a.Sumw2();
        hist_Signal4a =TH1D("hist_Signal4a","hist_Signal4a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4a.Sumw2();
        hist_WJets=TH1D("hist_WJets","hist_WJets"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_WJets.Sumw2();
        hist_TTbar=TH1D("hist_TTbar","hist_TTbar"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TTbar.Sumw2();
        hist_STop =TH1D("hist_STop","hist_STop"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_STop.Sumw2();
        hist_VV   =TH1D("hist_VV","hist_VV"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_VV.Sumw2();
        hist_TotalMC =TH1D("hist_TotalMC","hist_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TotalMC.Sumw2();


        hstack_TotalMC = THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle))
        if TTBarControl==0:
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);
            
            hstack_TotalMC.Add(hist_WJets); 
        else:
            hstack_TotalMC.Add(hist_WJets); 
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);


        hist_data.SetLineColor(self.color_palet["data"]); hist_data.SetFillColor(self.color_palet["data"]);
        hist_Signal.SetLineColor(self.color_palet["Signal"]); hist_Signal.SetFillColor(self.color_palet["Signal"]); hist_Signal.SetFillStyle(0);hist_Signal.SetLineWidth(2);
        hist_Signal2.SetLineColor(kBlue); hist_Signal2.SetFillColor(self.color_palet["Signal"]); hist_Signal2.SetFillStyle(0);hist_Signal2.SetLineWidth(2);
        hist_Signal3.SetLineColor(kGreen); hist_Signal3.SetFillColor(self.color_palet["Signal"]); hist_Signal3.SetFillStyle(0);hist_Signal3.SetLineWidth(2);
        hist_Signal4.SetLineColor(kOrange); hist_Signal4.SetFillColor(self.color_palet["Signal"]); hist_Signal4.SetFillStyle(0);hist_Signal4.SetLineWidth(2);
        hist_Signal3a.SetLineColor(kGreen); hist_Signal3a.SetFillColor(self.color_palet["Signal"]); hist_Signal3a.SetFillStyle(0);hist_Signal3a.SetLineWidth(2);hist_Signal3a.SetLineStyle(2);
        hist_Signal4a.SetLineColor(kOrange); hist_Signal4a.SetFillColor(self.color_palet["Signal"]); hist_Signal4a.SetFillStyle(0);hist_Signal4a.SetLineWidth(2);hist_Signal4a.SetLineStyle(2);


        hist_WJets.SetLineColor(kBlack); hist_WJets.SetFillColor(self.color_palet["WJets"]);
        hist_TTbar.SetLineColor(kBlack); hist_TTbar.SetFillColor(self.color_palet["TTbar"]);
        hist_STop.SetLineColor(kBlack); hist_STop.SetFillColor(self.color_palet["STop"]);
        hist_VV.SetLineColor(kBlack); hist_VV.SetFillColor(self.color_palet["VV"]);



        tree_data   =TChain("PKUTree");  tree_data.Add(self.file_Directory+self.file_data);
        tree_Signal =TChain("PKUTree");  tree_Signal.Add(self.file_Directory+self.file_signal);
        tree_Signal2 =TChain("PKUTree");  tree_Signal2.Add(self.file_Directory+self.file_signal2);
        tree_Signal3 =TChain("PKUTree");  tree_Signal3.Add(self.file_Directory+self.file_signal3);
        tree_Signal4 =TChain("PKUTree");  tree_Signal4.Add(self.file_Directory+self.file_signal4);

        tree_WJets  =TChain("PKUTree");tree_WJets.Add(self.file_Directory+self.file_WJets0_mc);
        tree_TTbar  =TChain("PKUTree");tree_TTbar.Add(self.file_Directory+self.file_TTbar_mc);
        tree_STop   =TChain("PKUTree");  tree_STop.Add(self.file_Directory+self.file_STop_mc);
        tree_VV     =TChain("PKUTree");      tree_VV.Add(self.file_Directory+self.file_VV_mc);

        
        tree_data.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_data"%(variable,variable1,variable2), weightcut_data);
        tree_Signal.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_Signal"%(variable,variable1,variable2), weightcut_mc_forSignal);
        tree_Signal2.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_Signal2"%(variable,variable1,variable2), weightcut_mc_forSignal1);
        tree_Signal3.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_Signal3"%(variable,variable1,variable2), weightcut_mc_forSignal);
        tree_Signal4.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_Signal4"%(variable,variable1,variable2), weightcut_mc_forSignal1);
        tree_Signal3.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_Signal3a"%(variable,variable1,variable2), weightcut_mc_forSignala);
        tree_Signal4.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_Signal4a"%(variable,variable1,variable2), weightcut_mc_forSignala1);
        #tree_WJets.Draw("%s >> hist_WJets"%(variable), weightcut_mc_forG);
        tree_WJets.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_WJets"%(variable,variable1,variable2), weightcut_mc_forWJets);
        #tree_TTbar.Draw("%s >> hist_TTbar"%(variable), weightcut_mc_forT);
        tree_TTbar.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_TTbar"%(variable,variable1,variable2), weightcut_mc_forTTBar);
        tree_STop.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_STop"%(variable,variable1,variable2), weightcut_mc_forT);
        tree_VV.Draw("pow(abs(%s)*abs(%s)*abs(%s),1./3) >> hist_VV"%(variable,variable1,variable2), weightcut_mc_forV);

 
       
        hist_WJets=UnderOverFlow1D(hist_WJets);
        hist_TTbar=UnderOverFlow1D(hist_TTbar);
        hist_STop=UnderOverFlow1D(hist_STop);
        hist_VV=UnderOverFlow1D(hist_VV);
        hist_Signal=UnderOverFlow1D(hist_Signal);
        hist_Signal2=UnderOverFlow1D(hist_Signal2);
        hist_Signal3=UnderOverFlow1D(hist_Signal3);
        hist_Signal4=UnderOverFlow1D(hist_Signal4);
        hist_Signal3a=UnderOverFlow1D(hist_Signal3a);
        hist_Signal4a=UnderOverFlow1D(hist_Signal4a);

        hist_TotalMC.Add(hist_WJets); hist_TotalMC.Add(hist_TTbar); hist_TotalMC.Add(hist_STop); hist_TotalMC.Add(hist_VV);


        canvas_controlplot = TCanvas("canvas_controlplot"+variable,"canvas_controlplot"+variable, 500,500);

        fPads1 = TPad("pad1", "Run2", 0.0, 0.0, 1.00, 1.00);
        #fPads2 = TPad("pad2", "", 0.00, 0.00, 1.00, 0.28);
        #canvas_controlplot.cd();
        fPads1.SetBottomMargin(0.07);
        fPads1.SetLeftMargin(0.10);
        fPads1.SetRightMargin(0.03);
##fPads2.SetTopMargin(0);
##      fPads2.SetBottomMargin(0.25);

        fPads1.Draw();
##       fPads2.Draw();
        fPads1.cd();

            
        
        hist_data.SetBinErrorOption(TH1D.kPoisson);
        #hist_data.Draw("e");
        #hstack_TotalMC.GetXaxis().SetTitle("%s"%(xtitle));
        #hist_data.Draw("same e");
        hist_Signal.Sumw2();
        hist_Signal2.Sumw2();
        hist_Signal3.Sumw2();
        hist_Signal4.Sumw2();
        hist_Signal3a.Sumw2();
        hist_Signal4a.Sumw2();

        hist_WJets.Sumw2();
        


        #hist_Signal.Scale(200);
        #hist_Signal2.Scale(200);
        #hist_Signal3a.Scale(200);
        #hist_Signal4a.Scale(200);
        #hist_Signal3.Scale(200);
        #hist_Signal4.Scale(200);


            #hist_Signal3a.Scale(hist_Signal3.GetMaximum()/hist_Signal3a.GetMaximum());
            #hist_Signal4a.Scale(hist_Signal4.GetMaximum()/hist_Signal4a.GetMaximum());

        histsigmax=TMath.Max(hist_Signal.GetMaximum(),hist_Signal2.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal3.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal4.GetMaximum());

        hist_Signal.GetYaxis().SetRangeUser(1,TMath.Max(histsigmax,hist_TotalMC.GetMaximum())*1.2);

        hist_Signal.Draw("HIST");
        hstack_TotalMC.Draw("same HIST");
        hist_Signal.Draw("same HIST");

        hist_Signal2.Draw("same HIST");
        hist_Signal3.Draw("same HIST");
        hist_Signal4.Draw("same HIST");
        hist_Signal3a.Draw("same HIST");
        hist_Signal4a.Draw("same HIST");



        hist_data.GetXaxis().SetTitleOffset(1.2);
        hist_data.GetYaxis().SetTitleOffset(1.3);
        hist_data.GetYaxis().SetTitleSize(0.07);
        hist_data.GetXaxis().SetTitleSize(0.08);
        hist_data.GetXaxis().SetLabelSize(0.06);
        hist_data.GetYaxis().SetLabelSize(0.06);


        banner = TLatex(0.95, 0.96, "35.9 fb^{-1} (13 TeV)");
        banner.SetNDC(); banner.SetTextSize(0.038); banner.SetTextFont(42); banner.SetTextAlign(31); banner.SetLineWidth(2); banner.Draw();
        CMStext = TLatex(0.10,0.96,"CMS");
        CMStext.SetNDC(); CMStext.SetTextSize(0.041); CMStext.SetTextFont(61); CMStext.SetTextAlign(11); CMStext.SetLineWidth(2); CMStext.Draw();
#    if self.channel=="el":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow e#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
#    elif self.channel=="mu":
        Extratext = TLatex(0.20, 0.96, "Preliminary");        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if(numak8jets=="1"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="1or2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1or2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
#    elif self.channel=="em":

#    elif self.channel=="em":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow l#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();

        theLeg = TLegend(0.51, 0.75, 0.76, 0.91, "", "NDC");
        theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);
        theLeg.SetTextSize(.025);
        # theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.025);
        theLeg.SetTextFont(42);

#theLeg.AddEntry(hist_data, "Data W#rightarrow #mu#nu","ep");
        theLeg.AddEntry(hist_WJets, "W+jets","F");
        theLeg.AddEntry(hist_VV, "WW/WZ","F");
        theLeg.AddEntry(hist_TTbar, "t#bar{t}","F");
        theLeg.AddEntry(hist_STop, "Single Top","F");
        #theLeg.AddEntry(gr_MCStat, "Sys.","F");
        theLeg.AddEntry(hist_Signal, "Signal(3000,1500) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal2, "Signal(1500,750) (#times %s)"%(tmp_signal_scale1),"L");
        theLeg.AddEntry(hist_Signal3, "Signal(3000,180) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal3a, "(3000,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal4, "Signal(1500,180) (#times %s)"%(tmp_signal_scale1),"L");

        theLeg.AddEntry(hist_Signal4a, "(1500,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale1),"L");


        theLeg.SetY1NDC(0.9 - 0.05*6 - 0.005);
        theLeg.SetY1(theLeg.GetY1NDC());
        theLeg.Draw();
        



        lineAtZero = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.0, hist_TotalMC.GetXaxis().GetXmax(), 1.0);
        lineAtZero.SetLineColor(2);
#lineAtZero.Draw();
        lineAtPlusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.5, hist_TotalMC.GetXaxis().GetXmax(), 1.5);
        lineAtPlusTwo.SetLineColor(2);
        lineAtPlusTwo.SetLineStyle(2);
#       lineAtPlusTwo.Draw();
        lineAtMinusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 0.5, hist_TotalMC.GetXaxis().GetXmax(), 0.5);
        lineAtMinusTwo.SetLineColor(2);
        lineAtMinusTwo.SetLineStyle(2);
#      lineAtMinusTwo.Draw();


        Directory=TString("ST_more_than_%sGeV/plots_nj%s/"%(options.stcut,numak8jets));


        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()): os.system("mkdir -p  "+Directory.Data());

        #only draw png
        if variable=="jet_tau2tau1_puppi":
            rlt_file=TString(Directory.Data()+"controlplot_"+"productof3tau"+"_"+tag+".pdf");
            rlt_file1=TString(Directory.Data()+"controlplot_"+"productof3tau"+"_"+tag+".png");
        if variable=="W_pt":
            rlt_file=TString(Directory.Data()+"controlplot_"+"productof3Wpt"+"_"+tag+".pdf");
            rlt_file1=TString(Directory.Data()+"controlplot_"+"productof3Wpt"+"_"+tag+".png");
        #fPads1.SetLogy();
        fPads1.Update();

        canvas_controlplot.SaveAs(rlt_file.Data());
        canvas_controlplot.SaveAs(rlt_file1.Data());

        if logy:
            #canvas_controlplot.SetLogy() ;
            #fPads1.SetLogy();
            fPads1.Update();
            #fPads2.Update();
            canvas_controlplot.Update();
            rlt_file.ReplaceAll(".pdf","_log.pdf"); 
            canvas_controlplot.SaveAs(rlt_file.Data());

    def make_controlplot_productof2Wpt(self,numak8jets,variable,variable1,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="", logy=0 ,TTBarControl=0,):
        tmp_lumi=self.GetLumi()
        tmp_signal_scale="%s"%(options.scale1)
        tmp_signal_scale1="%s"%(options.scale2)
        
        weight_mc_forSignal="weight*%s*%s"%(tmp_lumi, tmp_signal_scale);
        weight_mc_forSignal1="weight*%s*%s"%(tmp_lumi, tmp_signal_scale1);


        weight_mc_forV="weight*%s"%(tmp_lumi);
        weight_mc_forT="weight*%s"%(tmp_lumi);
        ##weight_mc_forV="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forV.getVal());#little error rrv_wtagger_eff_reweight_forV
        ##weight_mc_forT="weight*%s*%s"%(tmp_lumi,self.rrv_htagger_eff_reweight_forT.getVal());#little error rrv_wtagger_eff_reweight_forT
        weight_mc_forG="weight*%s"%(tmp_lumi); #General

        tmp_WJets_scale=0.8
        tmp_TTBar_scale=1.0

        ##if TTBarControl==0: 
        ##    if self.channel=="mu": tmp_WJets_scale=1.2.18
        ##    if self.channel=="el": tmp_WJets_scale=1.2.01 
        ##else: 
        ##    if self.channel=="mu": tmp_TTBar_scale=0.85
        ##    if self.channel=="el": tmp_TTBar_scale=0.70

        if self.channel=="mu": tmp_WJets_scale=1.2
        if self.channel=="el": tmp_WJets_scale=0.97
        if self.channel=="mu": tmp_TTBar_scale=1
        if self.channel=="el": tmp_TTBar_scale=0.81

        weight_mc_forWJets="weight*%s*%s"%(tmp_lumi, tmp_WJets_scale); #General
        weight_mc_forTTBar="weight*%s*%s"%(tmp_lumi, tmp_TTBar_scale); #General

        weightcut_mc_forSignal="(%s)*(%s)"%(weight_mc_forSignal,cut);
        weightcut_mc_forSignala="(%s)*(%s)"%(weight_mc_forSignal,cut1);
        weightcut_mc_forSignal1="(%s)*(%s)"%(weight_mc_forSignal1,cut);
        weightcut_mc_forSignala1="(%s)*(%s)"%(weight_mc_forSignal1,cut1);

        weightcut_mc_forV="(%s)*(%s)"%(weight_mc_forV,cut);
        weightcut_mc_forT="(%s)*(%s)"%(weight_mc_forT,cut);
        weightcut_mc_forG="(%s)*(%s)"%(weight_mc_forG,cut);
        weightcut_mc_forWJets="(%s)*(%s)"%(weight_mc_forWJets,cut);
        weightcut_mc_forTTBar="(%s)*(%s)"%(weight_mc_forTTBar,cut);
        weightcut_data="%s"%(cut);
        print "weightcut_mc_forV="+weightcut_mc_forV;
        print "weightcut_mc_forT="+weightcut_mc_forT;
        print "weightcut_mc_forG="+weightcut_mc_forG;
        print "weightcut_mc_forWJets="+weightcut_mc_forWJets;
        print "weightcut_mc_forTTBar="+weightcut_mc_forTTBar;
        hist_data =TH1D("hist_data","hist_data"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal =TH1D("hist_Signal","hist_Signal"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal.Sumw2();
        hist_Signal2 =TH1D("hist_Signal2","hist_Signal2"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal2.Sumw2();
        hist_Signal3 =TH1D("hist_Signal3","hist_Signal3"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3.Sumw2();
        hist_Signal4 =TH1D("hist_Signal4","hist_Signal4"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4.Sumw2();
        hist_Signal3a =TH1D("hist_Signal3a","hist_Signal3a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal3a.Sumw2();
        hist_Signal4a =TH1D("hist_Signal4a","hist_Signal4a"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_Signal4a.Sumw2();
        hist_WJets=TH1D("hist_WJets","hist_WJets"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_WJets.Sumw2();
        hist_TTbar=TH1D("hist_TTbar","hist_TTbar"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TTbar.Sumw2();
        hist_STop =TH1D("hist_STop","hist_STop"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_STop.Sumw2();
        hist_VV   =TH1D("hist_VV","hist_VV"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_VV.Sumw2();
        hist_TotalMC =TH1D("hist_TotalMC","hist_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);
        hist_TotalMC.Sumw2();


        hstack_TotalMC = THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle))
        if TTBarControl==0:
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);
            
            hstack_TotalMC.Add(hist_WJets); 
        else:
            hstack_TotalMC.Add(hist_WJets); 
            hstack_TotalMC.Add(hist_VV);
            hstack_TotalMC.Add(hist_STop);
            hstack_TotalMC.Add(hist_TTbar);


        hist_data.SetLineColor(self.color_palet["data"]); hist_data.SetFillColor(self.color_palet["data"]);
        hist_Signal.SetLineColor(self.color_palet["Signal"]); hist_Signal.SetFillColor(self.color_palet["Signal"]); hist_Signal.SetFillStyle(0);hist_Signal.SetLineWidth(2);
        hist_Signal2.SetLineColor(kBlue); hist_Signal2.SetFillColor(self.color_palet["Signal"]); hist_Signal2.SetFillStyle(0);hist_Signal2.SetLineWidth(2);
        hist_Signal3.SetLineColor(kGreen); hist_Signal3.SetFillColor(self.color_palet["Signal"]); hist_Signal3.SetFillStyle(0);hist_Signal3.SetLineWidth(2);
        hist_Signal4.SetLineColor(kOrange); hist_Signal4.SetFillColor(self.color_palet["Signal"]); hist_Signal4.SetFillStyle(0);hist_Signal4.SetLineWidth(2);
        hist_Signal3a.SetLineColor(kGreen); hist_Signal3a.SetFillColor(self.color_palet["Signal"]); hist_Signal3a.SetFillStyle(0);hist_Signal3a.SetLineWidth(2);hist_Signal3a.SetLineStyle(2);
        hist_Signal4a.SetLineColor(kOrange); hist_Signal4a.SetFillColor(self.color_palet["Signal"]); hist_Signal4a.SetFillStyle(0);hist_Signal4a.SetLineWidth(2);hist_Signal4a.SetLineStyle(2);


        hist_WJets.SetLineColor(kBlack); hist_WJets.SetFillColor(self.color_palet["WJets"]);
        hist_TTbar.SetLineColor(kBlack); hist_TTbar.SetFillColor(self.color_palet["TTbar"]);
        hist_STop.SetLineColor(kBlack); hist_STop.SetFillColor(self.color_palet["STop"]);
        hist_VV.SetLineColor(kBlack); hist_VV.SetFillColor(self.color_palet["VV"]);



        tree_data   =TChain("PKUTree");  tree_data.Add(self.file_Directory+self.file_data);
        tree_Signal =TChain("PKUTree");  tree_Signal.Add(self.file_Directory+self.file_signal);
        tree_Signal2 =TChain("PKUTree");  tree_Signal2.Add(self.file_Directory+self.file_signal2);
        tree_Signal3 =TChain("PKUTree");  tree_Signal3.Add(self.file_Directory+self.file_signal3);
        tree_Signal4 =TChain("PKUTree");  tree_Signal4.Add(self.file_Directory+self.file_signal4);

        tree_WJets  =TChain("PKUTree");tree_WJets.Add(self.file_Directory+self.file_WJets0_mc);
        tree_TTbar  =TChain("PKUTree");tree_TTbar.Add(self.file_Directory+self.file_TTbar_mc);
        tree_STop   =TChain("PKUTree");  tree_STop.Add(self.file_Directory+self.file_STop_mc);
        tree_VV     =TChain("PKUTree");      tree_VV.Add(self.file_Directory+self.file_VV_mc);


        if(numak8jets=="2" or numak8jets=="1" or numak8jets=="1or2" ):
            tree_data.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_data"%(variable,variable1), weightcut_data);
            tree_Signal.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_Signal"%(variable,variable1), weightcut_mc_forSignal);
            tree_Signal2.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_Signal2"%(variable,variable1), weightcut_mc_forSignal1);
            tree_Signal3.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_Signal3"%(variable,variable1), weightcut_mc_forSignal);
            tree_Signal4.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_Signal4"%(variable,variable1), weightcut_mc_forSignal1);
            tree_Signal3.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_Signal3a"%(variable,variable1), weightcut_mc_forSignala);
            tree_Signal4.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_Signal4a"%(variable,variable1), weightcut_mc_forSignala1);
            #tree_WJets.Draw("%s >> hist_WJets"%(variable), weightcut_mc_forG);
            tree_WJets.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_WJets"%(variable,variable1), weightcut_mc_forWJets);
            #tree_TTbar.Draw("%s >> hist_TTbar"%(variable), weightcut_mc_forT);
            tree_TTbar.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_TTbar"%(variable,variable1), weightcut_mc_forTTBar);
            tree_STop.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_STop"%(variable,variable1), weightcut_mc_forT);
            tree_VV.Draw("pow(abs(%s)*abs(%s),1./2) >> hist_VV"%(variable,variable1), weightcut_mc_forV);

 
       
        hist_WJets=UnderOverFlow1D(hist_WJets);
        hist_TTbar=UnderOverFlow1D(hist_TTbar);
        hist_STop=UnderOverFlow1D(hist_STop);
        hist_VV=UnderOverFlow1D(hist_VV);
        hist_Signal=UnderOverFlow1D(hist_Signal);
        hist_Signal2=UnderOverFlow1D(hist_Signal2);
        hist_Signal3=UnderOverFlow1D(hist_Signal3);
        hist_Signal4=UnderOverFlow1D(hist_Signal4);
        hist_Signal3a=UnderOverFlow1D(hist_Signal3a);
        hist_Signal4a=UnderOverFlow1D(hist_Signal4a);

        hist_TotalMC.Add(hist_WJets); hist_TotalMC.Add(hist_TTbar); hist_TotalMC.Add(hist_STop); hist_TotalMC.Add(hist_VV);


        canvas_controlplot = TCanvas("canvas_controlplot"+variable,"canvas_controlplot"+variable, 500,500);

        fPads1 = TPad("pad1", "Run2", 0.0, 0.0, 1.00, 1.00);
        #fPads2 = TPad("pad2", "", 0.00, 0.00, 1.00, 0.28);
        #canvas_controlplot.cd();
        fPads1.SetBottomMargin(0.07);
        fPads1.SetLeftMargin(0.10);
        fPads1.SetRightMargin(0.03);
##fPads2.SetTopMargin(0);
##      fPads2.SetBottomMargin(0.25);

        fPads1.Draw();
##       fPads2.Draw();
        fPads1.cd();

            
        
        hist_data.SetBinErrorOption(TH1D.kPoisson);
        #hist_data.Draw("e");
        #hstack_TotalMC.GetXaxis().SetTitle("%s"%(xtitle));
        #hist_data.Draw("same e");
        hist_Signal.Sumw2();
        hist_Signal2.Sumw2();
        hist_Signal3.Sumw2();
        hist_Signal4.Sumw2();
        hist_Signal3a.Sumw2();
        hist_Signal4a.Sumw2();

        hist_WJets.Sumw2();
        


        #hist_Signal.Scale(200);
        #hist_Signal2.Scale(200);
        #hist_Signal3a.Scale(200);
        #hist_Signal4a.Scale(200);
        #hist_Signal3.Scale(200);
        #hist_Signal4.Scale(200);


            #hist_Signal3a.Scale(hist_Signal3.GetMaximum()/hist_Signal3a.GetMaximum());
            #hist_Signal4a.Scale(hist_Signal4.GetMaximum()/hist_Signal4a.GetMaximum());

        histsigmax=TMath.Max(hist_Signal.GetMaximum(),hist_Signal2.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal3.GetMaximum());
        histsigmax=TMath.Max(histsigmax,hist_Signal4.GetMaximum());

        hist_Signal.GetYaxis().SetRangeUser(1,TMath.Max(histsigmax,hist_TotalMC.GetMaximum())*1.2);

        hist_Signal.Draw("HIST");
        hstack_TotalMC.Draw("same HIST");
        hist_Signal.Draw("same HIST");

        hist_Signal2.Draw("same HIST");
        hist_Signal3.Draw("same HIST");
        hist_Signal4.Draw("same HIST");
        hist_Signal3a.Draw("same HIST");
        hist_Signal4a.Draw("same HIST");



        hist_data.GetXaxis().SetTitleOffset(1.2);
        hist_data.GetYaxis().SetTitleOffset(1.3);
        hist_data.GetYaxis().SetTitleSize(0.07);
        hist_data.GetXaxis().SetTitleSize(0.08);
        hist_data.GetXaxis().SetLabelSize(0.06);
        hist_data.GetYaxis().SetLabelSize(0.06);


        banner = TLatex(0.95, 0.96, "35.9 fb^{-1} (13 TeV)");
        banner.SetNDC(); banner.SetTextSize(0.038); banner.SetTextFont(42); banner.SetTextAlign(31); banner.SetLineWidth(2); banner.Draw();
        CMStext = TLatex(0.10,0.96,"CMS");
        CMStext.SetNDC(); CMStext.SetTextSize(0.041); CMStext.SetTextFont(61); CMStext.SetTextAlign(11); CMStext.SetLineWidth(2); CMStext.Draw();
#    if self.channel=="el":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow e#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
#    elif self.channel=="mu":
        Extratext = TLatex(0.20, 0.96, "Preliminary");        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if(numak8jets=="1"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
        if(numak8jets=="1or2"):
            Extratext1 = TLatex(0.37, 0.96, "N_{j} = 1or2");
            Extratext1.SetNDC(); Extratext1.SetTextSize(0.032); Extratext1.SetTextFont(52); Extratext1.SetTextAlign(11); Extratext1.SetLineWidth(2); Extratext1.Draw();
#    elif self.channel=="em":

#    elif self.channel=="em":
#        Extratext = TLatex(0.241, 0.96, "Preliminary W#rightarrow l#nu");
#        Extratext.SetNDC(); Extratext.SetTextSize(0.032); Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();

        theLeg = TLegend(0.51, 0.75, 0.76, 0.91, "", "NDC");
        theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);
        theLeg.SetTextSize(.025);
        # theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.025);
        theLeg.SetTextFont(42);

#theLeg.AddEntry(hist_data, "Data W#rightarrow #mu#nu","ep");
        theLeg.AddEntry(hist_WJets, "W+jets","F");
        theLeg.AddEntry(hist_VV, "WW/WZ","F");
        theLeg.AddEntry(hist_TTbar, "t#bar{t}","F");
        theLeg.AddEntry(hist_STop, "Single Top","F");
        #theLeg.AddEntry(gr_MCStat, "Sys.","F");
        theLeg.AddEntry(hist_Signal, "Signal(3000,1500) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal2, "Signal(1500,750) (#times %s)"%(tmp_signal_scale1),"L");
        theLeg.AddEntry(hist_Signal3, "Signal(3000,180) (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal3a, "(3000,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale),"L");
        theLeg.AddEntry(hist_Signal4, "Signal(1500,180) (#times %s)"%(tmp_signal_scale1),"L");

        theLeg.AddEntry(hist_Signal4a, "(1500,180) MERGED  with W(lv) as PR (#times %s)"%(tmp_signal_scale1),"L");


        theLeg.SetY1NDC(0.9 - 0.05*6 - 0.005);
        theLeg.SetY1(theLeg.GetY1NDC());
        theLeg.Draw();
        



        lineAtZero = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.0, hist_TotalMC.GetXaxis().GetXmax(), 1.0);
        lineAtZero.SetLineColor(2);
#lineAtZero.Draw();
        lineAtPlusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 1.5, hist_TotalMC.GetXaxis().GetXmax(), 1.5);
        lineAtPlusTwo.SetLineColor(2);
        lineAtPlusTwo.SetLineStyle(2);
#       lineAtPlusTwo.Draw();
        lineAtMinusTwo = TLine(hist_TotalMC.GetXaxis().GetXmin(), 0.5, hist_TotalMC.GetXaxis().GetXmax(), 0.5);
        lineAtMinusTwo.SetLineColor(2);
        lineAtMinusTwo.SetLineStyle(2);
#      lineAtMinusTwo.Draw();


        Directory=TString("ST_more_than_%sGeV/plots_nj%s/"%(options.stcut,numak8jets));


        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()): os.system("mkdir -p  "+Directory.Data());

        #only draw png
        rlt_file=TString(Directory.Data()+"controlplot_"+"productof2Wpt"+"_"+tag+".pdf");
        rlt_file1=TString(Directory.Data()+"controlplot_"+"Productof2Wpt"+"_"+tag+".png");
        #fPads1.SetLogy();
        fPads1.Update();

        canvas_controlplot.SaveAs(rlt_file.Data());
        canvas_controlplot.SaveAs(rlt_file1.Data());
    
        if logy:
            #canvas_controlplot.SetLogy() ;
            #fPads1.SetLogy();
            fPads1.Update();
            #fPads2.Update();
            canvas_controlplot.Update();

    def GetLumi(self):
        if self.channel=="el":   return 35.9;
        elif self.channel=="mu": return 35.9;

def drawcontrolplot(channel):
    boostedW_fitter=doFit_wj_and_wlvj(channel);
    boostedW_fitter.ControlPlots();



#### Main Code
if __name__ == '__main__':
    channel=options.channel;

    if options.control:
        print 'control for %s sample'%(channel);
        drawcontrolplot(channel);

