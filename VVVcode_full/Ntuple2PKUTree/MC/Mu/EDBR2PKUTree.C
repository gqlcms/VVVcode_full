#define EDBR2PKUTree_cxx
#include "EDBR2PKUTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#define Pi 3.141593
//#include "BTagCalibrationStandalone.h"
//
vector<Double_t> generate_weights(TH1* data_npu_estimated, Int_t isForSynch){
	// see SimGeneral/MixingModule/python/mix_2015_25ns_Startup_PoissonOOTPU_cfi.pyy; copy and paste from there:
	const Double_t npu_probs[75] = {
/*
		4.8551E-07,
		1.74806E-06,
		3.30868E-06,
		1.62972E-05,
		4.95667E-05,
		0.000606966,
		0.003307249,
		0.010340741,
		0.022852296,
		0.041948781,
		0.058609363,
		0.067475755,
		0.072817826,
		0.075931405,
		0.076782504,
		0.076202319,
		0.074502547,
		0.072355135,
		0.069642102,
		0.064920999,
		0.05725576,
		0.047289348,
		0.036528446,
		0.026376131,
		0.017806872,
		0.011249422,
		0.006643385,
		0.003662904,
		0.001899681,
		0.00095614,
		0.00050028,
		0.000297353,
		0.000208717,
		0.000165856,
		0.000139974,
		0.000120481,
		0.000103826,
		8.88868E-05,
		7.53323E-05,
		6.30863E-05,
		5.21356E-05,
		4.24754E-05,
		3.40876E-05,
		2.69282E-05,
		2.09267E-05,
		1.5989E-05,
		4.8551E-06,
		2.42755E-06,
		4.8551E-07,
		2.42755E-07,
		1.21378E-07,
		4.8551E-08
*/
/*
                0.000108643,
                0.000388957,
                0.000332882,
                0.00038397,
                0.000549167,
                0.00105412,
                0.00459007,
                0.0210314,
                0.0573688,
                0.103986,
                0.142369,
                0.157729,
                0.147685,
                0.121027,
                0.08855,
                0.0582866,
                0.0348526,
                0.019457,
                0.0107907,
                0.00654313,
                0.00463195,
                0.00370927,
                0.0031137,
                0.00261141,
                0.00215499,
                0.00174491,
                0.00138268,
                0.00106731,
                0.000798828,
                0.00057785,
                0.00040336,
                0.00027161,
                0.000176535,
                0.00011092,
                6.75502e-05,
                4.00323e-05,
                2.32123e-05,
                1.32585e-05,
                7.51611e-06,
                4.25902e-06,
                2.42513e-06,
                1.39077e-06,
                8.02452e-07,
                4.64159e-07,
                2.67845e-07,
                1.5344e-07,
                8.68966e-08,
                4.84931e-08,
                2.6606e-08,
                1.433e-08	};
*/
/*
		0.000829312873542,
 		0.00124276120498,
 		0.00339329181587,
 		0.00408224735376,
 		0.00383036590008,
		0.00659159288946,
 		0.00816022734493,
 		0.00943640833116,
 		0.0137777376066,
 		0.017059392038,
 		0.0213193035468,
 		0.0247343174676,
 		0.0280848773878,
 		0.0323308476564,
 		0.0370394341409,
 		0.0456917721191,
 		0.0558762890594,
 		0.0576956187107,
 		0.0625325287017,
 		0.0591603758776,
 		0.0656650815128,
 		0.0678329011676,
 		0.0625142146389,
 		0.0548068448797,
 		0.0503893295063,
 		0.040209818868,
 		0.0374446988111,
 		0.0299661572042,
 		0.0272024759921,
 		0.0219328403791,
 		0.0179586571619,
 		0.0142926728247,
 		0.00839941654725,
 		0.00522366397213,
 		0.00224457976761,
 		0.000779274977993,
 		0.000197066585944,
 		7.16031761328e-05,
 		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
 		0.0,
 		0.0,
		0.0 };
*/
                1.78653e-05 ,
                2.56602e-05 ,
                5.27857e-05 ,
                8.88954e-05 ,
                0.000109362 ,
                0.000140973 ,
                0.000240998 ,
                0.00071209 ,
                0.00130121 ,
                0.00245255 , 
                0.00502589 ,
                0.00919534 ,
                0.0146697 ,
                0.0204126 , 
                0.0267586 ,
                0.0337697 ,
                0.0401478 ,
                0.0450159 ,
                0.0490577 ,
                0.0524855 ,
                0.0548159 ,
                0.0559937 ,
                0.0554468 ,
                0.0537687 ,
                0.0512055 ,
                0.0476713 ,
                0.0435312 ,
                0.0393107 ,
                0.0349812 ,
                0.0307413 ,
                0.0272425 ,
                0.0237115 ,
                0.0208329 ,
                0.0182459 ,
                0.0160712 ,
                0.0142498 ,
                0.012804 ,
                0.011571 ,
                0.010547 ,
                0.00959489 ,
                0.00891718 ,
                0.00829292 ,
                0.0076195 ,
                0.0069806 ,
                0.0062025 ,
                0.00546581 ,
                0.00484127 ,
                0.00407168 ,
                0.00337681 ,
                0.00269893 ,
                0.00212473 ,
                0.00160208 , 
                0.00117884 ,
                0.000859662 ,
                0.000569085 ,
                0.000365431 ,
                0.000243565 ,
                0.00015688 ,
                9.88128e-05 ,
                6.53783e-05 ,
                3.73924e-05 ,
                2.61382e-05 ,
                2.0307e-05 ,
                1.73032e-05 ,
                1.435e-05 ,
                1.36486e-05 ,
                1.35555e-05 ,
                1.37491e-05 ,
                1.34255e-05 ,
                1.33987e-05 ,
                1.34061e-05 ,
                1.34211e-05 ,
                1.34177e-05 ,
                1.32959e-05 ,
                1.33287e-05 };
	if (isForSynch==0) { //OFFICIAL RECIPE
		vector<Double_t> result(75);
		Double_t s = 0.0;
		for(Int_t npu=0; npu<75; ++npu){
			Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));              
			result[npu] = npu_estimated / npu_probs[npu];
			s += npu_estimated;
		}
		// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
		for(Int_t npu=0; npu<75; ++npu){
			result[npu] /= s;
		}
		return result;
	}
	else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
		vector<Double_t> result(60);
		for(Int_t npu=0; npu<60; ++npu){
			if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==NULL)
			  result[npu] = 0.;
			else {
				Double_t npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));            
				result[npu] = npu_estimated;
			}
		}
		return result;
	}

}

/*
Double_t bsv (Int_t cud, Double_t x ) // cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.892452;
  }
  else if (cud==2) { //up
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452+0.017849041149020195;}
   else if(x<70)  {result=0.892452+0.017849041149020195;}
   else if(x<100) {result=0.892452+0.017849041149020195;}
   else if(x<140) {result=0.892452+0.020885121077299118;}
   else if(x<200) {result=0.892452+0.025080939754843712;}
   else if(x<300) {result=0.892452+0.10671335458755493;}
   else if(x<670) {result=0.892452+0.16398745775222778;}
   else {result=1;}
  }
  else if (cud==3) {//down
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452-0.017849041149020195;}
   else if(x<70)  {result=0.892452-0.017849041149020195;}
   else if(x<100) {result=0.892452-0.017849041149020195;}
   else if(x<140) {result=0.892452-0.020885121077299118;}
   else if(x<200) {result=0.892452-0.025080939754843712;}
   else if(x<300) {result=0.892452-0.10671335458755493;}
   else if(x<670) {result=0.892452-0.16398745775222778;}
   else {result=1;}
  }
  return result;  
}



Double_t csv (Int_t cud, Double_t x ) // c-jet; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.892452;
  }
  else if (cud==2) { //up 
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452+0.03569808229804039;}
   else if(x<70)  {result=0.892452+0.03569808229804039;}
   else if(x<100) {result=0.892452+0.03569808229804039;}
   else if(x<140) {result=0.892452+0.041770242154598236;}
   else if(x<200) {result=0.892452+0.050161879509687424;}
   else if(x<300) {result=0.892452+0.21342670917510986;}
   else if(x<670) {result=0.892452+0.32797491550445557;}
   else {result=1;}
  }
  else if (cud==3) {//down
   if(x<30) {result=1;}
   else if(x<50)  {result=0.892452-0.03569808229804039 ;}
   else if(x<70)  {result=0.892452-0.03569808229804039 ;}
   else if(x<100) {result=0.892452-0.03569808229804039 ;}
   else if(x<140) {result=0.892452-0.041770242154598236 ;}
   else if(x<200) {result=0.892452-0.050161879509687424 ;}
   else if(x<300) {result=0.892452-0.21342670917510986 ;}
   else if(x<670) {result=0.892452-0.32797491550445557 ;}
   else {result=1;}
  }
  return result;
}

Double_t lsv (Int_t cud, Double_t x ) // light flavor; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.924144-0.000861952*x+3.46078e-06*x*x-2.4028e-09*x*x*x;
  }
  else if (cud==2) { //up 
   result=1.02632+-0.000806342*x+3.54292e-06*x*x+-2.50001e-09*x*x*x;
  }
  else if (cud==3) {//down
   result=0.821939-0.00091531*x+3.37305e-06*x*x-2.30372e-09*x*x*x;
  }
  return result;
}
*/

Double_t bsv (Int_t cud, Double_t x ) // cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.498094*((1.+(0.422991*x))/(1.+(0.210944*x)));
  }
  else if (cud==2) { //up
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.039442747831344604;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.012669667601585388;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.011243613436818123;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.010686126537621021;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.010994619689881802;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.011888998560607433;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.015397069044411182;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.020918292924761772;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.031905386596918106;}
   else {result=1;}
  }
 else if (cud==3) {//down
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.039442747831344604;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.012669667601585388;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.011243613436818123;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.010686126537621021;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.010994619689881802;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.011888998560607433;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.015397069044411182;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.020918292924761772;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.031905386596918106;}
   else {result=1;}
  }
  return result;
}


Double_t csv (Int_t cud, Double_t x ) // c-jet; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=0.498094*((1.+(0.422991*x))/(1.+(0.210944*x)));
  }
  else if (cud==2) { //up 
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.078885495662689209;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.025339335203170776;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.022487226873636246;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.021372253075242043;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.021989239379763603;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.023777997121214867;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.030794138088822365;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.041836585849523544;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))+0.063810773193836212;}
   else {result=1;}
  } 
  else if (cud==3) {//down
   if(x<20) {result=1;}
   else if(x<30)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.078885495662689209;}
   else if(x<50)  {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.025339335203170776;}
   else if(x<70) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.022487226873636246;}
   else if(x<100) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.021372253075242043;}
   else if(x<140) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.021989239379763603;}
   else if(x<200) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.023777997121214867;}
   else if(x<300) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.030794138088822365;}
   else if(x<600) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.041836585849523544;}
   else if(x<1000) {result=(0.498094*((1.+(0.422991*x))/(1.+(0.210944*x))))-0.063810773193836212;}
   else {result=1;}
  }
  return result;
}
  


Double_t lsv (Int_t cud, Double_t x ) // light flavor; cud=1,2,3 for central,up,down; x for pt
{
  double result=1.0;

  if (cud==1) {  //central
   result=1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x;
  }
  else if (cud==2) { //up 
   result=(1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1+(0.100485+3.95509e-05*x+-4.90326e-08*x*x));
  }
  else if (cud==3) {//down
   result=(1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x)*(1-(0.100485+3.95509e-05*x+-4.90326e-08*x*x));
  }
  return result;
}

void EDBR2PKUTree::Loop(TString channelname, Double_t XS, Double_t totaleventnumber, Int_t IsData) {

	std::vector<Double_t> weights_pu1; //these are made with our recipe
	std::vector<Double_t> weights_pu2; //these are made with the official recipe
	TFile* pileupFile1 = TFile::Open("MyDataPileupHistogram.root");//pileupDataRun2016BH_63mb_80X.root");  
	TH1F* pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
	weights_pu1 = generate_weights(pileupHisto1,0);
	pileupFile1->Close();

	//  TFile* pileupFile2 = TFile::Open("puweights.root");  
	TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
	TH1F *pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
	weights_pu2 = generate_weights(pileupHisto2,1);
	pileupFile2->Close();

	//TFile * input1 = new TFile ("puweights.root");
	//TH1F* hR1= (TH1F*)input1->Get("puweights");
	//zixu
	TFile * input1 = new TFile ("puweight.root");	
	TH1F* hR1= (TH1F*)input1->Get("h2");
	//TFile * input1 = new TFile ("test_mu.root");
	//TH1F* hR1= (TH1F*)input1->Get("hRatio"); //"pileup");//hRatio");


	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();

	Double_t n_deltaRlepjet = 0; 
	Double_t n_delPhijetlep = 0; 
	Double_t ntau = 0;
	Double_t number_qq = 0; 
	Double_t nmassVhad = 0; 
	Double_t nptVlepJEC = 0;
	Double_t nID_e = 0;
	Double_t npt_e = 0;
	Double_t nmet_e = 0; 
	Double_t nnum_bJet_e = 0; 
	Double_t n_delPhijetmet = 0; 

	Double_t nID_mu = 0;
	Double_t npt_mu = 0;
	Double_t nmet_mu = 0; 
	Double_t nnum_bJet_mu = 0; 
	//Double_t nbtb_mu = 0; 

	Double_t nptVhad = 0;
	Double_t yields = 0;
	//TLorentzVector jetV, genjetV;
	//some constants inside this analysis
	Double_t pi_2=1.57079632679;
	Long64_t npp = fChain->GetEntries("theWeight>0.");
	Long64_t nmm = fChain->GetEntries("theWeight<0.");
	cout<<"npp="<<npp<<" nmm="<<nmm<<" totaleventnumber="<<totaleventnumber<<endl;

	Double_t nn;
	Double_t eff_and_pu_Weight;
	Double_t eff_and_pu_Weight1;
	Float_t Identical_lumiWeight = XS;//All the events inside a sample are same lumiweight
	//Float_t Identical_lumiWeight = XS/totaleventnumber;//All the events inside a sample are same lumiweight

	Long64_t nbytes = 0, nb = 0;
	//for (Long64_t jentry=0; jentry<10;jentry++)
	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
                   if (jentry%40000==0)
                      {std::cout<<jentry<<std::endl;}
		nb = fChain->GetEntry(jentry); 
		nbytes += nb;
		pfMET             = Float_t(met);
		pfMETPhi          = Float_t(metPhi);
		l_pt              = Float_t(ptlep1);
		l_eta             = Float_t(etalep1);
		l_phi             = Float_t(philep1);
		ptVhad            = Float_t(ptVhad);
		jet_eta           = Float_t(yVhad);
		jet_phi           = Float_t(phiVhad);
		jet_mass_pruned   = Float_t(massVhadJEC);
                jet_mass_puppi    = Float_t(jetAK8puppi_sdJEC);
                jet_mass_puppi_un = Float_t(jetAK8puppi_sd);
                jet_tau2tau1_puppi    = Float_t(jetAK8puppi_tau21);
                jet_pt_puppi      = Float_t(jetAK8puppi_ptJEC);
		jetAK8_mass       = Float_t(jetAK8_mass);
//		jet_mass_softdrop = Float_t(sdropJEC);
		jet_tau2tau1      = Float_t(tau21);
		W_pt              = Float_t(ptVlepJEC);
		W_eta             = Float_t(yVlep);
		W_phi             = Float_t(phiVlep);
                m_lvj             = Float_t(candMasspuppiJEC);
		
		jet_mass_puppi_2    = Float_t(jetAK8puppi_sdJEC_2);
                jet_mass_puppi_un_2 = Float_t(jetAK8puppi_sd_2);
                jet_tau2tau1_puppi_2    = Float_t(jetAK8puppi_tau21_2);
                jet_pt_puppi_2      = Float_t(jetAK8puppi_ptJEC_2);
		jet_tau4tau2_puppi    = Float_t(jetAK8puppi_tau42);
		jet_tau4tau2_puppi_2    = Float_t(jetAK8puppi_tau42_2);
for(Int_t iii=0;iii<3;iii++) 
MassVV[iii] =massww[iii];
//cout<<jet_mass_puppi_un_2<<endl;
//cout<<jetAK8puppi_sd_2<<endl;
//		
//		fjet2_pt          = Float_t(jet2_pt);
//		fjet2_btag        = Float_t(jet2_btag);
//		fjet3_pt          = Float_t(jet3_pt);
//		fjet3_btag        = Float_t(jet3_btag);
                mtVlepnew         = Float_t(sqrt(2*ptlep1*met*(1.0-cos(philep1-metPhi))));
                MTVlep            = Float_t(sqrt(2*ptlep1*MET_et*(1.0-cos(philep1-MET_phi))));
                //puppi+softdrop recorrected by Thea's "JEC"
                Double_t gencorrect=1.0;
                Double_t recocorrect_0eta1p3=1.0;
                Double_t recocorrect_1p3eta2p5=1.0;

//                gencorrect=1.0-0.321*pow(jet_pt_puppi*0.0354,-1.1);
//                recocorrect_0eta1p3=1.09-1.69e-04*jet_pt_puppi+3.34e-07*pow(jet_pt_puppi,2)-2.47e-10*pow(jet_pt_puppi,3)+7.8e-14*pow(jet_pt_puppi,4)-8.83e-18*pow(jet_pt_puppi,5);
//                recocorrect_1p3eta2p5=1.3-7.76e-04*jet_pt_puppi+1.11e-06*pow(jet_pt_puppi,2)-6.79e-10*pow(jet_pt_puppi,3)+1.87e-13*pow(jet_pt_puppi,4)-1.9e-17*pow(jet_pt_puppi,5);
                gencorrect=1.006-1.062*pow(jet_pt_puppi*0.08,-1.2);
                recocorrect_0eta1p3=1.093-1.501e-04*jet_pt_puppi+3.449e-07*pow(jet_pt_puppi,2)-2.681e-10*pow(jet_pt_puppi,3)+8.674e-14*pow(jet_pt_puppi,4)-1.001e-17*pow(jet_pt_puppi,5);
                recocorrect_1p3eta2p5=1.272-5.72e-04*jet_pt_puppi+8.37e-07*pow(jet_pt_puppi,2)-5.204e-10*pow(jet_pt_puppi,3)+1.454e-13*pow(jet_pt_puppi,4)-1.504e-17*pow(jet_pt_puppi,5);
                if (fabs(jetAK8puppi_eta)<=1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_0eta1p3;}
                else if (fabs(jetAK8puppi_eta)<2.5 && fabs(jetAK8puppi_eta)>1.3){jet_mass_puppi_corr=jet_mass_puppi_un*gencorrect*recocorrect_1p3eta2p5;}


Double_t gencorrect_2=1.0;
                Double_t recocorrect_0eta1p3_2=1.0;
                Double_t recocorrect_1p3eta2p5_2=1.0;
gencorrect_2=1.006-1.062*pow(jet_pt_puppi_2*0.08,-1.2);
                recocorrect_0eta1p3_2=1.093-1.501e-04*jet_pt_puppi_2+3.449e-07*pow(jet_pt_puppi_2,2)-2.681e-10*pow(jet_pt_puppi_2,3)+8.674e-14*pow(jet_pt_puppi_2,4)-1.001e-17*pow(jet_pt_puppi_2,5);
                recocorrect_1p3eta2p5_2=1.272-5.72e-04*jet_pt_puppi_2+8.37e-07*pow(jet_pt_puppi_2,2)-5.204e-10*pow(jet_pt_puppi_2,3)+1.454e-13*pow(jet_pt_puppi_2,4)-1.504e-17*pow(jet_pt_puppi_2,5);
                if (fabs(jetAK8puppi_eta_2)<=1.3){jet_mass_puppi_corr_2=jet_mass_puppi_un_2*gencorrect_2*recocorrect_0eta1p3_2;}
                else if (fabs(jetAK8puppi_eta_2)<2.5 && fabs(jetAK8puppi_eta_2)>1.3){jet_mass_puppi_corr_2=jet_mass_puppi_un_2*gencorrect_2*recocorrect_1p3eta2p5_2;}

		// GEN-RECO match
		//deltaRleplep = deltaR(etalep1,philep1,etalep2,philep2);
		//deltaRWlepGen = deltaR(etaGenVlep, phiGenVlep, yVlep, phiVlep);
		//jetV.SetPtEtaPhiM(ptVhad, yVhad, phiVhad, massVhad);
		//genjetV.SetPtEtaPhiM(ptGenVhad, etaGenVhad, phiGenVhad, massGenVhad);
		//deltaRWhadGen = deltaR(etaGenVhad, phiGenVhad, yVhad, phiVhad);
		//Double_t deltaRWhadGen = sqrt(pow(etaGenVhad-yVhad,2) + pow(phiGenVhad-phiVhad,2));
		Double_t deltaRWhadGen = sqrt(pow(etaGenVhad-jetAK8puppi_eta,2) + pow(phiGenVhad-jetAK8puppi_phi,2));

		//Weight Calculation
		Int_t bin = hR1->FindBin(npT);
		pileupWeight = hR1->GetBinContent(bin);		

		eff_and_pu_Weight = 0;
		eff_and_pu_Weight1 = 0;
		if(IsData>0) {
			if(npT < weights_pu1.size()){
				eff_and_pu_Weight = weights_pu1[npT];
			}
			if(npT < weights_pu2.size()){
				eff_and_pu_Weight1 = weights_pu2[npT];
			}
		}
/*
      trigger_eff=1.0;

      TFile * input_trigger = new TFile ("MuonEff-TH2D.root");
      input_trigger->cd("");
      TH2D* HLTeff= (TH2D*) input_trigger->Get("2Dh");

      Double_t mu_pt=ptlep1;
      if (mu_pt>500) {mu_pt=499.0;}
      int mubin=HLTeff->FindBin(fabs(etalep1),mu_pt);
      trigger_eff=HLTeff->GetBinContent(mubin); 
      input_trigger->Close();

*/
//      std::cout<<trigger_eff<<std::endl;

            trigger_eff=1.0;
            if (MET_et>50. && MET_et<80.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.9985;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.998136;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.998161;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.995559;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.998492;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.998266;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.998287;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.995788;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.998405;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.998167;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.998242;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.996395;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.998144;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997348;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.998569;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.997289;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997899;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997447;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.997746;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.988402;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.998579;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.995749;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999496;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.996213;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.998111;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.996574;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996428;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.993633;}
                   }
                }
            else if (MET_et>80. && MET_et<100.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997644;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997270;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996798;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.992983;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997635;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997459;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996983;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.993342;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997509;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997320;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996941;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.994178;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997151;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.99619;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.997361;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.995318;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.996682;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.996094;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996094;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.983748;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997456;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.993761;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.998277;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.994784;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.996758;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.995032;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.994452;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.99123;}
                   }
                }
            else if (MET_et>100. && MET_et<150.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.996256;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.995956;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.99444;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.988785;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.996245;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.996235;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.994714;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.989353;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.99606;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.996038;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.994691;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.990507;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.995571;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.994475;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.99522;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.991925;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.994699;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.993945;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.993782;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.977099;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.995489;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.990738;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.995866;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.992801;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.998111;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.996574;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996428;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.993633;}
                   }
                }
            else if (MET_et>150. && MET_et<200.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997839;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997755;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.99665;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.993519;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997834;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997908;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996802;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.993845;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997732;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.997803;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996802;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.994456;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997473;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.996972;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.997064;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.995154;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.996933;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.996553;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.996326;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.987679;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.997255;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.994855;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.997192;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.996176;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.996664;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.995986;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.994996;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.99357;}
                   }
                }
            else if (MET_et>200. && MET_et<250.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999483;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999459;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999207;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.998453;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999482;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999496;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999244;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.998531;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999458;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.99947;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999243;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.99868;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999394;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999268;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999308;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.998853;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999267;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999174;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999127;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.997008;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999351;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.99876;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999351;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999068;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999209;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999031;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.998803;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.998433;}
                   }
                }
            else if (MET_et>250. && MET_et<300.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999802;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999773;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999728;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999411;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999801;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999801;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999743;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999441;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999791;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999777;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.99974;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.99951;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999761;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999683;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999774;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999603;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999721;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999673;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999685;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.998655;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999783;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.99948;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999846;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999569;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999725;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999587;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999535;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999275;}
                   }
                }
            else if (MET_et>300. && MET_et<500.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.99984;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999824;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999768;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999523;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.99984;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999836;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.99978;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999547;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999831;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999828;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999778;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999598;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.99981;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999758;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999802;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999663;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999774;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.99974;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999738;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.998991;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999812;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999598;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999839;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999681;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999767;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999682;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999627;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999627;}
                   }
                }
            else if (MET_et>500.){
                if (ptlep1>55. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999831;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999826;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999737;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999495;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999831;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999838;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999749;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.99952;}
                   }
                else if (ptlep1>=80. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999823;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999829;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999749;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999567;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999803;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999765;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999769;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999621;}
                   }
                else if (ptlep1>=200. && ptlep1<300.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999761;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999732;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999713;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999048;}
                   }
                else if (ptlep1>=300. && ptlep1<400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999785;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999601;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.999777;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999705;}
                   }
                else if (ptlep1>=400.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.999739;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.999689;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.99961;}
                    else if (fabs(etalep1)<2.4){trigger_eff=0.999504;}
                   }
                }


/*
                trigger_eff=1.0;
                if (ptlep1>50. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.913068;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.888577;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.842243;}
                   }
                else if (ptlep1>=60. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.912803;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.886555;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.83609;}
                   }
                else if (ptlep1>=120. && ptlep1<200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.909302;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.862569;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.757962;}
                   }
                else if (ptlep1>=200.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.887649;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.81576;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.626224;}
                   }
*/
/*
                if (ptlep1>50. && ptlep1<60.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.975;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.922;}
                    else if (fabs(etalep1)<1.44){trigger_eff=0.918;}
                    else if (fabs(etalep1)<1.56){trigger_eff=0.867;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.873;}
                   }
                else if (ptlep1>=60. && ptlep1<80.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.974;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.936;}
                    else if (fabs(etalep1)<1.44){trigger_eff=0.906;}
                    else if (fabs(etalep1)<1.56){trigger_eff=0.844;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.882;}
                   }
                else if (ptlep1>=80. && ptlep1<100.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.965;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.899;}
                    else if (fabs(etalep1)<1.44){trigger_eff=0.914;}
                    else if (fabs(etalep1)<1.56){trigger_eff=0.717;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.862;}
                   }
                else if (ptlep1>=100. && ptlep1<120.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.968;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.915;}
                    else if (fabs(etalep1)<1.44){trigger_eff=0.842;}
                    else if (fabs(etalep1)<1.56){trigger_eff=0.635;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.938;}
                   }
                else if (ptlep1>=120. && ptlep1<140.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.97;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.913;}
                    else if (fabs(etalep1)<1.44){trigger_eff=0.905;}
                    else if (fabs(etalep1)<1.56){trigger_eff=0.818;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.848;}
                   }
                else if (ptlep1>=140. && ptlep1<240.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.964;}
                    else if (fabs(etalep1)<1.2){trigger_eff=0.968;}
                    else if (fabs(etalep1)<1.44){trigger_eff=0.93;}
                    else if (fabs(etalep1)<1.56){trigger_eff=0.564;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.922;}
                   }
                else if (ptlep1>=240.)
                   {if (fabs(etalep1)<0.9){trigger_eff=0.974;}
                    else if (fabs(etalep1)<1.2){trigger_eff=1.0;}
                    else if (fabs(etalep1)<1.44){trigger_eff=1.0;}
                    else if (fabs(etalep1)<1.56){trigger_eff=1.0;}
                    else if (fabs(etalep1)<2.1){trigger_eff=0.8;}
                   }
*/

		//cout << "pileupWeight:"<<pileupWeight<< " eff_and_pu_Weight:" << eff_and_pu_Weight << " eff_and_pu_Weight1:" << eff_and_pu_Weight1 << endl;
		if(theWeight>0) nn=1;
		else nn= -1;
		if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
		else lumiWeight=Identical_lumiWeight/nentries;
		weight_nobtag=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff;
		//weight=lumiWeight*triggerWeight*pileupWeight*nn;
		if (IsData>1 ) weight_nobtag = weight_nobtag*1.21;

		//lumiWeight=Identical_lumiWeight;
		//if(npp>0) weight=lumiWeight*triggerWeight*pileupWeight/(npp-nmm)*nn*0.04024;//0.00559;
		//else weight=lumiWeight*triggerWeight*pileupWeight/nentries*0.04024;//0.00559;
		if ( IsData==0 ) weight_nobtag=1;
		//Weight Calculation Done


      IDweight=1.0;
                if (ptlep1>=55 && ptlep1<60)
                   {if (fabs(etalep1)<0.9){IDweight=0.982688;}
                    else if (fabs(etalep1)<1.2){IDweight=0.976116;}
                    else if (fabs(etalep1)<2.1){IDweight=0.986235;}
                    else if (fabs(etalep1)<2.4){IDweight=0.966154;}
                   }
                else if (ptlep1>=60)
                   {if (fabs(etalep1)<0.9){IDweight=0.99385;}
                    else if (fabs(etalep1)<1.2){IDweight=0.976965;}
                    else if (fabs(etalep1)<2.1){IDweight=0.990001;}
                    else if (fabs(etalep1)<2.4){IDweight=0.967312;}
                   }

      IDweightISO=1.0;
                if (ptlep1>=55 && ptlep1<60)
                   {if (fabs(etalep1)<0.9){IDweightISO=0.997968;}
                    else if (fabs(etalep1)<1.2){IDweightISO=0.998721;}
                    else if (fabs(etalep1)<2.1){IDweightISO=1.000081;}
                    else if (fabs(etalep1)<2.4){IDweightISO=1.00056;}
                   }
                else if (ptlep1>=60)
                   {if (fabs(etalep1)<0.9){IDweightISO=0.998726;}
                    else if (fabs(etalep1)<1.2){IDweightISO=0.999314;}
                    else if (fabs(etalep1)<2.1){IDweightISO=0.999228;}
                    else if (fabs(etalep1)<2.4){IDweightISO=1.002153;}
                   }

      IDweighttrk=1.0;
/*
        if (etalep1<-2.1){IDweighttrk=9.82399009186853522e-01;}
        else if(etalep1<-1.6){IDweighttrk=9.91746789037933008e-01;}
        else if(etalep1<-1.1){IDweighttrk=9.95944961092376846e-01;}
        else if(etalep1<-0.6){IDweighttrk=9.93413142541369476e-01;}
        else if(etalep1<0.0){IDweighttrk=9.91460688530866996e-01;}
        else if(etalep1<0.6){IDweighttrk=9.94680143661991423e-01;}
        else if(etalep1<1.1){IDweighttrk=9.96666389348924819e-01;}
        else if(etalep1<1.6){IDweighttrk=9.94933892427240618e-01;}
        else if(etalep1<2.1){IDweighttrk=9.91186607207322878e-01;}
        else if(etalep1<2.4){IDweighttrk=9.76811919457875155e-01;}
*/

        ToppTweight=1.0;
        double a_top=0.0615;
        double b_top=- 0.0005;
        if (gentop_pt<0 || genantitop_pt<0){ToppTweight=1.0;}
        else if(gentop_pt>0 && genantitop_pt>0){ToppTweight=sqrt(exp(a_top+b_top*gentop_pt)*exp(a_top+b_top*genantitop_pt));}



//--BSF----------------------------
      btagweight_center=1.0;
      btagweight_up=1.0;
      btagweight_down=1.0;
      

      double  bweight=1.0, dweight=1.0;
      double  bweightup=1.0;
      double  bweightdown=1.0;

      double  beff=1.0, ceff=1.0, leff=1.0;



      for(Int_t i=0; i<8; i++)  {
       deltaRAK4AK8_new[i]=0.;
       deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));

       if(ak4jet_pt[i]>30 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 ) {
         if(abs(ak4jet_hf[i])==5 && ak4jet_pt[i]<3000.) 
          { 
            double jet_sf    = bsv(1,ak4jet_pt[i]); 
            double jet_sfu    = bsv(2,ak4jet_pt[i]);
            double jet_sfd    = bsv(3,ak4jet_pt[i]);
//            int bbin=hb->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            beff=hb->GetBinContent(bbin);
           if (ak4jet_pt[i]>30. && ak4jet_pt[i]<50.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.647202;}
               else if (ak4jet_eta[i]<-0.8){beff=0.57439;}
               else if (ak4jet_eta[i]>0.8){beff=0.580744;}
              }
           else if (ak4jet_pt[i]>=50. && ak4jet_pt[i]<70.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.681678;}
               else if (ak4jet_eta[i]<-0.8){beff=0.608012;}
              else if (ak4jet_eta[i]>0.8){beff=0.615113;}
              }
           else if (ak4jet_pt[i]>=70. && ak4jet_pt[i]<100.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.695347;}
               else if (ak4jet_eta[i]<-0.8){beff=0.62226;}
               else if (ak4jet_eta[i]>0.8){beff=0.629468;}
              }
           else if (ak4jet_pt[i]>=100. && ak4jet_pt[i]<140.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.696113;}
               else if (ak4jet_eta[i]<-0.8){beff=0.619795;}
               else if (ak4jet_eta[i]>0.8){beff=0.627692;}
              }
           else if (ak4jet_pt[i]>=140. && ak4jet_pt[i]<200.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.680334;}
               else if (ak4jet_eta[i]<-0.8){beff=0.611209;}
               else if (ak4jet_eta[i]>0.8){beff=0.616186;}
              }
           else if (ak4jet_pt[i]>=200. && ak4jet_pt[i]<300.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.639435;}
               else if (ak4jet_eta[i]<-0.8){beff=0.563817;}
               else if (ak4jet_eta[i]>0.8){beff=0.568733;}
              }
           else if (ak4jet_pt[i]>=300. && ak4jet_pt[i]<600.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.560819;}
               else if (ak4jet_eta[i]<-0.8){beff=0.522108;}
               else if (ak4jet_eta[i]>0.8){beff=0.52554;}
              }
           else if (ak4jet_pt[i]>=600. && ak4jet_pt[i]<1000.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.470791;}
               else if (ak4jet_eta[i]<-0.8){beff=0.469051;}
               else if (ak4jet_eta[i]>0.8){beff=0.477012;}
              }
           else if (ak4jet_pt[i]>=1000. && ak4jet_pt[i]<3000.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){beff=0.415562;}
               else if (ak4jet_eta[i]<-0.8){beff=0.464442;}
               else if (ak4jet_eta[i]>0.8){beff=0.433394;}
              }

             if(ak4jet_icsv[i]>0.8484) {
                   bweight=bweight*beff*jet_sf; dweight=dweight*beff;
                   bweightup=bweightup*beff*jet_sfu; 
                   bweightdown=bweightdown*beff*jet_sfd; 
                   }
             else {
                   bweight=bweight*(1-beff*jet_sf); dweight=dweight*(1-beff);
                   bweightup=bweightup*(1-beff*jet_sfu); 
                   bweightdown=bweightdown*(1-beff*jet_sfd); 
                  }
          }
        else if(abs(ak4jet_hf[i])==4 && ak4jet_pt[i]<3000.)
          {   
            double jet_sf    = csv(1,ak4jet_pt[i]);
            double jet_sfu    = csv(2,ak4jet_pt[i]);
            double jet_sfd    = csv(3,ak4jet_pt[i]);
//            int cbin=hc->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            ceff=hc->GetBinContent(cbin);
           if (ak4jet_pt[i]>30. && ak4jet_pt[i]<50.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.13475;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.121195;}
               else if (ak4jet_eta[i]>0.8){ceff=0.121158;}
              }
           else if (ak4jet_pt[i]>=50. && ak4jet_pt[i]<70.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.131149;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.117452;}
              else if (ak4jet_eta[i]>0.8){ceff=0.118667;}
              }
           else if (ak4jet_pt[i]>=70. && ak4jet_pt[i]<100.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.136182;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.124056;}
               else if (ak4jet_eta[i]>0.8){ceff=0.125816;}
              }
           else if (ak4jet_pt[i]>=100. && ak4jet_pt[i]<140.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.143095;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.130459;}
               else if (ak4jet_eta[i]>0.8){ceff=0.132076;}
              }
           else if (ak4jet_pt[i]>=140. && ak4jet_pt[i]<200.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.141677;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.13829;}
               else if (ak4jet_eta[i]>0.8){ceff=0.13913;}
              }
           else if (ak4jet_pt[i]>=200. && ak4jet_pt[i]<300.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.134463;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.129813;}
               else if (ak4jet_eta[i]>0.8){ceff=0.130635;}
              }
           else if (ak4jet_pt[i]>=300. && ak4jet_pt[i]<600.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.121893;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.139806;}
               else if (ak4jet_eta[i]>0.8){ceff=0.138651;}
              }
           else if (ak4jet_pt[i]>=600. && ak4jet_pt[i]<1000.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.121765;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.149446;}
               else if (ak4jet_eta[i]>0.8){ceff=0.1143389;}
              }
           else if (ak4jet_pt[i]>=1000. && ak4jet_pt[i]<3000.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){ceff=0.13307;}
               else if (ak4jet_eta[i]<-0.8){ceff=0.178248;}
               else if (ak4jet_eta[i]>0.8){ceff=0.143266;}
              }

             if(ak4jet_icsv[i]>0.8484) {
                   bweight=bweight*ceff*jet_sf; dweight=dweight*ceff;
                   bweightup=bweightup*ceff*jet_sfu; 
                   bweightdown=bweightdown*ceff*jet_sfd; 
                   }
             else {
                   bweight=bweight*(1-ceff*jet_sf); dweight=dweight*(1-ceff);
                   bweightup=bweightup*(1-ceff*jet_sfu); 
                   bweightdown=bweightdown*(1-ceff*jet_sfd); 
                  }
          }
        else if (abs(ak4jet_hf[i])==0 && ak4jet_pt[i]<3000.)
          {  
            double jet_sf    = lsv(1,ak4jet_pt[i]);
            double jet_sfu    = lsv(2,ak4jet_pt[i]);
            double jet_sfd    = lsv(3,ak4jet_pt[i]);
//            int lbin=hl->FindBin(ak4jet_pt[i],ak4jet_eta[i]);
//            leff=hl->GetBinContent(lbin);
           if (ak4jet_pt[i]>=20. && ak4jet_pt[i]<200.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0088292;}
               else if (ak4jet_eta[i]<-0.8){leff=0.0126681;}
               else if (ak4jet_eta[i]>0.8){leff=0.0117989;}
              }
           else if (ak4jet_pt[i]>=200. && ak4jet_pt[i]<500.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0107827;}
               else if (ak4jet_eta[i]<-0.8){leff=0.0208594;}
               else if (ak4jet_eta[i]>0.8){leff=0.0200988;}
              }
           else if (ak4jet_pt[i]>=500. && ak4jet_pt[i]<1000.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.019908;}
               else if (ak4jet_eta[i]<-0.8){leff=0.0357475;}
               else if (ak4jet_eta[i]>0.8){leff=0.0326975;}
              }
           else if (ak4jet_pt[i]>=1000. && ak4jet_pt[i]<3000.)
              {if (ak4jet_eta[i]<0.8 && ak4jet_eta[i]>-0.8){leff=0.0353687;}
               else if (ak4jet_eta[i]<-0.8){leff=0.0530973;}
               else if (ak4jet_eta[i]>0.8){leff=0.0435917;}
              }

             if(ak4jet_icsv[i]>0.8484) {
                   bweight=bweight*leff*jet_sf; dweight=dweight*leff;
                   bweightup=bweightup*leff*jet_sfu; 
                   bweightdown=bweightdown*leff*jet_sfd; 
                   }

             else {
                   bweight=bweight*(1-leff*jet_sf); dweight=dweight*(1-leff);
                   bweightup=bweightup*(1-leff*jet_sfu); 
                   bweightdown=bweightdown*(1-leff*jet_sfd); 
                  }
          }
      //std::cout<<" yy1 "<<abs(ak4jet_hf[i])<<" "<<dweight<<" "<<bweight<<std::endl;
        } 
      //std::cout<<" yy2 "<<abs(ak4jet_hf[i])<<" "<<ak4jet_pt[i]<<" "<<fabs(ak4jet_eta[i])<<" "<<ak4jet_IDLoose[i]<<" "<<deltaRAK4AK8[i]<<std::endl;
     }
//     input->Close();
//     std::cout<<" xx "<<bweight/dweight<<" "<<bweightup/dweight<<" "<<bweightdown/dweight<<std::endl;

     btagweight_center=bweight/dweight;
     btagweight_up=bweightup/dweight;
     btagweight_down=bweightdown/dweight;

//     std::cout<<btagweight_center<<std::endl;

     
//--BSF----------------------------

                if(theWeight>0) nn=1;
                else nn= -1;
                if(npp>0) lumiWeight=Identical_lumiWeight/(npp-nmm);
                else lumiWeight=Identical_lumiWeight/nentries;
//                std::cout<<IDweight<<"        "<<btagweight_center<<"   "<<trigger_eff<<std::endl;
                weight=lumiWeight*triggerWeight*eff_and_pu_Weight*nn*trigger_eff*btagweight_center*IDweight*IDweightISO*IDweighttrk*ToppTweight;
                if (IsData==2 ) weight = weight*1.21*1.06684;
                if (IsData==3 ) weight = weight*1.21*1.046;
                if (IsData==4 ) weight = weight*1.21*0.99246;
                if (IsData==5 ) weight = weight*1.21*0.9155;
                if (IsData==6 ) weight = weight*1.21*0.8093;
                if (IsData==7 ) weight = weight*1.21*0.6498;
                if (IsData==8 ) weight = weight*1.21*0.4843;
                if ( IsData==0 ) weight=1;

		//number of bjet calculation
		num_bJet=0.;
		num_bJet_loose=0.;
		num_bJet_tight=0.;
		for(Int_t i=0; i<8; i++)  {
                        deltaRAK4AK8_new[i]=0.;
                        deltaRAK4AK8_new[i]=sqrt(pow(fabs(ak4jet_eta[i]-jetAK8puppi_eta),2)+pow(TMath::Min(fabs(ak4jet_phi[i]-jetAK8puppi_phi),2*Pi-fabs(ak4jet_phi[i]-jetAK8puppi_phi)),2));
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.8484 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 ) {num_bJet=num_bJet+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.5426 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 ) {num_bJet_loose=num_bJet_loose+1;}
			if(ak4jet_pt[i]>30 && ak4jet_icsv[i]>0.9535 && fabs(ak4jet_eta[i])<2.4 && ak4jet_IDLoose[i]>0 && deltaRAK4AK8_new[i]>=0.8 ) {num_bJet_tight=num_bJet_tight+1;}
		}
		nbtag=num_bJet;
		//number of bjet calculation Done

		Int_t nLooseLep=nLooseEle+nLooseMu;//the tight Lep included

		Double_t isAnaHP=1.;
		Double_t isAnaLP=1.;
		Double_t isAnaNP=1.;
		Double_t isTTBarControl=1.;
		Int_t tmp_categoryID_channel=0;
		if( channelname=="el" ){
		tmp_categoryID_channel=-1;// -1 for el; 1 for mu

			//HP: 0<jetAK8puppi_tau21<=0.5;
			if (isAnaHP>0 && lep==11 && nLooseLep==1){ nID_e = nID_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptlep1>55 && fabs(etalep1)<2.5){ npt_e = npt_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && MET_et>80) { nmet_e = nmet_e+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVlepJEC > 200.) { nptVlepJEC = nptVlepJEC +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 ){ nptVhad = nptVhad+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && num_bJet<1){ nnum_bJet_e = nnum_bJet_e +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && deltaRlepjet>pi_2) {n_deltaRlepjet = n_deltaRlepjet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetmet) >2.0)  {n_delPhijetmet = n_delPhijetmet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetlep)>2.0) {n_delPhijetlep = n_delPhijetlep +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && jetAK8puppi_tau21>0. && jetAK8puppi_tau21<=0.45) {ntau = ntau+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0)// && jetAK8puppi_sdcorr>0 && jetAK8puppi_sdcorr <150)// && m_lvj>200 && m_lvj<5000)
			{
				nmassVhad = nmassVhad +1;
				yields = yields + weight;
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<jetAK8puppi_tau21<=0.75;
			if ( lep==11 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200.  && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && jetAK8puppi_tau21>0.45 && jetAK8puppi_tau21<=0.75)// && (( jetAK8puppi_sdcorr >0 &&  jetAK8puppi_sdcorr< 150 )) )
			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }
			//NP: 0.75<jetAK8puppi_tau21
			if ( lep==11 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200.  && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && jetAK8puppi_tau21>0.75)// && (( jetAK8puppi_sdcorr >0 &&  jetAK8puppi_sdcorr< 150 )) )
			{ isAnaNP=1.; } 
			else{ isAnaNP=-1.; }


			//TTbar control
			if ( lep==11 && nLooseLep==1 && ptlep1>55 && fabs(etalep1)<2.5 && MET_et>80 && ptVlepJEC > 200. && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0  && num_bJet>0)// &&  jetAK8puppi_sdcorr>0 && jetAK8puppi_sdcorr <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }
		}
		else if( channelname=="mu" ){
		tmp_categoryID_channel=1;// -1 for el; 1 for mu
			//HP: 0<jetAK8puppi_tau21<=0.5;
			if (isAnaHP>0 && lep==13 && (HLT_Mu2>0 || HLT_Mu3>0 ) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && nLooseLep==1 ) { nID_mu = nID_mu+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptlep1>55){ npt_mu = npt_mu+1; }else{ isAnaHP=-1; }
			if (isAnaHP>0 && MET_et>40) { nmet_mu = nmet_mu+1; }else{ isAnaHP=-1; }
			if (isAnaHP>0 && ptVlepJEC>200) { nptVlepJEC = nptVlepJEC +1;} else{ isAnaHP=-1; }
			if (isAnaHP>0 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0){ nptVhad = nptVhad+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && num_bJet<1){ nnum_bJet_mu = nnum_bJet_mu +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && deltaRlepjet>pi_2) {n_deltaRlepjet = n_deltaRlepjet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetmet) >2.0)  {n_delPhijetmet = n_delPhijetmet+1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && fabs(delPhijetlep)>2.0) {n_delPhijetlep = n_delPhijetlep +1; } else{ isAnaHP=-1; }
			if (isAnaHP>0 && jetAK8puppi_tau21>0. && jetAK8puppi_tau21<=0.45) {ntau = ntau+1;} else{ isAnaHP=-1; }
			if (isAnaHP>0)// && (( jetAK8puppi_sdcorr >0&& jetAK8puppi_sdcorr< 150 )))// && m_lvj>100 && m_lvj<5000 )
			{ 
				nmassVhad = nmassVhad +1; 
				(*file_cutflow)<<"event:"<<event<<endl;
			} else{ isAnaHP=-1; }

			//LP: 0.5<jetAK8puppi_tau21<=0.75;
			if (lep==13 && (HLT_Mu2>0 || HLT_Mu3>0 ) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && nLooseLep==1 && ptlep1>55 && MET_et>40 && ptVlepJEC>200 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && jetAK8puppi_tau21>0.45 && jetAK8puppi_tau21<=0.75)// && (( jetAK8puppi_sdcorr >0&& jetAK8puppi_sdcorr< 150 )))
			{ isAnaLP=1.; } 
			else{ isAnaLP=-1.; }

			//NP: 0.75<jetAK8puppi_tau21;
			if (lep==13 && (HLT_Mu2>0 || HLT_Mu3>0 ) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && nLooseLep==1 && ptlep1>55 && MET_et>40 && ptVlepJEC>200 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 && num_bJet<1 && deltaRlepjet>pi_2 && fabs(delPhijetmet) >2.0 && fabs(delPhijetlep)>2.0 && jetAK8puppi_tau21>0.75)//  && (( jetAK8puppi_sdcorr >0&& jetAK8puppi_sdcorr< 150 )))
			{ isAnaNP=1.; } 
			else{ isAnaNP=-1.; }

			//TTbar control
			if (lep==13 && (HLT_Mu2>0 || HLT_Mu3>0 ) && trackIso/ptlep1<0.1 && muisolation<0.05 && fabs(etalep1)<2.4 && nLooseLep==1 && ptlep1>55 && MET_et>40 && ptVlepJEC>200 && jet_pt_puppi>200 && fabs(jetAK8puppi_eta)<2.4 && IDLoose>0 && num_bJet>0)// && jetAK8puppi_sdcorr>0 && jetAK8puppi_sdcorr <150)
			{ isTTBarControl=1.; } 
			else{ isTTBarControl=-1.; }

		}else{
			cout<<"We don't know channelname:"<<channelname<<endl;
		}

		Int_t tmp_categoryID_eventselection=0;
		if(isAnaHP>0)tmp_categoryID_eventselection=1;
		else if(isAnaLP>0)tmp_categoryID_eventselection=2;
		else if(isAnaNP>0)tmp_categoryID_eventselection=4;
		else if(isTTBarControl>0)tmp_categoryID_eventselection=3;
		else tmp_categoryID_eventselection=100;

		CategoryID=tmp_categoryID_channel* tmp_categoryID_eventselection;

		isMatch=1.;
		if(deltaRWhadGen >= 0.3) isMatch=-1;
		//cout << "massVhad" << massVhad << "jet_mass_pruned " << jet_mass_pruned << endl;
		if(jetAK8puppi_tau21<=0){vTagID=2;}
		else if(jetAK8puppi_tau21>0.45 && jetAK8puppi_tau21<=0.60){vTagID=1;}
		else if(jetAK8puppi_tau21>0.60 && jetAK8puppi_tau21<=0.75){vTagID=0;}
		else if(jetAK8puppi_tau21>0.75 && jetAK8puppi_tau21<=1){vTagID=-1;}
		else {vTagID=-2;}

		if(TMath::Abs(CategoryID)<10) ExTree->Fill();
	}

	if(channelname=="el"){ 
		std::cout << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << "; yields " << yields << std::endl;
		(*file_cutflow) << "nID_e" << nID_e << "; npt_e" << npt_e << "; nmet_e" << nmet_e << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<<"; nnum_bJet_e" << nnum_bJet_e <<"; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad << "; number_qq" << number_qq << std::endl;
	}
	if(channelname=="mu"){
		std::cout << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
		(*file_cutflow) << "nID_mu" << nID_mu << "; npt_mu" << npt_mu << "; nmet_mu" << nmet_mu << "; nptVlepJEC" << nptVlepJEC << "; nptVhad" << nptVhad<< "; nnum_bJet_mu" << nnum_bJet_mu << "; n_deltaRlepjet"<<n_deltaRlepjet<< "; n_delPhijetmet" << n_delPhijetmet <<"; n_delPhijetlep"<<n_delPhijetlep<<"; ntau"<<ntau<< "; nmassVhad" << nmassVhad<< "; number_qq" << number_qq << std::endl;
	}

}
