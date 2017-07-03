#include "RiceStyle.h"

using namespace std;

double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;

void makeTGraph_cme_pPb_differentials(){

	TFile* file = new TFile("../rootfiles/test.root");
	TH1D* Ntrk = (TH1D*) file->Get("Ntrkoffline");
	
	TH1D* QvsdEta[48][3];
	TH1D* PvsdEta[48][3];

	for(int sign = 0; sign < 3; sign++){
		for(int deta = 0; deta < NdEtaBins; deta++){
	  
			QvsdEta[deta][sign] = (TH1D*) file->Get( Form("c3_real_%d_%d",deta,sign) ); 
			PvsdEta[deta][sign] = (TH1D*) file->Get( Form("c2_real_%d_%d",deta,sign) );
		}
	}
	
	TH1D* hist1[3];
	TH1D* hist2[3];

	for(int sign = 0; sign < 3; sign++){

		hist1[sign] = new TH1D(Form("hist1_%d_%d",sign),"test1", NdEtaBins, dEtaBins);
		hist2[sign] = new TH1D(Form("hist2_%d_%d",sign),"test1", NdEtaBins, dEtaBins);
	
	}

	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){

			double P_total_real_dEta1 = PvsdEta[deta][sign]->GetMean();
			double P_total_real_dEta_error1 = PvsdEta[deta][sign]->GetMeanError();

			double value = P_total_real_dEta1;
			double error = P_total_real_dEta_error1;
			
			hist2[sign]->SetBinContent(deta+1, value );
			hist2[sign]->SetBinError(deta+1, error );

			double Q_total_real_dEta1 = QvsdEta[deta][sign]->GetMean();
			double Q_total_real_dEta_error1 = QvsdEta[deta][sign]->GetMeanError();

			double value = Q_total_real_dEta1;
			double error = Q_total_real_dEta_error1;

			hist1[sign]->SetBinContent(deta+1, value );
			hist1[sign]->SetBinError(deta+1, error );
		}
	}

//begin of deta

	TH1D* base1 = makeHist("base1", "", "|#Delta#eta|", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0012,0.0014);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);


	TH1D* temp1 = (TH1D*)hist1[0]->Clone("temp1");
	temp1->Add(hist1[1], +1);
	temp1->Scale(0.5);
	temp1->SetMarkerStyle(20);
	temp1->SetMarkerSize(1.4);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp2 = (TH1D*) hist1[2]->Clone("temp2");
	temp2->SetMarkerStyle(21);
	temp2->SetMarkerColor(kBlue);
	temp2->SetMarkerSize(1.4);
	temp2->SetLineColor(kBlue);

	TH1D* temp3 = (TH1D*)hist2[0]->Clone("temp3");
	temp3->Add(hist2[1], +1);
	temp3->Scale(0.5);
	temp3->SetMarkerStyle(24);
	temp3->SetMarkerSize(1.4);
	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);

	TH1D* temp4 = (TH1D*) hist2[2]->Clone("temp4");
	temp4->SetMarkerStyle(25);
	temp4->SetMarkerSize(1.3);
	temp4->SetMarkerColor(kBlue);
	temp4->SetLineColor(kBlue);

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	base1->Draw();
	temp1->Draw("Psame");
	temp2->Draw("Psame");
	
	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();

	TH1D* base2 = makeHist("base2", "", "|#Delta#eta|", "#LTcos(#phi_{#alpha}-#phi_{#beta})#GT", 48,0,4.8,kBlack);

	base2->GetYaxis()->SetRangeUser(-0.012,0.014);
	base2->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetXaxis()->SetTitleOffset(0.95);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);

	base2->Draw();
	temp3->Draw("Psame");
	temp4->Draw("Psame");




/*
Save in a root file
*/




}