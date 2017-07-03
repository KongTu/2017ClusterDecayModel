#include "../interface/EventGenerator.h"

using namespace std;

/*
ETA values are multiplied by 10
*/

int main(){

	cout << "----beginning of the Producer----" << endl;

/*
temp config
*/
	TH2D* map = (TH2D*) file1->Get("mother_Spectra");
	f2->SetParameter(0,0.07);//data v2 is about 7% at 185-250.
	f2->SetParameter(1,0.0);// event plane is at zero.

	f3->SetParameter(0, 6.59660e+04);
	f3->SetParameter(1, -1.93218e+00);
	f3->SetParameter(2, -6.91020e-01);
	f3->SetParameter(3, 1.75804e+01);

	double clusterMass = 0.0;//will generate realistically
	const int Nevents = 3000000;

/*
*/
	TCanvas* c1 = new TCanvas();
	TH1D* c2_tracker = new TH1D("c2_tracker","c2",100,-1,1);
	TH1D* Ntrkoffline = new TH1D("Ntrkoffline","Ntrk", 1000,0,1000);
	TH1D* Ntrkoffline_Count = new TH1D("Ntrkoffline_Count","Ntrk", 500,0,500);	
	TH1D* ptspectra = new TH1D("ptspectra","pt", 100,0,10);

	TH1D* c3_real[48][3];
	TH1D* c2_real[48][3];
	TH1D* delEta3p[3];

	for(int sign = 0; sign < 3; sign++){

		delEta3p[sign] = new TH1D(Form("delEta3p_%d",sign),Form("delEta3p_%d",sign), NdEtaBins, dEtaBinsArray);
		
		for(int deta = 0;deta < NdEtaBins; deta++){

			c3_real[deta][sign] = new TH1D(Form("c3_real_%d_%d", deta, sign),"c3_real", 1000, -1,1 );
			c2_real[deta][sign] = new TH1D(Form("c2_real_%d_%d", deta, sign),"c2_real", 1000, -1,1 );

		}
	}

/*
Generate events with multiplicity 185<Ntrkoffline<250
*/


	Event evt;
	evt.GenerateEvent(Nevents);
	vector<int> eventMult = evt.EventSize();
	
//Start the event loop:
	for(int i = 0; i < eventMult.size(); i++){

		//cout << "event: " << i << endl;

		TComplex Q_n3_trk, Q_0_trk;

		TComplex Q_n1_1[NetaBins][2], Q_n2_1[NetaBins][2];

		TComplex Q_n1n2_2[NetaBins][2];

		TComplex Q_0_1[NetaBins][2], Q_0_2[NetaBins][2];

		TComplex P_n1_1[NetaBins][2], P_n2_1[NetaBins][2];

		TComplex P_n1n2_2[NetaBins][2];

		TComplex P_0_1[NetaBins][2], P_0_2[NetaBins][2];
	
	/*
	Generate particles with rho mass
	*/
		int mult = eventMult[i];
		Ntrkoffline->Fill( mult );

		
		int mult_counting = 0;
		TVector3 v3(0,0,0);
		for(int j = 0; j < int(mult-105); j++){

			//step3: start to decay

			clusterMass = evt.GetMotherMass();//according to particle ratio;
		   
			evt.GenerateParticle( true, map, clusterMass);

			double pt = evt.GetMomPt();
			double eta = evt.GetMomEta();
			double phi = evt.GetMomPhi();

			vector<double> cluster4Momentum = get4Momentum(pt, eta, phi, clusterMass);

			double energy_total = cluster4Momentum[0];
		    double cluster_px = cluster4Momentum[1];
		    double cluster_py = cluster4Momentum[2];
		    double cluster_pz = cluster4Momentum[3];

		    TVector3 v1(cluster_px,cluster_py,cluster_pz);
		    
			vector<double> motherParticles;

			int k = 20;//local momentum conservation
			if( j % k != k-1 ){

				v3 += v1;
				motherParticles = getLightConeVar(v1(0), v1(1), v1(2));
			}
			else{

				motherParticles = getLightConeVar(-v3(0), -v3(1), -v3(2));//modify the k particle to conserve momentum of the previous k-1 particles
				v3(0) = 0.0; v3(1) = 0.0; v3(2) = 0.0;
			}

			vector< vector<double>> daugterParticles_total = evt.DecayParticle(motherParticles, clusterMass);//particle decay recursively. 

	   		ptspectra->Fill( motherParticles[0] );

		   	for(int i = 0; i < daugterParticles_total.size(); i++){

		   		mult_counting++;

		   		double dau1_pt = daugterParticles_total[i][0];
			   	double dau1_eta = daugterParticles_total[i][1]*10;//multiply all eta by 10.
			   	double dau1_phi = daugterParticles_total[i][2];
			   	double dau1_charge = daugterParticles_total[i][3];

			   	if( dau1_pt > 0.3 && dau1_pt < 3.0 && fabs(dau1_eta) < 24 ){
					for(int eta = 0; eta < NetaBins; eta++){
						if( dau1_eta > etaBins[eta] && dau1_eta < etaBins[eta+1] ){
							if( dau1_charge == 1 ){

								Q_n1_1[eta][0] += q_vector(1, 1, 1, dau1_phi);
								Q_n2_1[eta][0] += q_vector(1, 1, 1, dau1_phi);

								Q_n1n2_2[eta][0] += q_vector(2, 2, 1, dau1_phi);

								Q_0_1[eta][0] += q_vector(0, 1, 1, dau1_phi);
								Q_0_1[eta][0] += q_vector(0, 2, 1, dau1_phi);

								P_n1_1[eta][0] += q_vector(1, 1, 1, dau1_phi);
								P_n2_1[eta][0] += q_vector(-1, 1, 1, dau1_phi);//it is a minus n2_ because n2_ = 1

								P_n1n2_2[eta][0] += q_vector(0, 2, 1, dau1_phi);

								P_0_1[eta][0] += q_vector(0, 1, 1, dau1_phi);
								P_0_2[eta][0] += q_vector(0, 2, 1, dau1_phi);
							}
							if( dau1_charge == -1 ){

								Q_n1_1[eta][1] += q_vector(1, 1, 1, dau1_phi);
								Q_n2_1[eta][1] += q_vector(1, 1, 1, dau1_phi);

								Q_n1n2_2[eta][1] += q_vector(2, 2, 1, dau1_phi);

								Q_0_1[eta][1] += q_vector(0, 1, 1, dau1_phi);
								Q_0_1[eta][1] += q_vector(0, 2, 1, dau1_phi);

								P_n1_1[eta][1] += q_vector(1, 1, 1, dau1_phi);
								P_n2_1[eta][1] += q_vector(-1, 1, 1, dau1_phi);//it is a minus n2_ because n2_ = 1

								P_n1n2_2[eta][1] += q_vector(0, 2, 1, dau1_phi);

								P_0_1[eta][1] += q_vector(0, 1, 1, dau1_phi);
								P_0_2[eta][1] += q_vector(0, 2, 1, dau1_phi);
							}
						}
					}
				
					Q_n3_trk += q_vector(2, 1, 1, dau1_phi);
					Q_0_trk += q_vector(0, 1, 1, dau1_phi);
				}
		   	
		   	}//daughter loop
		}//track loop
	

		TComplex N_2_trk;
		double D_2_trk;

		N_2_trk = Q_n3_trk*TComplex::Conjugate(Q_n3_trk) - Q_0_trk;
		D_2_trk = Q_0_trk.Re()*(Q_0_trk.Re() - 1.);

		c2_tracker->Fill(N_2_trk.Re()/D_2_trk, D_2_trk);

		for(int ieta = 0; ieta < NetaBins; ieta++){
			for(int jeta = 0; jeta < NetaBins; jeta++){

	  			double deltaEta = fabs( etaBins[ieta] - etaBins[jeta] );
	  			
	  			for(int deta = 0; deta < NdEtaBins; deta++){
			        if( deltaEta > dEtaBins[deta]-0.001 && deltaEta < dEtaBins[deta+1]-0.001 ){
			          
						TComplex N_2;
						TComplex D_2;

						TComplex N_2_P;
						TComplex D_2_P;

						//same sign correlator:
						for(int sign = 0; sign < 2; sign++){
							if( ieta == jeta ){

							  delEta3p[sign]->Fill( deltaEta );

							  N_2 = Q_n1_1[ieta][sign]*Q_n2_1[ieta][sign] - Q_n1n2_2[ieta][sign];
							  D_2 = Q_0_1[ieta][sign]*Q_0_1[ieta][sign] - Q_0_2[ieta][sign];

							  N_2_P = P_n1_1[ieta][sign]*P_n2_1[ieta][sign] - P_n1n2_2[ieta][sign];
							  D_2_P = P_0_1[ieta][sign]*P_0_1[ieta][sign] - P_0_2[ieta][sign];

							}
							else{

							  delEta3p[sign]->Fill( deltaEta );

							  N_2 = Q_n1_1[ieta][sign]*Q_n2_1[jeta][sign];
							  D_2 = Q_0_1[ieta][sign]*Q_0_1[jeta][sign];

							  N_2_P = P_n1_1[ieta][sign]*P_n2_1[jeta][sign];
							  D_2_P = P_0_1[ieta][sign]*P_0_1[jeta][sign];

							}

							c3_real[deta][sign]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
							c2_real[deta][sign]->Fill(N_2_P.Re()/D_2_P.Re(), D_2_P.Re() );

						}

						delEta3p[2]->Fill( deltaEta );

						N_2 = Q_n1_1[ieta][0]*Q_n2_1[jeta][1];
						D_2 = Q_0_1[ieta][0]*Q_0_1[jeta][1];


						N_2_P = P_n1_1[ieta][0]*P_n2_1[jeta][1];
					  	D_2_P = P_0_1[ieta][0]*P_0_1[jeta][1];

						c3_real[deta][2]->Fill(N_2.Re()/D_2.Re(), D_2.Re() );
						c2_real[deta][2]->Fill(N_2_P.Re()/D_2_P.Re(), D_2_P.Re() );

					}
	    		}
			}
		}
		Ntrkoffline_Count->Fill( mult_counting );
	}//event loop

	cout << "v2: " << sqrt( c2_tracker->GetMean() ) << endl;

	TCanvas* c10 = new TCanvas();
	gPad->SetLogy(1);
	ptspectra->Draw();
	c10->Print("pt.pdf");

	TCanvas* c11 = new TCanvas();
	gPad->SetLogy(1);
	Ntrkoffline_Count->Draw();
	c11->Print("Ntrk.pdf");

	TFile f1("./rootfiles/job_2.root","RECREATE");
	Ntrkoffline->Write();
	c2_tracker->Write();

	for(int sign = 0; sign < 3; sign++){
		delEta3p[sign]->Write();
		for(int i = 0; i < NdEtaBins; i++){
			c3_real[i][sign]->Write();
			c2_real[i][sign]->Write();
		}
	}


}