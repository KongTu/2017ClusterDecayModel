#include "../interface/EventGenerator.h"
#include "../interface/ClusterGenerator.h"

using namespace std;

/*
ETA values are multiplied by 10
*/

int main(){

	cout << "----beginning of the Producer----" << endl;

/*
temp config
free parameters are:

1. v2 of the clusters, 
2. mass distribution of clusters,
3. Number of daughters, 
4. daughters' mass contributions,
5. mother cluster's momentum,
*/

TH2D* map = (TH2D*) file1->Get("mother_Spectra");

double pt_slope[] = {-0.1,-0.3,-0.5,-0.7,-0.9,-1.1,-1.3,-1.5,-1.7,-1.9,-2.2,-2.4,-2.6,-3.0};
double Nbody[] = {0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,0.79};
double mass_slope[] = {-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.8,-2.2,-3.0};
int momentumCons[] = {2,3,4,10,15,20,30,40,120,300};

for( int scan = 0; scan < 1; scan++){

	//f2 is v2 input
	f2->SetParameter(0,0.3);//data v2 is about 7% at 185-250.
	f2->SetParameter(1,0.0);// event plane is at zero.

	//f3 is pt spectra of input
	f3->SetParameter(0, 6.59660e+04);//6.5
	f3->SetParameter(1, -2.0);//-1.9
	f3->SetParameter(2, -6.91020e-01);
	f3->SetParameter(3, 1.75804e+01);

	//f5 is mass spectra of input
	f5->SetParameter(0, 1.0);
	f5->SetParameter(1, -1.2);
	f5->SetParameter(2, 0);

	//daughter composition input
	daugterSample->SetBinContent(1, 0.7);//PION
	daugterSample->SetBinContent(2, 0.1);//KAON
	daugterSample->SetBinContent(3, 0.05);//P

	NumberOfDau->SetBinContent(1, 0.99);//2 body
	NumberOfDau->SetBinContent(2, 0.01);//3 body
	NumberOfDau->SetBinContent(3, 0.0);//4 body
	NumberOfDau->SetBinContent(4, 0.0);//5 body

	bool doFlow_ = true;	
	int Ndaughters = 2;//reset later according to contribution
	const int NclusterToConserveP = 6;
	const int Nevents = 50000;
	const int Nclusters = 430;

/*
*/
	TCanvas* c1 = new TCanvas();
	TH1D* c2_tracker = new TH1D("c2_tracker","c2",100,-1,1);
	TH1D* Ntrkoffline = new TH1D("Ntrkoffline","Ntrk", 1000,0,1000);
	TH1D* Ntrkoffline_Count = new TH1D("Ntrkoffline_Count","Ntrk", 500,0,500);	
	TH1D* ptspectra = new TH1D("ptspectra","pt", 100,0.0,20);
	TH1D* massSpectra = new TH1D("massSpectra","mass", 100,0.0,20);

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
	
	Cluster clus;
	
//Start the event loop:
	for(int i = 0; i < Nevents; i++){

		TComplex Q_n3_trk, Q_0_trk;

		TComplex Q_n1_1[NetaBins][2], Q_n2_1[NetaBins][2];

		TComplex Q_n1n2_2[NetaBins][2];

		TComplex Q_0_1[NetaBins][2], Q_0_2[NetaBins][2];

		TComplex P_n1_1[NetaBins][2], P_n2_1[NetaBins][2];

		TComplex P_n1n2_2[NetaBins][2];

		TComplex P_0_1[NetaBins][2], P_0_2[NetaBins][2];
	
	/*
	Generate events
	*/

		int mult_counting = 0;
		TVector3 v3(0,0,0);
		for(int j = 0; j < Nclusters; j++){

	//step3: start to generate clusters
			clus.GenerateCluster( doFlow_, map );
			vector< double> motherCluster;
			
			motherCluster.push_back( clus.GetMomPt() );
			motherCluster.push_back( clus.GetMomEta() );
			motherCluster.push_back( clus.GetMomPhi() );
//			motherCluster.push_back( clus.GetMomMass() );
			motherCluster.push_back( 0.776 );//rho mass

			motherCluster.push_back( clus.GetMomCharge() );

			vector<double> cluster4Momentum = get4Momentum(motherCluster[0], motherCluster[1], motherCluster[2], motherCluster[3]);

			double energy_total = cluster4Momentum[0];
		    double cluster_px = cluster4Momentum[1];
		    double cluster_py = cluster4Momentum[2];
		    double cluster_pz = cluster4Momentum[3];

		    TVector3 v1(cluster_px,cluster_py,cluster_pz);
		    
			vector<double> motherParticles;
	
	//conserve N clusters' momentum
			int k = NclusterToConserveP;//local momentum conservation
			if( j % k != k-1 ){

				v3 += v1;
				motherParticles = getLightConeVar(v1(0), v1(1), v1(2));
			}
			else{

				motherParticles = getLightConeVar(-v3(0), -v3(1), -v3(2));//modify the k particle to conserve momentum of the previous k-1 particles
				v3(0) = 0.0; v3(1) = 0.0; v3(2) = 0.0;
			}

			motherParticles.push_back( motherCluster[3] );//mass
			motherParticles.push_back( motherCluster[4] );//charge

	//Cluster decay:
			Ndaughters = (int) NumberOfDau->GetRandom();
			vector< vector<double>> daugterParticles_total = clus.GetDecayProducts(motherParticles, Ndaughters);
			
		   	for(int i = 0; i < daugterParticles_total.size(); i++){

		   		double dau1_pt = daugterParticles_total[i][0];
			   	double dau1_eta = daugterParticles_total[i][1]*10;//multiply all eta by 10.
			   	double dau1_phi = daugterParticles_total[i][2];
			   	double dau1_charge = daugterParticles_total[i][4];//charge

   				massSpectra->Fill( motherParticles[3] );
			   	ptspectra->Fill( dau1_pt );
   				
			   	if( dau1_pt > 0.4 && fabs( dau1_eta ) < 24 ) {mult_counting++;}
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
	c10->Print(Form("./figures/pt_%d.pdf", scan));

	TCanvas* c11 = new TCanvas();
	gPad->SetLogy(1);
	Ntrkoffline_Count->Draw();
	c11->Print(Form("./figures/Ntrk_%d.pdf", scan));

	TCanvas* c12 = new TCanvas();
	gPad->SetLogy(1);
	massSpectra->Draw();
	c12->Print(Form("./figures/mass_%d.pdf", scan));

	TFile f1(Form("./rootfiles/decayModel_rho_matchedData.root",scan),"RECREATE");
	Ntrkoffline->Write();
	c2_tracker->Write();
	Ntrkoffline_Count->Write();
	ptspectra->Write(); 
	massSpectra->Write();

	for(int sign = 0; sign < 3; sign++){
		delEta3p[sign]->Write();
		for(int i = 0; i < NdEtaBins; i++){
			c3_real[i][sign]->Write();
			c2_real[i][sign]->Write();
		}
	}

	c2_tracker->Delete();  
	Ntrkoffline->Delete(); 
	Ntrkoffline_Count->Delete();
	ptspectra->Delete(); 
	massSpectra->Delete(); 
	for(int sign = 0; sign < 3; sign++){
		delEta3p[sign]->Delete();
		for(int i = 0; i < NdEtaBins; i++){
			c3_real[i][sign]->Delete();
			c2_real[i][sign]->Delete();
		}
	}

}

}