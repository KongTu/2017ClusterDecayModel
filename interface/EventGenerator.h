#ifndef EVENT_GENERATOR_H
#define EVENT_GENERATOR_H


#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>

#include <TMath.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TComplex.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TF1.h>

#include "ClusterProducer.h"

using namespace std;

#define PION 0.139
#define KAON 0.493
#define P    0.938
#define K0S  0.497
#define LAM  1.116
#define XI   1.321
#define OMEGA 1.672
#define RHO  0.775
#define DELTApp 1.232

double massP[] = {RHO, XI, LAM, K0S, DELTApp, PION, KAON, P, OMEGA};
double ratioP[] = {0.3125, 0.3125, 0.0625,0.125,0.0,0.1125,0.025,0.0625,0.0125};
double chargeP[] = {0, -1, 0, 0, 2, 1, 1, 1, -1};
bool IsDecay[] = {true, true, true, true, true, false, false, false, true};
double dau1[] = {PION, LAM, P, PION, P, 0.0, 0.0, 0.0, LAM};
double dau2[] = {PION, PION, PION, PION, PION, 0.0, 0.0, 0.0, KAON};
double charge1[] = {+1, 0, +1, +1, +1, 0, 0, 0, 0};
double charge2[] = {-1, -1, -1, -1, +1, 0, 0, 0, -1};
const int Nparticles = sizeof(massP) / sizeof(massP[0]) - 1;

TFile* file1 = new TFile("~/cernbox/2015RUN2work/2017ClusterDecayModel/data/rho_map.root"); 
TFile* file2 = new TFile("~/cernbox/2015RUN2work/2017ClusterDecayModel/data/Ntrkoffline.root");
TF1* f1 = new TF1("f1","1",-3.14,3.14);//flat phi distribution.
TF1* f2 = new TF1("f2","1+2*[0]*cos(2*x-2*[1])",-PI,PI);//flow
TF1* f3 = new TF1("f3","[0]*x*TMath::Exp([1]*x)+[3]*pow(x,[2])",0.3,100);
TH1D* prop = new TH1D("prop","prop",Nparticles,0,Nparticles);

class Event {

  private:
  	vector<int> eventMultiplicity;
  	double pt_m;
  	double eta_m;
  	double phi_m;

  public:
    void GenerateEvent (int);
    vector< vector<double>> DecayParticle( vector<double>, double );
    void GenerateParticle (bool, TH2D*, double);
    double GetMotherMass();
    vector<int> EventSize();
    double GetMomPt(){return pt_m;}
    double GetMomEta(){return eta_m;}
    double GetMomPhi(){return phi_m;}
} evt;


void Event::GenerateEvent(int N){

	TH1D* Ntrk = (TH1D*) file2->Get("Ntrk");

	for(int i = 0; i < N; i++){

		int mult = Ntrk->GetRandom();
		eventMultiplicity.push_back(mult);
	}
}
void Event::GenerateParticle(bool doFlow_, TH2D* map, double mass){
 
	double pt, eta, phi;
	map->GetRandom2(eta,pt);

	//pt = f3->GetRandom();
	
	if( doFlow_ ){
	 	phi = f2->GetRandom();
	 	pt = pt*(1+0.2*(mass/RHO));//radial flow
	}
    else{phi = f1->GetRandom();}

	pt_m = pt;
	eta_m = eta;
	phi_m = phi;

}
double Event::GetMotherMass(){
    
    for(int index = 0; index < Nparticles; index++){

        prop->SetBinContent(index+1, ratioP[index]);
    }
    double temp = prop->GetRandom();

    return massP[int(temp)];
}
vector< vector<double>> Event::DecayParticle( vector<double> kinematics, double mass ){

	double pt = kinematics[0];
	double eta = kinematics[1];
	double phi = kinematics[2];
	double clusterMass = mass;

    vector< vector<double>> total_daughter;
    vector<double> dau_mass, dau_charge;

    int counter = 0;
    for(int index = 0; index < Nparticles; index++){
        if( IsDecay[index] == true && mass == massP[index] ){

            dau_mass.push_back( dau1[index] );
            dau_mass.push_back( dau2[index] );

            if( f1->GetRandom() > 0 ){
                charge1[index] = -charge1[index];
                charge2[index] = -charge2[index];
            }

            dau_charge.push_back( charge1[index] );
            dau_charge.push_back( charge2[index] );

            counter++;
        }
        else if( IsDecay[index] == false && mass == massP[index] ){

            vector<double> temp;
            if( f1->GetRandom() > 0 ){
                chargeP[index] = -chargeP[index];
            }
            temp.push_back( pt ); temp.push_back( eta ); temp.push_back( phi ); temp.push_back( chargeP[index] );
            total_daughter.push_back( temp );
            temp.clear();

            counter++;
            return total_daughter;
        }
    }

    if( counter == 0 ) {cout << "something is wrong!" << endl; return total_daughter;}

	vector<double> cluster4Momentum = get4Momentum(pt, eta, phi, clusterMass);

    double energy_total = cluster4Momentum[0];
    double cluster_px = cluster4Momentum[1];
    double cluster_py = cluster4Momentum[2];
    double cluster_pz = cluster4Momentum[3];

    TLorentzVector p4(cluster_px,cluster_py,cluster_pz,energy_total);
    TGenPhaseSpace event;

    const int Ndau = dau_mass.size();
    double masses[Ndau];
    for(int i = 0; i< Ndau; i++){
    	masses[i] = dau_mass[i];
    }

    event.SetDecay(p4,Ndau,masses);
    double weight1 = event.Generate();

    for(int n = 0; n < Ndau; n++){

    	TLorentzVector* pPion1 = event.GetDecay(n);
    	double dau1_px = pPion1->Px();
	    double dau1_py = pPion1->Py();
	    double dau1_pz = pPion1->Pz();

        vector<double> dau1_LightConeVar = getLightConeVar(dau1_px, dau1_py, dau1_pz);

        bool decay_Again = false;

        for(int index = 0; index < Nparticles; index++){

            if(dau_mass[n] == massP[index]){

                decay_Again = IsDecay[index];
            }
        }

        if( decay_Again ){

            total_daughter = Event::DecayParticle(dau1_LightConeVar, dau_mass[n]);
        }
        else{

            dau1_LightConeVar.push_back( dau_charge[n] );
            total_daughter.push_back( dau1_LightConeVar );
            dau1_LightConeVar.clear();

        }

    }	

    return total_daughter;

}
vector<int> Event::EventSize(){

	return eventMultiplicity;
}




#endif












