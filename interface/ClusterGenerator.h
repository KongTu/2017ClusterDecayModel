#ifndef CLUSTER_GENERATOR_H
#define CLUSTER_GENERATOR_H

#include "ClusterProducer.h"

using namespace std;

TH1D* daugterSample = new TH1D("daugterSample","daugterSample",3,0,3);
TH1D* NumberOfDau = new TH1D("NumberOfDau","NumberOfDau",4,2,6);

double massFinal[] = {PION, KAON, P};
int clusterCharge = 10;

class Cluster {

  private:
    double pt_m;
    double eta_m;
    double phi_m;
    double mass_m;
    double charge_m;

  public:
    void GenerateCluster (bool, TH2D*);
    double GetMomPt(){return pt_m;}
    double GetMomEta(){return eta_m;}
    double GetMomPhi(){return phi_m;}
    double GetMomMass(){return mass_m;}
    double GetMomCharge(){return charge_m;}
    vector< vector<double>> GetDecayProducts( vector<double>, int );

} clus;


void Cluster::GenerateCluster(bool doFlow_, TH2D* map){

	double pt, eta, phi, mass;
    int charge;
    map->GetRandom2(eta, pt);
    //pt = f3->GetRandom();
    eta = f4->GetRandom();
    mass = f5->GetRandom();
    unsigned int prop = gRandom->Integer(2);
    if(prop == 0){ charge = (int) gRandom->Integer( clusterCharge ); }
    else if( prop == 1){ charge = (int) -(gRandom->Integer( clusterCharge ));} 
    
    if( doFlow_ ){
        phi = f2->GetRandom();
    }
    else{phi = f1->GetRandom();}

    pt_m = pt;
    eta_m = eta;
    phi_m = phi;
    mass_m = mass;
    charge_m = charge;
}

vector< vector<double>> Cluster::GetDecayProducts( vector<double> mother, int N ){

    const int Ndau = N;
    double masses[Ndau];
    double charges[Ndau];
    vector< vector<double>> daugterParticles;

    double pt = mother[0];
    double eta = mother[1];
    double phi = mother[2];
    double mass = mother[3];
    int charge = mother[4];//mother charge not used yet. 
    
    unsigned int j = gRandom->Integer(2);
    TGenPhaseSpace event;

    if( Ndau == 2 ){
        charges[0] = 1; charges[1] = -1;
    }
    else if( Ndau == 3 ){
        charges[0] = 1; charges[1] = -1; charges[2] = 0;
    }
    else if( Ndau == 4 ){
        charges[0] = 1; charges[1] = -1; charges[2] = 1; charges[3] = -1;
    }
    else if( Ndau == 5 ){
        charges[0] = 1; charges[1] = -1; charges[2] = 1; charges[3] = -1; charges[4] = 0; 
    }
            
    for(int i = 0; i < Ndau; i++){

        double temp = daugterSample->GetRandom();
        masses[i] = massFinal[int(temp)];
        //charges[i] = pow(-1,i+j);
    }

    vector<double> cluster4Momentum = get4Momentum(pt, eta, phi, mass);

    double energy_total = cluster4Momentum[0];
    double cluster_px = cluster4Momentum[1];
    double cluster_py = cluster4Momentum[2];
    double cluster_pz = cluster4Momentum[3];

    TLorentzVector p4(cluster_px,cluster_py,cluster_pz,energy_total);

    bool DECAY_ = event.SetDecay(p4,Ndau,masses);

    if( !DECAY_ ) return daugterParticles;
        
    double weight1 = event.Generate();

    for(int i = 0; i < Ndau; i++){

        TLorentzVector* pPion1 = event.GetDecay(i);
        double dau1_px = pPion1->Px();
        double dau1_py = pPion1->Py();
        double dau1_pz = pPion1->Pz();

        vector<double> dau1_LightConeVar = getLightConeVar(dau1_px, dau1_py, dau1_pz);
        dau1_LightConeVar.push_back( masses[i] );
        dau1_LightConeVar.push_back( charges[i] );

        daugterParticles.push_back( dau1_LightConeVar );

        dau1_LightConeVar.clear();

    }


    return daugterParticles;
}

#endif