#ifndef CLUSTER_PRODUCER_H
#define CLUSTER_PRODUCER_H

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

using namespace std;

#define PI 3.1415926

int etaBins[] = {-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
const int NetaBins = sizeof(etaBins) / sizeof(etaBins[0]) - 1;
int dEtaBins[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,34,38,42,48};
double dEtaBinsArray[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,34,38,42,48};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;

TFile* file1 = new TFile("~/cernbox/2015RUN2work/2017ClusterDecayModel/data/rho_map.root"); 
TFile* file2 = new TFile("~/cernbox/2015RUN2work/2017ClusterDecayModel/data/Ntrkoffline.root");
TF1* f1 = new TF1("f1","1",-3.14,3.14);//flat phi distribution.
TF1* f2 = new TF1("f2","1+2*[0]*cos(2*x-2*[1])",-PI,PI);//flow
TF1* f3 = new TF1("f3","[0]*x*TMath::Exp([1]*x)+[3]*pow(x,[2])",0.01,3);
TF1* f4 = new TF1("f4","1",-5,5);//flat eta distribution.
TF1* f5 = new TF1("f5","[0]*TMath::Exp([1]*x)+[2]",0.3,2.0);

vector<double> get4Momentum(double pt, double eta, double phi, double mass)
{
  double polar_angle = 2*TMath::ATan( TMath::Exp(-eta) );
  double pz = pt/TMath::Tan( polar_angle );
  double px = sqrt(pt*pt/( 1+TMath::Tan(phi)*TMath::Tan(phi) ) );
  double py = sqrt(pt*pt - px*px);
  double E = sqrt(px*px+py*py+pz*pz + mass*mass);

  vector<double> temp;

  if( phi > 0 && phi < PI){
    py = py; 
  }
  else if( phi < 0 && phi > -PI){
    py = -py;
  }

  if( phi > -PI/2.0 && phi < PI/2.0 ){
    px = px;
  }
  else if( phi < -PI/2.0 || phi > PI/2.0 ){
    px = -px;
  }
  
  temp.push_back( E );
  temp.push_back( px );
  temp.push_back( py );
  temp.push_back( pz );
  temp.push_back( polar_angle ); 

  return temp;

}
vector<double> getLightConeVar(double px, double py, double pz){

  double pt = sqrt(px*px + py*py);
  double phi = TMath::ATan(py/px);
  double three_momentum = sqrt(px*px+py*py+pz*pz);
  double eta = TMath::ATanH( pz/three_momentum );

  vector<double> temp;

  if( px > 0 && py > 0){
    phi = phi;
  }
  else if( px > 0 && py < 0){
    phi = phi;
  }
  else if( px < 0 && py > 0){
    phi = PI + phi;
  }
  else if( px < 0 && py < 0){
    phi = -PI + phi;
  }

  temp.push_back( pt );
  temp.push_back( eta );
  temp.push_back( phi );

  return temp; 
}
TComplex q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}


#endif


