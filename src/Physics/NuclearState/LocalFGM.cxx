//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Author: Joe Johnston, University of Pittsburgh (Advisor Steven Dytman)

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>

#include "TFile.h"
#include "TRandom.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/LocalFGM.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include "AV18.hh"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;
using namespace genie::utils;


//____________________________________________________________________________
LocalFGM::LocalFGM() :
NuclearModelI("genie::LocalFGM")
{

}
//____________________________________________________________________________
LocalFGM::LocalFGM(string config) :
NuclearModelI("genie::LocalFGM", config)
{

}
//____________________________________________________________________________
LocalFGM::~LocalFGM()
{

}
//____________________________________________________________________________
bool LocalFGM::GenerateNucleon(const Target & target,
				      double hitNucleonRadius) const
{
  assert(target.HitNucIsSet());

  fCurrRemovalEnergy = 0;
  fCurrMomentum.SetXYZ(-1,-1,-1);


  //-- set fermi momentum vector
  //
  TH1D * prob = this->ProbDistro(target,hitNucleonRadius);
  if(!prob) {
    LOG("LocalFGM", pNOTICE)
              << "Null nucleon momentum probability distribution";
    exit(1);
  }
  double p = prob->GetRandom();
  delete prob;
  LOG("LocalFGM", pINFO) << "|p,nucleon| = " << p;

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px, py, pz;  

  int nucleon_pdgc = target.HitNucPdg();
  assert(pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc));
  int A = target.A();

  assert(target.HitNucIsSet());
  bool is_p = pdg::IsProton(nucleon_pdgc);
  double numNuc = (is_p) ? (double)target.Z():(double)target.N();

  // Calculate Fermi Momentum using Local FG equations
  double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi;
  double KF= TMath::Power(3*kPi2*numNuc*genie::utils::nuclear::Density(hitNucleonRadius,A),
			    1.0/3.0) *hbarc;



  if(p < KF){
        px = p*sintheta*cosfi;
        py = p*sintheta*sinfi;
        pz = p*costheta;
  } else {
        

        px = p*sintheta*cosfi + .5*fCOMCurrMomentum.X();
        py = p*sintheta*sinfi + .5*fCOMCurrMomentum.Y();
        pz = p*costheta + .5*fCOMCurrMomentum.Z();



  }



  fCurrMomentum.SetXYZ(px,py,pz);
  fFermiMomentum = KF;
  fRecoilPDG = recoil_pdg;


  //-- set removal energy
  //
  int Z = target.Z();
  map<int,double>::const_iterator it = fNucRmvE.find(Z);
  if(it != fNucRmvE.end()) fCurrRemovalEnergy = it->second;
  else fCurrRemovalEnergy = nuclear::BindEnergyPerNucleon(target);


  return true;
}
//____________________________________________________________________________
double LocalFGM::Prob(double p, double w, const Target & target,
			     double hitNucleonRadius) const
{
  if(w<0) {
    TH1D * prob = this->ProbDistro(target, hitNucleonRadius);
    int bin = prob->FindBin(p);
    double y  = prob->GetBinContent(bin);
    double dx = prob->GetBinWidth(bin);
    double pr  = y*dx;
    delete prob;
    return pr;
  }
  return 1;
}
//____________________________________________________________________________
// *** The TH1D object must be deleted after it is used ***
TH1D * LocalFGM::ProbDistro(const Target & target, double r) const
{
  LOG("LocalFGM", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString()
	     << ", Nucleon Radius = " << r;
  LOG("LocalFGM", pNOTICE)
             << ", P(max) = " << fPMax;

  assert(target.HitNucIsSet());

  //-- get information for the nuclear target
  int nucleon_pdgc = target.HitNucPdg();
  assert(pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc));

  // bool is_p = pdg::IsProton(nucleon_pdgc);
  // double numNuc = (is_p) ? (double)target.Z():(double)target.N();

  // Calculate Fermi Momentum using Local FG equations
  double KF = LocalFermiMomentum( target, nucleon_pdgc, r ) ; 


  double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi;
  

  LOG("LocalFGM",pNOTICE) << "KF = " << KF;

  double a  = 2.0;
  double C  = 4. * kPi * TMath::Power(KF,3) / 3.;

  // Do not include nucleon correlation tail
  //double R  = 1. / (1.- KF/4.);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("LocalFGM", pDEBUG) << "a  = " << a;
  LOG("LocalFGM", pDEBUG) << "C  = " << C;
  //LOG("LocalFGM", pDEBUG) << "R  = " << R;
#endif

  //-- create the probability distribution

  int npbins = (int) (1000*fPMax);
  TH1D *prob = new TH1D("", "", npbins, 0, fPMax);

  //apapadop

  fCOMCurrMomentum.SetXYZ(-1,-1,-1);
  TF1 COMGaus("COMGaus","gaus",-.5,.5);
  COMGaus.SetParameters(1,0,0.15);
  fCOMCurrMomentum.SetXYZ(COMGaus.GetRandom(),COMGaus.GetRandom(),COMGaus.GetRandom());
 
 // fCOMCurrMomentum.SetXYZ(0,0,0);



  double CFpMin = KF; // [GeV/c]                                                                                                                      
  double CFpMax = .7; // [GeV/c]


  TF1 RL_func("RL_func","[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4" /*+ [5]*x^5 + [6]*x^6 + [7]*x^7 + [8]*x^8"*/,CFpMin,CFpMax);                                                   
  RL_func.SetParameter(0,213.143);                                                                                                                   
  RL_func.SetParameter(1,-1581.17);                                                                                                                    
  RL_func.SetParameter(2,4398.29);                                                                                                                   
  RL_func.SetParameter(3,-5403.76); 
  RL_func.SetParameter(4,2467.83);
  


  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double fi = 2 * kPi * rnd->RndGen().Rndm();

  double RLMag = 0;
  RLMag = RL_func.GetRandom();
	    
  TVector3 RLMotion(-1,-1,-1);
  RLMotion.SetMagThetaPhi(RLMag,TMath::ACos(costheta),fi);

  //

  double dp = fPMax / (npbins-1);
//  double iC = (C>0) ? 1./C : 0.; // unused variables
//  double kfa_pi_2 = TMath::Power(KF*a/kPi,2); // unused variables

  bool hit_nuc_p = pdg::IsProton(target.HitNucPdg());

  double contact_sum = fC12_Cpn0 + fC12_Cpn1 + 2*fC12_Cpp0;
  double pair_prob = rnd->RndGen().Rndm() * contact_sum;
    
  if(pair_prob < fC12_Cpn0 + fC12_Cpn1) {
    if(hit_nuc_p) recoil_pdg = kPdgNeutron;
    if(!hit_nuc_p) recoil_pdg = kPdgProton;
  }
  else {
    if(hit_nuc_p) recoil_pdg = kPdgProton;
    if(!hit_nuc_p) recoil_pdg = kPdgNeutron;
  }


  vector<double> AV18_univ_function;
  if(pair_prob < fC12_Cpn1) AV18_univ_function = fAV18_pn1;
  else if(pair_prob < fC12_Cpn1 + fC12_Cpn0) AV18_univ_function = fAV18_pn0;
  else AV18_univ_function = fAV18_pp0;

 // KF = .25;
  for(int i = 0; i < npbins; i++) {
     double p  = i * dp;
     double p2 = TMath::Power(p,2);


     // use expression with fSRC_Fraction to allow the possibility of 
     // using the Correlated Fermi Gas Model with a high momentum tail

     // calculate |phi(p)|^2

     

     double phi2 = -1;
        if (p <= KF){

            phi2 = (1./(4*kPi)) * (3/TMath::Power(KF,3.)) * ( 1 - fSRC_Fraction );


        }else if( p > KF && p < fPCutOff ){            

          
            double bin = p / hbarc / 0.1;
        
            if (bin < 0.)
              phi2 = 0.;
            if (bin < 1.)
              phi2 = bin * AV18_univ_function[0];
            if (bin > 100.)
              phi2 = 0.;
                
            int b = bin;
            double x = bin - b;
            if(phi2 == -1){
              phi2 = x*AV18_univ_function[b] + (1.-x)*AV18_univ_function[b-1];
            }   
            
            phi2 *= (1./(4*kPi));

        
            p2 = p*p; 


        } else if (p >=fPCutOff){phi2 = 0;}        

     // calculate probability density : dProbability/dp
     double dP_dp = 4*kPi * p2 * phi2;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("LocalFGM", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
#endif
     if(p < 0 || dP_dp < 0){
       std::cout << "p < 0: " << p << " dp_dp < 0: "<< dP_dp << std::endl;
     }
     prob->Fill(p, dP_dp);
  }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  return prob;
}
//____________________________________________________________________________
double LocalFGM::LocalFermiMomentum( const Target & t, int nucleon_pdg, double radius ) const {

  assert(pdg::IsProton(nucleon_pdg) || pdg::IsNeutron(nucleon_pdg)) ;

  bool is_p = pdg::IsProton(nucleon_pdg);
  double numNuc = (double) ( (is_p) ? t.Z() : t.N() );

  //  double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi;
  
  double kF = TMath::Power( 3*kPi2*numNuc*genie::utils::nuclear::Density( radius, t.A() ),
			   1.0/3.0 ) 
    / genie::units::fermi ;

  return kF ;
}
//____________________________________________________________________________
void LocalFGM::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void LocalFGM::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void LocalFGM::LoadConfig(void)
{
  

  this->GetParamDef("LFG-MomentumMax", fPMax, 1.0);
  assert(fPMax > 0);

  this->GetParamDef("SRC-Fraction", fSRC_Fraction, 0.0);

  this->GetParamVect("AV18_pn1", fAV18_pn1);
  this->GetParamVect("AV18_pn0", fAV18_pn0);
  this->GetParamVect("AV18_pp0", fAV18_pp0);

  this->GetParam("C12_Cpn1", fC12_Cpn1);
  this->GetParam("C12_Cpn0", fC12_Cpn0);
  this->GetParam("C12_Cpp0", fC12_Cpp0);

  this->GetParam("LFG-MomentumCutOff", fPCutOff);


  if (fPCutOff > fPMax) {
        LOG("LocalFGM", pFATAL) << "Momentum CutOff greater than Momentum Max";
        exit(78);
  }


  // Load removal energy for specific nuclei from either the algorithm's
  // configuration file or the UserPhysicsOptions file.
  // If none is used use Wapstra's semi-empirical formula.
  //

  for(int Z=1; Z<140; Z++) {
    for(int A=Z; A<3*Z; A++) {
      ostringstream key ;
      int pdgc = pdg::IonPdgCode(A,Z);
      key << "RFG-NucRemovalE@Pdg=" << pdgc;
      RgKey rgkey = key.str();
      double eb ;
      if ( GetParam( rgkey, eb, false ) ) {
    	eb = TMath::Max(eb, 0.);
        LOG("LocalFGM", pINFO) << "Nucleus: " << pdgc << " -> using Eb =  " << eb << " GeV";
        fNucRmvE.insert(map<int,double>::value_type(Z,eb));
      }
    }
  }
}
//____________________________________________________________________________
