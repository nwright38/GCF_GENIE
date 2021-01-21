//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 
*/
//____________________________________________________________________________

#include "Physics/QuasiElastic/EventGen/QELPrimaryLeptonGenerator.h"

#include <TVector3.h>
#include <TLorentzVector.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/Common/PrimaryLeptonGenerator.h"
#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Numerical/MathUtils.h"

using namespace genie;

using namespace genie::constants;


//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::~QELPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void QELPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in QEL events


  //return(PrimaryLeptonGenerator::ProcessEventRecord(evrec));


  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB
  double pnuc_beta = pnuc4.Beta();
/*  if(pnuc_beta < 1.){ 
    // no modification is required to the std implementation
    return(PrimaryLeptonGenerator::ProcessEventRecord(evrec));
  }
 */ 


  TLorentzVector * p4v = init_state.GetProbeP4(kRfLab); // v 4p @ Nucleon rest frame
  TVector3 beta = (*p4v + pnuc4).BoostVector();

  double s = (pnuc4 + *p4v).M2();
  

  p4v->Boost(-beta);  


  // Look-up selected kinematics & other needed kinematical params
  double Q2  = interaction->Kine().Q2(true);
  double Ev  = p4v->E();
  double ml  = interaction->FSPrimLepton()->Mass();
  double ml2 = TMath::Power(ml,2);

  std::cout << "Q2: " << Q2 << std::endl;
  std::cout << "ml2: " << ml2 << std::endl;

  double W = interaction->Kine().W(true);
  double W2 = W*W;
  double mv = p4v->M();
  double mv2 = mv*mv;

  std::cout << "W: " << W << std::endl;
  std::cout << "mv2: " << mv2 << std::endl;
  std::cout << "s: " << s << std::endl;


  LOG("LeptonicVertex", pNOTICE)
             << "Ev = " << Ev << ", Q2 = " << Q2;

  // Compute the final state primary lepton energy and momentum components
  // in COM frame
  double El  = (s + ml2 - W2)/(2*sqrt(s));
  double plMag = sqrt(El*El - ml2);
  double cthl = (2*Ev*El - mv2 - ml2 - Q2)/(2*p4v->Vect().Mag()*plMag);

  std::cout << "El: " << El << std::endl;
  std::cout << "pl Mag: " << plMag << std::endl;
  std::cout << "cthl: " << cthl << std::endl;

  // Randomize phi components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();

  std::cout << "phi: " << phi << std::endl;



//Just in case
/*  double plp = El - 0.5*(Q2+ml2)/Ev;                          // p(//)
  double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)

  LOG("LeptonicVertex", pNOTICE)
        << "fsl: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2*kPi * rnd->RndLep().Rndm();
  double pltx = plt * TMath::Cos(phi);
  double plty = plt * TMath::Sin(phi);

  */

  // Take a unit vector along the neutrino direction @ the nucleon rest frame
  TVector3 unit_nudir = p4v->Vect().Unit();

  // Rotate lepton momentum vector from the reference frame (x'y'z') where
  // {z':(neutrino direction), z'x':(theta plane)} to the nucleon rest frame
  //TVector3 p3l(pltx,plty,plp);


  TVector3 p3l(-1,-1,-1);  
  p3l.SetMagThetaPhi(plMag,TMath::ACos(cthl),phi);
  p3l.RotateUz(unit_nudir);


  // Lepton 4-momentum in the nucleon rest frame
  TLorentzVector p4l(p3l,El);

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ NRF: " << utils::print::P4AsString(&p4l);

  std::cout << "p4l x BEFORE boost: " << p4l.X() << std::endl;
  std::cout << "p4l y BEFORE boost: " << p4l.Y() << std::endl;
  std::cout << "p4l z BEFORE boost: " << p4l.Z() << std::endl;
  std::cout << "p4l E BEFORE boost: " << p4l.E() << std::endl;
  std::cout << "p4l M BEFORE boost: " << p4l.M() << std::endl;

  // Boost final state primary lepton to the lab frame
  p4l.Boost(beta); // active Lorentz transform

  LOG("LeptonicVertex", pNOTICE)
       << "fsl @ LAB: " << utils::print::P4AsString(&p4l);

  std::cout << "p4l x AFTER boost: " << p4l.X() << std::endl;
  std::cout << "p4l y AFTER boost: " << p4l.Y() << std::endl;
  std::cout << "p4l z AFTER boost: " << p4l.Z() << std::endl;
  std::cout << "p4l E AFTER boost: " << p4l.E() << std::endl;
  std::cout << "p4l M AFTER boost: " << p4l.M() << std::endl;

  // Figure out the Final State Lepton PDG Code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Create a GHepParticle and add it to the event record
  this->AddToEventRecord(evrec, pdgc, p4l);

  // Set final state lepton polarization
  this->SetPolarization(evrec);

  delete p4v;



  //Get kinematic variables from interaction kinematics

  //Solve for pl and pT 

  //Randomize transverse component

  //

}//___________________________________________________________________________
