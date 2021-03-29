///____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Afroditi Papadopoulou <apapadop \at mit.edu>
         Massachusetts Institute of Technology - October 04, 2019

 @ October 4, 2019 - Afroditi Papadopoulou (AP)
   Created this new module that controls the addition of the recoil nucleon in the event record 
   and extracts its kinematics 
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/SRCNuclearRecoil.h"

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
SRCNuclearRecoil::SRCNuclearRecoil() :
SecondNucleonEmissionI("genie::SRCNuclearRecoil")
{

}
//___________________________________________________________________________
SRCNuclearRecoil::SRCNuclearRecoil(string config) :
  SecondNucleonEmissionI("genie::SRCNuclearRecoil", config )
{

}

//___________________________________________________________________________

SRCNuclearRecoil::~SRCNuclearRecoil()
{

}

//___________________________________________________________________________

void SRCNuclearRecoil::ProcessEventRecord(GHepRecord * evrec) const
{

  Interaction *  interaction = evrec       -> Summary();
  InitialState * init_state  = interaction -> InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  // do nothing for non-nuclear targets
  if(!tgt->IsNucleus()) return;

  // access the hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = evrec->HitNucleon();
  GHepParticle * nucleus = evrec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  double pN2 = TMath::Power(nucleon->P4()->Rho(),2.); // (momentum of struck nucleon)^2

  // Set this to either a proton or neutron to eject a secondary particle
  int eject_nucleon_pdg = this->SRCRecoilPDG( *nucleon, *tgt );
  // Ejection of secondary particle
  if (eject_nucleon_pdg != 0) {
	if (fGaussianEmission || fGeneralizedContactFormalismEmission) { EmitSecondNucleon(evrec,eject_nucleon_pdg); } 
        else { SecondNucleonEmissionI::EmitSecondNucleon(evrec,eject_nucleon_pdg); }
  }

}

//___________________________________________________________________________

bool SRCNuclearRecoil::EmitSecondNucleon(GHepRecord * evrec, const int eject_nucleon_pdg) const {

  LOG("SRCNuclearRecoil", pINFO) << "Adding a recoil nucleon with PDG " << eject_nucleon_pdg ;


 //  EventRecord & event = *(evrec->event);
                
  const Interaction *summary = evrec->Summary();
  const ProcessInfo &proc = summary->ProcInfo();
  const Kinematics &kine = summary->Kine();
  const InitialState & init_state = summary->InitState();


  GHepParticle * nucleon = evrec->HitNucleon();
  GHepParticle * nucleus = evrec->TargetNucleus();



  int nucleon_pdgc = nucleon->Pdg();
  bool is_hit_p  = pdg::IsProton(nucleon_pdgc);
  bool is_rec_p  = pdg::IsProton(eject_nucleon_pdg);
  int Z;
  if((is_hit_p && !is_rec_p) || (!is_hit_p && is_rec_p)){ Z = nucleus->Z()-1; }
  if(is_hit_p && is_rec_p){ Z = nucleus->Z()-2; }
  if(!is_hit_p && !is_rec_p){ Z = nucleus->Z(); }
  int A = nucleus->A() - 2;


  TParticlePDG * remn_nucleus = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  remn_nucleus = PDGLibrary::Instance()->Find(ipdgc);
  
  GHepStatus_t status = kIStHadronInTheNucleus;
  int imom = evrec->TargetNucleusPosition();

  //-- Has opposite momentum from the struck nucleon
  double vx = nucleon->Vx();
  double vy = nucleon->Vy();
  double vz = nucleon->Vz();
  // recoil nucleon not exactly in the opposite direction of the hit nucleon 
  // Use gaussian distribution for center-of-mass motion
  double px = fGaussianEmission ? gRandom->Gaus(0,fGaussianSigma) - nucleon->Px() : fNuclModel->COMMomentum3().X() - fNuclModel->Momentum3().X();
  double py = fGaussianEmission ? gRandom->Gaus(0,fGaussianSigma) - nucleon->Py() : fNuclModel->COMMomentum3().Y() - fNuclModel->Momentum3().Y();
  double pz = fGaussianEmission ? gRandom->Gaus(0,fGaussianSigma) - nucleon->Pz() : fNuclModel->COMMomentum3().Z() - fNuclModel->Momentum3().Z();

  double M  = PDGLibrary::Instance()->Find(eject_nucleon_pdg)->Mass();
  double E  = TMath::Sqrt(px*px+py*py+pz*pz+M*M);


  
  evrec->AddParticle( eject_nucleon_pdg, status, imom, -1, -1, -1, px, py, pz, E, vx, vy, vz, 0 );

  return true ;

}
//___________________________________________________________________________

int SRCNuclearRecoil::SRCRecoilPDG( const GHepParticle & nucleon, const Target & tgt) const {

      int eject_nucleon_pdg = 0;

      // const int nucleus_pdgc = tgt->Pdg();
      const int nucleon_pdgc = nucleon.Pdg();

      //double pN2 = TMath::Power(nucleon.P4()->Rho(),2.); // (momentum of struck nucleon)^2

      double kF = fNuclModel -> LocalFermiMomentum( tgt, 
						    nucleon_pdgc, 
						    nucleon.X4()->Vect().Mag() );


      
      if (fNuclModel->Momentum() > kF) {
        if (fNuclModel->RecoilPDG() != 0) eject_nucleon_pdg = fNuclModel->RecoilPDG();
      /*  else{
          double Pp = (nucleon.Pdg() == kPdgProton) ? fPPPairPercentage : fPNPairPercentage;
          RandomGen * rnd = RandomGen::Instance();
          double prob = rnd->RndGen().Rndm();
          eject_nucleon_pdg = (prob > Pp) ? kPdgNeutron : kPdgProton;
        }
        */
      }

      

      return eject_nucleon_pdg;
}

//___________________________________________________________________________
void SRCNuclearRecoil::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SRCNuclearRecoil::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SRCNuclearRecoil::LoadConfig(void)
{

  SecondNucleonEmissionI::LoadConfig() ;

  this->GetParamDef("PNPairPercentage",       fPNPairPercentage,    0.95);

  if (fPNPairPercentage < 0. || fPNPairPercentage > 1.) { 

	LOG("SRCNuclearRecoil", pFATAL)
	<< "PNPairPercentage either less than 0 or greater than 1: Exiting" ;

	exit(78); 
  }

  fPPPairPercentage = 1. - fPNPairPercentage;

  // Get the Fermi momentum table for relativistic Fermi gas
  GetParam( "FermiMomentumTable", fKFTableName ) ;
  fKFTable = 0;
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  fKFTable = kftp->GetTable(fKFTableName);
  assert(fKFTable);

  this->GetParam("SRC-GaussianEmission",fGaussianEmission);
  this->GetParamDef("SRC-GaussianSigma",fGaussianSigma,0.);

  this->GetParam("SRC-GeneralizedContactFormalism-GaussianEmission",fGeneralizedContactFormalismEmission);
  this->GetParamDef("SRC-GenaralizedContactFormalism-GaussianSigma",fGeneralizedContactFormalismGaussianSigma,0.);

}
//____________________________________________________________________________
