//____________________________________________________________________________
/*!

\class    genie::LocalFGM

\brief    local Fermi gas model. Implements the NuclearModelI 
          interface.

\ref      

\author   Joe Johnston, Steven Dytman

\created  December 2015

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LOCAL_FGM_H_
#define _LOCAL_FGM_H_

#include <map>

#include <TH1D.h>

#include "Physics/NuclearState/NuclearModelI.h"

using std::map;

namespace genie {

class LocalFGM : public NuclearModelI {

public:
  LocalFGM();
  LocalFGM(string config);
  virtual ~LocalFGM();

  using NuclearModelI::GenerateNucleon;  // inherit versions not overridden here
  using NuclearModelI::Prob;

  //-- allow methods to be called with a radius
  bool   GenerateNucleon (const Target & t, double hitNucleonRadius) const;
  double Prob            (double p, double w, const Target & t,
			  double hitNucleonRadius) const;

  //-- implement the NuclearModelI interface
  bool GenerateNucleon (const Target & t) const {
    return GenerateNucleon(t,0.0);
  }
  double Prob (double p, double w, const Target & t) const {
    return Prob(p,w,t,0.0);
  }
  NuclearModel_t ModelType       (const Target &) const 
  { 
    return kNucmLocalFermiGas; 
  }

  virtual double LocalFermiMomentum( const Target & t, int nucleon_pdg, double radius ) const ;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set) ;

  mutable int recoil_pdg;

 protected:
  void   LoadConfig (void);


private:
  TH1D * ProbDistro (const Target & t, double r) const;
  vector<double> fAV18_pn1;
  vector<double> fAV18_pn0;
  vector<double> fAV18_pp0;


  double fC12_Cpn1;
  double fC12_Cpn0;
  double fC12_Cpp0;


  map<int, double> fNucRmvE;

  double fPMax;
   
  const NuclearModelI *  fNuclModel;   ///< nuclear model
  // options related to SRC pairs
  double fSRC_Fraction;
  double fPCutOff;
};

}         // genie namespace
#endif    // _LOCAL_FGM_H_

