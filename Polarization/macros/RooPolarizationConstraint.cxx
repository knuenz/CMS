#include "RooPolarizationConstraint.h"
#include <string>
#include <math.h>

//ClassImp(RooPolarizationConstraint)

RooPolarizationConstraint::RooPolarizationConstraint() {
}


RooPolarizationConstraint::RooPolarizationConstraint(const char* name,
						     const char* title,
						     RooAbsReal& lambda_theta,
						     RooAbsReal& lambda_phi,
						     RooAbsReal& lambda_thetaphi):
  RooAbsPdf(name,title),
  lTheta("lTheta","Lambda Theta",this,lambda_theta),
  lPhi("lPhi","Lambda Phi",this,lambda_phi),
  lThetaPhi("lThetaPhi","Lambda Theta-Phi",this,lambda_thetaphi)
{}

RooPolarizationConstraint::RooPolarizationConstraint(const RooPolarizationConstraint& other,
						     const char* name):
  RooAbsPdf(other,name),
  lTheta("lTheta",this,other.lTheta),
  lPhi("lPhi",this,other.lPhi),
  lThetaPhi("lThetaPhi",this,other.lThetaPhi)
{}

RooPolarizationConstraint::~RooPolarizationConstraint() {
}

Int_t RooPolarizationConstraint::getAnalyticalIntegral(RooArgSet& allVars,
						       RooArgSet& analVars,
						       const char*) const {
  return 1;
}

Double_t RooPolarizationConstraint::analyticalIntegral(Int_t code,
						       const char* rN) const {
  assert(code == 1);
  return 1.0;
}

Double_t RooPolarizationConstraint::evaluate() const {
  if(fabs(lPhi) > .5*(1+lTheta)) return 1e-120;
  if(fabs(lThetaPhi) > .5*(1-lPhi)) return 1e-120;
  if(lTheta*lTheta + 2*lThetaPhi*lThetaPhi > 1.) return 1e-120;
  if(lPhi < -1./3. && (1.+2*lPhi)*(1.+2*lPhi) + 2.*lThetaPhi*lThetaPhi > 1) return 1e-120;

  return 1.0;
}


