#ifndef _JPsiPolarization_PolarizationModel_h_
#define _JPsiPolarization_PolarizationModel_h_
#include "RooPolarizationPdf.h"

#include <vector>
#include <string>
#include <map>

class RooRealVar;
class RooGenericPdf;
class TDirectory;

namespace JPsiPolarization {
  class PolarizationModel {
  public:
    PolarizationModel(const std::string componentname,
		      const std::string& framename);
    ~PolarizationModel();
    
    void loadParameters(TDirectory&);
    void saveParameters(TDirectory&);

    void fix(const std::string& s = "");
    void unfix(const std::string& s = "");
    void setVal(const std::string&,double newval);
    void setErr(const std::string&,double newerr);

    void initModel(RooRealVar& costh,
		   RooRealVar& phi,
		   RooAbsReal& map_uniform,
		   RooAbsReal& map_theta,
		   RooAbsReal& map_phi,
		   RooAbsReal& map_thetaphi);
		   
//    RooGenericPdf* model() { return model_; }
    RooPolarizationPdf* model() { return model_; }

    RooRealVar* lambdathetaCS() { return vars_["promptlambda_theta_CS"];}
    RooRealVar* lambdaphiCS() { return vars_["promptlambda_phi_CS"];}
    RooRealVar* lambdathetaphiCS() { return vars_["promptlambda_thetaphi_CS"];}

    RooRealVar* lambdathetaHX() { return vars_["promptlambda_theta_HX"];}
    RooRealVar* lambdaphiHX() { return vars_["promptlambda_phi_HX"];}
    RooRealVar* lambdathetaphiHX() { return vars_["promptlambda_thetaphi_HX"];}

    RooRealVar* promptNorm() { return vars_["nPrompt"]; }

  private:
    PolarizationModel(const PolarizationModel& p) {}
    PolarizationModel& operator=(const PolarizationModel& p) {}
        
    std::string compName_,fName_;

    std::map<std::string,RooRealVar*> vars_;    
    
    bool varsLoaded;

//    RooGenericPdf* model_;
    RooPolarizationPdf* model_;

  };
}

#endif
