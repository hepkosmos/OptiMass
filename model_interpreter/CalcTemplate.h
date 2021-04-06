#ifndef OptiMass_CalcTemplate__
#define OptiMass_CalcTemplate__

#include "CalcTemplate.h"
#include <vector>
#include "MassMinimizer.h"

namespace OptiMass {

class CalcTemplate : public MassMinimizer {

    public: 

    CalcTemplate();
    ~CalcTemplate();

    virtual void InitContainersProlog();
    virtual void InitContainersEpilog();

    virtual void CalcProlog();
    virtual void CalcStrategy(ROOT::Minuit2::MnUserParameters& params);
    virtual void CalcEpilog();

    virtual void CalcConstraints(std::vector<double>& vec_constraints, std::vector<bool>& vec_constraints_using);
	
	void RunOptiMass(
		{run_optimass_args_});

};

} //end namespace OptiMass

#endif // OptiMass_CalcTemplate__
