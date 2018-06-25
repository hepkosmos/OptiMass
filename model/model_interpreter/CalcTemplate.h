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

    void InitContainersProlog();
    void InitContainersEpilog();

    void CalcProlog();
    void CalcStrategy(ROOT::Minuit2::MnUserParameters& params);
    void CalcEpilog();

    void CalcConstraints(std::vector<double>& vec_constraints, std::vector<bool>& vec_constraints_using);
};

}//end namespace optimass

#endif // OptiMass_CalcTemplate__
