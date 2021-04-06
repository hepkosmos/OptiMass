#ifndef ALM_BASE_SRC_MASSFUNCTION_H_
#define ALM_BASE_SRC_MASSFUNCTION_H_

#include <vector>
#include "ALMController.h"
#include "MassFunctionInterface.h"
#include "Minuit2/FCNBase.h"
#include "ProcessTree.h"
#include "TLorentzVector.h"

namespace OptiMass {
// class FCNBase for Minuit minimizer
class MassFunction : public ROOT::Minuit2::FCNBase {
private:
    ProcessTree *process_tree_;      // Tree
    ALMController *alm_controller_;  // ALM
    MassFunctionInterface *mass_interface_;
    // mutable std::vector< double > vec_masses;
    // mutable std::vector< TVector2 > vec_ptl_last_momenta_;

public:
    explicit MassFunction(ProcessTree *process_tree)
        : process_tree_(process_tree) {}
    MassFunction() {} 
    ~MassFunction() {}

    void SetFactory(ProcessTree *process_tree, ALMController *alm_controller) {
        process_tree_ = process_tree;
        alm_controller_ = alm_controller;
    }
    void Init(MassFunctionInterface *interface) { mass_interface_ = interface; }

    double operator()(const std::vector<double> &par) const;

    double Up() const { return 1; }
};
}  // namespace OptiMass

#endif  // ALM_BASE_SRC_MASSFUNCTION_H_
