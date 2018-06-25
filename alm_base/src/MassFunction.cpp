#include "MassFunction.h"
#include <vector>
#include "ALMController.h"
#include "MassFunctionInterface.h"
#include "ProcessTree.h"
#include "TLorentzVector.h"

namespace OptiMass {
// ========================================================================
// MassFunction Class
// ========================================================================
double MassFunction::operator()(const std::vector<double> &par) const {
    // Update Momentum
    process_tree_->UpdateInvisibleMomenta(par);

    // Mass Value
    // if(!(process_tree_->vec_ptl_optimize_not_final_label_.empty())) {
    process_tree_->RefreshMomentum();
    // }
    const double massFtn = mass_interface_->Calc();

    // If constraints are supplied:
    if (alm_controller_->NeedALM()) {
        // Update momentum
        // if(!(process_tree_->vec_ptl_optimize_not_final_label_.empty())) {
        process_tree_->RefreshMomentum();
        // }
        alm_controller_->CalcConstraints();

        double lagrangeMultiplier = 0.;
        double penaltyTerm = 0.;

        alm_controller_->CalcALM(lagrangeMultiplier, penaltyTerm);

        // return the sum of above three
        return massFtn + lagrangeMultiplier + penaltyTerm;
    } else {  // else constraints are not supplied:
        // return mass
        return massFtn;
    }
}
}  // namespace OptiMass
