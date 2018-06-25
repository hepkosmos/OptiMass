#include "ALMController.h"
#include <cmath>
#include "ConstraintBase.h"

namespace OptiMass {
void ALMController::SetMinimizer(ConstraintBase *constraint_base) {
    constraint_base_ = constraint_base;
}

void ALMController::CalcConstraints() {
    constraint_base_->CalcConstraints(vec_constraints_, vec_constraints_using_);
}

void ALMController::InitContainers() {
    vec_alm_lagrange_multiplier_.assign(num_constraints_, 0.);
    vec_alm_lagrange_multiplier_init_.assign(num_constraints_, 0.);
    vec_constraints_using_.assign(num_constraints_, false);
    vec_constraints_.assign(num_constraints_, 0.);
    initialized_ = true;
}

void ALMController::Init() {
    // If ALMController is first run, or number of constraints modified,
    // initialize container once.
    if (!initialized_) { InitContainers(); }

    // ALM initialization
    if (need_alm_) {
        // ALM parameter initialization
        for (unsigned int i = 0; i < num_constraints_; ++i) {
            vec_alm_lagrange_multiplier_[i] =
                vec_alm_lagrange_multiplier_init_[i];
        }
        alm_penalty_ = alm_penalty_init_;

        // ALM Convergence / Tolerance Parameter intialization
        map_alm_control_param_["alpha"] =
            fmin(alm_penalty_, map_alm_control_param_.at("gamma"));
        map_alm_control_param_["eta"] =
            map_alm_control_param_.at("eta_ratio") *
            map_alm_control_param_.at("eta_s") *
            std::pow(map_alm_control_param_.at("alpha"),
                     map_alm_control_param_.at("b_eta0"));
    }
    n_iter_ = 0;
    n_alm_phase_1_ = 0;
    n_alm_phase_2_ = 0;
}

double ALMController::CalcALM(double &lagrange_multiplier,
                              double &penalty_term) {
    double lagrangeMultiplier = 0.;
    double mod_c_square = 0.;

    // Lagrange Multiplier
    auto itParALM = vec_alm_lagrange_multiplier_.cbegin();
    auto itConstraints = vec_constraints_.cbegin();
    do {
        lagrangeMultiplier += -(*itParALM) * (*itConstraints);
        mod_c_square += std::pow(*itConstraints, 2);
        ++itParALM;
        ++itConstraints;
    } while (itConstraints != vec_constraints_.end());

    lagrange_multiplier = lagrangeMultiplier;
    penalty_term = 1. / (2. * alm_penalty_) * mod_c_square;

    mod_c_ = std::sqrt(mod_c_square);
    return mod_c_;
}

bool ALMController::CheckConverged(double mod_c) {
    //#######################################################//
    // Phase (2,3) : Augmented-Lagrange multipler driven mode//
    //#######################################################//
    // if( mod_c < fmax(alm_controller_.map_alm_control_param_.at("eta_s"),
    // map_alm_control_param_.at("eta"))) {
    if (mod_c < map_alm_control_param_.at("eta")) {
        //####################################//
        // Phase-3 ( |C_k| < eta_s ) : STOP   //
        //####################################//
        if (mod_c <= map_alm_control_param_.at("eta_s")) {
            return true;
        }  // End-if(phase-3)

        //###################################//
        // Phase-2 ( eta_s < |C_k| < eta_k ) //
        //###################################//
        else if (mod_c > map_alm_control_param_.at("eta_s")) {
            ++n_alm_phase_2_;

            for (unsigned int i = 0; i < num_constraints_; ++i) {
                vec_alm_lagrange_multiplier_[i] +=
                    -vec_constraints_[i] / alm_penalty_;
            }

            //// Update the convergence tolerance/parameters
            map_alm_control_param_["alpha"] = alm_penalty_;
            map_alm_control_param_["eta"] =
                map_alm_control_param_.at("eta") *
                std::pow(map_alm_control_param_.at("alpha"),
                         map_alm_control_param_.at("b_eta"));
        }
    }
    //######################################################//
    // Phase-1 ( |C_k| > eta_k ) : Penalty-term driven mode //
    //######################################################//
    else {
        ++n_alm_phase_1_;

        //// Update the penalty parameter in phase-1 : Decrease mu_k+1
        alm_penalty_ = map_alm_control_param_.at("tau_mu") * alm_penalty_;

        //// Update the convergence tolerance/parameters
        map_alm_control_param_["alpha"] =
            alm_penalty_ * map_alm_control_param_.at("gamma");
        map_alm_control_param_["eta"] =
            map_alm_control_param_.at("eta_ratio") *
            map_alm_control_param_.at("eta_s") *
            std::pow(map_alm_control_param_.at("alpha"),
                     map_alm_control_param_.at("b_eta"));

        // cout << "phase 1" << endl;
    }
    return false;
}
}  // namespace OptiMass
