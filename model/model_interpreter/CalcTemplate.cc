#include "CalcTemplate.h"

#include "ProcessTree.h"
#include "MassMinimizer.h"
#include "MassFunctionInterface.h"
#include "MathUtils.h"

#include "Minuit2/VariableMetricMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/FunctionMinimum.h"

#include <iostream>

OptiMass::CalcTemplate::CalcTemplate() : MassMinimizer({debug}) {
    //std::cout << "welcome from Template" << std::endl;
}

OptiMass::CalcTemplate::~CalcTemplate() {
    //std::cout << "bye from Template" << std::endl;
    delete mass_interface_;
}

void OptiMass::CalcTemplate::InitContainersProlog() {
    {process_syntax_}

    {particle_invisibles_}

    {map_particle_mass_}

    {vec_ptl_target_label_}

    {vec_ptl_optimize_label_}

    {vec_ptl_invisible_subsystem_label_}

    // Minuit Parameters
    {map_minuit_param_}

    // ALM Control Parameters (fixed)
    {map_alm_control_param_}

    // set number of constraints
    {num_constraints_}
}
void OptiMass::CalcTemplate::InitContainersEpilog() {
    //
    {vec_ptl_invisible_subsystem_calc_auto_}
    // ALM lagrange multiplier and penalty initial values
    {alm_param_initial_values_}
    // Transverse Projection ( use this for calculate MT2 type variable
    {transverse_projection_}
    //std::cout << "configuration done." << std::endl;
}

void OptiMass::CalcTemplate::CalcProlog() {
    {calc_prolog_}
}
void OptiMass::CalcTemplate::CalcStrategy(ROOT::Minuit2::MnUserParameters& params) {
    {calc_strategy_}
}
void OptiMass::CalcTemplate::CalcEpilog() {
    {calc_epilog_}
}
    
void OptiMass::CalcTemplate::CalcConstraints(std::vector<double>& vec_constraints, std::vector<bool>& vec_constraints_using) {
    {constraints_}
}
