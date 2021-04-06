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

}


OptiMass::CalcTemplate::~CalcTemplate() {
    delete mass_interface_;
}


void OptiMass::CalcTemplate::InitContainersProlog() {

	// Process Topology Syntax
    {process_syntax_}

	// Set All Invisibles
    {particle_invisibles_}

	// Set All of the Particle Masses
    {map_particle_mass_}

	// Objective Mass Function
    {objective_mass_function_}

	// Set Invisible particles to be Optimized
    {ptls_to_optimize_}

	// Set a Group of Mother Particles for each PT-conserving system
    {vec_ptl_invisible_subsystem_label_}

    // Minuit Parameters
    {map_minuit_param_}

    // ALM Control Parameters (fixed)
    {map_alm_control_param_}

	// ALM Control Parameters from User's Customization Code (CDATA)
    {map_alm_control_param_code_}

    // Set the number of constraints
    {num_constraints_}

}


void OptiMass::CalcTemplate::InitContainersEpilog() {

	// True/False Switch for Calculating the MET for each PT conservation system  
	// set True/Fasle if <Subsystem set_value="manual"/"automatic"> 
    {vec_ptl_invisible_subsystem_calc_auto_}

    // ALM lagrange multiplier and penalty initial values
    {alm_param_initial_values_}

    // Transverse Projection (use this for calculate MT2 type variable)
    {transverse_projection_}

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


void OptiMass::CalcTemplate::RunOptiMass(
		{run_optimass_args_})
{

	// Visible Momentum Input 
	{visible_momentum_input_}

	// Set MissingET of each PT-conservation system
	{met_input_}
	
	// Constraint Function Switch
	{constraint_function_switch_}

    // Call OptiMass Calculation
    this->Calc();

	// Assign OptiMass and SumSquaredConstraints
	om = sqrt(this->GetOptiMass());
    ssc = this->alm_controller_.GetSumSquaredConstraints();

}

