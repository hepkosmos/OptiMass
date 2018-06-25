#include "MassMinimizer.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>

#include "TLorentzVector.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/VariableMetricMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnPrint.h"

#include "MassFunction.h"
#include "ProcessTree.h"
#include "Types.h"
#include "ConstraintBase.h"

// ========================================================================
// MassMinimizer Constructor
// ========================================================================

OptiMass::MassMinimizer::~MassMinimizer() {
    //delete minimizer_;
    if(debug_)
        std::cout << "bye from MassMinimizer" << std::endl;
}

// ========================================================================
// MassMinimizer Member functions
// ========================================================================

void OptiMass::MassMinimizer::InitContainers(){
    // Run InitContainter Prolog
    InitContainersProlog();
    // Initialize ProcessTree
    process_tree_.ParseProcess();

    // Initialize ALM Controller
    alm_controller_.InitContainers();
    alm_controller_.SetMinimizer(dynamic_cast<ConstraintBase*>(this));

    // Initialize mass function interface
    mass_interface_->Init();

    // Initialize Minuit FCNBase
    ftnn_.SetFactory(&process_tree_, &alm_controller_);
    ftnn_.Init(mass_interface_);

    // Get fitting parameters as MnUserParamters
    user_param_init_ = process_tree_.GetMnUserParameters();

    // Run InitContainter Epilog
    InitContainersEpilog();
}

void OptiMass::MassMinimizer::Calc(){
    // Refresh missing ET
    process_tree_.CalcMissingET();

    process_tree_.InitializeVisibleMomenta();

    // load preconfigured initial condition
    CalcProlog();

    // Calculating Effective Missing ET
    process_tree_.CalcMissingETEffective();

    // Initialize ALM controller
    alm_controller_.Init();

    // Entering Minuit Minimization
    ROOT::Minuit2::MnUserParameters output_params = user_param_init_;
    minimum_valid_ = false;
    //output_params.SetPrecision(1.E-7);
    
    // Mimization
    for(unsigned int i = 0 ; i < 4*process_tree_.GetNumInvisibles(); ++i){
        output_params.SetError(i,init_step_size_);
    }
    CalcStrategy(output_params);

    // Entering ALM routine
    if(alm_controller_.NeedALM()){

        //double outputMigrad, outputSimplex;
        double lagrangeMultiplier, penaltyTerm, mod_c;

        bool converged = false;

	    do{
            // Update Momentum
            process_tree_.UpdateInvisibleMomenta(output_params.Params());
            // refresh momentum
            process_tree_.RefreshMomentum();

            // Evaluate constraint formula
            alm_controller_.CalcConstraints();

            // Calculate Lagrange Multiplier
            lagrangeMultiplier = 0.;
            penaltyTerm = 0.;
            mod_c = alm_controller_.CalcALM(lagrangeMultiplier,penaltyTerm);

            output_ -= penaltyTerm;
            output_ -= lagrangeMultiplier;

            // check ALM convergence 
            converged = alm_controller_.CheckConverged(mod_c);

            if(!converged){
                alm_controller_.Next();
                //++nIterOutput;
            }
            // If it does not meet convergence criterion, do ALM minimization
            if(alm_controller_.GetNumberIteration() < alm_controller_.GetALMIterMax() && !converged){

                for(unsigned int i = 0 ; i < 4*process_tree_.GetNumInvisibles(); ++i){
                    output_params.SetError(i,init_step_size_);
                }                
                CalcStrategy(output_params);
            }
        }while(alm_controller_.GetNumberIteration() <alm_controller_.GetALMIterMax() && !converged);

        // if number of iteration saturated, update momenta of last iteration and exit.
        if( !converged ) {
            // Update Momentum
            process_tree_.UpdateInvisibleMomenta(output_params.Params());
            // refresh momentum
            process_tree_.RefreshMomentum();

            // Evaluate constraint formula
            alm_controller_.CalcConstraints();

            // Calculate Lagrange Multiplier
            lagrangeMultiplier = 0.;
            penaltyTerm = 0.;
            mod_c = alm_controller_.CalcALM(lagrangeMultiplier,penaltyTerm);

            output_ -= penaltyTerm;
            output_ -= lagrangeMultiplier;
        }

    }// end ALM

    CalcEpilog();
}

void OptiMass::MassMinimizer::TransverseProjection(bool chk){
    if(chk){
        for(unsigned int i = 0; i < process_tree_.GetNumInvisibles();++i) {
             user_param_init_.Fix(i+2);
        } 
    }else{
        for(unsigned int i = 0; i < process_tree_.GetNumInvisibles();++i) {
             user_param_init_.Release(i+2);
        } 
    }
}

