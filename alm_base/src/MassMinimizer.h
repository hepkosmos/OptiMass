#ifndef ALM_BASE_SRC_MASSMINIMIZER_H_
#define ALM_BASE_SRC_MASSMINIMIZER_H_

#include <map>
#include <string>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnStrategy.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/VariableMetricMinimizer.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include "ALMController.h"
#include "ConstraintBase.h"
#include "MassFunction.h"
#include "MathUtils.h"
#include "ProcessTree.h"
#include "Types.h"

namespace OptiMass {
class MassMinimizer : public ConstraintBase {
protected:
    // Minuit control parameters
    int maxfcn_;             // maximum minuit FCNbase call
    double tolerance_;       // minuit tolerance
    double init_step_size_;  // minuit err param for scanning variable

    ProcessTree process_tree_;
    ALMController alm_controller_;

    bool debug_;

    // buffer objects for evaluation routine
    TLorentzVector buf_vector_;

    // minimizer and minimizing function
    MassFunction ftnn_;

    ROOT::Minuit2::MnUserParameters user_param_init_;

    ROOT::Minuit2::VariableMetricMinimizer migrad_minimizer_;
    ROOT::Minuit2::SimplexMinimizer simplex_minimizer_;

    double output_;
    bool minimum_valid_;
    double edm_;
    // ROOT::Minuit2::FunctionMinimum function_minimum_;

    MassFunctionInterface *mass_interface_;

public:
    // ============================================
    // Constructor and Destructors
    // ============================================
    explicit MassMinimizer(bool debugInput = false)
        : maxfcn_(100),
          tolerance_(0.1),
          init_step_size_(1.),
          process_tree_(debugInput),
          debug_(debugInput),
          buf_vector_(0, 0, 0, 0),
          output_(0),
          minimum_valid_(false),
          edm_(0.) {}
    virtual ~MassMinimizer();

    // ============================================
    // Pre-configuration functions
    // ============================================
    inline void AddProcess(const std::string &syntax) {
        process_tree_.AddProcess(syntax);
    }

    // Set Minuit param functions
    inline void SetMinuitIterMax(int nIter) { maxfcn_ = nIter; }
    inline void SetMinuitTolerance(double val) { tolerance_ = val; }
    inline void SetMinuitParamErr(double val) { init_step_size_ = val; }

    // Set Minuit Fitting param configuration
    inline void SetInitInvisibleMomentum(const std::string &str, double val) {
        user_param_init_.SetValue(str, val);
    }
    inline double GetInitInvisibleMomentum(const std::string &str) {
        return user_param_init_.Value(str);
    }
    inline void Fix(const std::string &str) { user_param_init_.Fix(str); }
    inline void Release(const std::string &str) {
        user_param_init_.Release(str);
    }

    // Set momentum inputs
    inline void SetMomentumValue(const std::string &str, double val) {
        process_tree_.SetVisibleMomentumComponent(str, val);
    }
    inline void SetMissingET(double val_x, double val_y) {
        process_tree_.SetMissingET(val_x, val_y);
    }
    inline void SetInvisibleSubsystemMomenta(int index, double val_x,
                                             double val_y) {
        process_tree_.SetInvisibleSubsystemMomenta(index, val_x, val_y);
    }

    // Transverse projection
    void TransverseProjection(bool chk);

    // void ParseCard(const char* filename);
    // inline void ParseCard(std::string& str){
    //    readCard(str.c_str());
    //}

    // initialize containers
    virtual void InitContainersProlog() = 0;
    void InitContainers();
    virtual void InitContainersEpilog() = 0;
    /*
        inline void SetInvisibleSubsystemCalcAuto ( int index, bool chk ) {
            process_tree_.vec_ptl_invisible_subsystem_calc_auto_[index] = chk;
        }*/
    // ============================================
    // Main calculation routines
    // ============================================

    //  void generateMomentumDict();

    virtual void CalcProlog() = 0;
    void Calc();
    virtual void CalcStrategy(ROOT::Minuit2::MnUserParameters &) = 0;
    virtual void CalcEpilog() = 0;

    inline void UseConstraint(int index, bool use);

    // ============================================
    // Output retrival functions
    // ============================================
    // Get output mass
    inline double GetOptiMass() { return output_; }
    inline bool IsValid() { return minimum_valid_; }
    inline bool GetEdm() { return minimum_valid_; }

    // ============================================
    // Module Retirival Functions
    // ============================================
    inline ProcessTree &GetProcessTree();
    inline ALMController &GetALMController();

protected:
    inline void MinimizeCombined(ROOT::Minuit2::MnUserParameters &params);
    inline void MinimizeMigrad(ROOT::Minuit2::MnUserParameters &params);
    inline void MinimizeSimplex(ROOT::Minuit2::MnUserParameters &params);
};
}  // namespace OptiMass

inline void OptiMass::MassMinimizer::UseConstraint(int index, bool use) {
    alm_controller_.UseConstraint(index, use);
}

inline OptiMass::ProcessTree &OptiMass::MassMinimizer::GetProcessTree() {
    return process_tree_;
}

inline OptiMass::ALMController &OptiMass::MassMinimizer::GetALMController() {
    return alm_controller_;
}

inline void OptiMass::MassMinimizer::MinimizeCombined(
    ROOT::Minuit2::MnUserParameters &params) {
    // Migrad
    ROOT::Minuit2::MnStrategy strategy(2);
    ROOT::Minuit2::FunctionMinimum min_M = migrad_minimizer_.Minimize(
        ftnn_, params, strategy, maxfcn_, tolerance_);

    // Simplex -> Migrad
    ROOT::Minuit2::FunctionMinimum min_S = simplex_minimizer_.Minimize(
        ftnn_, params, strategy, maxfcn_, tolerance_);
    ROOT::Minuit2::MnUserParameters paramsBuf = min_S.UserParameters();
    for (unsigned int i = 0; i < 4 * process_tree_.GetNumInvisibles(); ++i) {
        paramsBuf.SetError(i, init_step_size_);
    }
    min_S = migrad_minimizer_.Minimize(ftnn_, paramsBuf, strategy, maxfcn_,
                                       tolerance_);

    if (min_M.Fval() <= min_S.Fval()) {
        output_ = min_M.Fval();
        params = min_M.UserParameters();
        minimum_valid_ = min_M.IsValid();
        edm_ = min_M.Edm();
    } else {
        output_ = min_S.Fval();
        params = min_S.UserParameters();
        minimum_valid_ = min_S.IsValid();
        edm_ = min_S.Edm();
    }
    /*
        if(min_M.Fval() <= min_S.Fval()){
            function_minimum_ = min_M;
        }else{
            function_minimum_ = min_S;
        }
        output = function_minimum_.Fval();
        params = function_minimum_.UserParameters();
        minimum_valid_ = function_minimum_.IsValid();
    */
}
inline void OptiMass::MassMinimizer::MinimizeMigrad(
    ROOT::Minuit2::MnUserParameters &params) {
    // Migrad
    ROOT::Minuit2::MnStrategy strategy(2);
    ROOT::Minuit2::FunctionMinimum min_M = migrad_minimizer_.Minimize(
        ftnn_, params, strategy, maxfcn_, tolerance_);

    output_ = min_M.Fval();
    params = min_M.UserParameters();
    minimum_valid_ = min_M.IsValid();
    edm_ = min_M.Edm();
}
inline void OptiMass::MassMinimizer::MinimizeSimplex(
    ROOT::Minuit2::MnUserParameters &params) {
    // Simplex
    ROOT::Minuit2::MnStrategy strategy(2);
    ROOT::Minuit2::FunctionMinimum min_S = simplex_minimizer_.Minimize(
        ftnn_, params, strategy, maxfcn_, tolerance_);

    output_ = min_S.Fval();
    params = min_S.UserParameters();
    minimum_valid_ = min_S.IsValid();
    edm_ = min_S.Edm();
}
#endif  // ALM_BASE_SRC_MASSMINIMIZER_H_
