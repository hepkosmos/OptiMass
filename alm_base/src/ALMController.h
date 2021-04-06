#ifndef ALM_BASE_SRC_ALMCONTROLLER_H_
#define ALM_BASE_SRC_ALMCONTROLLER_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "ConstraintBase.h"

namespace OptiMass {

class ALMController {
private:
    // is ALMController initialized?
    bool initialized_;
    // ALM control params container
    std::unordered_map<std::string, double> map_alm_control_param_;

    unsigned int num_constraints_;
    std::vector<double> vec_constraints_;
    std::vector<bool> vec_constraints_using_;

    // ALM Lagrange Multiplier
    std::vector<double> vec_alm_lagrange_multiplier_,
        vec_alm_lagrange_multiplier_init_;
    // ALM Penalty Param
    double alm_penalty_;
    double alm_penalty_init_;
    double mod_c_;

    // ALM Phase Counter;
    unsigned int n_iter_, n_alm_phase_1_, n_alm_phase_2_;
    // bool for run ALM or not
    bool need_alm_;

    ConstraintBase *constraint_base_;

public:
    ALMController()
        : initialized_(false),
          num_constraints_(0),
          alm_penalty_(0.1),
          alm_penalty_init_(0.1),
          mod_c_(0.),
          n_iter_(0),
          n_alm_phase_1_(0),
          n_alm_phase_2_(0),
          need_alm_(false) {}
	~ALMController() {} //rev

    //
    void SetMinimizer(ConstraintBase *constraint_base);
    //inline void UseConstraint(int index, bool use);//rev:inline
    void UseConstraint(int index, bool use);
    void CalcConstraints();

    // Set Ininitial state of ALM parameters
    void InitContainers();
    void Init();

    // Calc lagrange multiplier and penalty term. returns sqrt(sum of const
    // squared)
    inline bool NeedALM();//rev:inline
    //bool NeedALM();
    // Increase itertion number
    inline void Next() { ++n_iter_; }//rev:inline
    //void Next() { ++n_iter_; }
    // Calculate ALM
    double CalcALM(double &lagrange_multiplier, double &penalty_term);
    // Check ALM convergence test
    bool CheckConverged(double mod_c);

    // Set ALM param functions //rev:inline
    inline void SetNumberConstraints(unsigned int num);
    inline void SetInitialLagrangeMultiplier(int index, double val);
    inline void SetInitialPenaltyParam(double val);
    inline void SetALMControlParam(const std::string &str, double val);
	inline double GetALMControlParam(const std::string &str); //rev:new
    inline void SetALMIterMax(int nIter);

    inline unsigned int GetALMIterMax();
    inline unsigned int GetNumberIteration();
    inline unsigned int GetNumberPhase1();
    inline unsigned int GetNumberPhase2();
    inline double GetSumSquaredConstraints();
//    void SetNumberConstraints(unsigned int num);
//    void SetInitialLagrangeMultiplier(int index, double val);
//    void SetInitialPenaltyParam(double val);
//    void SetALMControlParam(const std::string &str, double val);
//    void SetALMIterMax(int nIter);
//
//    unsigned int GetALMIterMax();
//    unsigned int GetNumberIteration();
//    unsigned int GetNumberPhase1();
//    unsigned int GetNumberPhase2();
//    double GetSumSquaredConstraints();
};


//rev:inline
inline bool ALMController::NeedALM() { return need_alm_; }
//bool ALMController::NeedALM() { return need_alm_; }

inline unsigned int ALMController::GetALMIterMax() {
//unsigned int ALMController::GetALMIterMax() {
    return (unsigned int)map_alm_control_param_["nIterMax"];
}
inline unsigned int ALMController::GetNumberIteration() { return n_iter_; }
//unsigned int ALMController::GetNumberIteration() { return n_iter_; }
inline unsigned int ALMController::GetNumberPhase1() { return n_alm_phase_1_; }
//unsigned int ALMController::GetNumberPhase1() { return n_alm_phase_1_; }
inline unsigned int ALMController::GetNumberPhase2() { return n_alm_phase_2_; }
//unsigned int ALMController::GetNumberPhase2() { return n_alm_phase_2_; }
inline double ALMController::GetSumSquaredConstraints() { return mod_c_; }
//double ALMController::GetSumSquaredConstraints() { return mod_c_; }

// Set ALM param functions
inline void ALMController::SetNumberConstraints(unsigned int num) { 
//void ALMController::SetNumberConstraints(unsigned int num) {
    num_constraints_ = num;
    initialized_ = false;
}

//rev
inline void ALMController::SetInitialLagrangeMultiplier(int index, double val) { 
//void ALMController::SetInitialLagrangeMultiplier(int index, double val) {
    vec_alm_lagrange_multiplier_init_[index] = val;
}
inline void ALMController::SetInitialPenaltyParam(double val) {
//void ALMController::SetInitialPenaltyParam(double val) {
    alm_penalty_init_ = val;
}

//rev:remove inline
inline void ALMController::SetALMControlParam(const std::string &str,
//void ALMController::SetALMControlParam(const std::string &str,
                                              double val) {
    if ((str.compare("tau_mu") == 0) || (str.compare("eta_s") == 0) ||
        (str.compare("eta_ratio") == 0) || (str.compare("gamma") == 0) ||
        (str.compare("b_eta0") == 0) || (str.compare("b_eta") == 0) ||
        (str.compare("nIterMax") == 0)) {
        map_alm_control_param_[str] = val;
    } else {
        std::cout << " (!) Invalid ALM control parameter name : " << str << '\n';
    }
}

//rev:new
inline double ALMController::GetALMControlParam(const std::string &str) {
    if ((str.compare("tau_mu") == 0) || (str.compare("eta_s") == 0) ||
        (str.compare("eta_ratio") == 0) || (str.compare("gamma") == 0) ||
        (str.compare("b_eta0") == 0) || (str.compare("b_eta") == 0) ||
        (str.compare("nIterMax") == 0)) {
        return map_alm_control_param_[str];
    } else {
        std::cout << " (!) Invalid ALM control parameter name : " << str << '\n';
		return 1;
    }
}

inline void ALMController::SetALMIterMax(int nIter) {
//void ALMController::SetALMIterMax(int nIter) {
    map_alm_control_param_["nIterMax"] = nIter;
}

}  // end namespace OptiMass

#endif  // ALM_BASE_SRC_ALMCONTROLLER_H_
