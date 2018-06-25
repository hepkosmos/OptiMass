#ifndef ALM_BASE_SRC_PROCESSTREE_H_
#define ALM_BASE_SRC_PROCESSTREE_H_

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "MathUtils.h"
#include "Minuit2/MnUserParameters.h"
#include "Types.h"

namespace OptiMass {
// ======================================================================================
// class ParticleNode: represents a particle node in a process tree.
// ======================================================================================
class ParticleNode {
private:
    std::string label;  // label for this node
    std::vector<ParticleNode *>
        siblingPtls;  // vector of sibling node pointers.
    bool debug;

public:
    // Construct a particle node with label taken
    explicit ParticleNode(const std::string &labelInput,
                          bool debugInput = false);
    // Destruct this node and its sibling nodes
    ~ParticleNode();

    // Add a sibling node to this particle node
    void AddSibling(ParticleNode *ptlNode);
    // Get final state of this particle state
    std::vector<std::string> getFinalStates();
    // Get particle label
    std::string GetLabel();
    // Get string which represents tree structure from this node
    std::string toString();
};
// ======================================================================================
// class ProcessTree: reads card and returns vectors required by main calculator
// ======================================================================================
class ProcessTree {
private:
    std::vector<std::string> processCmd;
    std::vector<ParticleNode *> processChain;
    std::vector<ParticleNode *> vec_state_node_;

    std::vector<std::string> vec_ptl_final_state_;
    std::vector<std::string> vec_ptl_invisible_state_;
    unsigned int num_ptl_invisible_;

    std::map<std::string, double> parMomentum;

    PtlIndexDict map_ptl_component_index_, map_ptl_component_effective_index_;
    PtlMomentumDict map_ptl_component_momentum_,
        map_ptl_component_effective_momentum_;

    std::unordered_map<std::string, double> map_ptl_mass_;

    // vectors containing component name vector of each particle
    std::vector<std::vector<std::string> > vec_ptl_visible_param_name_,
        vec_ptl_invisible_param_name_, vec_ptl_optimize_param_name_;
    // vectors containing particle name in given category
    std::vector<std::string> vec_ptl_visible_label_, vec_ptl_invisible_label_,
        vec_ptl_optimize_label_, vec_ptl_optimize_not_final_label_;
    // vectors containing momentum of particle in given cateogory, index
    // informations are in map_ptl_component_index_
    std::vector<TLorentzVector> vec_ptl_visible_momenta_,
        vec_ptl_invisible_momenta_, vec_ptl_optimize_momenta_;

    std::vector<std::vector<std::string> > vec_ptl_invisible_subsystem_label_;
    std::vector<bool> vec_ptl_invisible_subsystem_calc_auto_;
    std::vector<std::vector<int> > vec_ptl_invisible_subsystem_index_;
    std::vector<TVector2> vec_ptl_invisible_subsystem_momenta_;
    std::vector<TVector2> vec_ptl_invisible_subsystem_momenta_eff_;
    // buffer objects for evaluation routine
    std::vector<TVector2> buf_vec_ptl_last_momenta_;

    bool debug;

public:
    explicit ProcessTree(bool debugInput = false);
    ~ProcessTree();

    // ============================================
    // Input functions
    // ============================================
    inline void AddProcess(const std::string &cmd);
    inline void AddInvisibleSubsystem(
        const std::vector<std::string> &vec_label);
    inline void SetFinalStatesVisible(std::vector<std::string> &vec_labels);
    inline void SetFinalStatesInvisible(std::vector<std::string> &vec_labels);
    inline void SetInvisibles(std::vector<std::string> &vec_labels);
    inline void SetPtlOptimize(std::vector<std::string> &vec_labels);

    // ============================================
    // Input Parsers
    // ============================================
    void ParseProcess();

private:
    // Process interpreting parsers
    ParticleNode *singleProcessParser(const std::string &arg);
    ParticleNode *processParser(const std::string &arg);

public:
    // ============================================
    // Container Initializers
    // ============================================
    void GeneratePtlIndexDict();
    void GenerateEffectiveDict();
    void GenerateMomentumContainers();
    void GenerateMomentumDict();
    void InitializeVisibleMomenta();
    inline std::vector<std::string> getFinalStatesVisible();
    inline std::vector<std::string> getFinalStatesInvisible();

    ROOT::Minuit2::MnUserParameters GetMnUserParameters();

    void CalcMissingET();
    void CalcMissingETEffective();
    void UpdateInvisibleMomenta(const std::vector<double> &par);
    void RefreshMomentum();

    inline void SetMass(const std::string &label, double mass);
    inline void SetVisibleMomentumComponent(const std::string &str, double val);
    inline void SetMissingET(double val_x, double val_y);
    inline void SetInvisibleSubsystemMomenta(int index, double val_x,
                                             double val_y);
    inline void SetInvisibleSubsystemMomentaAutoCalc(int index, bool chk);

    // ============================================
    // Output retriving functions
    // ============================================
    double GetEffectiveScale();
    inline unsigned int GetNumInvisibles();
    // Invaraint mass of given labeled subsystem
    inline double GetSubsystemMass(const std::string &ptl_name);
    // Invaraint mass of given labeled subsystem
    inline double GetSubsystemMassSquare(const std::string &ptl_name);
    // Px of given labeled subsystem
    inline double GetSubsystemPx(const std::string &ptl_name);
    // Py of given labeled subsystem
    inline double GetSubsystemPy(const std::string &ptl_name);
    // Pz of given labeled subsystem
    inline double GetSubsystemPz(const std::string &ptl_name);
    // E of given labeled subsystem
    inline double GetSubsystemE(const std::string &ptl_name);
    // transverse momentum of given labeled subsystem
    inline double GetSubsystemPt(const std::string &ptl_name);
    // momentum of given labeled subsystem
    inline TLorentzVector GetSubsystemMomentum(const std::string &ptl_name);
    inline ParticleSubelements<TLorentzVector *> &GetSubsystemSubelements(
        const std::string &ptl_name);
    inline ParticleSubelements<TLorentzVector *>
        &GetSubsystemSubelementsEffective(const std::string &ptl_name);

    template <typename T>
    inline T GetSubsystemObservable(
        const std::string &ptl_name,
        T (*ftn)(const std::vector<TLorentzVector *> &,
                 const std::vector<TLorentzVector *> &));

    // Maintenance ftn
    void PrintDict();
};

// ======================================================================================
// class Utility functions for ParticleNode
// ======================================================================================
bool ptlIn(ParticleSubelements<int> vec_visible_index,
           ParticleSubelements<int> vec_invisible_index);
ParticleSubelements<int> ptlRemove(ParticleSubelements<int> vec_visible_index,
                                   ParticleSubelements<int> vec_invisible_index,
                                   int ptl_invisible_index);
void printPtlDict(PtlIndexDict dict);

}  // end namespace OptiMass

inline void OptiMass::ProcessTree::AddProcess(const std::string &cmd) {
    processCmd.push_back(cmd);
}
inline void OptiMass::ProcessTree::AddInvisibleSubsystem(
    const std::vector<std::string> &vec_label) {
    vec_ptl_invisible_subsystem_label_.push_back(vec_label);
}

inline void OptiMass::ProcessTree::SetFinalStatesVisible(
    std::vector<std::string> &vec_labels) {
    vec_ptl_visible_label_ = vec_labels;
}
inline void OptiMass::ProcessTree::SetFinalStatesInvisible(
    std::vector<std::string> &vec_labels) {
    vec_ptl_invisible_label_ = vec_labels;
}
inline void OptiMass::ProcessTree::SetInvisibles(
    std::vector<std::string> &vec_labels) {
    num_ptl_invisible_ = vec_labels.size();
    vec_ptl_invisible_state_ = vec_labels;
}
inline void OptiMass::ProcessTree::SetPtlOptimize(
    std::vector<std::string> &vec_labels) {
    vec_ptl_optimize_label_ = vec_labels;
}
inline void OptiMass::ProcessTree::SetMass(const std::string &label,
                                           double mass) {
    map_ptl_mass_[label] = mass;
}

// Set momentum inputs
inline void OptiMass::ProcessTree::SetVisibleMomentumComponent(
    const std::string &str, double val) {
    parMomentum[str] = val;
}
inline void OptiMass::ProcessTree::SetMissingET(double val_x, double val_y) {
    vec_ptl_invisible_subsystem_momenta_[0].Set(val_x, val_y);
}
inline void OptiMass::ProcessTree::SetInvisibleSubsystemMomenta(int index,
                                                                double val_x,
                                                                double val_y) {
    vec_ptl_invisible_subsystem_momenta_[index].Set(val_x, val_y);
}
inline void OptiMass::ProcessTree::SetInvisibleSubsystemMomentaAutoCalc(
    int index, bool chk) {
    vec_ptl_invisible_subsystem_calc_auto_[index] = chk;
}

inline std::vector<std::string> OptiMass::ProcessTree::getFinalStatesVisible() {
    return vec_ptl_visible_label_;
}
inline std::vector<std::string>
OptiMass::ProcessTree::getFinalStatesInvisible() {
    return vec_ptl_invisible_label_;
}

inline unsigned int OptiMass::ProcessTree::GetNumInvisibles() {
    return num_ptl_invisible_;
}

// Invaraint mass of given labeled subsystem
inline double OptiMass::ProcessTree::GetSubsystemMass(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcMassFromIndice);
}
// Invaraint mass of given labeled subsystem
inline double OptiMass::ProcessTree::GetSubsystemMassSquare(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcMassSquareFromIndice);
}
// Px of given labeled subsystem
inline double OptiMass::ProcessTree::GetSubsystemPx(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcPxFromIndice);
}
// Py of given labeled subsystem
inline double OptiMass::ProcessTree::GetSubsystemPy(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcPyFromIndice);
}
// Pz of given labeled subsystem
inline double OptiMass::ProcessTree::GetSubsystemPz(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcPzFromIndice);
}
// E of given labeled subsystem
inline double OptiMass::ProcessTree::GetSubsystemE(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcEFromIndice);
}
// transverse momentum of given labeled subsystem
inline double OptiMass::ProcessTree::GetSubsystemPt(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcPtFromIndice);
}
// momentum of given labeled subsystem
inline TLorentzVector OptiMass::ProcessTree::GetSubsystemMomentum(
    const std::string &ptl_name) {
    return GetSubsystemObservable(ptl_name, &CalcMomentumFromIndice);
}
// momentum of given labeled subsystem
inline OptiMass::ParticleSubelements<TLorentzVector *> &
OptiMass::ProcessTree::GetSubsystemSubelements(const std::string &ptl_name) {
    return map_ptl_component_momentum_.at(ptl_name);
}
// momentum of given labeled subsystem
inline OptiMass::ParticleSubelements<TLorentzVector *>
    &OptiMass::ProcessTree::GetSubsystemSubelementsEffective(
        const std::string &ptl_name) {
    return map_ptl_component_effective_momentum_.at(ptl_name);
}

template <typename T>
inline T OptiMass::ProcessTree::GetSubsystemObservable(
    const std::string &ptl_name,
    T (*ftn)(const std::vector<TLorentzVector *> &,
             const std::vector<TLorentzVector *> &)) {
    auto &ptl_component_effective_momentum =
        map_ptl_component_effective_momentum_.at(ptl_name);
    if (ptl_component_effective_momentum.invisible.size() == 0) {
        // cout << "using optimized system" << endl;
        return (*ftn)(ptl_component_effective_momentum.visible,
                      ptl_component_effective_momentum.optimize);
    } else {
        auto &ptl_component_momentum = map_ptl_component_momentum_.at(ptl_name);
        // cout << "using reconstructed system" << endl;
        return (*ftn)(ptl_component_momentum.visible,
                      ptl_component_momentum.invisible);
    }
}

#endif  // ALM_BASE_SRC_PROCESSTREE_H_
