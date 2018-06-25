#include "MassFunctionInterface.h"
#include "MassFunction.h"

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>

#include "TLorentzVector.h"

#include "ProcessTree.h"

using std::vector;
using std::string;

namespace OptiMass {
//======================================================================================
// Class MassFunctionInterface
//======================================================================================
MassFunctionInterface::MassFunctionInterface(ProcessTree *process_tree)
    : process_tree_(process_tree) {}
MassFunctionInterface::~MassFunctionInterface() {}
//======================================================================================

//======================================================================================
// Class MassFunctionParticle
//======================================================================================
MassFunctionParticle::MassFunctionParticle(
    ProcessTree *process_tree, string &label,
    double (*ftn)(const vector<TLorentzVector *> &,
                  const vector<TLorentzVector *> &))
    : MassFunctionInterface(process_tree), ptl_label_(label), ftn_(ftn) {}
MassFunctionParticle::MassFunctionParticle(
    ProcessTree *process_tree, const char *label,
    double (*ftn)(const vector<TLorentzVector *> &,
                  const vector<TLorentzVector *> &))
    : MassFunctionInterface(process_tree), ptl_label_(label), ftn_(ftn) {
    // std::cout << "create label " << ptl_label_ << std::endl;
}
MassFunctionParticle::~MassFunctionParticle() {
    // std::cout << "delete label " << ptl_label_ << std::endl;
}
//======================================================================================
void MassFunctionParticle::Init() {
    auto &ptl_component_effective_momentum =
        process_tree_->GetSubsystemSubelementsEffective(ptl_label_);
    vec_visibles_ = ptl_component_effective_momentum.visible;
    vec_invisibles_ = ptl_component_effective_momentum.optimize;
}
double MassFunctionParticle::Calc() {
    // std::cout << "label " << ptl_label_ << " | " <<
    // process_tree_->GetSubsystemObservable(ptl_label_,ftn_) << std::endl;
    // return process_tree_->GetSubsystemObservable(ptl_label_,ftn_);
    return ftn_(vec_visibles_, vec_invisibles_);
}
//======================================================================================
// Class MassFunctionGroup
//======================================================================================
MassFunctionGroup::MassFunctionGroup(ProcessTree *process_tree,
                                     double (*ftn)(const vector<double> &))
    : MassFunctionInterface(process_tree), num_entries_(0), ftn_(ftn) {}
MassFunctionGroup::MassFunctionGroup(
    ProcessTree *process_tree, double (*ftn)(const vector<double> &),
    vector<MassFunctionInterface *> &vec_mass_interface)
    : MassFunctionInterface(process_tree),
      num_entries_(0),
      ftn_(ftn),
      vec_mass_interface_(vec_mass_interface),
      vec_buffer_(vec_mass_interface.size(), 0.) {}
MassFunctionGroup::MassFunctionGroup(
    ProcessTree *process_tree, double (*ftn)(const vector<double> &),
    std::initializer_list<MassFunctionInterface *> vec_mass_interface)
    : MassFunctionInterface(process_tree),
      num_entries_(0),
      ftn_(ftn),
      vec_mass_interface_(vec_mass_interface),
      vec_buffer_(vec_mass_interface.size(), 0.) {
    // std::cout << "create group" << std::endl;
}
MassFunctionGroup::~MassFunctionGroup() {
    // std::cout << "delete group" << std::endl;
    for (unsigned int i = 0, size = vec_mass_interface_.size(); i < size; ++i) {
        delete vec_mass_interface_[i];
    }
}
//======================================================================================
void MassFunctionGroup::AddSubInterface(MassFunctionInterface *interface) {
    vec_mass_interface_.push_back(interface);
    vec_buffer_.push_back(0.);
}

MassFunctionParticle *MassFunctionGroup::MakeSubParticle(
    string &label, double (*ftn)(const vector<TLorentzVector *> &,
                                 const vector<TLorentzVector *> &)) {
    MassFunctionParticle *p_mass_particle =
        new MassFunctionParticle(process_tree_, label, ftn);
    vec_mass_interface_.push_back(
        dynamic_cast<MassFunctionInterface *>(p_mass_particle));
    vec_buffer_.push_back(0.);
    return p_mass_particle;
}
MassFunctionParticle *MassFunctionGroup::MakeSubParticle(
    const char *label, double (*ftn)(const vector<TLorentzVector *> &,
                                     const vector<TLorentzVector *> &)) {
    MassFunctionParticle *p_mass_particle =
        new MassFunctionParticle(process_tree_, label, ftn);
    vec_mass_interface_.push_back(
        dynamic_cast<MassFunctionInterface *>(p_mass_particle));
    vec_buffer_.push_back(0.);
    return p_mass_particle;
}
MassFunctionGroup *MassFunctionGroup::MakeSubGroup(
    double (*ftn)(const vector<double> &)) {
    MassFunctionGroup *p_mass_group = new MassFunctionGroup(process_tree_, ftn);
    vec_mass_interface_.push_back(
        dynamic_cast<MassFunctionInterface *>(p_mass_group));
    vec_buffer_.push_back(0.);
    return p_mass_group;
}
void MassFunctionGroup::MassFunctionGroup::Init() {
    num_entries_ = vec_mass_interface_.size();
    for (unsigned int i = 0; i < num_entries_; ++i) {
        vec_mass_interface_[i]->Init();
    }
}
double MassFunctionGroup::MassFunctionGroup::Calc() {
    unsigned int i = num_entries_;
    while (--i) { vec_buffer_[i] = vec_mass_interface_[i]->Calc(); }
    vec_buffer_[0] = vec_mass_interface_[0]->Calc();
    return ftn_(vec_buffer_);
}
}  // namespace OptiMass
