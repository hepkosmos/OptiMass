#ifndef OptiMass_MassFunctionInterface_
#define OptiMass_MassFunctionInterface_

#include <initializer_list>
#include <string>
#include <vector>

#include "TLorentzVector.h"

#include "ProcessTree.h"

namespace OptiMass {

class MassFunctionInterface {
protected: 
    ProcessTree* process_tree_;
public:
    MassFunctionInterface(ProcessTree* process_tree);
    virtual ~MassFunctionInterface();

    virtual void Init() = 0;
    virtual double Calc() = 0;
};

class MassFunctionParticle : public MassFunctionInterface {
private: 
    std::string ptl_label_;
    double (*ftn_)(const std::vector< TLorentzVector* >& vec_visibles, const std::vector< TLorentzVector* >& vec_invisibles);
    std::vector<TLorentzVector* > vec_visibles_;
    std::vector<TLorentzVector* > vec_invisibles_;
public:
    MassFunctionParticle(ProcessTree* process_tree, std::string& label, double (*ftn)(const std::vector< TLorentzVector* >&, const std::vector< TLorentzVector* >& ));
    MassFunctionParticle(ProcessTree* process_tree, const char* label, double (*ftn)(const std::vector< TLorentzVector* >&, const std::vector< TLorentzVector* >& ));
    virtual ~MassFunctionParticle();

    void Init();
    double Calc();
};

class MassFunctionGroup : public MassFunctionInterface {
private:
    unsigned int num_entries_;
    double (*ftn_)(const std::vector<double>& masses);
    std::vector<MassFunctionInterface*> vec_mass_interface_;
    std::vector<double> vec_buffer_;
public:
    MassFunctionGroup(ProcessTree* process_tree, double (*ftn)(const std::vector<double>&));
    MassFunctionGroup(ProcessTree* process_tree, double (*ftn)(const std::vector<double>&), std::vector<MassFunctionInterface*>& vec_mass_interface);
    MassFunctionGroup(ProcessTree* process_tree, double (*ftn)(const std::vector<double>&), std::initializer_list<MassFunctionInterface*> vec_mass_interface);
    virtual ~MassFunctionGroup();

    void AddSubInterface(MassFunctionInterface* interface);
    MassFunctionParticle* MakeSubParticle(std::string& label, double (*ftn)(const std::vector< TLorentzVector* >&, const std::vector< TLorentzVector* >& ));
    MassFunctionParticle* MakeSubParticle(const char* label, double (*ftn)(const std::vector< TLorentzVector* >&, const std::vector< TLorentzVector* >& ));
    MassFunctionGroup*  MakeSubGroup(double (*ftn)(const std::vector<double>&));
    void Init();
    double Calc();
};

}//end namespace OptiMass

#endif // OptiMass_MassFunctionInterface_
