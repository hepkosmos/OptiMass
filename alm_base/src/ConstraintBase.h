#ifndef ALM_BASE_SRC_CONSTRAINTBASE_H_
#define ALM_BASE_SRC_CONSTRAINTBASE_H_

#include <vector>

namespace OptiMass {
class ConstraintBase {
public:
    ConstraintBase() {}
    virtual ~ConstraintBase() {}
    virtual void CalcConstraints(std::vector<double> &vec_constraints,
                                 std::vector<bool> &vec_constraints_using) = 0;
};
}  // namespace OptiMass

#endif  // ALM_BASE_SRC_CONSTRAINTBASE_H_
