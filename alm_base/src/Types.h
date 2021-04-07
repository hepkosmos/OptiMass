#ifndef ALM_BASE_SRC_TYPES_H_
#define ALM_BASE_SRC_TYPES_H_

#include <string>
#include <unordered_map>
#include <vector>
#include "TLorentzVector.h"

namespace OptiMass {
template <typename T>
struct ParticleSubelements {
    std::vector<T> visible;
    std::vector<T> invisible;
    std::vector<T> optimize;
};

using PtlIndexDict = std::unordered_map<std::string, ParticleSubelements<int>>;
using PtlMomentumDict =
    std::unordered_map<std::string, ParticleSubelements<TLorentzVector *>>;
}  // end namespace OptiMass

#endif  // ALM_BASE_SRC_TYPES_H_
