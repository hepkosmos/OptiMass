#include "MathUtils.h"
#include <algorithm>
#include <numeric>
#include <vector>
#include "TLorentzVector.h"

using std::vector;

namespace OptiMass {
double CalcMassSquareFromIndice(
    const vector<TLorentzVector *> &vec_visibles,
    const vector<TLorentzVector *> &vec_invisibles) {
    const auto v = CalcMomentumFromIndice(vec_visibles, vec_invisibles);
    return v.M2();
}

TLorentzVector sumVectors(const vector<TLorentzVector *> &vs) {
    TLorentzVector sumV{0, 0, 0, 0};
    for (const auto &v : vs) { sumV += *v; }
    return sumV;
}

double CalcMTSquareFromIndice(const vector<TLorentzVector *> &vec_visibles,
                              const vector<TLorentzVector *> &vec_invisibles) {
    const auto visSum = sumVectors(vec_visibles);
    double momenta_et = visSum.Mt();  // Mt = sqrt(E*E - Z*Z)
    TVector2 momenta_xy{visSum.Px(), visSum.Py()};
    for (const auto &inv : vec_invisibles) {
        momenta_et += inv->Mt();
        momenta_xy += TVector2{inv->Px(), inv->Py()};
    }
    return momenta_et * momenta_et - momenta_xy.Mod2();  // Mod2 = X*X + Y*Y
}

TLorentzVector CalcMomentumFromIndice(
    const vector<TLorentzVector *> &vec_visibles,
    const vector<TLorentzVector *> &vec_invisibles) {
    return sumVectors(vec_visibles) + sumVectors(vec_invisibles);
}

double CalcMaximum(const vector<double> &masses) {
    return *std::max_element(masses.cbegin(), masses.cend());
}

double CalcMean(const vector<double> &masses) {
    const double buf = std::accumulate(masses.cbegin(), masses.cend(), 0.);
    return buf / masses.size();
}
}  // namespace OptiMass
