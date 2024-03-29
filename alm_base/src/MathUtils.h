#ifndef ALM_BASE_SRC_MATHUTILS_H_
#define ALM_BASE_SRC_MATHUTILS_H_

#include <cmath>
#include <vector>
#include "TLorentzVector.h"

namespace OptiMass {
double CalcMassSquareFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

inline double CalcMassFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    const double m2 = CalcMassSquareFromIndice(vec_visibles, vec_invisibles);
    return m2 < 0 ? -std::sqrt(-m2) : std::sqrt(m2);
}

double CalcMTSquareFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

inline double CalcPxFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double px = 0;
    for (const auto &v : vec_visibles) { px += v->Px(); }
    for (const auto &v : vec_invisibles) { px += v->Px(); }
    return px;
}

inline double CalcPyFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double py = 0;
    for (const auto &v : vec_visibles) { py += v->Py(); }
    for (const auto &v : vec_invisibles) { py += v->Py(); }
    return py;
}

inline double CalcPzFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double pz = 0;
    for (const auto &v : vec_visibles) { pz += v->Pz(); }
    for (const auto &v : vec_invisibles) { pz += v->Pz(); }
    return pz;
}

inline double CalcEFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double e = 0;
    for (const auto &v : vec_visibles) { e += v->E(); }
    for (const auto &v : vec_invisibles) { e += v->E(); }
    return e;
}

inline double CalcPtFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    TVector2 pt{0, 0};
    for (const auto &v : vec_visibles) { pt += TVector2{v->Px(), v->Py()}; }
    for (const auto &v : vec_invisibles) { pt += TVector2{v->Px(), v->Py()}; }
    return pt.Mod();  // Mod = sqrt(X*X + Y*Y)
}

TLorentzVector CalcMomentumFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcMaximum(const std::vector<double> &masses);

double CalcMean(const std::vector<double> &masses);
}  // namespace OptiMass

#endif  // ALM_BASE_SRC_MATHUTILS_H_
