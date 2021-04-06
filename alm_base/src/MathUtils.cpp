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


double CalcMassFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    const double m2 = CalcMassSquareFromIndice(vec_visibles, vec_invisibles);
    return m2 < 0 ? -std::sqrt(-m2) : std::sqrt(m2);
}

double CalcPxFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double px = 0;
    for (const auto &v : vec_visibles) { px += v->Px(); }
    for (const auto &v : vec_invisibles) { px += v->Px(); }
    return px;
}

double CalcPyFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double py = 0;
    for (const auto &v : vec_visibles) { py += v->Py(); }
    for (const auto &v : vec_invisibles) { py += v->Py(); }
    return py;
}

double CalcPzFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double pz = 0;
    for (const auto &v : vec_visibles) { pz += v->Pz(); }
    for (const auto &v : vec_invisibles) { pz += v->Pz(); }
    return pz;
}

double CalcEFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    double e = 0;
    for (const auto &v : vec_visibles) { e += v->E(); }
    for (const auto &v : vec_invisibles) { e += v->E(); }
    return e;
}

double CalcPtFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles) {
    TVector2 pt{0, 0};
    for (const auto &v : vec_visibles) { pt += TVector2{v->Px(), v->Py()}; }
    for (const auto &v : vec_invisibles) { pt += TVector2{v->Px(), v->Py()}; }
    return pt.Mod();  // Mod = sqrt(X*X + Y*Y)
}

}  // namespace OptiMass
