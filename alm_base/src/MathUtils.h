#ifndef ALM_BASE_SRC_MATHUTILS_H_
#define ALM_BASE_SRC_MATHUTILS_H_

#include <cmath>
#include <vector>
#include "TLorentzVector.h"

namespace OptiMass {
double CalcMassSquareFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcMTSquareFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

TLorentzVector CalcMomentumFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcMaximum(const std::vector<double> &masses);

double CalcMean(const std::vector<double> &masses);


double CalcMassFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcPxFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcPyFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcPzFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcEFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

double CalcPtFromIndice(
    const std::vector<TLorentzVector *> &vec_visibles,
    const std::vector<TLorentzVector *> &vec_invisibles);

}  // namespace OptiMass

#endif  // ALM_BASE_SRC_MATHUTILS_H_
