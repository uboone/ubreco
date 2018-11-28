#include "LinearModel.h"

#include <Eigen/Dense>
using namespace Eigen;

wcopreco::LinearModel::LinearModel()
{}

wcopreco::LinearModel::~LinearModel()
{}

VectorXd wcopreco::LinearModel::Predict()
{
    return _X * _beta;
}

double wcopreco::LinearModel::chi2_base()
{
    return ( _y - Predict() ).squaredNorm();
}


double wcopreco::LinearModel::MeanResidual()
{
    return ( _y - Predict() ).norm() / _y.size();
}
