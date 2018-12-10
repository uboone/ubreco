#ifndef WIRECELLRESS_LASSOMODEL_H
#define WIRECELLRESS_LASSOMODEL_H

#include "ElasticNetModel.h"

namespace wcopreco {

class LassoModel: public ElasticNetModel {
public:
    LassoModel(double lambda=1., int max_iter=100000, double TOL=1e-3, bool non_negtive=true);
    ~LassoModel();

    void Fit();
    void Set_init_values(std::vector<double> values);

    double chi2_l1();

 private:
    bool flag_initial_values;
    std::vector<double> init_betas;

};

}

#endif
