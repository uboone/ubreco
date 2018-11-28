#include "HitFinder_beam.h"

namespace wcopreco {

  wcopreco::HitFinder_beam::HitFinder_beam(OpWaveformCollection &deconvolved_beam, std::vector<kernel_fourier_container> &kernel_container_v, const Config_Hitfinder_Beam &cfg_HB, const Config_Deconvolver &cfg_DC)
  :_cfg(cfg_HB)
  {
    
    //main function for beam hit finder
    //perform deconvolution with kernals add and with filters (true,true)
    wcopreco::Deconvolver filtered_wfm(deconvolved_beam, true, kernel_container_v, cfg_DC);

    OpWaveformCollection filtered_collection = filtered_wfm.Deconvolve_Collection(deconvolved_beam);

    totPE_v.resize(_cfg._nbins_beam/_cfg._rebin_frac);
    mult_v.resize(_cfg._nbins_beam/_cfg._rebin_frac);
    l1_totPE_v.resize(_cfg._nbins_beam/_cfg._rebin_frac);
    l1_mult_v.resize(_cfg._nbins_beam/_cfg._rebin_frac);
    decon_vv.resize(_cfg._num_channels);

    //loop through each channel and perform the l1 fit
    for (int ch=0; ch<_cfg._num_channels; ch++){
      //totPE mult, and their l1 versions are additive (each element is always +=). Each iteration of ch will add to these values.
      decon_vv.at(ch).reserve(300);
      Perform_L1( filtered_collection.at(ch),
                  decon_vv,
                  totPE_v,
                  mult_v,
                  l1_totPE_v,
                  l1_mult_v,
                  ch);
    }
  }

  void HitFinder_beam::Perform_L1(std::vector<double> inverse_res1,
                               std::vector< std::vector<double> > &decon_vv,
                               std::vector<double> &totPE_v,
                               std::vector<double> &mult_v,
                               std::vector<double> &l1_totPE_v,
                               std::vector<double> &l1_mult_v,
                               int ch)
  {
    // prepare L1 fit ...
    std::vector<float> rebin_v;
    rebin_v.resize(_cfg._nbins_beam/_cfg._rebin_frac);

    double new_bin_total;
    for (int i=0;i!=_cfg._nbins_beam/_cfg._rebin_frac;i++){
      new_bin_total =0;
      for (int p=0;p!=_cfg._rebin_frac;p++){
         new_bin_total += inverse_res1.at(_cfg._rebin_frac*i+p) ;
      }
      rebin_v[i] = new_bin_total ;
    }

    decon_vv[ch].resize(_cfg._nbins_beam/_cfg._rebin_frac);
    for (int i=0;i!=_cfg._nbins_beam/_cfg._rebin_frac;i++){
      decon_vv[ch].at(i) = rebin_v[i];

    }

    // work on the L1 ...
    std::vector<double> vals_y;
    std::vector<double> vals_x;
    std::vector<int> vals_bin;

    //save bin content if > threshold of .3
    for (int i=0;i!=_cfg._nbins_beam/_cfg._rebin_frac;i++){
      double content = rebin_v[i];
      if (content > _cfg._l1_content_thresh){
       vals_y.push_back(content);
       vals_x.push_back(i+0.5);
       vals_bin.push_back(i);
      }
    }

    int nbin_fit = vals_x.size();
    Eigen::VectorXd W = Eigen::VectorXd::Zero(nbin_fit);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(nbin_fit,nbin_fit);
    for (int i=0;i!=nbin_fit;i++){
        W(i) = vals_y.at(i) / sqrt(vals_y.at(i));
        double t1 = vals_x.at(i); // measured time
        for (int k=0;k!=nbin_fit;k++){
           double t2 = vals_x.at(k); // real time
           if (t1>t2) {
               G(i,k) = (_cfg._frac_G_t2_first * (exp(-((t1-t2)*_cfg._G_p0*_cfg._tick_width_us-_cfg._G_p1*_cfg._tick_width_us)/_cfg._G_p2)-exp(-((t1-t2)*_cfg._G_p0*_cfg._tick_width_us+3*_cfg._tick_width_us)/_cfg._G_p2))) / sqrt(vals_y.at(i));
           }
            else if (t1==t2){
               G(i,k) = (_cfg._frac_G_sametime + _cfg._frac_G_t2_first *(1-exp(-_cfg._G_p1*_cfg._tick_width_us/_cfg._G_p2))) / sqrt(vals_y.at(i));
           }
            else{
               continue;
           }
        }
    }

    wcopreco::LassoModel m2(_cfg._Lasso_p0, _cfg._Lasso_p1, _cfg._Lasso_p2);
    m2.SetData(G, W);
    m2.Fit();
    Eigen::VectorXd beta = m2.Getbeta();

    //Make vector to hold L1 fit values
    std::vector<double> l1_v;
    l1_v.resize(_cfg._nbins_beam/_cfg._rebin_frac);
    for (int i=0;i!=nbin_fit;i++){
        l1_v[vals_bin.at(i)] = beta(i);
    }

    for (int j=0;j!=_cfg._nbins_beam/_cfg._rebin_frac;j++){
      double content = decon_vv[ch].at(j);

      if (content > _cfg._totPE_v_thresh ) {;
          totPE_v.at(j)= totPE_v.at(j) + content;
        }
      if (content > _cfg._mult_v_thresh) {// ~2 PE threshold ...
          mult_v.at(j)= mult_v.at(j) + 1 ;
        }

      content = l1_v.at(j);

      l1_totPE_v.at(j) = l1_totPE_v.at(j) + content;
      if (content > _cfg._l1_mult_v_thresh) {// 1 PE threshold
          l1_mult_v.at(j) = l1_mult_v.at(j) + 1;
        }
      }
  }

}
