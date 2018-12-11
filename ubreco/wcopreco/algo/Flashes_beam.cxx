#include "Flashes_beam.h"

namespace wcopreco {

  wcopreco::Flashes_beam::Flashes_beam(std::vector<double> *totPE_v,
                                        std::vector<double> *mult_v,
                                        std::vector<double> *l1_totPE_v,
                                        std::vector<double> *l1_mult_v,
                                        std::vector< std::vector<double> > decon_vv,
                                        double beam_start_time,
                                        const Config_FlashesBeam &configFB,
                                        const Config_Opflash &configOpF)
    : _cfgOpF(configOpF) ,  _cfgFB(configFB)
  {
    //Module for flash finding for beams

    std::vector<int> flash_time;
    flash_time.reserve(_cfgFB._nbins_beam / _cfgFB._rebin_frac);
    std::vector<double> flash_pe;

    double *prev_pe_a = new double[_cfgFB._num_channels];
    double *curr_pe_a = new double[_cfgFB._num_channels];

    for (int i=0;i!=_cfgFB._nbins_beam / _cfgFB._rebin_frac;i++){
      double pe = totPE_v->at(i);
      double mult = mult_v->at(i);

      if (pe >= _cfgFB._bflash_pe_thresh && mult >= _cfgFB._bflash_mult_thresh){
        bool flag_save = false;
        //first time through:
        if (flash_time.size()==0){
  	       flag_save = true;
  	       for (int j=0;j!=_cfgFB._num_channels;j++){
  	          prev_pe_a[j] = decon_vv.at(j).at(i);
  	       }
        }
        // all other times through:
        else{
  	       for (int j=0;j!=_cfgFB._num_channels;j++){
              curr_pe_a[j] = decon_vv[j].at(i);
  	       }
  	       if (i - flash_time.back() >= _cfgFB._bflash_bin_diff_p0){
  	          flag_save = true;
           }
           else if (i-flash_time.back() > _cfgFB._bflash_bin_diff_p1 && pe > flash_pe.back()){
  	          if (i-flash_time.back()> _cfgFB._bflash_bin_diff_p2){
  	             flag_save = true;
  	          }
              else{
                 if (KS_maxdiff(_cfgFB._num_channels,prev_pe_a,curr_pe_a) > _cfgFB._KS_test_thresh){
  	               flag_save = true;
  	             }
  	          }
  	       }

           for (int j=0;j!=_cfgFB._num_channels;j++){
              prev_pe_a[j] = decon_vv[j].at(i);
  	       }
        }

        if (flag_save){
          	flash_time.push_back(i);
          	flash_pe.push_back(pe);
        }
        else{
  	       if (i - flash_time.back()<= _cfgFB._bflash_pe_thresh && pe > flash_pe.back())
             {
             flash_pe.back()=pe;
           }
        }
      }
    }
    delete[] prev_pe_a;
    delete[] curr_pe_a;

    for (size_t i=0; i!=flash_time.size(); i++){
      //dertermine start and end bin of flash
      int start_bin = flash_time.at(i)-_cfgFB._bflash_bin_start_cushion;
      if (start_bin <0) start_bin = 0;

      int end_bin = start_bin + _cfgFB._bflash_bin_diff_p0;
      if (end_bin > _cfgFB._nbins_beam / _cfgFB._rebin_frac) end_bin = _cfgFB._nbins_beam / _cfgFB._rebin_frac;
      if (i+1<flash_time.size()){
        if (end_bin > flash_time.at(i+1)-_cfgFB._bflash_bin_start_cushion) {
          end_bin = flash_time.at(i+1)-_cfgFB._bflash_bin_start_cushion;
        }
      }

      //check with the next bin content ...
      //create Opflash

      Opflash *flash = new Opflash(decon_vv, beam_start_time, start_bin, end_bin, _cfgOpF);
      for (size_t p=0; p< l1_totPE_v->size(); p++){
      }
      flash->Add_l1info(l1_totPE_v, l1_mult_v, beam_start_time, start_bin, end_bin, _cfgOpF);
      beam_flashes.push_back(flash);

      }
    }


  double wcopreco::Flashes_beam::KS_maxdiff(int n, double *array1, double *array2){
      //function to calculate KS_maxdiff

     for (int index =1; index<n; index++){
       array1[index] += array1[index-1];
       array2[index] += array2[index-1];
     }
     double maxdiff=0;
     double thisdiff=0;
     for (int index=0; index<n; index++){
       array1[index] = array1[index]/array1[n-1];
       array2[index] = array2[index]/array2[n-1];
       thisdiff =fabs(array1[index]-array2[index]);
       if (thisdiff > maxdiff) {maxdiff = thisdiff;}
     }
     return maxdiff;
   }

}
