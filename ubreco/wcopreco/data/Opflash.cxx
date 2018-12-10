#include "Opflash.h"
#include <algorithm>
#include <iostream>

using namespace wcopreco;

wcopreco::Opflash::Opflash(COphitSelection &ophits, const Config_Opflash &configOpF) //cosmic
  : type(1)
  , flash_id (-1)
  , _cfgOpF(configOpF)
{
  PE.resize(_cfgOpF._num_channels);
  PE_err.resize(_cfgOpF._num_channels);
  for (int i=0;i!=_cfgOpF._num_channels;i++){
    PE[i] = 0;
    PE_err[i] = _cfgOpF._PE_err_cosmic;
  }
  time = 0;
  total_PE = 0;
  for (size_t i=0; i!= ophits.size(); i++){
    fired_channels.push_back(ophits.at(i)->get_ch_no());
    PE[ophits.at(i)->get_ch_no()] = ophits.at(i)->get_PE() - _cfgOpF._PE_subtract*2;
    PE_err[ophits.at(i)->get_ch_no()] = ophits.at(i)->get_PE_err();
    time += ophits.at(i)->get_PE() * ophits.at(i)->get_time();
    total_PE += ophits.at(i)->get_PE() ;
  }

  if (total_PE !=0){
    time /= total_PE;
  }else{
    time = ophits.at(0)->get_time();
  }
  float bin_width = _cfgOpF._tick_width_us; //this  should be derived from det constants in larlite/larsoft
  low_time = time - _cfgOpF._flash_low_time_cushion * bin_width;
  high_time = time + _cfgOpF._flash_high_time_cushion * bin_width;

}

wcopreco::Opflash::Opflash(const std::vector<std::vector<double>> &vec_v, double start_time, int start_bin, int end_bin , const Config_Opflash &configOpF) //beam
  : type(2)
  , flash_id (-1)
  , _cfgOpF (configOpF)
{
  PE.resize(_cfgOpF._num_channels);
  PE_err.resize(_cfgOpF._num_channels);
  float bin_width = _cfgOpF._rebin_frac*_cfgOpF._tick_width_us;
  low_time = start_time + (start_bin+0.5)*bin_width;
  high_time = start_time + (end_bin-0.5)*bin_width;

  for (int i=0;i!=_cfgOpF._num_channels;i++){
    PE[i] = 0;
    PE_err[i] = _cfgOpF._PE_err_beam;
  }

  double max = 0;
  int max_bin = 0;

  for (int i=start_bin; i!=end_bin;i++){
    double peak = 0;
    double mult = 0;
    for (int j=0;j!=_cfgOpF._num_channels;j++){
      double content = vec_v[j].at(i);
      if (content < _cfgOpF._PE_err_beam) content = 0;
      peak += content;
      PE[j] += content;
      if (content>_cfgOpF._mult_content_thresh){
	       mult++;
      }
    }

    if (peak > max && mult >=_cfgOpF._mult_required){
      max = peak;
      max_bin = i;
    }
  }

  time = start_time + (max_bin+0.5) * bin_width;


  total_PE = 0;
  for (int i=0;i!=_cfgOpF._num_channels;i++){
    PE[i] -= _cfgOpF._PE_err_stat_beam ;
    if (PE[i]<0) PE[i] = 0;
    total_PE += PE[i];

    if (PE[i]>=2){
      fired_channels.push_back(i);
    }

    PE_err[i] = sqrt(pow(PE_err[i],2) // standard error
		     + (PE[i] + _cfgOpF._PE_err_stat_beam*2) //statistical term
		     + pow(PE[i]*_cfgOpF._PE_err_unc_beam,2) // 2% base systematic uncertainties
		     );

  }

}

wcopreco::Opflash::~Opflash(){

}

void wcopreco::Opflash::swap_channels(){
  //harded coded for now ...
  for (auto it = fired_channels.begin(); it!=fired_channels.end(); it++){
    if (*it == 26) *it = 27;
    if (*it == 27) *it = 28;
    if (*it == 28) *it = 29;
    if (*it == 29) *it = 30;
    if (*it == 30) *it = 31;
    if (*it == 31) *it = 26;
  }

  double temp;
  temp = PE[31];
  PE[31] = PE[30];
  PE[30] = PE[29];
  PE[29] = PE[28];
  PE[28] = PE[27];
  PE[27] = PE[26];
  PE[26] = temp;

  temp = PE_err[31];
  PE_err[31] = PE_err[30];
  PE_err[30] = PE_err[29];
  PE_err[29] = PE_err[28];
  PE_err[28] = PE_err[27];
  PE_err[27] = PE_err[26];
  PE_err[26] = temp;
}

void wcopreco::Opflash::Add_l1info(std::vector<double> *totPE_v, std::vector<double> *mult_v, double start_time, int start_bin, int end_bin, const Config_Opflash &configOpF){
  float bin_width = _cfgOpF._rebin_frac*_cfgOpF._tick_width_us;
  std::vector<int> fired_bin;
  std::vector<double> fired_pe;
  for (int i=start_bin; i!=end_bin;i++){
    double pe = totPE_v->at(i);
    double mult = mult_v->at(i);
    if (pe >=_cfgOpF._addl1_pe_thresh && mult>=_cfgOpF._addl1_mult_thresh){
      fired_bin.push_back(i);
      fired_pe.push_back(pe);
    }
  }

  for (size_t i=0; i!=fired_bin.size();i++){
    if (i==0){
      double time = start_time + (fired_bin.at(0)+0.5)*bin_width;
      l1_fired_time.push_back(time);
      l1_fired_pe.push_back(fired_pe.at(0));
    }else{
      if (fired_bin.at(i)-fired_bin.at(i-1)==1){
	l1_fired_pe.back() += fired_pe.at(i);
      }else{
	double time = start_time + (fired_bin.at(i)+0.5)*bin_width;
	l1_fired_time.push_back(time);
	l1_fired_pe.push_back(fired_pe.at(i));
      }
    }
  }

}


bool wcopreco::Opflash::get_fired(int ch){
  if (std::find(fired_channels.begin(),fired_channels.end(),ch)==fired_channels.end()){
    return false;
  }else{
    return true;
  }
}
