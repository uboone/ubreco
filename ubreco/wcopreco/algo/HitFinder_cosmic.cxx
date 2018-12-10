#include "HitFinder_cosmic.h"

namespace wcopreco {

  wcopreco::HitFinder_cosmic::HitFinder_cosmic(OpWaveformCollection* merged_cosmic,
                                              std::vector<float> *op_gain,
                                              std::vector<float> *op_gainerror,
                                              const Config_Hitfinder_Cosmic &configHC,
                                              const Config_COpHit &configCOpH)
    : _cfgCOpH(configCOpH), _cfgHC(configHC)
  {
    //Module for hit finding for cosmics
    //Much of this code can be left the way it is in WC
      int count =0;
      for (size_t i=0; i!=merged_cosmic->size(); i++){
        OpWaveform wfm_cosmic = merged_cosmic->at(i);
        int channel = wfm_cosmic.get_ChannelNum();
        double timestamp = wfm_cosmic.get_time_from_trigger();
        COphit *op_hit = new COphit(channel, &wfm_cosmic, timestamp, op_gain->at(channel), op_gainerror->at(channel), _cfgCOpH);

        op_hits.push_back(op_hit);

        if (op_hit->get_type()){
          //get_type returns flag for good baseline
          bool flag_used = false;

          //first time through, make ophits_group
          if (ophits_group.size()==0){
          	COphitSelection ophits;
          	ophits.push_back(op_hit);
            count++;
            // std::cout << i <<  " pushed back an ophit, now totalling: " << count << "\n";
          	ophits_group.push_back(ophits);
          	flag_used = true;
          }

          else{
          	for (size_t j=0; j!=ophits_group.size();j++){
          	  for (size_t k=0; k!= ophits_group.at(j).size(); k++){
          	    if (fabs(op_hit->get_time() - ophits_group.at(j).at(k)->get_time()) < _cfgHC._ophit_group_t_diff_max ){
          	      ophits_group.at(j).push_back(op_hit);
          	      flag_used = true;
          	      break;
          	    }
          	  }
          	  if (flag_used)
          	    break;
          	}
          }

          if (!flag_used){
          	COphitSelection ophits;
          	ophits.push_back(op_hit);
            count++;
            // std::cout << i <<  " pushed back an ophit, now totalling: " << count << "\n";
          	ophits_group.push_back(ophits);
          }
        }

        //if not good baseline
        else{
          left_ophits.push_back(op_hit);
        }
      }

      for (size_t i=0;i!=left_ophits.size();i++){
        bool flag_used = false;
        for (size_t j=0; j!=ophits_group.size();j++){
          for (size_t k=0; k!= ophits_group.at(j).size(); k++){
          	if (fabs(left_ophits.at(i)->get_time() - ophits_group.at(j).at(k)->get_time())< _cfgHC._ophit_group_t_diff_max ){
          	  ophits_group.at(j).push_back(left_ophits.at(i));
          	  flag_used = true;
          	  break;
          	}
          }
          if (flag_used)
      	     break;
        }
      }

  }

  void HitFinder_cosmic::clear_ophits(){
    for (auto it = op_hits.begin(); it!=op_hits.end(); it++){
      delete (*it);
    }
    op_hits.clear();
  }

}
