#include "FlashFiltering.h"

namespace wcopreco {

  wcopreco::FlashFiltering::FlashFiltering(OpflashSelection *cosmic_flashes, OpflashSelection *beam_flashes, const Config_FlashFiltering &configFF)
    : _cfg(configFF)
  {
    //code to perform flash filtering
    //inputs come from Flashes_cosmic and Flashes_beam
    Opflash *prev_cflash = 0;

    //loop to decide if each cosmic flash should be saved
    for (size_t i=0; i!=cosmic_flashes->size();i++){
      Opflash *cflash = cosmic_flashes->at(i);
      bool save = true;
      for (size_t j=0; j!=beam_flashes->size();j++){
        Opflash *bflash = beam_flashes->at(j);
        if (cflash->get_time() >= bflash->get_low_time() && cflash->get_time() <= bflash->get_high_time()){
        	save = false;
        	break;
        }
      }

      if (save){
        flashes.push_back(cflash);
        if (prev_cflash==0 ){
        	if (cflash->get_time()<0)
        	  prev_cflash = cflash;
        }
        else{
        	if (cflash->get_time() < 0 && cflash->get_time() > prev_cflash->get_time())
        	  prev_cflash = cflash;
        }
      }
    }

    //if the cosmic isn't saved, save beam instead
    for (size_t j=0; j!=beam_flashes->size();j++){
      Opflash *bflash = beam_flashes->at(j);
      if (prev_cflash!=0){;
        if (bflash->get_time() - prev_cflash->get_time() < _cfg._flash_filter_time_thresh && // veto for 3 us
  	    bflash->get_total_PE() < _cfg._flash_filter_pe_thresh * prev_cflash->get_total_PE())
  	       continue;
      }
      flashes.push_back(bflash);
    }

    sort_flashes();
    update_pmt_map();

  }

  void FlashFiltering::sort_flashes(){
    //sort flashes by time for beam, cosmic, and total

    for (auto it= cosmic_flashes.begin(); it!= cosmic_flashes.end(); it++){
      cosmic_set.insert(*it);
    }
    cosmic_flashes.clear();
    std::copy(cosmic_set.begin(), cosmic_set.end(), std::back_inserter(cosmic_flashes));

    for (auto it=beam_flashes.begin(); it!=beam_flashes.end(); it++){
      beam_set.insert(*it);
    }
    beam_flashes.clear();
    std::copy(beam_set.begin(), beam_set.end(), std::back_inserter(beam_flashes));

    for (auto it=flashes.begin(); it!=flashes.end(); it++){
      all_set.insert(*it);
    }
    flashes.clear();
    std::copy(all_set.begin(), all_set.end(), std::back_inserter(flashes));
  }


  void wcopreco::FlashFiltering::update_pmt_map(){
    for (auto it=flashes.begin(); it!=flashes.end(); it++){
      Opflash *flash = *it;
      if (_cfg._do_swap_channels )flash->swap_channels();
    }
  }


  // void wcopreco::FlashFiltering::clear_flashes(){
  //   for (auto it = beam_flashes.begin(); it!=beam_flashes.end(); it++){
  //       delete (*it);
  //     }
  //     for (auto it=cosmic_flashes.begin(); it!= cosmic_flashes.end(); it++){
  //       delete (*it);
  //     }
  //     beam_flashes.clear();
  //     cosmic_flashes.clear();
  //     flashes.clear();
  // }


}
