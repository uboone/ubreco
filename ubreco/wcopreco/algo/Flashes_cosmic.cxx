#include "Flashes_cosmic.h"

namespace wcopreco {

  wcopreco::Flashes_cosmic::Flashes_cosmic(std::vector<COphitSelection> *ophits_group,  const Config_Opflash &configOpF)
    : _cfgOpF(configOpF)
  {
    //Module for flash finding for cosmics

    int count =0;
    int count_deleted=0;

    for (size_t j=0; j!=ophits_group->size();j++){
      Opflash *flash = new Opflash(ophits_group->at(j), _cfgOpF);
      if (flash->get_total_PE()!=0){
        count++;
        cosmic_flashes.push_back(flash);
      }
      else{
        delete flash;
        count_deleted++;

      }
    }

  }

}
