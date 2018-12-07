#ifndef HITFINDER_COSMIC_H
#define HITFINDER_COSMIC_H

//data
#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "EventOpWaveforms.h"
#include "COphit.h"

#include "Config_COpHit.h"
#include "Config_Hitfinder_Cosmic.h"

namespace wcopreco{


  class HitFinder_cosmic {
  public:
    HitFinder_cosmic(OpWaveformCollection* merged_cosmic,
                    std::vector<float> *op_gain,
                    std::vector<float> *op_gainerror ,
                    const Config_Hitfinder_Cosmic &configHC,
                    const Config_COpHit &configCOpH);
    ~HitFinder_cosmic() {
      for (auto it = op_hits.begin(); it!=op_hits.end(); it++){
        delete (*it);
      }
      op_hits.clear();
    };

    void clear_ophits();
    std::vector<COphitSelection>  get_ophits_group() {return ophits_group;}
    COphitSelection               get_left_ophits() {return left_ophits;}
    COphitSelection               get_op_hits(){return op_hits;}

  protected:

    std::vector<COphitSelection> ophits_group;
    COphitSelection left_ophits;
    COphitSelection op_hits;

    Config_COpHit _cfgCOpH;
    Config_Hitfinder_Cosmic _cfgHC;

  };

}

#endif
