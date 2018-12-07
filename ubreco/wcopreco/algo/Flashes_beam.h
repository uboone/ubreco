#ifndef FLASHES_BEAM_H
#define FLASHES_BEAM_H

//data
#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "EventOpWaveforms.h"
#include "Opflash.h"
#include "HitFinder_beam.h"

#include "Config_Opflash.h"
#include "Config_FlashesBeam.h"

namespace wcopreco{
  // class for finding flashes in beam data
  // takes in totalPE, totalMult, totalPE_l1, totalMult_l1, decon_vv, and the  start time of the beam
  // all inputs are from HitFinder_beam
  // decon_vv is a vector of vectors

  class Flashes_beam {
  public:
    Flashes_beam(std::vector<double> *totPE_v,
                std::vector<double> *mult_v,
                std::vector<double> *l1_totPE_v,
                std::vector<double> *l1_mult_v,
                std::vector< std::vector<double> > decon_vv,
                double beam_start_time,
                const Config_FlashesBeam &configFB,
                const Config_Opflash &configOpF
              );
    ~Flashes_beam() {};

    OpflashSelection get_beam_flashes(){return beam_flashes;};

  protected:
    double KS_maxdiff(int n, double *array1, double *array2);
    OpflashSelection beam_flashes;

    Config_Opflash _cfgOpF;
    Config_FlashesBeam _cfgFB;

  };

}

#endif
