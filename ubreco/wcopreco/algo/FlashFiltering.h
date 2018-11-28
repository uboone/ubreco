#ifndef FLASHFILTERING_H
#define FLASHFILTERING_H

#include "Flashes_beam.h"
#include "Flashes_cosmic.h"
#include "Opflash.h"
#include "Config_FlashFiltering.h"

namespace wcopreco{
  //class to filter Flashes
  //don't save cosmic flashes that occur at the same time as a beam flash

  class FlashFiltering {
  public:
    FlashFiltering(OpflashSelection *cosmic_flashes, OpflashSelection *beam_flashes, const Config_FlashFiltering &configFF);
    ~FlashFiltering() {};

    OpflashSelection& get_flashes(){return flashes;};
    OpflashSelection& get_beam_flashes(){return beam_flashes;};
    OpflashSelection& get_cosmic_flashes(){return cosmic_flashes;};

    OpFlashSet& get_all_set(){return all_set;};
    OpFlashSet& get_beam_set(){return beam_set;};
    OpFlashSet& get_cosmic_set(){return cosmic_set;};

    void sort_flashes();
    void update_pmt_map();
    // void clear_flashes();

  protected:

    OpflashSelection flashes;
    OpflashSelection beam_flashes;
    OpflashSelection cosmic_flashes;
    OpFlashSet all_set;
    OpFlashSet beam_set;
    OpFlashSet cosmic_set;

    Config_FlashFiltering _cfg;

  };

}

#endif
