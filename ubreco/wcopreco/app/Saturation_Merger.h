#ifndef SATURATION_MERGER_H
#define SATURATION_MERGER_H

#include "DataReader.h"
#include "UBEventWaveform.h"
#include "Config_Saturation_Merger.h"
namespace wcopreco{

  // This class is microboone specific. It is designed to merge the high and
  // low gain channels for beam, and also for cosmic into a merged_beam and
  // merged_cosmic OpWaveformCollection datamembers, replacing the high gain
  // waveform regions that get saturated if there is a low gain corresponding
  // to it.
  class Saturation_Merger {
  public:
    Saturation_Merger(UBEventWaveform, const Config_Saturation_Merger &);
    ~Saturation_Merger() {};

    OpWaveformCollection get_merged_beam() {return merged_beam;}
    OpWaveformCollection get_merged_cosmic() {return merged_cosmic;}
    UBEventWaveform get_merged_UB_Ev() {return UB_Ev_Merged;}


  protected:
    Config_Saturation_Merger _cfg;
    float findScaling(size_t channel);
    void scale_lowgains(OpWaveformCollection *BLG, OpWaveformCollection *CLG);
    double findBaselineLg(OpWaveform *wfm, int nbin);
    std::vector<std::pair<short,short> > findSaturationTick(OpWaveform *wfm, short saturation_threshold );
    OpWaveform replaceSaturatedBin(OpWaveform &high, OpWaveform &low, std::vector<std::pair<short,short>> saturation_ranges);

    OpWaveformCollection* cosmic_merger(OpWaveformCollection *CHG, OpWaveformCollection *CLG);
    //This function is designed to merge a OpWaveformCollection of cosmic High Gain Waveforms
    //with cosmic Low Gain Waveforms. High gains are the first argument, and low gains are the
    //second argument. A third optional argument allows for customized saturation threshold
    OpWaveformCollection* beam_merger(OpWaveformCollection *BHG, OpWaveformCollection *BLG);
    //This function is designed to merge a OpWaveformCollection of Beam High Gain Waveforms
    //with Beam Low Gain Waveforms. High gains are the first argument, and low gains are the
    //second argument. A third optional argument allows for customized saturation threshold
    OpWaveformCollection merged_cosmic;
    OpWaveformCollection merged_beam;
    UBEventWaveform UB_Ev_Merged;


  };

}

#endif
