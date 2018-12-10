#include "OpWaveform.h"
#include <vector>

using namespace wcopreco;

wcopreco::OpWaveform::OpWaveform(int ChNum, double time_from_trig, int typ, const std::vector<double>& wfm)
 : std::vector<double> (wfm)
  {
    ChannelNum = ChNum;
    time_from_trigger = time_from_trig;
    type = typ;
  }


wcopreco::OpWaveform::OpWaveform(int ChNum, double time_from_trig, int typ, const int nelems)
 : std::vector<double> (nelems,0)
  {
    ChannelNum = ChNum;
    time_from_trigger = time_from_trig;
    type = typ;
  }
