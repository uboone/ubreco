#ifndef OPWAVEFORM_H
#define OPWAVEFORM_H

#include <vector>

namespace wcopreco {

  class OpWaveform : public std::vector<double> {
  public:
    OpWaveform(int ChNum, double time_from_trig, int typ, const std::vector<double>& wfm );

    OpWaveform(int ChNum, double time_from_trig, int typ, const int nelems );
    virtual ~OpWaveform() {};

    int get_ChannelNum() const {return ChannelNum;};
    double get_time_from_trigger() const {return time_from_trigger;};
    int get_type() const {return type;};
    void set_ChannelNum(int ch) {ChannelNum = ch;};
    void set_time_from_trigger(double time_diff) {time_from_trigger = time_diff;};
    void set_type(int typ) {type = typ;};



  protected:
    int ChannelNum;
    double time_from_trigger;
    int type;
  };

}

#endif
