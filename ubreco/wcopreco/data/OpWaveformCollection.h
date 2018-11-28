#ifndef OPWAVEFORMCOLLECTION_H
#define OPWAVEFORMCOLLECTION_H

#include "OpWaveform.h"
#include <vector>
#include <map>

#include <iostream>
#include <sstream>

namespace wcopreco {

  class OpWaveformCollection : public std::vector<OpWaveform> {

  public:
    OpWaveformCollection();
    virtual ~OpWaveformCollection() {};

    void set_channel2index(std::map <int,std::vector<int>>);
    void set_index2channel(std::map <int,int>);
    void insert_channel2index(int channel, int index);
    void insert_index2channel(int index, int channel) {index2channel.insert(std::pair<int,int>(index,channel));};
    int get_index2channel(int index) {return index2channel[index];};
    std::vector<int> get_channel2index(int channel) {return channel2index[channel];};
    void add_waveform(OpWaveform wfm);
    
    std::vector<float> get_op_gain() const {return op_gain;}
    void set_op_gain(std::vector<float> gains_v) {op_gain = gains_v;}

    std::vector<float> get_op_gainerror() const {return op_gainerror;}
    void set_op_gainerror(std::vector<float> gainserror_v) {op_gainerror = gainserror_v;}

  protected:

    std::map <int,std::vector<int>> channel2index;
    std::map <int,int> index2channel;
    std::vector<float> op_gain;
    std::vector<float> op_gainerror;

  };

}

#endif
