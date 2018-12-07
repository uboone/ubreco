#include "OpWaveformCollection.h"
#include "OpWaveform.h"

#include <vector>
#include <map>

namespace wcopreco{

  OpWaveformCollection::OpWaveformCollection()
   : std::vector<OpWaveform> (0,OpWaveform(0, 0.0, 0, (0)))
    {
      reserve(500);
      // std::vector<int> v;
      // for (int i =0; i<36; i++){
      //   channel2index[i] = v;
      // }
    }

  void OpWaveformCollection::set_channel2index(std::map <int,std::vector<int>> input_map)
  {channel2index = input_map;}

  void OpWaveformCollection::insert_channel2index(int channel, int index) {
    // for (int i=0; i<channel2index[channel].size(); i++ ) {std::cout << i << " ";}
    // std::cout << std::endl;
    std::vector<int> v = channel2index[channel];
    // std::cout << channel2index[channel].size() << " Map size before insertion ";
    v.emplace_back(std::move(index));

    // for (int i=0; i<v.size(); i++ ) {std::cout << i << ".";}

    channel2index[channel] = v;
    // channel2index.insert(std::pair<int,std::vector<int>>(channel,v));
    // std::cout << channel2index[channel].size() << " Map size after insertion\n";
  }

  void OpWaveformCollection::set_index2channel(std::map <int,int> input_map)
  {index2channel = input_map;}

  void OpWaveformCollection::add_waveform(OpWaveform wfm) {
    emplace_back(std::move(wfm));
    insert_channel2index(wfm.get_ChannelNum(), size()-1);
    insert_index2channel(size()-1, wfm.get_ChannelNum());
  }

}
