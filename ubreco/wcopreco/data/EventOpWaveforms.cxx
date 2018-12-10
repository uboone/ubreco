#include "EventOpWaveforms.h"
#include "OpWaveformCollection.h"
#include "OpWaveform.h"
#include <vector>
#include <map>


namespace wcopreco {

  void EventOpWaveforms::set_type2index(std::map <int,int> input_map)
  {type2index = input_map;}

  void EventOpWaveforms::set_index2type(std::map <int,int> input_map)
  {index2type = input_map;}

  void EventOpWaveforms::set_wfm_v(std::vector<OpWaveformCollection> input_vector)
  {_wfm_v = input_vector;}

  void EventOpWaveforms::set_wfm_v(OpWaveformCollection input_collection)
  {
    _wfm_v.resize(0);
    _wfm_v.emplace_back(std::move(input_collection));
  }

  void EventOpWaveforms::emplace_back_wfm_v(OpWaveformCollection input_collection)
  {_wfm_v.emplace_back(std::move(input_collection));}

  void EventOpWaveforms::emplace_back_wfm_v(std::vector<OpWaveformCollection> input_vector_collection)
  {

    _wfm_v.resize( _wfm_v.size() + input_vector_collection.size() );

    for (size_t i = 0; i<input_vector_collection.size(); i++) {

      _wfm_v.emplace_back(std::move(input_vector_collection[i]));

    }
  }
  void EventOpWaveforms::add_entry(OpWaveformCollection input_collection, int type ) {
    _wfm_v.emplace_back(std::move(input_collection));
    type2index.insert(std::pair<int,int>(type,_wfm_v.size()-1));
    index2type.insert(std::pair<int,int>(_wfm_v.size()-1,type));
  }

  void EventOpWaveforms::push_back_wfm( int type, const OpWaveform& wfm ) {

    auto it = type2index.find(type);
    if ( it==type2index.end() ) {
      throw std::runtime_error("filling for unrecognized type");
    }

    _wfm_v[ it->second ].push_back( wfm );

  }


}
