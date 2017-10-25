////////////////////////////////////////////////////////////////////////
// Class:       VPlaneHandler
// Module Type: producer
// File:        VPlaneHandler_module.cc
//
// Generated at Wed May 10 08:51:17 2017 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
//#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Wire.h"

#include <memory>
#include <vector>

class VPlaneHandler;

class VPlaneHandler : public art::EDProducer {
public:
  explicit VPlaneHandler(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VPlaneHandler(VPlaneHandler const &) = delete;
  VPlaneHandler(VPlaneHandler &&) = delete;
  VPlaneHandler & operator = (VPlaneHandler const &) = delete;
  VPlaneHandler & operator = (VPlaneHandler &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::string _producer;
  bool _handle_plateau;
  bool _smooth;
  bool _shift;

  float truncated_mean(const std::vector<float>& data,
		       size_t start_idx, size_t end_idx);

  void rolling_mean(std::vector<float>& data,
		    size_t start_idx, size_t end_idx);
};


VPlaneHandler::VPlaneHandler(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  produces< std::vector<recob::Wire> >();
  _producer = p.get<std::string>("ProducerLabel");
  _handle_plateau = p.get<bool>("HandlePlateau",true);
  _smooth = p.get<bool>("Smooth",true);
  _shift  = p.get<bool>("Shift",true);
}

void VPlaneHandler::produce(art::Event & e)
{
  std::unique_ptr< std::vector<recob::Wire> > wire_v_ptr(new std::vector<recob::Wire>);

  art::Handle<std::vector<recob::Wire> > data_h;
  e.getByLabel(_producer,data_h);
  if(!data_h.isValid()) 
    std::cout<<"\033[93m" << "Could not locate data!" << "\033[00m" << std::endl;

  for(size_t wire_idx=0; wire_idx<data_h->size(); ++wire_idx) {

    auto const& w = (*data_h)[wire_idx];
      
    if( w.View() != 1 ) {
      wire_v_ptr->push_back(w);
      continue;
    }

    auto const& roi_v = w.SignalROI().get_ranges();
    std::vector<bool> status_v(roi_v.size(),true);
    for(size_t roi_idx = 0; roi_idx < roi_v.size(); ++roi_idx) {

      auto const& wf = roi_v[roi_idx].data();
      bool status = true;
      size_t cont_negative = 0;
      for(size_t adc_idx=1; adc_idx<wf.size(); ++adc_idx) {
	if(wf[adc_idx]>-5) { cont_negative = 0; continue; }
	++cont_negative;
	if(cont_negative > 5) {
	  status = false;
	  break;
	}
      }
      status_v[roi_idx] = status;
    }
    bool fix=false;
    for(auto const& status : status_v) {
      if(!status) {fix = true; break;}
    }
    if(!fix) {
      wire_v_ptr->push_back(w);
      continue;
    }

    recob::Wire::RegionsOfInterest_t fixed_roi_v;
    for(size_t roi_idx=0; roi_idx<status_v.size(); ++roi_idx) {
      auto const& orig_roi = roi_v[roi_idx];
      if(status_v[roi_idx]) {
	fixed_roi_v.add_range(orig_roi.begin_index(),orig_roi.data());
	continue;
      }
      std::vector<float> copy_data = orig_roi.data();
      bool negative_state=false;
      bool positive_state=false;
      size_t start_idx = 0;
      size_t peak_idx  = 0;
      size_t positive_idx = 0;
      size_t negative_ctr=0; // consecutive negative counter
      float  negative_max=0;
      std::vector<std::array<size_t,4> > flip_range_v;
      for(size_t adc_idx=1; adc_idx<copy_data.size(); ++adc_idx) {

	if(copy_data[adc_idx]<-2 && copy_data[adc_idx-1]<-2) negative_ctr++;
	else negative_ctr = 0;

	if(negative_ctr>5 && !negative_state) {
	  negative_state = true;
	  start_idx = adc_idx-5;
	}

	if(negative_state && !positive_state) {
	  if(copy_data[adc_idx] < negative_max) {
	    negative_max = copy_data[adc_idx];
	    peak_idx = adc_idx;
	  }
	}

	if(negative_state && copy_data[adc_idx] > 2.0) {
	  if(!positive_state) positive_idx = adc_idx-1;
	  positive_state = true;
	}

	if(positive_state && copy_data[adc_idx] < 1.5) {
	  std::array<size_t,4> idxarr;
	  idxarr[0] = start_idx;
	  idxarr[1] = peak_idx;
	  idxarr[2] = positive_idx;
	  idxarr[3] = adc_idx;
	  flip_range_v.push_back(idxarr);
	  negative_state = false;
	  positive_state = false;
	  start_idx = 0;
	  negative_ctr = 0;
	  positive_idx = 0;
	  negative_max = 0;
	  peak_idx = 0;
	}
      }
      if(negative_state) {
	std::array<size_t,4> idxarr;
	idxarr[0] = start_idx;
	idxarr[1] = peak_idx;
	idxarr[2] = positive_idx;
	idxarr[3] = copy_data.size()-1;
	flip_range_v.push_back(idxarr);
      }

      for(auto const& idxarr : flip_range_v) {
	auto const& neg_start = idxarr[0];
	auto const& peak = idxarr[1];
	auto const& pos_start = idxarr[2];
	auto const& end   = idxarr[3];
	for(size_t adc_idx=neg_start; adc_idx<pos_start; ++adc_idx)
	  copy_data[adc_idx] *= -1;
	if(_handle_plateau) {
	  auto mean = truncated_mean(copy_data,pos_start,end);
	  for(size_t adc_idx=pos_start; adc_idx<=end; ++adc_idx)
	    copy_data[adc_idx] -= mean;
	}
	if(_smooth) rolling_mean(copy_data,neg_start,end);
	if(_shift) {
	  size_t shift = pos_start - peak;
	  auto copy_copy_data = copy_data;
	  for(size_t adc_idx=neg_start; (adc_idx+shift)<copy_data.size(); ++adc_idx) {
	    copy_copy_data[adc_idx+shift] = copy_data[adc_idx];
	    if(adc_idx==pos_start) break;
	  }
	  for(size_t adc_idx=neg_start; adc_idx<copy_data.size(); ++adc_idx) {
	    if(adc_idx>=(neg_start+shift)) break;
	    copy_copy_data[adc_idx] = 0;
	  }
	  copy_data = copy_copy_data;
	}
      }
      fixed_roi_v.add_range(orig_roi.begin_index(),copy_data);
    }
    wire_v_ptr->emplace_back(recob::Wire(std::move(fixed_roi_v),w.Channel(), w.View()));
  }

  e.put(std::move(wire_v_ptr));
  e.removeCachedProduct(data_h);
}

float VPlaneHandler::truncated_mean(const std::vector<float>& data,
				    size_t start_idx, size_t end_idx)
{
  float last_mean = -1.;
  float mean = 0.;
  float std = 0.;
  std::vector<bool> skip_v(data.size(),false);
  while( last_mean<0 || std::fabs(last_mean-mean) > std*1.5 ) {
    last_mean = mean;
    // Compute mean
    float sum=0.;
    float ctr=0.;
    for(size_t i=start_idx; i<end_idx; ++i) {
      if(skip_v[i]) continue;
      sum += data[i];
      ctr += 1.;
    }
    mean = sum / ctr;
    // Compute std
    sum = 0.;
    ctr = 0.;
    for(size_t i=start_idx; i<end_idx; ++i) {
      if(skip_v[i]) continue;
      sum += pow(data[i]-mean,2);
      ctr += 1.;
    }
    std = sqrt(sum/ctr);
    // Exclude samples > 2*std from mean
    for(size_t i=start_idx; i<end_idx; ++i) {
      float diff = std::fabs(data[i] - mean);
      if(diff > 2*std) skip_v[i] = true;
    }
  }
  return mean;
}

void VPlaneHandler::rolling_mean(std::vector<float>& data,
				 size_t start_idx, size_t end_idx)
{
  size_t start, end;
  if(data.size() < 9) return;
  for(size_t idx=start_idx; idx<=end_idx; ++idx) {
      
    if(idx>=(start_idx+4)) {
      if(end_idx>=(idx+4)) {
	start = idx-4;
	end   = idx+4;
      }else{
	end   = end_idx;
	start = end_idx-8;
      }
    }else{
      start = start_idx;
      end   = start_idx+8;
    }
    float mean=0;
    for(size_t subidx=start; subidx<=end; ++subidx) mean += data[subidx];
    mean /= 9.;
    data[idx] = mean;
  }
}

DEFINE_ART_MODULE(VPlaneHandler)