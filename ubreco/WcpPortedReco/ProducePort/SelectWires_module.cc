#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <memory>

class SelectWires;


class SelectWires : public art::EDProducer {
public:
  explicit SelectWires(fhicl::ParameterSet const& p);

  SelectWires(SelectWires const&) = delete;
  SelectWires(SelectWires&&) = delete;
  SelectWires& operator=(SelectWires const&) = delete;
  SelectWires& operator=(SelectWires&&) = delete;

  void produce(art::Event& e) override;

private:

  std::string fHitProducer;
  std::string fWireProducer;
  //bool fChargeSupplement;
};


SelectWires::SelectWires(fhicl::ParameterSet const& p)
  : EDProducer{p}  
{
  fHitProducer = p.get<std::string>("HitProducer");
  fWireProducer = p.get<std::string>("WireProducer");
  //fChargeSupplement = p.get<bool>("ChargeSupplement");

  produces<std::vector<raw::RawDigit> >();
  produces<std::vector<recob::Wire> >();
  produces<art::Assns<raw::RawDigit,recob::Wire> >();
}

void SelectWires::produce(art::Event& e)
{
  auto outputRawdigitVec = std::make_unique<std::vector<raw::RawDigit> >();
  auto outputWireVec = std::make_unique<std::vector<recob::Wire> >();
  auto outputRawdigitWireAssoc = std::make_unique<art::Assns<raw::RawDigit, recob::Wire> >();

  art::PtrMaker<raw::RawDigit> rawPtr(e,"");
  art::PtrMaker<recob::Wire> wirePtr(e,"");

  auto const& channelMap = art::ServiceHandle<geo::WireReadout const>()->Get();

  art::Handle<std::vector<recob::Hit> > hit_handle;
  e.getByLabel(fHitProducer, hit_handle);
  if(!hit_handle.isValid())
    {
      mf::LogDebug("SelectWires") << "Failed to find Hit label"<<std::endl;
    }
  lar_pandora::HitVector hit_vec;
  art::fill_ptr_vector(hit_vec, hit_handle);

  art::Handle<std::vector<recob::Wire> > wire_handle;
  e.getByLabel(fWireProducer, wire_handle);
  if(!wire_handle.isValid())
    {
      mf::LogDebug("SelectWires") << "Failed to find Wire label"<<std::endl;
    }
  lar_pandora::WireVector wire_vec;
  art::fill_ptr_vector(wire_vec, wire_handle);

  std::map<int,std::vector<std::pair<int,int> > > hit_map;
  std::pair<int,int> hit_roi;

  for(size_t i_h=0; i_h<hit_vec.size(); ++i_h)
    {
      art::Ptr<recob::Hit> hit = hit_vec.at(i_h);
      hit_roi = std::make_pair( (int)hit->StartTick(), (int)hit->EndTick() );
      hit_map[ (int)hit->Channel() ].emplace_back( hit_roi );
    }

  for(size_t i_w=0; i_w<wire_vec.size(); ++i_w)
    {
      art::Ptr<recob::Wire> wire = wire_vec.at(i_w);
      auto wire_channel = wire->Channel();
      std::vector<float> wf = wire->Signal();
      recob::Wire::RegionsOfInterest_t roi( (int)wf.size() );
      bool found = false;

      for(auto const& m : hit_map)
	{
	  if( m.first != (int)wire_channel) continue;
	  found = true;

	  for(auto & e : m.second)
	    {
	      std::vector<float> wf_temp(wf.size(), 0.);
	      for(int a=e.first; a<=e.second; a++)
		{
		  wf_temp[a] = wf[a];
		}
	      
	      auto first = wf_temp.begin();
	      auto done = wf_temp.end();
	      auto beg = first;
	      while(true)
		{
		  beg = std::find_if(beg, done, [](float v){ return v!=0.0; } );
		  if(beg==done)
		    {
		      break;
		    }
		  auto end = std::find_if(beg, done, [](float v){ return v==0.0;} );
		  
		  std::vector<float> c(beg, end);
		  roi.add_range(beg-first, c.begin(), c.end());
		  beg=end;
		}
	    }
	}

      if(found==true)
	{
          outputWireVec->emplace_back(recob::Wire(roi,wire_channel,channelMap.View(wire_channel)));

	  raw::RawDigit::ADCvector_t outadc;
	  outadc.resize((int)wf.size(), 1);
	  outputRawdigitVec->emplace_back(raw::RawDigit(wire_channel,(int)wf.size(),outadc,raw::kNone));
	  const size_t outind = outputWireVec->size();
	  auto const rawptr = rawPtr(outind-1);
	  auto const sigptr = wirePtr(outind-1);
	  outputRawdigitWireAssoc->addSingle(rawptr,sigptr);

	}
    }
  std::cout<<"[select wire] size: "<<outputWireVec->size()<<std::endl;
  e.put(std::move(outputWireVec));
  e.put(std::move(outputRawdigitVec));
  e.put(std::move(outputRawdigitWireAssoc));

  hit_handle.clear();
  hit_vec.clear();
  wire_handle.clear();
  wire_vec.clear();
}

DEFINE_ART_MODULE(SelectWires)
