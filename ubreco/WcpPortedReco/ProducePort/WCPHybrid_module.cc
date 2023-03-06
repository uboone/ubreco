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

#include "TF1.h"
#include "TH1.h"
#include "TMath.h"

#include <memory>

class WCPHybrid;


class WCPHybrid : public art::EDProducer {
public:
  explicit WCPHybrid(fhicl::ParameterSet const& p);

  WCPHybrid(WCPHybrid const&) = delete;
  WCPHybrid(WCPHybrid&&) = delete;
  WCPHybrid& operator=(WCPHybrid const&) = delete;
  WCPHybrid& operator=(WCPHybrid&&) = delete;

  void produce(art::Event& e) override;

  recob::Wire::RegionsOfInterest_t unbinROI(recob::Wire::RegionsOfInterest_t ROI, int wfSize);

private:

  std::string fHitProducer;
  std::string fWireProducer;
  bool fChargeSupplement;
  bool fUnbin;
  
};


WCPHybrid::WCPHybrid(fhicl::ParameterSet const& p)
  : EDProducer{p}
{
  fHitProducer = p.get<std::string>("HitProducer");
  fWireProducer = p.get<std::string>("WireProducer");
  fChargeSupplement = p.get<bool>("ChargeSupplement");
  fUnbin = p.get<bool>("Unbin");

  produces<std::vector<raw::RawDigit> >();
  produces<std::vector<recob::Wire> >();
  produces<art::Assns<raw::RawDigit,recob::Wire> >();

}

void WCPHybrid::produce(art::Event& e)
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
      mf::LogDebug("WCPHybrid") << "Failed to find Hit label"<<std::endl;
    }
  lar_pandora::HitVector hit_vec;
  art::fill_ptr_vector(hit_vec, hit_handle);

  art::Handle<std::vector<recob::Wire> > wire_handle;
  e.getByLabel(fWireProducer, wire_handle);
  if(!wire_handle.isValid())
    {
      mf::LogDebug("WCPHybrid") << "Failed to find Wire label"<<std::endl;
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

  for(size_t i_w=0; i_w<wire_vec.size(); ++i_w) // loop over all 8256 wires
    {
      art::Ptr<recob::Wire> wire = wire_vec.at(i_w);
      auto wire_channel = wire->Channel();
      std::vector<float> wf = wire->Signal();
      recob::Wire::RegionsOfInterest_t roi( (int)wf.size() );
      //recob::Wire::RegionsOfInterest_t unbinnedROI( (int)wf.size() );

      bool found = false;
      bool dead = false;

      for(auto const& m : hit_map) // loop over (WCP-ported) hit channel-tick pairs
	{
	  if( m.first != (int)wire_channel) {
	    continue;
	  }
	  found=true;

	  double totalsignal=0;
	  for(int i=0; i<(int)wf.size(); i++){
	    totalsignal+=wf.at(i);
	  }

	  if(totalsignal==0.0){
	    dead=true;
	  }

	  for(auto & e : m.second) // loop over pairs of ticks of for a given <hit channel-time tick> pair
	    {
	      std::vector<float> wf_temp(wf.size(), 0.);
	      for(int a=e.first; a<=e.second; a++) // e.first = start tick, e.second = end tick
		{
		  if(a<0 || a>(int)wf.size()-1) continue;
		  wf_temp[a] = wf[a];
		}

	      auto first = wf_temp.begin();
	      auto done = wf_temp.end();
	      auto beg = first;

	      while(true)
		{
		  beg = std::find_if(beg, done, [](float v){ return v!=0.0; } ); // find the first "non-zero" element's iterator 
		  if(beg==done) // if there's no element that satisfies the condition "element non-zero", then quit
		    {
		      break;
		    }
		  auto end = std::find_if(beg, done, [](float v){ return v==0.0;} ); // find the first "zero" element "after" beg iterator found above
		  
		  std::vector<float> c(beg, end);
		  roi.add_range(beg-first, c.begin(), c.end()); //(offset:where to add the elements, first, last)
		  beg=end;
		}
	    }
	}

      if(found==true && dead==false) // live channels
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
      else if(found==true && dead==true) // charge supplement for dead channels                  
	{
	  if(fChargeSupplement){

	    for(size_t i_h=0; i_h<hit_vec.size(); ++i_h)
	      {
		art::Ptr<recob::Hit> hit = hit_vec.at(i_h);
		if((int)hit->Channel()==(int)wire_channel){
		  for(int i=(int)hit->StartTick(); i<(int)hit->EndTick()+1; ++i){
		    roi.set_at(i, hit->PeakAmplitude());
		  }
		}
	      }
	  }
	  if(fUnbin==true) // unbinned WCP hit
            outputWireVec->emplace_back(recob::Wire(unbinROI(roi,(int)wf.size()),wire_channel,channelMap.View(wire_channel)));
	  else if(fUnbin==false) // raw WCP hit
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

  std::cout<<"[wcp hybrid wire] size: "<<outputWireVec->size()<<std::endl;

  e.put(std::move(outputWireVec));
  e.put(std::move(outputRawdigitVec));
  e.put(std::move(outputRawdigitWireAssoc));

  hit_handle.clear();
  hit_vec.clear();
  wire_handle.clear();
  wire_vec.clear();
}

recob::Wire::RegionsOfInterest_t WCPHybrid::unbinROI(const recob::Wire::RegionsOfInterest_t ROI, const int wfSize)
{
  recob::Wire::RegionsOfInterest_t unbinnedROI(wfSize);

  //Find peak
  Int_t startPeak[100]={0};
  Int_t endPeak[100]={0};
  Float_t area[100]={0.0};
  Int_t nPeaks = 0;
  Bool_t inPeak = false;

  //std::cout<<"wcp original"<<std::endl;
  //for(int i=0; i<wfSize; i++){
  //if(ROI[i]!=0){
  //std::cout<<i<<" "<<ROI[i]<<std::endl;
  //}
  //}
  
  for(int i=0;i<wfSize;i++){
    if(ROI[i]>0 && ROI[i+1]>0 && ROI[i+2]>0 && ROI[i+3]>0){
      if(inPeak==false){
	startPeak[nPeaks]=i;
	area[nPeaks]=0;
      }
      area[nPeaks] += ROI[i]+ROI[i+1]+ROI[i+2]+ROI[i+3];

      inPeak = true;
      i += 3;
    }else{
      if(inPeak==true){
	endPeak[nPeaks]=i-1;
	nPeaks += 1;
	inPeak = false;
      }
    }
  }
  
  if(inPeak==true){
    endPeak[nPeaks]=6400-1;
    nPeaks += 1;
    inPeak = false;
  }

  //////////////////////////////////////
  // Define Gaussian functions for peaks

  //std::cout<<"number of peaks: "<<nPeaks<<std::endl;
  for(int p=0; p<nPeaks; p++){
    //std::cout<<"peak "<<p<<std::endl;
    //std::cout<<"start peak, end peak: "<<startPeak[p]<<" "<<endPeak[p]<<std::endl;
    float u=0;
    float sigma=0;
    float AA=0;

    u=startPeak[p]+ (endPeak[p]-startPeak[p])/2.0;
    sigma=(endPeak[p]-startPeak[p])/5.0;
    AA=area[p]/(sigma*TMath::Sqrt(2*TMath::Pi()));
    //std::cout<<"u, sigma, AA: "<<u<<" "<<sigma<<" "<<AA<<std::endl;

    TF1 *gauss = new TF1("gauss", "gaus(0)",startPeak[p],endPeak[p]);
    gauss->FixParameter(0,AA);
    gauss->FixParameter(1,u);
    gauss->FixParameter(2,sigma);

    //std::cout<<"unbinned"<<std::endl;
    for(int t=startPeak[p]; t<endPeak[p]+1; ++t){
      unbinnedROI.set_at(t, gauss->Eval(t));
      //if(unbinnedROI[t]!=0)
	//std::cout<<t<<" "<<unbinnedROI[t]<<std::endl;
    }
  }

  return unbinnedROI;

}

DEFINE_ART_MODULE(WCPHybrid)
