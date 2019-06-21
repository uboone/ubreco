////////////////////////////////////////////////////////////////////////
// Class:       OpNoiseCreateMask
// Plugin Type: producer (art v3_01_02)
// File:        OpNoiseCreateMask_module.cc
//
// Generated at Thu Feb 28 14:45:40 2019 by Sophie Berkman using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"

#include "lardata/ArtDataHelper/HitCreator.h"
#include "larreco/HitFinder/HitFilterAlg.h"


#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

class OpNoiseCreateMask;

class OpNoiseCreateMask : public art::EDProducer {
public:
  explicit OpNoiseCreateMask(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpNoiseCreateMask(OpNoiseCreateMask const&) = delete;
  OpNoiseCreateMask(OpNoiseCreateMask&&) = delete;
  OpNoiseCreateMask& operator=(OpNoiseCreateMask const&) = delete;
  OpNoiseCreateMask& operator=(OpNoiseCreateMask&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  float fPEThresh;
  float fFlagFrac;
  float fFlagWireHigh;
  float fFlagWireLow;
  float fFlagTimeHigh;
  float fFlagTimeLow;
  float fNoiseWireHigh;
  float fNoiseWireLow;
  float fNoiseTimeHigh;
  float fNoiseTimeLow;
  bool fApplyOpNoiseMask;
  std::string fHitProducer;

};


OpNoiseCreateMask::OpNoiseCreateMask(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  fPEThresh = p.get<float> ("PEThresh");
  fFlagFrac = p.get<float> ("FlagFrac");

  fFlagWireHigh = p.get<float> ("FlagWireHigh");
  fFlagWireLow = p.get<float> ("FlagWireLow");
  fFlagTimeHigh = p.get<float> ("FlagTimeHigh");
  fFlagTimeLow = p.get<float> ("FlagTimeLow");

  fNoiseWireHigh =p.get<float> ("NoiseWireHigh");
  fNoiseWireLow = p.get<float> ("NoiseWireLow");
  fNoiseTimeHigh =p.get<float> ("NoiseTimeHigh");
  fNoiseTimeLow = p.get<float> ("NoiseTimeLow");

  fApplyOpNoiseMask = p.get<bool> ("ApplyOpNoiseMask");

  fHitProducer = p.get<std::string>("HitProducer");


  // Call appropriate produces<>() functions here.
 
  produces<std::vector<bool> > ();

  if(fApplyOpNoiseMask==true){
    produces<std::vector<recob::Hit> >(); 
  }
 // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void OpNoiseCreateMask::produce(art::Event& e)
{
  // Implementation of required member function here.
  art::InputTag HitInputTag(fHitProducer);
  art::ValidHandle<std::vector<recob::Hit> > inputHits = e.getValidHandle<std::vector<recob::Hit> >(HitInputTag);

  art::InputTag OpHitCosmicInputTag("ophitCosmic");
  art::ValidHandle<std::vector<recob::OpHit> > inputOpHits = e.getValidHandle<std::vector<recob::OpHit> >(OpHitCosmicInputTag);

  art::ServiceHandle<geo::Geometry> geom;

  std::vector<double> pe;
  std::vector<double> peaktime;
  std::vector< std::vector<double> > ophit_wire;
  std::vector<double> frac_flag_to_noise_vec;
  
  bool ispl1_flag=false;
  bool ispl1_noise=false;

  int n_flag_pl1=0;
  int n_noise_pl1=0;

//initialize mask
  const int mask_size=inputHits->size();
  bool default_value=true;
  std::vector< bool > hit_mask(mask_size, default_value);

  //loop over optical hits
  for ( unsigned int ioh=0; ioh<inputOpHits->size(); ioh++) {
    const auto& ophit = inputOpHits->at(ioh);

    std::vector<double> opwiretmp;

    //save if "high" p.e. PMT Pulse and within TPC truncated time window
    //look for noise if over some p.e. threshold
    if(ophit.PE()>fPEThresh && ophit.PeakTime()>-400 && ophit.PeakTime()<3200){

      double PMTxyz[3];
      unsigned int opch;

      pe.push_back(ophit.PE());
      peaktime.push_back(ophit.PeakTime());

      opch=ophit.OpChannel();
      
      geom->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);

      for(int ipl=0; ipl<3; ipl++){
	auto plid = geo::PlaneID(0,0,ipl);
	auto wire = geom->WireCoordinate(PMTxyz[1],PMTxyz[2],plid);
	opwiretmp.push_back(wire);
      }
      ophit_wire.push_back(opwiretmp);
    }
  }//finish ophit loop


  //only look at TPC if there are high pe pmt pulses
  if(pe.size()!=0){
    
    //loop over high p.e. PMT hits
    for(unsigned int i=0; i<pe.size(); i++){
      //loop over TPC hits
      for ( unsigned int ih=0; ih<inputHits->size(); ih++) {
	const auto& tpchit = inputHits->at(ih);     
	auto hit_w=tpchit.WireID().Wire;
	auto hit_pl=tpchit.WireID().Plane;
	
	//convert hit time to coordinates of TPC time, relative  to trigger
	float hit_time_pmtcoord  = (tpchit.PeakTime()/2)-400;
	
	//only try to mask hits on plane 1
	if(hit_pl==1){
	  
	  //cut on position in flag region
	  if(ophit_wire[i][hit_pl]-hit_w>fFlagWireLow && ophit_wire[i][hit_pl]-hit_w<fFlagWireHigh && ispl1_flag==false){
	    //cut on time
	    if(peaktime[i]-hit_time_pmtcoord>fFlagTimeLow && peaktime[i]-hit_time_pmtcoord<fFlagTimeHigh){
	      n_flag_pl1=n_flag_pl1+1;//count number of hits in flag region
	      ispl1_flag=true;
	    }
	  }
	  //cut on noise box
	  if(ophit_wire[i][hit_pl]-hit_w>fFlagWireLow && ophit_wire[i][hit_pl]-hit_w<fFlagWireHigh && ispl1_noise==false){
	    //cut on time                                                                                            
	    if(peaktime[i]-hit_time_pmtcoord>fNoiseTimeLow && peaktime[i]-hit_time_pmtcoord<fNoiseTimeHigh){
	      n_noise_pl1=n_noise_pl1+1;//count number of hits in flag region   
	      ispl1_noise=true;
	    }
	  }
	  
	}//plane 1 only
      }//tpc hits
      double frac_flag_to_noise=(double) n_flag_pl1/n_noise_pl1;
      frac_flag_to_noise_vec.push_back(frac_flag_to_noise);
      ispl1_flag=false; 
      ispl1_noise=false;
    }//pmt loop
  }//pe.size>0

  //loop over high p.e. PMT hits to mask noise
  for(unsigned int i=0; i<pe.size(); i++){
    if(frac_flag_to_noise_vec[i]>fFlagFrac){//require at least 2% in the flag box
      //loop over TPC hits                                                                                                                
      for ( unsigned int ih=0; ih<inputHits->size(); ih++) {
	const auto& tpchit = inputHits->at(ih);
	auto hit_w=tpchit.WireID().Wire;
	auto hit_pl=tpchit.WireID().Plane;
	//convert hit time to coordinates of TPC time, relative  to trigger                                                          
	float hit_time_pmtcoord  = (tpchit.PeakTime()/2)-400;
	//only try to mask hits on plane 1                                                                                         
	if(hit_pl==1){
	  //cut on noise box                                                                            
	  if(ophit_wire[i][hit_pl]-hit_w>fNoiseWireLow && ophit_wire[i][hit_pl]-hit_w<fNoiseWireHigh && ispl1_noise==false){
	    //cut on time                                                                                                              
	    if(peaktime[i]-hit_time_pmtcoord>fNoiseTimeLow && peaktime[i]-hit_time_pmtcoord<fNoiseTimeHigh){
	      hit_mask[ih]=false;
	    }
	  }
	}//plane 1 only                                                                                                               
      }//tpc hits                                                          
    }//flag box
  }//PMT loop      

  //put mask in reco2 file                                                                                                                   
  std::unique_ptr< std::vector<bool> > opnoise_hitmask(new std::vector<bool>);
  for(int i=0; i<mask_size; i++){
    opnoise_hitmask->push_back(hit_mask[i]);
  }

  e.put(std::move(opnoise_hitmask));

  //only make a new hit vector if fcl parameter is set
  if(fApplyOpNoiseMask==true){ 
    std::unique_ptr< std::vector<recob::Hit> > opfilterhit(new std::vector<recob::Hit>);
    for( unsigned int ih=0; ih<inputHits->size(); ih++){
      const auto &tpchit = inputHits->at(ih);
      if(hit_mask[ih]==true){
	opfilterhit->emplace_back(std::move(tpchit));
      }
    }
    e.put(std::move(opfilterhit));
  }

}

DEFINE_ART_MODULE(OpNoiseCreateMask)
