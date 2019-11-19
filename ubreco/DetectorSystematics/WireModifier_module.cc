////////////////////////////////////////////////////////////////////////
// Class:       WireModifier
// Plugin Type: producer (art v3_01_02)
// File:        WireModifier_module.cc
//
// Generated at Thu Aug 22 12:19:37 2019 by Wesley Ketchum using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "TFile.h"
#include "TSpline.h"

namespace sys {
  class WireModifier;
}


class sys::WireModifier : public art::EDProducer {
public:
  explicit WireModifier(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireModifier(WireModifier const&) = delete;
  WireModifier(WireModifier&&) = delete;
  WireModifier& operator=(WireModifier const&) = delete;
  WireModifier& operator=(WireModifier&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fWireInputTag;
  art::InputTag fSimEdepShiftedInputTag;
  art::InputTag fSimEdepOrigInputTag;
  art::InputTag fHitInputTag;
  bool          fMakeRawDigitAssn;
  bool          fCopyRawDigitAssn;

  std::string fSplinesFileName;
  std::vector<std::string> fSplineNames_Charge_X;
  std::vector<std::string> fSplineNames_Sigma_X;
  
  std::vector<TSpline3*> fTSplines_Charge_X;
  std::vector<TSpline3*> fTSplines_Sigma_X;

  //useful math things
  //static constexpr double ONE_OVER_SQRT_2PI = 1./std::sqrt(2*util::pi());
  double GAUSSIAN(double t, double mean,double sigma,double a=1.0){
    return ( (a/sigma /std::sqrt(2*util::pi()) * std::exp( -1.*(t-mean)*(t-mean)*0.5/sigma/sigma) ));
  }

  static constexpr double A_w = 3.33328;
  static constexpr double C_U = 338.140;
  static constexpr double C_V = 2732.53;
  static constexpr double C_Y = 4799.19;
  static constexpr double A_t = 18.2148;
  static constexpr double C_t = 818.351;
  static constexpr double SIN_SIXTY = std::sqrt(3)/2;
  static constexpr double COS_SIXTY = 0.5;


  typedef std::pair<unsigned int,unsigned int> ROI_Key_t;
  std::map< ROI_Key_t,std::vector<size_t> > fROIMatchedEdepMap;
  std::map< ROI_Key_t,std::vector<size_t> > fROIMatchedHitMap;
  
  typedef struct ROIProperties{
    ROI_Key_t key;
    unsigned int plane;
    float begin;
    float end;
    float total_q;
    float center;   //charge weighted center of ROI
    float sigma;    //charge weighted RMS of ROI
  } ROIProperties_t;

  typedef std::pair<ROI_Key_t, unsigned int> SubROI_Key_t;
  
  typedef struct SubROIProperties{
    SubROI_Key_t key;
    unsigned int plane;
    float total_q;
    float center;
    float sigma;
  } SubROIProperties_t;

  typedef struct TruthProperties{
    float x;
    float x_rms;
    float x_rms_noWeight;
    float tick;
    float tick_rms;
    float tick_rms_noWeight;
    float total_energy;
    float x_min;
    float x_max;
    float tick_min;
    float tick_max;
  } TruthProperties_t;

  typedef struct ScaleValues{
    float r_Q;
    float r_sigma;
  } ScaleValues_t;

  ROIProperties_t CalcROIProperties(recob::Wire::RegionsOfInterest_t::datarange_t const&);

  std::vector< std::pair<unsigned int, unsigned int> > GetTargetROIs(sim::SimEnergyDeposit const&);
  std::vector< std::pair<unsigned int, unsigned int> > GetHitTargetROIs(recob::Hit const&);
  
  void FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const&, std::vector<recob::Wire> const&);
  void FillROIMatchedHitMap(std::vector<recob::Hit> const&, std::vector<recob::Wire> const&);

  std::vector<SubROIProperties_t> CalcSubROIProperties(ROIProperties_t const&, std::vector<const recob::Hit*> const&);

  std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> MatchEdepsToSubROIs(std::vector<SubROIProperties_t> const&, std::vector<const sim::SimEnergyDeposit*> const&);

  TruthProperties_t CalcPropertiesFromEdeps(std::vector<const sim::SimEnergyDeposit*> const&);

  ScaleValues_t GetScaleValues(TruthProperties_t const&,ROIProperties_t const&);

  void ModifyROI(std::vector<float> &,
		 ROIProperties_t const &, 
		 std::vector<SubROIProperties_t> const&,
		 std::map<SubROI_Key_t, ScaleValues_t> const&);

};


sys::WireModifier::ROIProperties_t 
sys::WireModifier::CalcROIProperties(recob::Wire::RegionsOfInterest_t::datarange_t const& roi)
{
  ROIProperties_t roi_vals;
  roi_vals.begin = roi.begin_index();
  roi_vals.end = roi.end_index();

  roi_vals.center=0;
  roi_vals.total_q=0;
  roi_vals.sigma=0;

  auto const& roi_data = roi.data();
  for(size_t i_t = 0; i_t<roi_data.size(); ++i_t){
    roi_vals.center += roi_data[i_t]*(i_t+roi_vals.begin);
    roi_vals.total_q += roi_data[i_t];
  }
  roi_vals.center = roi_vals.center/roi_vals.total_q;

  for(size_t i_t = 0; i_t<roi_data.size(); ++i_t)
    roi_vals.sigma += roi_data[i_t]*(i_t+roi_vals.begin-roi_vals.center)*(i_t-+roi_vals.begin-roi_vals.center);
  roi_vals.sigma = std::sqrt(roi_vals.sigma/roi_vals.total_q);

  return roi_vals;
}

std::vector< std::pair<unsigned int, unsigned int> > 
sys::WireModifier::GetTargetROIs(sim::SimEnergyDeposit const& shifted_edep)
{
  //channel number, time tick
  std::vector< std::pair<unsigned int,unsigned int> > target_roi_vec;

  int edep_U_wire = std::round( A_w*(-SIN_SIXTY*shifted_edep.Y() + COS_SIXTY*shifted_edep.Z()) + C_U );
  int edep_V_wire = std::round( A_w*( SIN_SIXTY*shifted_edep.Y() + COS_SIXTY*shifted_edep.Z()) + C_V );
  int edep_Y_wire = std::round( A_w*shifted_edep.Z() + C_Y );
  int edep_tick   = std::round( A_t*shifted_edep.X() + C_t );

  if (edep_tick<0 || edep_tick>=6400)
    return target_roi_vec;

  if(edep_U_wire>=0 && edep_U_wire<2400)
    target_roi_vec.emplace_back((unsigned int)edep_U_wire,(unsigned int)edep_tick);

  if(edep_V_wire>=2400 && edep_V_wire<4800)
    target_roi_vec.emplace_back((unsigned int)edep_V_wire,(unsigned int)edep_tick);

  if(edep_Y_wire>=4800 && edep_Y_wire<8256)
    target_roi_vec.emplace_back((unsigned int)edep_Y_wire,(unsigned int)edep_tick);

  return target_roi_vec;
}

std::vector< std::pair<unsigned int, unsigned int> > 
sys::WireModifier::GetHitTargetROIs(recob::Hit const& hit)
{
  //vector of pairs of channel number, time tick
  std::vector< std::pair<unsigned int,unsigned int> > target_roi_vec;
  
  int hit_wire = hit.Channel();
  int hit_tick = int(round(hit.PeakTime()));
  
  if ( hit_tick<0 || hit_tick>=6400 )
    return target_roi_vec;

  target_roi_vec.emplace_back((unsigned int)hit_wire, (unsigned int)hit_tick);
  
  return target_roi_vec;
}

void sys::WireModifier::FillROIMatchedEdepMap(std::vector<sim::SimEnergyDeposit> const& edepVec,
					      std::vector<recob::Wire> const& wireVec)
{
  fROIMatchedEdepMap.clear();

  std::unordered_map<unsigned int,unsigned int> wireChannelMap;
  for(size_t i_w=0; i_w<wireVec.size(); ++i_w)
    wireChannelMap[wireVec[i_w].Channel()] = i_w;

  for(size_t i_e=0; i_e<edepVec.size(); ++i_e){

    std::vector< std::pair<unsigned int, unsigned int> > target_rois = GetTargetROIs(edepVec[i_e]); 

    for( auto const& target_roi : target_rois){

      //std::cout << "Matched edeps to channel " << target_roi.first << " time tick " << target_roi.second << std::endl;

      auto const& target_wire = wireVec.at(wireChannelMap[target_roi.first]);

      //if(target_roi.first!=target_wire.Channel())
      //throw std::runtime_error("ERROR! Channel ordering doesn't match wire ordering.");

      //std::cout << "\tGot wire " << target_wire.Channel() << std::endl;

      //std::cout << "\tWire has " << target_wire.SignalROI().n_ranges() << std::endl;

      if(target_wire.SignalROI().n_ranges()==0) continue;
      if(target_wire.SignalROI().is_void(target_roi.second)) continue;

      auto range_number = target_wire.SignalROI().find_range_iterator(target_roi.second) - target_wire.SignalROI().begin_range();

      fROIMatchedEdepMap[std::make_pair(target_wire.Channel(),range_number)].push_back(i_e);

    }//end loop over target rois

  }//end loop over all edeps

}

void sys::WireModifier::FillROIMatchedHitMap(std::vector<recob::Hit> const& hitVec,
					     std::vector<recob::Wire> const& wireVec)
{
  fROIMatchedHitMap.clear();
  
  std::unordered_map<unsigned int,unsigned int> wireChannelMap;
  for(size_t i_w=0; i_w<wireVec.size(); ++i_w)
    wireChannelMap[wireVec[i_w].Channel()] = i_w;
  
  for(size_t i_h=0; i_h<hitVec.size(); ++i_h){
    
    std::vector< std::pair<unsigned int, unsigned int> > target_rois = GetHitTargetROIs(hitVec[i_h]);
    
    for( auto const& target_roi : target_rois){
      
      auto const& target_wire = wireVec.at(wireChannelMap[target_roi.first]);
      
      if(target_wire.SignalROI().n_ranges()==0) continue;
      if(target_wire.SignalROI().is_void(target_roi.second)) continue;

      auto range_number = target_wire.SignalROI().find_range_iterator(target_roi.second) - target_wire.SignalROI().begin_range();
      
      fROIMatchedHitMap[std::make_pair(target_wire.Channel(),range_number)].push_back(i_h);
      
    }//end loop over target rois
    
  }//end loop over all hits
 
}

std::vector<sys::WireModifier::SubROIProperties_t>
sys::WireModifier::CalcSubROIProperties(sys::WireModifier::ROIProperties_t const& roi_properties, std::vector<const recob::Hit*> const& hitPtrVec) {

  std::vector<SubROIProperties_t> subroi_properties_vec;
  SubROIProperties_t subroi_properties;
  subroi_properties.plane = roi_properties.plane;

  // if this ROI doesn't contain any hits, define SubROI based on ROI properities
  if ( hitPtrVec.size() == 0 ) {
    subroi_properties.key     = std::make_pair( roi_properties.key, 0 );
    subroi_properties.total_q = roi_properties.total_q;
    subroi_properties.center  = roi_properties.center;
    subroi_properties.sigma   = roi_properties.sigma;
    subroi_properties_vec.push_back(subroi_properties);
  }

  // otherwise, define SubROIs based on hits
  else {
    for ( unsigned int i_h=0; i_h < hitPtrVec.size(); i_h++ ) {
      auto hit_ptr = hitPtrVec[i_h];
      subroi_properties.key     = std::make_pair( roi_properties.key, i_h );
      subroi_properties.total_q = hit_ptr->Integral();
      subroi_properties.center  = hit_ptr->PeakTime();
      subroi_properties.sigma   = hit_ptr->RMS();
      subroi_properties_vec.push_back(subroi_properties);
    }
  }

  return subroi_properties_vec;

}

std::map<sys::WireModifier::SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>>
sys::WireModifier::MatchEdepsToSubROIs(std::vector<sys::WireModifier::SubROIProperties_t> const& subROIPropVec,
				       std::vector<const sim::SimEnergyDeposit*> const& edepPtrVec) {

  std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> SubROIMatchedEdepMap;

  // dumbest possible thing for now...
  for ( auto sub_roi_prop : subROIPropVec ) {
    auto key = sub_roi_prop.key;
    for ( auto edep_ptr : edepPtrVec ) {
      auto edep_tick = A_t * edep_ptr->X() + C_t;
      if ( edep_tick > sub_roi_prop.center-sub_roi_prop.sigma && edep_tick < sub_roi_prop.center+sub_roi_prop.sigma ) 
	SubROIMatchedEdepMap[key].push_back(edep_ptr);
    }
  }

  return SubROIMatchedEdepMap;

}

/*
std::vector<const sim::SimEnergyDeposit*>
sys::WireModifier::DummyFunction(std::vector<const sim::SimEnergyDeposit*> const& edepPtrVec, std::vector<const recob::Hit*> const& hitPtrVec, 
				 sys::WireModifier::ROIProperties_t const& roi_properties) {

  std::vector<const sim::SimEnergyDeposit*> returnVec;

  int flag = 0;

  // for each TrackID, which EDeps are associated with it?
  std::map<int, std::vector<const sim::SimEnergyDeposit*>> TrackIDMatchedEDepMap;
  // total energy of EDeps matched to the ROI
  double total_energy = 0.;
  for ( auto edep_ptr : edepPtrVec ) {
    TrackIDMatchedEDepMap[edep_ptr->TrackID()].push_back(edep_ptr);
    total_energy += edep_ptr->E();
  }
  // calculate EDep properties by TrackID
  std::map<int, sys::WireModifier::TruthProperties_t> TrackIDMatchedPropertyMap;
  for ( auto const& track_edeps : TrackIDMatchedEDepMap ) { TrackIDMatchedPropertyMap[track_edeps.first] = CalcPropertiesFromEdeps(track_edeps.second); }  
  // TODO -- clean the above up after we decide what we need in production
  
  // TODO --  handle case where hitPtrVec.size() == 0
  if ( hitPtrVec.size() == 0 ) return returnVec;
  
  // do only ROIs with multiple hits for now
  if ( hitPtrVec.size() == 1 ) return returnVec;
  
  // for each EDep, which hit(s) (if any) is it plausibly matched to? based on whether the EDep's projected tick is within +/-1 RMS of the hit peak
  std::map<unsigned int, std::vector<unsigned int>> EDepMatchedHitMap;   // keys are indexes of edepPtrVec, values are vectors of indexes of hitPtrVec
  // for each TrackID, which hit(s) (if any) are its EDeps matched to? based on EDeps matched to exactly one hit
  std::map<int, std::unordered_set<unsigned int>> TrackIDMatchedHitMap;  // keys are TrackIDs, values are sets of indexes of hitPtrVec
  std::map<int, std::unordered_set<unsigned int>> TrackIDLooseMatchedHitMap; // based on all hits that all EDeps get matched to
  // for each hit, which EDep(s) (if any) are matched to it? based on multi-pass matching
  std::map<unsigned int, std::vector<unsigned int>> HitMatchedEDepMap;   // keys are indexes of hitPtrVec, values are vectors of indexes of edepPtrVec
  // for each hit, which TrackID(s) (if any) are matched? for each TrackID, how much energy is matched? based on multi-pass matching
  std::map<unsigned int, std::map<int, double>> HitMatchedTrackEnergyMap; // keys are indexes of hitPtrVec, values are maps of TrackIDs to matched energy (in MeV)

  // which EDeps have not yet been assigned?
  std::unordered_set<unsigned int> UnassignedEDeps;
  // for each pass, which EDeps were assigned in this pass?
  std::unordered_set<unsigned int> PassAssignedEDeps;
  

  // first pass: fill EDepMatchedHitMap, assign EDeps matched to exactly one hit
  //   also fill TrackIDMatchedHitMap based on EDeps matched to exactly one hit
  // loop over EDeps
  for ( unsigned int i_e=0; i_e < edepPtrVec.size(); i_e++ ) {
    
    // get EDep properties
    auto edep_ptr  = edepPtrVec[i_e];
    auto edep_tick = A_t * edep_ptr->X() + C_t;

    // loop over hits
    for ( unsigned int i_h=0; i_h < hitPtrVec.size(); i_h++ ) {
      auto hit_ptr = hitPtrVec[i_h];
      if ( edep_tick > hit_ptr->PeakTime()-hit_ptr->RMS() && edep_tick < hit_ptr->PeakTime()+hit_ptr->RMS() ) {
	EDepMatchedHitMap[i_e].push_back(i_h);
	TrackIDLooseMatchedHitMap[edep_ptr->TrackID()].emplace(i_h);
      }
    } // end loop over hits

    // if EDep is matched to exactly one hit, assign it to that hit => update HitMatchedEDepMap, HitMatchedTrackEnergyMap
    //   also fill TrackIDMatchedHitMap
    if ( EDepMatchedHitMap[i_e].size() == 1 ) {
      auto i_h = EDepMatchedHitMap[i_e][0];
      TrackIDMatchedHitMap[edep_ptr->TrackID()].emplace(i_h);
      HitMatchedEDepMap[i_h].push_back(i_e);
      HitMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()] += edep_ptr->E();
      //PassAssignedEDeps.emplace(i_e);
      returnVec.push_back(edep_ptr);
    }
    // else, add this EDep to the set of unassigned EDeps
    else UnassignedEDeps.emplace(i_e);

  } // end loop over EDeps
  // end first pass
  
  // if first pass assigned every EDep, return here
  //if ( UnassignedEDeps.size() == 0 ) return;

  return returnVec;

  //std::cout << "EDeps not assigned: ";
  //for ( auto i_e : UnassignedEDeps ) std::cout << i_e << " ";
  //std::cout << std::endl;
  
  
  // second pass: for EDeps that are matched to no hits...
  //                if EDep's TrackID is not matched to any hits, do not match it
  //                if EDep's TrackID is matched to exactly one hit, match it to that hit
  // reset PassAssignedEDeps
  PassAssignedEDeps.clear();
  // loop over unassigned EDeps
  for ( auto i_e :  UnassignedEDeps ) {
    //std::cout << "Doing second pass for EDep #" << i_e << "..." << std::endl;
    
    // if EDep is matched to multiple hits, skip it for this pass
    if ( EDepMatchedHitMap[i_e].size() != 0 ) continue;

    // get EDep pointer
    auto edep_ptr = edepPtrVec[i_e];

    // if EDep's TrackID is not matched to any hits, do not match it
    if ( TrackIDMatchedHitMap[edep_ptr->TrackID()].size() == 0 ) {
      //std::cout << "  Resloved by case 0-0" << std::endl;
      PassAssignedEDeps.emplace(i_e);
      continue;
    }

    // if EDep's TrackID is matched to exactly one hit, assign it to that hit
    if ( TrackIDMatchedHitMap[edep_ptr->TrackID()].size() == 1 ) { 
      //std::cout << "  Resolved by case 0-1" << std::endl;
      auto i_h = *TrackIDMatchedHitMap[edep_ptr->TrackID()].begin();

      // if this approach sems boarderline, print out some info
      if ( HitMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()] < 0.75*TrackIDMatchedPropertyMap[edep_ptr->TrackID()].total_energy ) {
	flag = 1;
	std::cout << "WARNING: assigning EDep #" << i_e << " to hit #" << i_h << ", but that hit contains "
		  << HitMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()] << " / " << TrackIDMatchedPropertyMap[edep_ptr->TrackID()].total_energy
		  << " = " << HitMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()]/TrackIDMatchedPropertyMap[edep_ptr->TrackID()].total_energy
		  << " < 0.75 of the TrackID's total energy (TrackID " << edep_ptr->TrackID() << ")" << std::endl;
      }

      HitMatchedEDepMap[i_h].push_back(i_e);
      HitMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()] += edep_ptr->E();
      PassAssignedEDeps.emplace(i_e);
    }
    
  } // end loop over EDeps
  // update UnassignedEDeps
  for ( auto i_e : PassAssignedEDeps ) UnassignedEDeps.erase(i_e);
  // end second pass

  // if second pass assigned every EDep, return here
  //if ( UnassignedEDeps.size() == 0 ) return;
  // if second pass was not sufficient, set flag = 1
  if ( UnassignedEDeps.size() != 0 ) flag = 1;

  // third pass: assign remaining EDeps to the closest hit, based on the hit RMS
  // reset PassAssignedEDeps
  PassAssignedEDeps.clear();
  // loop over unassigned EDeps
  for ( auto i_e :  UnassignedEDeps ) {
    //std::cout << "Doing third pass for EDep #" << i_e << "..." << std::endl;
    
    // get EDep properties
    auto edep_ptr = edepPtrVec[i_e];
    auto edep_tick = A_t * edep_ptr->X() + C_t;

    // get list of hits to consider assigning the EDep to
    std::vector<unsigned int> hit_list;
    if ( EDepMatchedHitMap[i_e].size() != 0 ) hit_list = EDepMatchedHitMap[i_e];
    else  { for ( auto i_h : TrackIDMatchedHitMap[edep_ptr->TrackID()] ) hit_list.push_back(i_h); }

    // get the hit that is closest to this EDep
    unsigned int closest_hit;
    float min_dist = 6400.;
    for ( auto i_h : hit_list ) {
      auto hit_ptr = hitPtrVec[i_h];
      float hit_dist = std::abs( edep_tick - hit_ptr->PeakTime() ) / hit_ptr->RMS();
      if ( hit_dist < min_dist ) {
	closest_hit = i_h;
	min_dist    = hit_dist;
      }
    }

    // assign the EDep to the closest hit
    HitMatchedEDepMap[closest_hit].push_back(i_e);
    HitMatchedTrackEnergyMap[closest_hit][edep_ptr->TrackID()] += edep_ptr->E();
    PassAssignedEDeps.emplace(i_e);
    
  } // end loop over EDeps
  // update UnassignedEDeps
  for ( auto i_e : PassAssignedEDeps ) UnassignedEDeps.erase(i_e);
  // end third pass

  //if ( UnassignedEDeps.size() == 0 ) return;
  if ( flag == 0 ) return returnVec;
  
  // ** INITIAL INFORMATIONAL COUTS **
  // print out ROI properties
  std::cout << "Channel " << roi_properties.key.first << ", signal ROI " << roi_properties.key.second << " (plane " << roi_properties.plane << ")" << std::endl;
  std::cout << "  ROI bounds: (" << roi_properties.begin << ", " << roi_properties.end << ")" << std::endl;
  std::cout << "  ROI center: " << roi_properties.center<< std::endl;
  std::cout << "  ROI sigma:  " << roi_properties.sigma << std::endl;
  std::cout << "  ROI charge: " << roi_properties.total_q << std::endl;
  
  // print out hit properties
  std::cout << "  Number of hits: " << hitPtrVec.size() << std::endl;
  float total_hit_charge = 0.;
  for ( auto hit_ptr : hitPtrVec ) total_hit_charge += hit_ptr->Integral();
  std::cout << "  Total hit charge: " << total_hit_charge << std::endl;
  for ( unsigned int i_h=0; i_h < hitPtrVec.size(); i_h++ ) { 
    auto hit_ptr = hitPtrVec[i_h];
    std::cout << "    For hit #" << i_h << ":" << std::endl;
    std::cout << "      Hit center: " << hit_ptr->PeakTime() << std::endl;
    std::cout << "      Hit width:  " << hit_ptr->RMS() << std::endl;
    std::cout << "      Hit charge: " << hit_ptr->Integral() << std::endl;
  }

  // print out global edep properties
  auto edep_col_properties = CalcPropertiesFromEdeps(edepPtrVec);
  //std::cout << "  Global EDep center: " << edep_col_properties.tick << std::endl;
  //std::cout << "  Global EDep width:  " << edep_col_properties.tick_rms << " (with charge-weighted RMS)" << std::endl;
  //std::cout << "  Global EDep width:  " << edep_col_properties.tick_rms_noWeight << " (with non-charge-weighted RMS)" << std::endl;
  std::cout << "  Number of contributing EDeps: " << edepPtrVec.size() << std::endl;
  std::cout << "  Number of contributing TrackIDs: " << TrackIDMatchedEDepMap.size() << std::endl;
  std::cout << "  Total energy: " << total_energy << "  MeV" << std::endl;

  // print out EDep properties by TrackID...
  // for each TrackID, print center, width, and total energy
  int trk_ctr = 0;
  for ( auto trk_edep_pair : TrackIDMatchedEDepMap ) {
    std::cout << "    For track #" << trk_ctr << ":" << std::endl;
    std::cout << "      TrackID:   " << trk_edep_pair.first << std::endl;
    std::cout << "      Num EDeps: " << trk_edep_pair.second.size() << std::endl;
    std::cout << "      PDG code:  " << trk_edep_pair.second[0]->PdgCode() << std::endl;
    edep_col_properties = CalcPropertiesFromEdeps(trk_edep_pair.second);
    std::cout << "      TrackID EDep center:    " << edep_col_properties.tick << std::endl;
    //std::cout << "      TrackID EDep width:     " << edep_col_properties.tick_rms << " (with charge-weighted RMS)" << std::endl;
    //std::cout << "      TrackID EDep width:     " << edep_col_properties.tick_rms_noWeight << " (with non-charge-weighted RMS)" << std::endl; 
    std::cout << "      TrackID EDep min tick:  " << edep_col_properties.tick_min << std::endl;
    std::cout << "      TrackID EDep max tick:  " << edep_col_properties.tick_max << std::endl;
    std::cout << "      TrackID total energy:   " << edep_col_properties.total_energy << " MeV" << std::endl;
    trk_ctr++;
  }

  // print out information on assigned EDeps
  std::cout << "  Info for assigned EDeps..." << std::endl;
  //for ( auto hit_trk_pair : HitMatchedTrackEnergyMap ) {
  for ( unsigned int i_h=0; i_h < hitPtrVec.size(); i_h++ ) {
    std::cout << "  Hit #" << std::setw(4) << i_h << ":";
    for ( auto trk_edep_pair : TrackIDMatchedEDepMap ) {
      std::cout << std::setw(10);
      auto trk = trk_edep_pair.first;
      auto it_trk = HitMatchedTrackEnergyMap[i_h].find(trk);
      if ( it_trk == HitMatchedTrackEnergyMap[i_h].end() ) std::cout << "--";
      else std::cout << it_trk->second;
    }
    std::cout << std::endl;
  }


  // print out information about unassigned EDeps
  std::cout << "  Number of unassigned EDeps: " << UnassignedEDeps.size() << std::endl;
  for ( auto i_e : UnassignedEDeps ) {

    // get EDep properties
    auto edep_ptr  = edepPtrVec[i_e];
    auto edep_tick = A_t * edep_ptr->X() + C_t;

    std::cout << "    For EDep #" << i_e << ":" << std::endl;
    std::cout << "      TrackID:         " << edep_ptr->TrackID() << std::endl;
    std::cout << "      EDep tick:       " << edep_tick << std::endl;
    std::cout << "      EDep energy:     " << edep_ptr->E() << " MeV" << std::endl;
    std::cout << "      Matched to hits: ";
    if ( EDepMatchedHitMap[i_e].size() == 0 ) std::cout << "none" << std::endl;
    else { 
      for ( auto i_h : EDepMatchedHitMap[i_e] ) std::cout << i_h << " ";
      std::cout << std::endl;
    }
    std::cout << "      TrackID matched to " << TrackIDMatchedHitMap[edep_ptr->TrackID()].size() << " hits..." << std::endl;
    for ( auto i_h : TrackIDMatchedHitMap[edep_ptr->TrackID()] ) {
      std::cout << "        Hit #" << i_h << " has " << HitMatchedTrackEnergyMap[i_h][edep_ptr->TrackID()] << " MeV of the TrackID's energy" << std::endl;
    }
  }

  


  // try out different distance metrics; if they disagree for any edep, print some info about that edep and all the ROI info
  int flag = 0;
  for ( unsigned int i_e=0; i_e < edepPtrVec.size(); i_e++ ) {

    // get the edep ptr and its tick
    auto edep_ptr  = edepPtrVec[i_e];
    auto edep_tick = A_t * edep_ptr->X() + C_t;

    // reset values
    float min_dist1 = 6400.;
    float min_dist2 = 6400.;
    float max_weight3 = -1.;
    int hit1 = -1;
    int hit2 = -1; 
    int hit3 = -1;
    std::vector<float> dist1_vec, dist2_vec, dist3_vec;

    // loop over hits
    for ( unsigned int i_h=0; i_h < hitPtrVec.size(); i_h++ ) {
      
      // get the hit ptr
      auto hit_ptr = hitPtrVec[i_h];

      // calculate hit using distance metric #1 -> distance in ticks
      float dist1 = std::abs( edep_tick - hit_ptr->PeakTime() );
      dist1_vec.push_back(dist1);
      if ( dist1 < min_dist1 ) {
	min_dist1 = dist1;
	hit1 = i_h;
      }

      // calculate hit using distance metric #2 -> distance in hit RMS
      float dist2 = std::abs( edep_tick - hit_ptr->PeakTime() ) / hit_ptr->RMS();
      dist2_vec.push_back(dist2);
      if ( dist2 < min_dist2 ){
	min_dist2 = dist2;
	hit2 = i_h;
      }
      
      // calculate hit using distance metric #3 -> "pull"
      float weight3 = GAUSSIAN( edep_tick, hit_ptr->PeakTime(), hit_ptr->RMS(), hit_ptr->PeakAmplitude() );
      dist3_vec.push_back(weight3);
      if ( weight3 > max_weight3 ) {
	max_weight3 = weight3;
	hit3 = i_h;
      }
      
    } // end loop over hits

    // if all metrics agree, be quiet
    if ( hit1 == hit2 && hit2 == hit3 && hit1 != -1 ) continue;

    // otherwise, flag this as interesting
    flag = 1;

    std::cout << "EDep #" << i_e << " (tick " << edep_tick << ", energy " << edep_ptr->E() << " MeV, trackID " << edep_ptr->TrackID() << "):" << std::endl;
    std::cout << "  dist1    ";
    for ( auto dist : dist1_vec) std::cout << "\t" << dist;
    std::cout << "\t -> hit #" << hit1 << std::endl;
    std::cout << "  dist2    ";
    for ( auto dist : dist2_vec) std::cout << "\t" << dist;
    std::cout << "\t -> hit #" << hit2 << std::endl;
    std::cout << "  wght3  ";
    for ( auto dist : dist3_vec) std::cout << "\t" << dist;
    std::cout << "\t -> hit #" << hit3 << std::endl;

  } // end loop over edeps

  // TODO -- merge elements of showers such that their energy depositions get counted together, recalculate these maps? not in this function...

}
*/

sys::WireModifier::TruthProperties_t 
sys::WireModifier::CalcPropertiesFromEdeps(std::vector<const sim::SimEnergyDeposit*> const& edepPtrVec){
  TruthProperties_t edep_col_properties;
  
  //do calculations here
  edep_col_properties.x = 0.;
  edep_col_properties.x_rms = 0.;
  edep_col_properties.x_rms_noWeight = 0.;
  edep_col_properties.x_min = 300.;
  edep_col_properties.x_max = -50.;

  double total_energy = 0.0;
  for(auto edep_ptr : edepPtrVec){
    edep_col_properties.x += edep_ptr->X()*edep_ptr->E();
    if ( edep_ptr->X() < edep_col_properties.x_min ) edep_col_properties.x_min = edep_ptr->X();
    if ( edep_ptr->X() > edep_col_properties.x_max ) edep_col_properties.x_max = edep_ptr->X();
    total_energy += edep_ptr->E();
  }

  if(total_energy>0.0)
    edep_col_properties.x = edep_col_properties.x/total_energy;

  for(auto edep_ptr : edepPtrVec) {
   edep_col_properties.x_rms += (edep_ptr->X()-edep_col_properties.x)*(edep_ptr->X()-edep_col_properties.x)*edep_ptr->E();
   edep_col_properties.x_rms_noWeight += (edep_ptr->X()-edep_col_properties.x)*(edep_ptr->X()-edep_col_properties.x);
  }

  edep_col_properties.x_rms_noWeight = std::sqrt(edep_col_properties.x_rms_noWeight);

  if(total_energy>0.0)
    edep_col_properties.x_rms = std::sqrt(edep_col_properties.x_rms/total_energy);

  // convert x-related proerties to ticks
  edep_col_properties.tick = A_t*edep_col_properties.x + C_t;
  edep_col_properties.tick_rms = A_t*edep_col_properties.x_rms;
  edep_col_properties.tick_rms_noWeight = A_t*edep_col_properties.x_rms_noWeight;
  edep_col_properties.tick_min = A_t*edep_col_properties.x_min + C_t;
  edep_col_properties.tick_max = A_t*edep_col_properties.x_max + C_t;

  edep_col_properties.total_energy = total_energy;

  return edep_col_properties;
}



sys::WireModifier::ScaleValues_t
sys::WireModifier::GetScaleValues(sys::WireModifier::TruthProperties_t const& truth_props, sys::WireModifier::ROIProperties_t const& roi_vals)
{
  ScaleValues_t scales;

  //get scales here
  if(roi_vals.plane==0){
    scales.r_Q = fTSplines_Charge_X[0]->Eval(truth_props.x);
    scales.r_sigma = fTSplines_Sigma_X[0]->Eval(truth_props.x);
  }
  else if(roi_vals.plane==1){
    scales.r_Q = fTSplines_Charge_X[1]->Eval(truth_props.x);
    scales.r_sigma = fTSplines_Sigma_X[1]->Eval(truth_props.x);
  }
  else if(roi_vals.plane==2){
    scales.r_Q = fTSplines_Charge_X[2]->Eval(truth_props.x);
    scales.r_sigma = fTSplines_Sigma_X[2]->Eval(truth_props.x);
  }

  //std::cout << "For plane=" << roi_vals.plane << " and x=" << truth_props.x
  //	    << " scales are " << scales.r_Q << " and " << scales.r_sigma << std::endl;

  return scales;
}

void sys::WireModifier::ModifyROI(std::vector<float> & roi_data,
				  sys::WireModifier::ROIProperties_t const& roi_prop,
				  std::vector<sys::WireModifier::SubROIProperties_t> const& subROIPropVec, 
				  std::map<sys::WireModifier::SubROI_Key_t, sys::WireModifier::ScaleValues_t> const& subROIScaleMap)
{

  std::cout << "Modifying wire " << roi_prop.key.first << ", " << roi_prop.key.second << std::endl;

  bool verbose=true;
  //if(roi_data.size()>100) verbose=true;

  double q_orig = 0.;
  double q_mod = 0.;
  double scale_ratio = 1.;

  for(size_t i_t=0; i_t < roi_data.size(); i_t++) {

    std::cout << "  Doing tick #" << i_t << "..." << std::endl;

    // reset q_orig, q_mod, scale_ratio
    q_orig = 0.;
    q_mod = 0.;
    scale_ratio = 1.;
    
    // calculate q_orig, q_mod for this tick
    for ( auto const& subroi_prop : subROIPropVec ) {
      auto scale_vals = subROIScaleMap.find(subroi_prop.key)->second;
      q_orig += GAUSSIAN( i_t+roi_prop.begin,
			  subroi_prop.center,
			  subroi_prop.sigma,
			  subroi_prop.total_q );
      std::cout << "    Incrementing q_mod by GAUSSIAN( " << i_t+roi_prop.begin << ", " << subroi_prop.center << ", " << scale_vals.r_sigma 
		<< " * " << subroi_prop.sigma << ", " << scale_vals.r_Q << " * " << subroi_prop.total_q << ") = " << GAUSSIAN( i_t+roi_prop.begin,
															       subroi_prop.center,
															       scale_vals.r_sigma * subroi_prop.sigma,
															       scale_vals.r_Q     * subroi_prop.total_q )
		<< std::endl;
      q_mod  += GAUSSIAN( i_t+roi_prop.begin,
			  subroi_prop.center,
			  scale_vals.r_sigma * subroi_prop.sigma, 
			  scale_vals.r_Q     * subroi_prop.total_q );
    }

    std::cout << "  q_orig = " << q_orig << ", q_mod = " << q_mod << std::endl;

    if(isnan(q_orig)) {
      std::cout << "WARNING: obtained q_orig = NaN... setting to zero" << std::endl;
      q_orig = 0.;
    }
    if(isnan(q_mod)) {
      std::cout << "WARNING: obtained q_mod = NaN... settign to zero" << std::endl;
      q_mod = 0.;
    }

    scale_ratio = q_mod / q_orig;
    std::cout << "  scale_ratio = " << scale_ratio << std::endl;
    
    if(isnan(scale_ratio)) {
      std::cout << "WARNING: obtained scale_ratio = " << q_mod << " / " << q_orig << " = NaN... setting to 1" << std::endl;
      scale_ratio = 1.;
    }
    if(isinf(scale_ratio)) {
      std::cout << "WARNING: obtained scale_ratio = " << q_mod << " / " << q_orig << " = inf... setting to 1" << std::endl;
      scale_ratio = 1.0;
    }
    
    roi_data[i_t] = scale_ratio * roi_data[i_t];

    if(verbose)
      std::cout << "\t\t\t tick " << i_t << ":"
		<< " data=" << roi_data[i_t]
		<< ", den. = " << q_orig
		<< ", num. = " << q_mod
		<< ", ratio=" << scale_ratio
		<< std::endl;
  }
  
  return;

}

sys::WireModifier::WireModifier(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fWireInputTag(p.get<art::InputTag>("WireInputTag")),
  fSimEdepShiftedInputTag(p.get<art::InputTag>("SimEdepShiftedInputTag")),
  fSimEdepOrigInputTag(p.get<art::InputTag>("SimEdepOrigInputTag")),
  fHitInputTag(p.get<art::InputTag>("HitInputTag")),
  fMakeRawDigitAssn(p.get<bool>("MakeRawDigitAssn",true)),
  fCopyRawDigitAssn(p.get<bool>("CopyRawDigitAssn",false)),
  fSplinesFileName(p.get<std::string>("SplinesFileName")),
  fSplineNames_Charge_X(p.get< std::vector<std::string> >("SplineNames_Charge_X")),
  fSplineNames_Sigma_X(p.get< std::vector<std::string> >("SplineNames_Sigma_X"))
  //  ONE_OVER_SQRT_2PI(1./std::sqrt(2*util::pi()))
{
  produces< std::vector< recob::Wire > >();

  if(fMakeRawDigitAssn)
    produces< art::Assns<raw::RawDigit,recob::Wire> >();


  fSplinesFileName = std::string(std::getenv("UBOONEDATA_DIR"))+"/systematics/det_sys/"+fSplinesFileName;
  std::cout << "Spline file is " << fSplinesFileName;

  TFile f_splines(fSplinesFileName.c_str(),"r");

  fTSplines_Charge_X.resize(fSplineNames_Charge_X.size());
  for(size_t i_s=0; i_s<fSplineNames_Charge_X.size(); ++i_s)
    f_splines.GetObject(fSplineNames_Charge_X[i_s].c_str(),fTSplines_Charge_X[i_s]);
  
  fTSplines_Sigma_X.resize(fSplineNames_Sigma_X.size());
  for(size_t i_s=0; i_s<fSplineNames_Sigma_X.size(); ++i_s)
    f_splines.GetObject(fSplineNames_Sigma_X[i_s].c_str(),fTSplines_Sigma_X[i_s]);
}

void sys::WireModifier::produce(art::Event& e)
{

  //get wires
  art::Handle< std::vector<recob::Wire> > wireHandle;
  e.getByLabel(fWireInputTag, wireHandle);
  auto const& wireVec(*wireHandle);

  //get association to rawdigit
  art::FindManyP<raw::RawDigit> digit_assn(wireHandle,e,fWireInputTag);

  //get sim edeps
  art::Handle< std::vector<sim::SimEnergyDeposit> > edepShiftedHandle;
  e.getByLabel(fSimEdepShiftedInputTag,edepShiftedHandle);
  auto const& edepShiftedVec(*edepShiftedHandle);
  
  art::Handle< std::vector<sim::SimEnergyDeposit> > edepOrigHandle;
  e.getByLabel(fSimEdepOrigInputTag,edepOrigHandle);
  auto const& edepOrigVec(*edepOrigHandle);

  // get hits
  art::Handle< std::vector<recob::Hit> > hitHandle;
  e.getByLabel(fHitInputTag, hitHandle);
  auto const& hitVec(*hitHandle);

  //output new wires and new associations
  std::unique_ptr< std::vector<recob::Wire> > new_wires(new std::vector<recob::Wire>());
  std::unique_ptr< art::Assns<raw::RawDigit,recob::Wire> > new_digit_assn(new art::Assns<raw::RawDigit,recob::Wire>());

  //first fill our roi to edep map
  FillROIMatchedEdepMap(edepShiftedVec,wireVec);
  // and fill our roi to hit map
  FillROIMatchedHitMap(hitVec,wireVec);

  //loop through the wires and rois per wire...
  for(size_t i_w=0; i_w<wireVec.size(); ++i_w){

    auto const& wire = wireVec[i_w];

    //make a new roi list
    recob::Wire::RegionsOfInterest_t new_rois;
    new_rois.resize(wire.SignalROI().size());

    unsigned int my_plane=wire.View();

    for(size_t i_r=0; i_r<wire.SignalROI().get_ranges().size(); ++i_r){


      auto const& range = wire.SignalROI().get_ranges()[i_r];
      ROI_Key_t roi_key(wire.Channel(),i_r);

      //get the matching edeps
      auto it_map = fROIMatchedEdepMap.find(roi_key);
      if(it_map==fROIMatchedEdepMap.end()){
	new_rois.add_range(range.begin_index(),range.data());
	continue;
      }
      std::vector<size_t> matchedEdepIdxVec = it_map->second;      
      if(matchedEdepIdxVec.size()==0){
	new_rois.add_range(range.begin_index(),range.data());
	continue;
      }
      std::vector<const sim::SimEnergyDeposit*> matchedEdepPtrVec;
      std::vector<const sim::SimEnergyDeposit*> matchedShiftedEdepPtrVec;
      for(auto i_e : matchedEdepIdxVec) {
	matchedEdepPtrVec.push_back(&edepOrigVec[i_e]);
	matchedShiftedEdepPtrVec.push_back(&edepShiftedVec[i_e]);
      }      

      // get the matching hits
      std::vector<const recob::Hit*> matchedHitPtrVec;
      auto it_hit_map = fROIMatchedHitMap.find(roi_key);
      if( it_hit_map != fROIMatchedHitMap.end() ) {
	for( auto i_h : it_hit_map->second ) {
	  matchedHitPtrVec.push_back(&hitVec[i_h]);
	}
      }

      //calc roi properties
      auto roi_properties = CalcROIProperties(range);
      roi_properties.key   = roi_key;
      roi_properties.plane = my_plane;

      // get the subROIs
      auto subROIPropVec = CalcSubROIProperties(roi_properties, matchedHitPtrVec);

      // get the edeps per subROI
      auto SubROIMatchedShifedEdepMap = MatchEdepsToSubROIs(subROIPropVec, matchedShiftedEdepPtrVec);
      // convert from shifted edep pointers to edep indexes
      std::map<SubROI_Key_t, std::vector<size_t>> SubROIMatchedEdepIdxMap;
      for ( auto const& key_edepPtrVec_pair : SubROIMatchedShifedEdepMap ) {
	auto key = key_edepPtrVec_pair.first;
	for ( auto i_e : matchedEdepIdxVec ) {
	  for ( auto const& edep_ptr : key_edepPtrVec_pair.second ) {
	    if ( matchedShiftedEdepPtrVec[i_e] == edep_ptr ) 
	      SubROIMatchedEdepIdxMap[key].push_back(i_e);
	  }
	}
      }
      // convert from edep indexes to original edep pointers
      //   maybe collapse with above and skip over edep indexes?
      std::map<SubROI_Key_t, std::vector<const sim::SimEnergyDeposit*>> SubROIMatchedEdepMap;
      for ( auto const& key_edepIdxVec_pair : SubROIMatchedEdepIdxMap ) {
	auto key = key_edepIdxVec_pair.first;
	for ( auto i_e : key_edepIdxVec_pair.second ) {
	  SubROIMatchedEdepMap[key].push_back(matchedEdepPtrVec[i_e]);
	}
      }

      //get the scaling values
      //auto edep_col_properties = CalcPropertiesFromEdeps(matchedEdepPtrVec);
      //auto scales = GetScaleValues(edep_col_properties,roi_properties);
      std::map<SubROI_Key_t, TruthProperties_t> SubROIMatchedTruthMap;
      std::map<SubROI_Key_t, ScaleValues_t>     SubROIMatchedScalesMap;
      for( auto const& key_edepPtrVec_pair : SubROIMatchedEdepMap ) {
	auto key = key_edepPtrVec_pair.first;
	ScaleValues_t scale_vals;
	if ( key_edepPtrVec_pair.second.size() == 0 ){
	  scale_vals.r_Q     = 1.;
	  scale_vals.r_sigma = 1.;
	}
	else {
	  auto truth_vals = CalcPropertiesFromEdeps(key_edepPtrVec_pair.second);
	  SubROIMatchedTruthMap[key]  = truth_vals;
	  scale_vals = GetScaleValues(truth_vals, roi_properties);
	}
	SubROIMatchedScalesMap[key] = scale_vals;
      }

      //get modified ROI given scales
      std::vector<float> modified_data(range.data());
      ModifyROI(modified_data, roi_properties, subROIPropVec, SubROIMatchedScalesMap);
      
      new_rois.add_range(roi_properties.begin,modified_data);
      
    }//end loop over rois
    
    //make our new wire object
    new_wires->emplace_back(new_rois,wire.Channel(),wire.View());

    
    //get the associated rawdigit
    if(fCopyRawDigitAssn){
      auto const& rd_ptrs = digit_assn.at(i_w);
      for(auto const& rd_ptr : rd_ptrs)
	util::CreateAssn(*this,e,*new_wires,rd_ptr,*new_digit_assn,new_wires->size()-1);
    }

  }//end loop over wires

  e.put(std::move(new_wires));

  if(fMakeRawDigitAssn)
    e.put(std::move(new_digit_assn));
  
  //get a list of ptrs to the matched edeps for a wire
}
DEFINE_ART_MODULE(sys::WireModifier)

