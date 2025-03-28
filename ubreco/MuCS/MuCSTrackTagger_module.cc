////////////////////////////////////////////////////////////////////////
/// \file  MuCSTagger_module.cc
/// \brief EDProducer for tagging tracks as MuCS.
///
/// \version $Id: MuCSTagger_module.cxx
/// \author  Matthew.Bass@physics.ox.ac.uk && 
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "larcorealg/Geometry/geo.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/WireReadout.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"

#include "TVector3.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"


class MuCSTrackTagger : public art::EDProducer {
public:
  explicit MuCSTrackTagger(fhicl::ParameterSet const & p);
  virtual ~MuCSTrackTagger();

  void produce(art::Event & e) override;

  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;


private:

  
  bool intersectsBoxes(const TVector3 & start, const TVector3& dir);

  std::string fTrackModuleLabel; //< Track label to find MuCS tags in
  std::string fFlashModuleLabel; //< flash producer to be used
  std::vector<float> fMuCSTopBox, fMuCSBottomBox; //< Box Edge Positions (x1,x2,y1,y2,z1,z2)
  float fBoxExtension; //< Amount to extend acceptance for box interception [cm]
  unsigned int fDirFromNPoints; //< Number of points to use to determine track direction (0=use track end direction)
  float fMinTrackLength; //< Minimum length of track to consider [cm]
  bool fUsePMTCalib; //< whether to use and apply PMT calibrations
  bool fDebug; //< debug flag
  
  //hists
  TH2F* fTopBoxPosHist;
  TH2F* fBottomBoxPosHist;

  TTree* _tree;
  
  int _run, _sub, _evt;
  int _ntag;
  float _trk_start_x, _trk_start_y, _trk_start_z;
  float _trk_end_x, _trk_end_y, _trk_end_z;
  float _trk_dir_x, _trk_dir_y, _trk_dir_z;
  float _trk_len;
  std::vector<float> _trk_x_v, _trk_y_v, _trk_z_v;
  int _nflash;
  float _flash_t, _flash_pe;
  float _flash_z, _flash_zw, _flash_y, _flash_yw;
  std::vector<double> _gain_v;
  std::vector<double> _flash_pe_v, _flash_pe_v_corr;
  

};

MuCSTrackTagger::MuCSTrackTagger(fhicl::ParameterSet const & p)
: EDProducer(p) {
  this->reconfigure(p);
  // Call appropriate Produces<>() functions here.
  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();
  //produces< art::Assns<recob::Track, recob::OpFlash> >();
}

MuCSTrackTagger::~MuCSTrackTagger() {}

void MuCSTrackTagger::beginJob(){
  art::ServiceHandle<art::TFileService> tfs;
  
  fTopBoxPosHist     = tfs->make<TH2F>("topboxpos","TopBoxPositions;X;Z",100, 
                                  0.9*(fMuCSTopBox[0]-fBoxExtension), 1.1*(fMuCSTopBox[1]+fBoxExtension)
                                  ,100, 0.9*(fMuCSTopBox[4]-fBoxExtension), 1.1*(fMuCSTopBox[5]+fBoxExtension));

  fBottomBoxPosHist  = tfs->make<TH2F>("bottomboxpos","BottomBoxPositions;X;Z",100, 
                                  0.9*(fMuCSBottomBox[0]-fBoxExtension), 1.1*(fMuCSBottomBox[1]+fBoxExtension)
                                  ,100, 0.9*(fMuCSBottomBox[4]-fBoxExtension), 1.1*(fMuCSBottomBox[5]+fBoxExtension));

  _tree = tfs->make<TTree>("tree","MuCS tagged tracks");
  _tree->Branch("_run",&_run,"run/I");
  _tree->Branch("_sub",&_sub,"sub/I");
  _tree->Branch("_evt",&_evt,"evt/I");
  _tree->Branch("_ntag",&_ntag,"ntag/I");
  _tree->Branch("_trk_len",&_trk_len,"_trk_len/F");
  _tree->Branch("_trk_start_x",&_trk_start_x,"_trk_start_x/F");
  _tree->Branch("_trk_start_y",&_trk_start_y,"_trk_start_y/F");
  _tree->Branch("_trk_start_z",&_trk_start_z,"_trk_start_z/F");
  _tree->Branch("_trk_end_x",&_trk_end_x,"_trk_end_x/F");
  _tree->Branch("_trk_end_y",&_trk_end_y,"_trk_end_y/F");
  _tree->Branch("_trk_end_z",&_trk_end_z,"_trk_end_z/F");
  _tree->Branch("_trk_dir_x",&_trk_dir_x,"_trk_dir_x/F");
  _tree->Branch("_trk_dir_y",&_trk_dir_y,"_trk_dir_y/F");
  _tree->Branch("_trk_dir_z",&_trk_dir_z,"_trk_dir_z/F");
  _tree->Branch("_trk_x_v","std::vector<float>",&_trk_x_v);
  _tree->Branch("_trk_y_v","std::vector<float>",&_trk_y_v);
  _tree->Branch("_trk_z_v","std::vector<float>",&_trk_z_v);
  _tree->Branch("_nflash",&_nflash,"nflash/I");
  _tree->Branch("_flash_pe",&_flash_pe,"_flash_pe/F");
  _tree->Branch("_flash_t",&_flash_t,"_flash_t/F");
  _tree->Branch("_flash_z",&_flash_z,"_flash_z/F");
  _tree->Branch("_flash_y",&_flash_y,"_flash_y/F");
  _tree->Branch("_flash_zw",&_flash_zw,"_flash_zw/F");
  _tree->Branch("_flash_yw",&_flash_yw,"_flash_yw/F");
  _tree->Branch("_gain_v","std::vector<double>",&_gain_v);
  _tree->Branch("_flash_pe_v","std::vector<double>",&_flash_pe_v);
  _tree->Branch("_flash_pe_v_corr","std::vector<double>",&_flash_pe_v_corr);

                                                                    
}


bool MuCSTrackTagger::intersectsBoxes(const TVector3 & start, const TVector3& dir){
  //return true if this trajector will intersect both boxes
  
  TVector3 newpTop, newpBottom;
  newpTop.SetXYZ(start.X() + (fMuCSTopBox[2]-start.Y())*dir.X()/dir.Y(),
                 fMuCSTopBox[2],
                 start.Z() + (fMuCSTopBox[2]-start.Y())*dir.Z()/dir.Y());
  newpBottom.SetXYZ(start.X() + (fMuCSBottomBox[2]-start.Y())*dir.X()/dir.Y(),
                 fMuCSBottomBox[2],
                 start.Z() + (fMuCSBottomBox[2]-start.Y())*dir.Z()/dir.Y());  
                 
  //populate hists regardless of intersection
  fTopBoxPosHist->Fill(newpTop.X(),newpTop.Z());
  fBottomBoxPosHist->Fill(newpBottom.X(),newpBottom.Z());

  
  if(newpTop.X() > fMuCSTopBox[0]-fBoxExtension && newpTop.X() < fMuCSTopBox[1]+fBoxExtension
     && newpTop.Z() > fMuCSTopBox[4]-fBoxExtension && newpTop.Z() < fMuCSTopBox[5]+fBoxExtension
     && newpBottom.X() > fMuCSBottomBox[0]-fBoxExtension && newpBottom.X() < fMuCSBottomBox[1]+fBoxExtension
     && newpBottom.Z() > fMuCSBottomBox[4]-fBoxExtension && newpBottom.Z() < fMuCSBottomBox[5]+fBoxExtension)
   return true;
   else
    return false;
}

void MuCSTrackTagger::produce(art::Event & e) {
  // Implementation of required member function here.

  std::unique_ptr< std::vector< anab::CosmicTag > >              cosmicTagTrackVector ( new std::vector<anab::CosmicTag> );
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > >  assnOutCosmicTagTrack( new art::Assns<recob::Track, anab::CosmicTag> );
  //std::unique_ptr< art::Assns<recob::Track, recob::OpFlash  > >  assnOutOpFlashTrack  ( new art::Assns<recob::Track, recob::OpFlash>  );

  _run = e.run();
  _sub = e.subRun();
  _evt = e.event();
  _ntag = 0;

  _trk_x_v.clear();
  _trk_y_v.clear();
  _trk_z_v.clear();

  art::Handle<std::vector<recob::Track> > Trk_h;
  e.getByLabel( fTrackModuleLabel, Trk_h );
  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, Trk_h);

  art::Handle<std::vector<recob::OpFlash> > Flash_h;
  e.getByLabel( fFlashModuleLabel, Flash_h );
  std::vector<art::Ptr<recob::OpFlash> > FlashVec;
  art::fill_ptr_vector(FlashVec, Flash_h);

  
  // grab PMT calibrations if needed
  auto const& channelMapAlg = art::ServiceHandle<geo::WireReadout const>()->Get();
  const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
  /*
  for (unsigned int i=0; i!= geo->NOpDets(); ++i) {
    if (geo->IsValidOpChannel(i) && i<32) {
      std::cout<<"Channel "<<i <<" "<<gain_provider.Gain(i)<<" "
               <<gain_provider.GainErr(i) <<" "
               <<gain_provider.ExtraInfo(i).GetStringData("gain_fit_status")<<" "
               <<gain_provider.ExtraInfo(i).GetFloatData("amplitude_gain")<<" "
               <<gain_provider.ExtraInfo(i).GetFloatData("amplitude_gain_err")<<" "
               <<gain_provider.ExtraInfo(i).GetStringData("amplitude_gain_fit_status")<<" "
               <<gain_provider.ExtraInfo(i).GetFloatData("pedestal_mean")<<" "
               <<gain_provider.ExtraInfo(i).GetFloatData("pedestal_mean_err")<<" "
               <<gain_provider.ExtraInfo(i).GetFloatData("pedestal_rms")<<" "
               <<gain_provider.ExtraInfo(i).GetFloatData("pedestal_rms_err")<<std::endl;
    }
    else if (geo->IsValidOpChannel(i)) {
      std::cout<<"Channel "<<i<<std::endl;
    }
   }
  */
  
  // find flash in time with MuCS
  
  _nflash = 0;
  //size_t whichflash = 0; // which flash is in time? this is needed to store association
  //size_t flash_ctr = 0;

  for (auto flash: FlashVec) {

    if (fDebug) { std::cout << "Found flash @ time " << flash->Time() << " with " << flash->TotalPE() << " PEs" << std::endl; }
    
    if ( (flash->Time() > -1.5) && (flash->Time() < -0.5) ) {

      if (fDebug) { std::cout << "\t Tagged!" << std::endl; }
      
      _nflash += 1;
      _flash_t    = flash->Time();
      _flash_pe_v = flash->PEs();
      std::cout << "there are " << _flash_pe_v.size() << " entries in flash vector " << std::endl;
      for (size_t p=0; p < _flash_pe_v.size(); p++)
	std::cout << "entry " << p << " has " << _flash_pe_v[p] << " PE" << std::endl;
      _flash_pe_v_corr = _flash_pe_v;
      _gain_v = std::vector<double>(_flash_pe_v.size(),0.);
      // apply gain
      if (fUsePMTCalib) {
	for (size_t ch=0; ch < _flash_pe_v.size(); ch++) {
	  if (ch%100 < 32) {
            if (channelMapAlg.IsValidOpChannel(ch%100) == true) {
	      auto chgain = gain_provider.ExtraInfo(ch%100).GetFloatData("amplitude_gain");
	      std::cout << "Gain for channel " << ch%100 << " is " << chgain << std::endl;
	      _gain_v.at(ch) = chgain;
	      _flash_pe_v_corr.at(ch) *= ( 20. / chgain ); 
	    }
	  }
	}// if channel is valid
      }// if we are using the PMT gains
      _flash_pe   = flash->TotalPE();
      _flash_z    = flash->ZCenter();
      _flash_y    = flash->YCenter();
      _flash_zw   = flash->ZWidth();
      _flash_yw   = flash->YWidth();
      
      //whichflash = flash_ctr;

    }// if flash in time

    //flash_ctr+= 1;

  }// for all flashes

  for (auto trk: TrkVec){
    if((trk->Length())<fMinTrackLength) continue;

    std::cout << std::endl << "New track" << std::endl;
    //choose highest edge as track start
    TVector3 start, end, startDir, endDir;
    if(trk->Vertex().Y()>trk->End().Y()){
      start=trk->Vertex<TVector3>();
      end=trk->End<TVector3>();
      startDir=trk->VertexDirection<TVector3>();
      endDir=trk->EndDirection<TVector3>();
    }else{
      start=trk->End<TVector3>();
      end=trk->Vertex<TVector3>();
      startDir=trk->EndDirection<TVector3>();
      endDir=trk->VertexDirection<TVector3>();
    }
    
    //find which end of the trajectory to use to get direction
    unsigned int pStart;
    TVector3 dir;
    int pSign;
    if(trk->LocationAtPoint<TVector3>(trk->FirstValidPoint())==start){
      pStart=0;
      pSign=1; //go forward for track direction
    }else if(trk->LocationAtPoint<TVector3>(trk->LastValidPoint())==start){
      pStart=trk->LastValidPoint();
      pSign=-1; //go backward for track direction
    }else{
      throw cet::exception("MuCSTrackTagger") << "Start seems to be in wrong position!\n";
    }
    
    if(fDirFromNPoints==0){
      //use reversed track start direction
      dir=-startDir;
    }else{
      //use diff between pstart and pstart+psign*(fDirFromNPoints-1)
      if(fDirFromNPoints>trk->CountValidPoints())
	mf::LogInfo("MuCSTrackTagger") << "Track has too few trajectory points ("<<trk->CountValidPoints()<<"), skipping it.\n";
      dir=(trk->LocationAtPoint<TVector3>(pStart) - trk->LocationAtPoint<TVector3>(trk->NextValidPoint(pStart+pSign*(fDirFromNPoints-1) ) ) ).Unit();
    }

    std::cout << "start @ [" << start.X() << ", " << start.Y() << ", " << start.Z() << "]" 
	      << "\t w dir   [" << dir.X()   << ", " << dir.Y()   << ", " << dir.Z()   << "]"  << std::endl;


    //find interesections and generate tags if appropriate
    bool btag=intersectsBoxes(start,dir);

    if (btag){
      std::cout << "\t ***** this track intersects boxes!   *****" << std::endl;
      cosmicTagTrackVector->emplace_back(-999.);
      util::CreateAssn(*this, e, *cosmicTagTrackVector, trk, *assnOutCosmicTagTrack );
      //util::CreateAssn(*this, e, *cosmicTagTrackVector, trk, *assnOutCosmicTagTrack );

      // save info to TTree
      _ntag += 1;

      _trk_len = trk->Length();
      
      _trk_start_x = start.X();
      _trk_start_y = start.Y();
      _trk_start_z = start.Z();

      _trk_end_x = end.X();
      _trk_end_y = end.Y();
      _trk_end_z = end.Z();

      _trk_dir_x = dir.X();
      _trk_dir_y = dir.Y();
      _trk_dir_z = dir.Z();

      _trk_x_v.clear();
      _trk_y_v.clear();
      _trk_z_v.clear();

      for (size_t p=0; p < trk->CountValidPoints(); p++) {
	auto pt = trk->LocationAtPoint( trk->NextValidPoint(pStart+pSign*p) );
	_trk_x_v.push_back( pt.X() );
	_trk_y_v.push_back( pt.Y() );
	_trk_z_v.push_back( pt.Z() );
      }// for all trajectory points

    }// if tagged
    
  }// for all tracks

  _tree->Fill();
 
  // e.put( std::move(outTracksForTags) );
  e.put( std::move(cosmicTagTrackVector) );
  e.put( std::move(assnOutCosmicTagTrack) );


} // end of produce


void MuCSTrackTagger::reconfigure(fhicl::ParameterSet const & p) {

  fTrackModuleLabel = p.get< std::string >("TrackModuleLabel", "track");
  fFlashModuleLabel = p.get< std::string >("FlashModuleLabel", "simpleFlashCosmic");

  fMuCSBottomBox=p.get< std::vector< float > >("MuCSBottomBox");
  fMuCSTopBox=p.get< std::  vector< float > >("MuCSTopBox");
  if(fMuCSBottomBox.size()!=6 || fMuCSTopBox.size()!=6)
    throw cet::exception("MuCSTrackTagger") << "MuCSBottomBox or MuCSTopBox has wrong size!\n";
  
  fBoxExtension=p.get<float>("BoxExtension",0.);
  fDirFromNPoints=p.get<unsigned int>("DirFromNPoints",0);
  fMinTrackLength=p.get<float>("MinTrackLength",0.);

  fUsePMTCalib=p.get<bool>("UsePMTCalib",false);

  fDebug=p.get<bool>("Debug",false);
  
}

DEFINE_ART_MODULE(MuCSTrackTagger)
