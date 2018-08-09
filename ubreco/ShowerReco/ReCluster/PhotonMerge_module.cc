////////////////////////////////////////////////////////////////////////
// Class:       PhotonMerge
// Plugin Type: producer (art v2_09_06)
// File:        PhotonMerge_module.cc
//
// Generated at Wed Mar 21 07:50:29 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
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

// larsoft data-products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardata/Utilities/FindManyInChainP.h"

#include "ubreco/ShowerReco/TwoDimTools/Linearity.h"
#include "ubreco/ShowerReco/TwoDimTools/Poly2D.h"

#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/Cluster.h"
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/ClusterMaker.h"

#include "art/Persistency/Common/PtrMaker.h"

class PhotonMerge;


class PhotonMerge : public art::EDProducer {
public:
  explicit PhotonMerge(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonMerge(PhotonMerge const &) = delete;
  PhotonMerge(PhotonMerge &&) = delete;
  PhotonMerge & operator = (PhotonMerge const &) = delete;
  PhotonMerge & operator = (PhotonMerge &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string fShrProducer, fVtxProducer, fPhotonProducer;

  bool fDebug;

  double _wire2cm, _time2cm, _trigoff;
  
  double _vtxW, _vtxT;

  // shower width (opening angle, in degrees)
  double fWidth; 
  // shower length (cm)
  double fShrLen;
  // maximum fraction of shower charge to be added by re-clustering
  double fFracShrQ;
  // maximum slope difference between shower cluster and to-be merged cluster
  double fMaxSlopeAngle;

  ::cluster::ClusterMaker _clusterMaker;

  // vector of all collection-plane hits for all showers in the event
  std::vector< art::Ptr<recob::Hit> > _allshr_hit_v;

  // map connecting photon cluster index to linearity object
  std::map< size_t, twodimtools::Linearity > _photon_lin_map; 
  // map connecting photon cluster index to poly2d object
  std::map< size_t, twodimtools::Poly2D > _photon_poly_map; 

  twodimtools::Poly2D projectShower(const art::Ptr<recob::Cluster> clus);

  double PhotonShowerAngle(const twodimtools::Poly2D& shr,
			   const twodimtools::Linearity& photon);

  double slopeCompat(const twodimtools::Poly2D& shr,
		     const twodimtools::Linearity& photon);

  bool photonCrossesShower(const twodimtools::Poly2D& shr,
			   const twodimtools::Poly2D& photon);

  // return number of hits actually merged
  int MergeHits(std::vector< art::Ptr<recob::Hit> >& shrhits, const std::vector< art::Ptr<recob::Hit> >& gammahits);
  
};


PhotonMerge::PhotonMerge(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces<std::vector<recob::PFParticle> >();
  produces<std::vector<recob::Cluster>    >();
  produces<std::vector<recob::Hit>        >();
  produces<art::Assns <recob::PFParticle, recob::Cluster>    >();
  produces<art::Assns <recob::Cluster   , recob::Hit>        >();

  fShrProducer    = p.get<std::string>("ShrProducer");
  fVtxProducer    = p.get<std::string>("VtxProducer");
  fPhotonProducer = p.get<std::string>("PhotonProducer");
  fWidth          = p.get<double>     ("Width");
  fShrLen         = p.get<double>     ("ShrLen");
  fFracShrQ       = p.get<double>     ("FracShrQ");
  fMaxSlopeAngle  = p.get<double>     ("MaxSlopeAngle");
  fDebug          = p.get<bool>       ("Debug");
  
  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
  _trigoff = detp->TriggerOffset();
}

void PhotonMerge::produce(art::Event & e)
{

  // produce recob::PFParticles
  std::unique_ptr< std::vector<recob::PFParticle> > PFParticle_v(new std::vector<recob::PFParticle> );
  std::unique_ptr< std::vector<recob::Cluster>    > Cluster_v   (new std::vector<recob::Cluster>    );
  std::unique_ptr< std::vector<recob::Hit>        > Hit_v       (new std::vector<recob::Hit>        );
  std::unique_ptr< art::Assns <recob::PFParticle, recob::Cluster> > PFParticle_Cluster_assn_v( new art::Assns<recob::PFParticle, recob::Cluster>   );
  std::unique_ptr< art::Assns <recob::Cluster  , recob::Hit>      > Cluster_Hit_assn_v       ( new art::Assns<recob::Cluster   , recob::Hit>       );

  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>(fShrProducer);
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>(fVtxProducer);
  // load input photon clusters
  auto const& photon_h = e.getValidHandle<std::vector<recob::Cluster>>(fPhotonProducer);

  // grab pfparticles associated with shower
  art::FindManyP<recob::PFParticle> shr_pfp_assn_v(shr_h, e, fShrProducer);
  
  // grab clusters associated with shower
  art::FindManyP<recob::Cluster> shr_clus_assn_v(shr_h, e, fShrProducer);

  // grab hits associated with shower
  art::FindManyP<recob::Hit> shr_hit_assn_v(shr_h, e, fShrProducer);

  // grab hits associated associated with photon clusters
  art::FindManyP<recob::Hit> photon_hit_assn_v(photon_h, e, fPhotonProducer);


  // grab the hits associated to the showers
  //auto shr_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(shr_h, e, fShrProducer);

  // pfparticle pointer maker for later to create associations
  art::PtrMaker<recob::PFParticle> PFPartPtrMaker(e, *this);
  // cluster pointer maker for later to create associations
  art::PtrMaker<recob::Cluster>    ClusPtrMaker  (e, *this);
  // hit pointer maker for later to create associations
  art::PtrMaker<recob::Hit>        HitPtrMaker   (e, *this);

  _photon_poly_map.clear();
  _photon_lin_map.clear();


  if (vtx_h->size() != 1)
    std::cout << "no vertex" << std::endl;
  if (photon_h->size() == 0)
    std::cout << "no photons" << std::endl;

  _clusterMaker.loadVertex(vtx_h);

  // load vertex and project on collection-plane
  auto const vtx = vtx_h->at(0);
  Double_t xyz[3] = {};
  vtx.XYZ(xyz);
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  _vtxW = geom->WireCoordinate(xyz[1],xyz[2],geo::PlaneID(0,0,2)) * _wire2cm;
  _vtxT = xyz[0];
  if (fDebug) std::cout << "\n\t Vertex @ Pl 2 : [w,t] -> [" << _vtxW << ", " << _vtxT << "]" << std::endl;

  // create polygon objects for each photon cluster.
  for (size_t p=0; p < photon_h->size(); p++){

    // get assciated hits
    std::vector< art::Ptr<recob::Hit> > photon_hit_v = photon_hit_assn_v.at(p);

    if (photon_hit_v.size() == 0) continue;

    auto plane = photon_hit_v.at(0)->WireID().Plane;

    // we are re-clustering only on the collection-plane
    if (plane != 2) continue;

    if (fDebug) { std::cout << "... new photon" << std::endl; }
    
    // create linearity objects   
    std::vector<double> hit_w_v;
    std::vector<double> hit_t_v;

    // get polygon   
    twodimtools::Poly2D clusPoly(photon_hit_v);
    if (clusPoly.Size() == 0) continue;

    for (auto hitptr : photon_hit_v){
      hit_w_v.push_back( hitptr->WireID().Wire * _wire2cm );
      hit_t_v.push_back( (hitptr->PeakTime() - _trigoff) * _time2cm + 0.6);
    }// for all hits in photon cluster                          

    twodimtools::Linearity clusLin(hit_w_v, hit_t_v);

    _photon_lin_map [ p ] = clusLin;
    _photon_poly_map[ p ] = clusPoly;

  }// for all clusters

  // get full list of all hits associated to all showers on the collection plane
  // these cannot be added to any cluster
  _allshr_hit_v.clear();
  for (size_t s=0; s < shr_h->size(); s++) {
    std::vector< art::Ptr<recob::Hit> > shr_hit_v = shr_hit_assn_v.at(s);
    for (auto hitPtr : shr_hit_v) {
      if (hitPtr->WireID().Plane == 2) { _allshr_hit_v.push_back( hitPtr ); }
    }// for all hits
  }// for all showers

  // save photons to be added to each shower if found compatible
  // first index is photon cluster index
  // pair is < shower index, angle compatibility >
  // if more then one shower the one with the best angle compatibility is chosen
  std::map< size_t, std::vector< std::pair<size_t, double> > > Photon_Shower_Map;
  
  // loop through reconstructed showers.                                         
  for (size_t s=0; s < shr_h->size(); s++) {
    
    //auto const& shr = shr_h->at(s);    
    if (fDebug) { std::cout << "new shower " << s << " of " << shr_h->size()-1 << std::endl; }
    
    // grab collection-plane hits and cluster associated with this shower
    std::vector< art::Ptr<recob::Cluster> > shr_clus_v = shr_clus_assn_v.at(s);
    size_t collidx = 0; // cluster index associated to collection-plane
    for (size_t c=0; c < shr_clus_v.size(); c++) {
      if (shr_clus_v.at(c)->Plane().Plane == 2) {
	collidx = c;
	break;
      }
    }
    std::vector< art::Ptr<recob::Hit> > shr_hit_v = shr_hit_assn_v.at(s);
    std::vector< std::vector< art::Ptr<recob::Hit> > > shr_hit_plv_v;
    shr_hit_plv_v.resize(3);
    for (auto const& shr_hit : shr_hit_v)
      shr_hit_plv_v.at(shr_hit->WireID().Plane).push_back( shr_hit ); 
    
    auto shrPoly = projectShower(shr_clus_v.at(collidx));
    
    if (fDebug) { std::cout << "hits on Y plane before merging : " << shr_hit_plv_v[2].size() << std::endl; }
    
    if (fDebug) { std::cout << "polygon has size : " << shrPoly.Size() << std::endl; }

      // if no hits associated on this plane -> skip
    if (shr_hit_plv_v[2].size() == 0) continue;
    
    // loop over photons for this plane
    for (std::map<size_t,twodimtools::Poly2D>::iterator it = _photon_poly_map.begin(); it != _photon_poly_map.end(); ++it) {
      
      auto const& photonPoly = it->second;
      auto const& photonIdx  = it->first;
      // get linearity
      auto const& photonLin = _photon_lin_map[ photonIdx ];

      if (fDebug) { std::cout << "\n\n new photon "; }

      // check the number of hits in the cluster
      auto nhits = photon_hit_assn_v.at(photonIdx).size();
      
      if (fDebug) { std::cout << "with " << nhits << " hits and COM [" << photonLin._meanx << ", " << photonLin._meany << "]. index " << photonIdx << " of " << photon_h->size()-1 << std::endl; }
      
      // apply cut on fraction of photon charge (# of hits here) that can be added to existing shower
      if (nhits > (shr_clus_v.at(collidx)->NHits() * fFracShrQ) ) {
	if (fDebug) std::cout << "\t hit num too large... continue" << std::endl;
	continue;
      }
      

      
      // do the polygons overlap?
      bool overlap = ( shrPoly.Overlap(photonPoly) || shrPoly.Contained(photonPoly) );
      
      if (overlap == false) {
	if (fDebug) std::cout << "\t no overlap w/ shower cone... continue" << std::endl;
	continue;
      }
      
      // for large clusters, add extra checks
      if (nhits > 8) {
	
	// apply cut on slope agreement
	if (slopeCompat( shrPoly, photonLin) > fMaxSlopeAngle) {
	  if (fDebug) std::cout << "\t slope not compatible with shower... continue" << std::endl;
	  continue;
	}
	
	// make sure the photon doesn't extend on both sides of the shower cone
	if (photonCrossesShower(shrPoly, photonPoly) == true) {
	  if (fDebug) std::cout << "\t photon fully crosses shower... continue" << std::endl;
	  continue;
	}

      }// if cluster is large
	
      // made it this far -> merge the photon with the shower
      // get set of hits to add (removing potential duplicates)
      if (fDebug) std::cout << "Merge photon. Photon hits : " << photon_hit_assn_v.at(photonIdx).size() 
			    << "\t shr hits : " << shr_hit_plv_v[2].size() << std::endl;

      // if photon not yet found compatible with any shower
      if (Photon_Shower_Map.find( photonIdx ) == Photon_Shower_Map.end() ) {
	Photon_Shower_Map[ photonIdx ] = { std::make_pair(s, PhotonShowerAngle(shrPoly,photonLin)) };
      }
      // else, if already compatible
      else {
      	Photon_Shower_Map[ photonIdx ].push_back( std::make_pair(s, PhotonShowerAngle(shrPoly,photonLin) ) );
      }

    }// for all gammas
  }// for all showers

  // now loop through showers, create new PFParticles and save new clusters which have merged hits from Y plane
  for (size_t s=0; s < shr_h->size(); s++) {

    if (fDebug) { std::cout << std::endl << std::endl << "New shower index " << s << std::endl << std::endl; }

    // split shower hits by plane
    std::vector< art::Ptr<recob::Hit> > shr_hit_v = shr_hit_assn_v.at(s);
    std::vector< std::vector< art::Ptr<recob::Hit> > > shr_hit_plv_v;
    shr_hit_plv_v.resize(3);
    for (auto const& shr_hit : shr_hit_v)
      shr_hit_plv_v.at(shr_hit->WireID().Plane).push_back( shr_hit ); 

    // loop through all merged photons, if merged to this shower, and this shower is "the best merge"
    // then add hits
    for (std::map<size_t,std::vector< std::pair<size_t,double> > >::iterator it = Photon_Shower_Map.begin(); it != Photon_Shower_Map.end(); ++it) {    

      // photon index
      auto photonIdx = it->first;
      // vector of showers to merge with?
      auto showerV = it->second;

      std::vector< art::Ptr<recob::Hit> > photon_hit_v = photon_hit_assn_v.at( photonIdx );

      if (fDebug && (showerV.size() > 0) ) { 
	auto photonlin = _photon_lin_map[photonIdx];
	std::cout << "cluster with COM [w, t] -> ["  << photonlin._meanx << ", " << photonlin._meany << "] and with " << photon_hit_v.size() << " hits" << std::endl; 
      }
      
      // merged with this shower?
      bool   mergeshower = false;
      double bestangle   = 1e6; // angle with best match shower
      double thisangle   = 1e6; // angle with the shower we are interested in
      for (auto const& shrinfo : showerV) {
	if (fDebug) { std::cout << "photon idx : " << photonIdx << " associated to shr " << shrinfo.first << " with angle " << shrinfo.second << std::endl; }
	if (shrinfo.first == s) { 
	  mergeshower = true; 
	  if (shrinfo.second < thisangle) { thisangle = shrinfo.second; }
	}
	if (shrinfo.second < bestangle) { bestangle = shrinfo.second; }
      }

      if (mergeshower == false) continue;
      if (thisangle > bestangle) continue;

      if (fDebug) { std::cout << "shower index " << s << " is the best match! add!" << std::endl; }
      

      // ok -> made it this far! merge
      auto nhitsmerged = MergeHits(shr_hit_plv_v[2],photon_hit_v);
      
      if (fDebug) { std::cout << "\t after merging -> " << shr_hit_plv_v[2].size() << " hits. Merged " << nhitsmerged << "/" << photon_hit_v.size() << " hits." << std::endl; }
      
    }// for all gamma polygons
    
    if (fDebug) std::cout << "\n\n Creating new PFParticles..." << std::endl;
    
    // create a new PFParticle
    recob::PFParticle newpfpart(*(shr_pfp_assn_v.at(s).at(0)));
    PFParticle_v->emplace_back(newpfpart);
    art::Ptr<recob::PFParticle> const PFPPtrAssn = PFPartPtrMaker(PFParticle_v->size()-1);

    // save clusters/hits
    for (size_t pl=0; pl < 3; pl++) {

      if (fDebug) std::cout << "\t plane " << pl << std::endl;

      auto const& plane_hits = shr_hit_plv_v.at(pl);

      // zero hits? skip...
      if (plane_hits.size() == 0) continue;

      if (fDebug) std::cout << "\t has valid cluster! " << std::endl;

      // make cluster from which to extract start/end wire/tick
      ::cluster::Cluster CMCluster;
      _clusterMaker.MakeCluster(plane_hits,CMCluster);

      float startW = CMCluster._start_pt._w / _wire2cm;
      float startT = CMCluster._start_pt._t / _time2cm;// + detp->TriggerOffset();
      
      float endW   = CMCluster._end_pt._w / _wire2cm;
      float endT   = CMCluster._end_pt._t / _time2cm;// + detp->TriggerOffset();

      if (pl ==2 && fDebug)
	std::cout << "\t new cluster: Shower Collection plane start wire/tick are @ " << startW << ", " << startT << std::endl;
      
      auto planeid = geo::PlaneID(0,0,pl);
      
      recob::Cluster clus(startW, 0., startT, 0., 0., CMCluster._angle, 0., 
			  endW,   0., endT,   0., 0., 0., 0., 
			  CMCluster._sum_charge, 0., CMCluster._sum_charge, 0., 
			  CMCluster.GetHits().size(), 0., 0., (s*3)+pl,
			  geom->View(planeid),
			  planeid);
      
      if (fDebug) std::cout << "\t creating cluster and hit in event record! " << std::endl;

      Cluster_v->emplace_back(clus);
      art::Ptr<recob::Cluster> const ClusPtrAssn = ClusPtrMaker(Cluster_v->size()-1);
      PFParticle_Cluster_assn_v->addSingle(PFPPtrAssn,ClusPtrAssn);
      
      for (auto hitPtr : plane_hits) {
	Hit_v->emplace_back( *hitPtr );
	art::Ptr<recob::Hit> const HitPtrAssn = HitPtrMaker(Hit_v->size()-1);
	Cluster_Hit_assn_v->addSingle(ClusPtrAssn,HitPtrAssn);
      }
    }
    
  }// for all showers
  
  e.put(std::move(PFParticle_v));
  e.put(std::move(Cluster_v));
  e.put(std::move(Hit_v));
  e.put(std::move(PFParticle_Cluster_assn_v));
  e.put(std::move(Cluster_Hit_assn_v));

  return;
  }

twodimtools::Poly2D PhotonMerge::projectShower(const art::Ptr<recob::Cluster> clus) {
  
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  double sW = clus->StartWire() * _wire2cm;
  double sT = ( clus->StartTick() - detp->TriggerOffset() ) * _time2cm;

  double dW = (sW-_vtxW);
  double dT = (sT-_vtxT);
  
  double dmag = sqrt((dW*dW)+(dT*dT));
  dW /= dmag;
  dT /= dmag;
  
  double oangle = 0.;
  if (fWidth > 0) oangle = fWidth * 3.14 / 180.;
  
  if (fDebug) { std::cout << "direction : dW " << dW << "\t dT : " << dT << std::endl; }
  
  double eW = sW + fShrLen * dW;
  double eT = sT + fShrLen * dT;
  
  // vector of coordiantes with which to consutrct triangle polygon
  std::vector<std::pair<float,float> > triangle_coordinates;
  
  // create polygon vertex for triangle
  triangle_coordinates.push_back( std::pair<float,float>(sW, sT) );
  
  // figure out how far to go on each side.
  double shrWidth = fShrLen * fabs(tan(oangle)) * 0.5;
  
  if (fDebug) { std::cout << "width : " << shrWidth << std::endl; }
  
  // unlike length, width is not stretched or compressed on projection.
  // extend end-point "left" and "right" by one width to complete triangle
  // extend in direction perpendicular to line connecting start_pl and end_pl
  double slope = -1. / ( (eT - sT) / (eW - sW) );
  if (fDebug)
    {
      std::cout << "\t End pt  = [ "  << eW << ", " << eT << " ]" << std::endl;
      std::cout << "\t SLOPE = " << slope << std::endl;
      std::cout << "\t WIDTH = " << shrWidth << std::endl;
    }
  
  // find triangle base coordinate points 
  double pt1_w = eW + shrWidth * sqrt( 1. / (1 + slope*slope) );
  double pt1_t = eT + slope * shrWidth * sqrt( 1. / (1 + slope*slope) );
  triangle_coordinates.push_back( std::pair<float,float>( pt1_w, pt1_t) );
  
  double pt2_w = eW - shrWidth * sqrt( 1. / (1 + slope*slope) );
  double pt2_t = eT - slope * shrWidth * sqrt( 1. / (1 + slope*slope) );
  triangle_coordinates.push_back( std::pair<float,float>( pt2_w, pt2_t) );
  
  if (fDebug){
    std::cout << "SHOWER COORDINATES : " << std::endl
	      << "\t Vertex @ [ " << sW << ", " << sT << " ]" << std::endl
	      << "\t Pt1    @ [ " << pt1_w      << ", " << pt1_t      << " ]" << std::endl
	      << "\t Pt2    @ [ " << pt2_w      << ", " << pt2_t      << " ]" << std::endl;
  }
  
  return twodimtools::Poly2D(triangle_coordinates);
}

double PhotonMerge::PhotonShowerAngle(const twodimtools::Poly2D& shr,
				      const twodimtools::Linearity& photon) {
  
  // get slope of photon (w.r.t. vertex)
  double photonYdiff = ( photon._meany - _vtxT );
  double photonXdiff = ( photon._meanx - _vtxW );
  double photon_slope_val = photonYdiff / photonXdiff;
  
  // get shower slope
  double shr_slope_val;
  
  double ydiff = ( (shr.Point(2).second + shr.Point(1).second) / 2. ) - shr.Point(0).second;
  double xdiff = ( (shr.Point(2).first + shr.Point(1).first) / 2. ) - shr.Point(0).first;
  
  shr_slope_val = ydiff/xdiff;

  double angle = fabs ( atan( ( shr_slope_val - photon_slope_val ) / ( 1 + photon_slope_val * shr_slope_val ) ) );
  
  return angle;
}

double PhotonMerge::slopeCompat(const twodimtools::Poly2D& shr,
				const twodimtools::Linearity& photon) {
  
  
  // get slope of photon
  auto const& photon_slope_val = photon._slope;
  
  // get shower slope
  double shr_slope_val;

  double ydiff = ( (shr.Point(2).second + shr.Point(1).second) / 2. ) - shr.Point(0).second;
  double xdiff = ( (shr.Point(2).first + shr.Point(1).first) / 2. ) - shr.Point(0).first;

  shr_slope_val = ydiff/xdiff;

  double angle = fabs ( atan( ( shr_slope_val - photon_slope_val ) / ( 1 + photon_slope_val * shr_slope_val ) ) );

  angle *= (180./3.14);

  if (fDebug) { std::cout << "\t photon slope angle = " << angle << std::endl; }

  return angle;

}

bool PhotonMerge::photonCrossesShower(const twodimtools::Poly2D& shr,
				      const twodimtools::Poly2D& photon) {
  
  // get the triangle points
  double Ox = shr.Point(0).first;
  double Oy = shr.Point(0).second;
  double Ax = shr.Point(1).first;
  double Ay = shr.Point(1).second;
  double Bx = shr.Point(2).first;
  double By = shr.Point(2).second;

  // check where the shower end-point lies
  double Ex = (Ax+Bx)/2.;
  double Ey = (Ay+By)/2.;

  double signEA = (Ex-Ox)*(Ay-Oy) - (Ey-Oy)*(Ax-Ox);
  double signEB = (Ex-Ox)*(By-Oy) - (Ey-Oy)*(Bx-Ox);

  // how many photon edges are on the left or right of the cone?
  int slope_left  = 0;
  int slope_right = 0;

  // being on one side or another of the segments OA and OB depends on the sign
  // of the equations used to compute signEA and signEB
  // to check if points are on either side of the cone compare their sign
  // with the signs of EA and EB knowing that these two are in the center.
  for (size_t i=0; i < photon.Size(); i++) {

    auto const& pt = photon.Point(i);

    double signA = (pt.first-Ox)*(Ay-Oy) - (pt.second-Oy)*(Ax-Ox);
    double signB = (pt.first-Ox)*(By-Oy) - (pt.second-Oy)*(Bx-Ox);

    if ( ((signA*signEA) < 0) && ((signB*signEB) > 0) ) slope_left  += 1;
    if ( ((signA*signEA) > 0) && ((signB*signEB) < 0) ) slope_right += 1;

  }// for all photon points
  
  if ( (slope_left > 0) and (slope_right > 0) ){

    if (fDebug) std::cout << "\t photon crosses shower" << std::endl;
    return true;
  }

  return false;
}


int PhotonMerge::MergeHits(std::vector< art::Ptr<recob::Hit> >& shrhits, const std::vector< art::Ptr<recob::Hit> >& gammahits) {

  int nhitsmerged = 0;

  for (auto const& gammahit : gammahits) {

    bool duplicate = false;

    for (auto const& shrhit : _allshr_hit_v) {
      if ( (shrhit->WireID().Wire == gammahit->WireID().Wire) && (shrhit->PeakTime() == gammahit->PeakTime() ) ) {
	duplicate = true;
	break; 
      }// if duplicate
    }// for hits in gamma hits
    
    if (duplicate == false) {
      nhitsmerged +=1;
      shrhits.push_back(gammahit);
    }
    
  }// for hits in shower hit
	
  return nhitsmerged;
}

void PhotonMerge::beginJob()
{
  // Implementation of optional member function here.
}

void PhotonMerge::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PhotonMerge)
