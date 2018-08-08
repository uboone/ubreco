#ifndef PROTOSHOWER_PROTOSHOWERALGBASE_CXX
#define PROTOSHOWER_PROTOSHOWERALGBASE_CXX

#include "ProtoShowerCMTool.h"

namespace protoshower {

  void ProtoShowerCMTool::GenerateProtoShowers(::art::Event & e,
					       const std::string& fPFPproducer,
					       const std::string& fClusterproducer,
					       const std::string& fVtxproducer,
					       std::vector<protoshower::ProtoShower> & proto_shower_v) {

    // grab PFParticles in event
    auto const& pfp_h = e.getValidHandle<std::vector<recob::PFParticle> >(fPFPproducer);

    // grab clusters associated with PFParticles
    art::FindManyP<recob::Cluster> pfp_clus_assn_v(pfp_h, e, fPFPproducer);

    // ADDITION FROM PETRILLO
    e.getValidHandle<std::vector<recob::Cluster>>(fClusterproducer);

    // grab the hits associated to the PFParticles
    auto pfp_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(pfp_h, e, fPFPproducer);

    // load event vertex associated to tagged neutrino interaction
    auto const& vertex_h = e.getValidHandle<std::vector<recob::Vertex> >(fVtxproducer);

    // loop through PFParticles
    for (size_t p=0; p < pfp_h->size(); p++) {

      ::protoshower::ProtoShower proto_shower;
      proto_shower.Reset();

      const recob::PFParticle pfp = pfp_h->at(p);

      // associated clusters
      const std::vector< art::Ptr<recob::Cluster> >& clus_v = pfp_clus_assn_v.at(p);

      // associated hits
      const std::vector< art::Ptr<recob::Hit> >& hit_v = pfp_hit_assn_v.at(p);

      // set number of clusters for protoshower
      proto_shower._clusters.resize(clus_v.size());

      // loop through clusters
      for (size_t c=0; c < clus_v.size(); c++) {

	auto const clus = clus_v.at(c);

	// find the hits associated to this cluster
	std::vector< art::Ptr<recob::Hit> > clusterhits;
	for (auto const& hit : hit_v) {
	  // if hit in same plane as cluster -> add to vector
	  if (hit->WireID().Plane == clus->Plane().Plane )
	    clusterhits.push_back( hit );
	}//for all hits associated to PFParticle

	proto_shower._clusters.at(c) = MakeCluster2D( clus, clusterhits );
	
	proto_shower.hasCluster2D(true);

      }// for all clusters
      
      // require a single vertex!
      if (vertex_h->size() == 1) {
	
	auto const vtx = vertex_h->at(0);
	Double_t xyz[3] = {};
	vtx.XYZ(xyz);
	proto_shower._vertex = TVector3(xyz[0],xyz[1],xyz[2]);
	proto_shower.hasVertex(true);

	std::cout << "\t\t DD VTX @ [ " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " ]" << std::endl;
	
      }// if there is only one vertex
      
      proto_shower_v.push_back( proto_shower );
      
    }// for all PFParticles
    
  }// GenerateProtoShower end

}// namespace

#endif
