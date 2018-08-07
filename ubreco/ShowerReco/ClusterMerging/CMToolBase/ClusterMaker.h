/**
 * \file ClusterMaker.h
 *
 * \ingroup CMToolBase
 * 
 * \brief Class def header for a class Cluster
 *
 * @author david caratelli
 */

/** \addtogroup CMToolBase

    @{*/
#ifndef CLUSTER_CLUSTERMAKER_H
#define CLUSTER_CLUSTERMAKER_H

#include <iostream>
#include <vector>

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "art/Framework/Principal/Event.h"

// Data Products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "Cluster.h"

/**
   \class Cluster
   User defined class Cluster ... these comments are used to generate
   doxygen documentation!
 */

namespace cluster {


  
  class ClusterMaker {
    
  public:
    
    /// Default constructor
    ClusterMaker();
    
    /// Default destructor
    ~ClusterMaker(){}

    void MakeClusters(const art::ValidHandle<std::vector<recob::Cluster> >& clus_h,
		      const art::FindManyP<recob::Hit>&  clus_hit_assn_v,
		      const art::ValidHandle<std::vector<recob::Vertex> >& vtx_h,
		      std::vector<::cluster::Cluster>& cluster);

    
    void MakeCluster(const std::vector<art::Ptr<recob::Hit> >& hit_v,
		     ::cluster::Cluster& cluster);

    bool loadVertex(const art::ValidHandle<std::vector<::recob::Vertex> > vtx_h);

  private:

    void GetClusterPts(const std::vector<art::Ptr<recob::Hit> >& hit_v,
		       std::vector<::cluster::pt>& pt_v);

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;
    
    /// conversion factors for hits
    double _wire2cm, _time2cm;

  };

}
#endif
/** @} */ // end of doxygen group 

