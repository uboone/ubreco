#ifndef SEAVIEWER_H
#define SEAVIEWER_H

#include <iostream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include <string>

#include <memory>

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPrincipal.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TEllipse.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <sys/stat.h>


#include "DBSCAN.h"
namespace seaview {

class cluster {

    public:
    
    cluster(int ID, int plane, std::vector<std::vector<double>> &pts, std::vector<art::Ptr<recob::Hit>> &hits) :f_ID(ID), f_plane(plane), f_pts(pts), f_hits(hits) {

        f_npts = f_pts.size();
        if(pts.size() != hits.size()){
            std::cerr<<"seaviewer::cluster, input hits and pts not alligned"<<std::endl;
        }
        std::vector<double> wires(f_npts);
        std::vector<double> ticks(f_npts);
        for(int p =0; p< f_npts; ++p){
            wires[p]=f_pts[p][0];
            ticks[p]=f_pts[p][1];
        }
        TGraph af_graph(f_npts,&wires[0],&ticks[0]);
        f_graph = af_graph;

    };

    int getID() {return f_ID;}
    int getN() {return f_npts;}
    int getPlane(){ return f_plane;}
    TGraph * getGraph(){ return &f_graph;}

    private:
    int f_ID;
    int f_npts;
    int f_plane;
    std::vector<std::vector<double>> f_pts;
    std::vector<art::Ptr<recob::Hit>> f_hits;
    TGraph f_graph;
};



class SEAviewer {

public:

  /// Default constructor
  SEAviewer(std::string tag,geo::GeometryCore const * geom,detinfo::DetectorProperties const * theDetector );

  void configure(const fhicl::ParameterSet& pset){};

  int loadVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z);
  int addTrueVertex(double x, double y,double z);
  int addSliceHits(std::vector<art::Ptr<recob::Hit>>& hits);
  int addAllHits(std::vector<art::Ptr<recob::Hit>>& hits);
  int addPFParticleHits(std::vector<art::Ptr<recob::Hit>>& hits, std::string leg );
  int setBadChannelList(std::vector<std::pair<int,int>> &in);
  int addShower(art::Ptr<recob::Shower>&shr);
  int calcUnassociatedHits();
  int setHitThreshold(double);
  int Print();
  int runDBSCAN(double min_pts, double eps);

  double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
        double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
        return wire;
    }

  double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorProperties const& detprop){
        double time = detprop.ConvertXToTicks(X, plane, fTPC,fCryostat);
        return time;
    }


    std::vector<std::vector<double>> to2D(std::vector<double> & threeD);




  protected:
  int n_pfps;
  int n_showers;
  std::string tag;
  double hit_threshold;
  bool has_been_clustered;  
  std::vector<std::vector<TGraph>> vec_graphs;

  std::vector<std::string> vec_pfp_legend;
  // PFP, Plane: index
  std::vector<std::vector<std::vector<double>>> vec_ticks;
  std::vector<std::vector<std::vector<double>>> vec_chans;

  geo::GeometryCore const * geom;
  detinfo::DetectorProperties const * theDetector ;

    double tick_shift;
    double chan_shift;

  double tick_max;
  double tick_min;
  std::vector<double> chan_max;
  std::vector<double> chan_min;

  std::vector<std::pair<int,int>> m_bad_channel_list;

  //Vertex
  std::vector<double> vertex_tick; 
  std::vector<double> vertex_chan; 
  std::vector<TGraph> vertex_graph;

  bool plot_true_vertex;
  std::vector<double> true_vertex_tick; 
  std::vector<double> true_vertex_chan; 
  std::vector<TGraph> true_vertex_graph;

  std::vector<art::Ptr<recob::Hit>> slice_hits;
  std::vector<art::Ptr<recob::Hit>> all_hits;
  std::map<art::Ptr<recob::Hit>,bool> map_unassociated_hits;
  std::map<art::Ptr<recob::Hit>, bool> map_slice_hits;
  
  std::vector<TGraph> vec_unass_graphs;
  std::vector<std::vector<double>> vec_unass_ticks;
  std::vector<std::vector<double>> vec_unass_chans;
  std::vector<std::vector<std::vector<double>>> vec_unass_pts;
  std::vector<std::vector<art::Ptr<recob::Hit>>> vec_unass_hits;

  std::vector<TGraph> vec_all_graphs;
  std::vector<std::vector<double>> vec_all_ticks;
  std::vector<std::vector<double>> vec_all_chans;

  std::vector<int> num_clusters;
  std::vector<std::vector<int>> cluster_labels;
  TRandom3 *rangen;

  std::vector<seaview::cluster> vec_clusters;  
  std::vector<art::Ptr<recob::Shower>> vec_showers;

};

}// namespace

#endif

