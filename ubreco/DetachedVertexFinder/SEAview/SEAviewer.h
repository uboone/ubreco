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

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <sys/stat.h>


namespace seaview {

class SEAviewer {

public:

  /// Default constructor
  SEAviewer(std::string tag,geo::GeometryCore const * geom,detinfo::DetectorProperties const * theDetector );

  void configure(const fhicl::ParameterSet& pset){};

  int loadVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z);
  int addPFParticleHits(std::vector<art::Ptr<recob::Hit>>& hits);
  int setBadChannelList(std::vector<std::pair<int,int>> &in);
  int Print();
   
  double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
        double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
        return wire;
    }

  double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorProperties const& detprop){
        double time = detprop.ConvertXToTicks(X, plane, fTPC,fCryostat);
        return time;
    }



  protected:
  int n_pfps;
  std::string tag;
  std::vector<std::vector<TGraph>> vec_graphs;

  // PFP, Plane: index
  std::vector<std::vector<std::vector<double>>> vec_ticks;
  std::vector<std::vector<std::vector<double>>> vec_chans;

  geo::GeometryCore const * geom;
  detinfo::DetectorProperties const * theDetector ;

  double tick_max;
  double tick_min;
  std::vector<double> chan_max;
  std::vector<double> chan_min;

  std::vector<std::pair<int,int>> m_bad_channel_list;

  //Vertex
  std::vector<double> vertex_tick; 
  std::vector<double> vertex_chan; 
  std::vector<TGraph> vertex_graph;


 
};

}// namespace

#endif

