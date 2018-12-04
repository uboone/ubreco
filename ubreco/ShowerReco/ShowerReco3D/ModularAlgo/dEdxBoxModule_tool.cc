#ifndef DEDXBOXMODULE_CXX
#define DEDXBOXMODULE_CXX

#include <iostream>
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"

//#include "ubreco/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
//#include "ubreco/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"

#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h"

/**
\class dedxModule : ShowerRecoModuleBase
This is meant to compute the 2D dedx along the start of the shower.
*/

namespace showerreco
{

  class dEdxBoxModule : public ShowerRecoModuleBase
  {

    public:
      /// Default constructor
      dEdxBoxModule(const fhicl::ParameterSet& pset);

      /// Default destructor
      ~dEdxBoxModule(){};

      void configure(const fhicl::ParameterSet& pset);

      void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);

      void initialize();

      double median(std::vector<double> values);

      bool isInside( std::vector<double> P, std::vector< std::vector<double> > V);
      double getPitch(const TVector3 &direction, const int &pl);
      void buildRectangle(double length, double width,
                                          std::vector<double> &start,
                                          std::vector<double> &axis,
                                          std::vector<std::vector<double>> &points);

    protected:

      double _rectangle_length;
      double _rectangle_width;
      double _rectangle_tolerance;
      bool _verbose;

      // debugging tree
      TTree* _dedx_tree;
      double _dedx0, _dedx1, _dedx2;
      std::vector<double> _dedx0_v, _dedx1_v, _dedx2_v;
      double _pitch0, _pitch1, _pitch2;
      int    _nhits0, _nhits1, _nhits2;
      int    _ntot0, _ntot1, _ntot2;
      double _px, _py, _pz;
  };

  dEdxBoxModule::dEdxBoxModule(const fhicl::ParameterSet& pset)
  {
    _name = "dEdxBoxModule";
    configure(pset);

    art::ServiceHandle<art::TFileService> tfs;
    _dedx_tree = tfs->make<TTree>("_dedxBox_tree", "dE/dx Box TTree");

    _dedx_tree->Branch("_px",&_px,"px/D");
    _dedx_tree->Branch("_py",&_py,"py/D");
    _dedx_tree->Branch("_pz",&_pz,"pz/D");

    _dedx_tree->Branch("_pitch0",&_pitch0,"pitch0/D");
    _dedx_tree->Branch("_ntot0",&_ntot0,"ntot0/I");
    _dedx_tree->Branch("_nhits0",&_nhits0,"nhits0/I");
    _dedx_tree->Branch("_dedx0",&_dedx0,"dedx0/D");
    _dedx_tree->Branch("_dedx0_v","std::vector<double>",&_dedx0_v);

    _dedx_tree->Branch("_pitch1",&_pitch1,"pitch1/D");
    _dedx_tree->Branch("_ntot1",&_ntot1,"ntot1/I");
    _dedx_tree->Branch("_nhits1",&_nhits1,"nhits1/I");
    _dedx_tree->Branch("_dedx1",&_dedx1,"dedx1/D");
    _dedx_tree->Branch("_dedx1_v","std::vector<double>",&_dedx1_v);

    _dedx_tree->Branch("_pitch2",&_pitch2,"pitch2/D");
    _dedx_tree->Branch("_ntot2",&_ntot2,"ntot2/I");
    _dedx_tree->Branch("_nhits2",&_nhits2,"nhits2/I");
    _dedx_tree->Branch("_dedx2",&_dedx2,"dedx2/D");
    _dedx_tree->Branch("_dedx2_v","std::vector<double>",&_dedx2_v);

  }

  void dEdxBoxModule::configure(const fhicl::ParameterSet& pset)
  {
    _rectangle_length = pset.get<double>("rectangleLength", 4);
    _rectangle_width = pset.get<double>("rectangleWidth", 1);
    _rectangle_tolerance = pset.get<double>("rectangleTolerance", 0.001);
    _verbose   = pset.get<bool>("verbose", false);
  }

  void dEdxBoxModule::initialize()
  {
    return;
  }

  void dEdxBoxModule::do_reconstruction(const ::protoshower::ProtoShower & proto_shower, Shower_t & resultShower)
  {
    //handle to tpc energy calibration provider
    //const lariov::TPCEnergyCalibProvider& energyCalibProvider  = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();

    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D())
    {
      std::stringstream ss;
      ss << "aaa";
      throw ShowerRecoException(ss.str());
    }

    auto& clusters = proto_shower.clusters();

    // grab shower direction
    auto const& dir3D = resultShower.fDCosStart;

    if (_verbose)
    {
      std::cout << "\n3D shower direction : " << dir3D[0] << ", " << dir3D[1] << ", " << dir3D[2] << std::endl;
    }

    _px = dir3D[0];
    _py = dir3D[1];
    _pz = dir3D[2];

    double pitches[3];
    for (int i=0; i<3; i++)
    {
      pitches[i] = getPitch(dir3D, i);
    }

    _pitch0 = pitches[0];
    _pitch1 = pitches[1];
    _pitch2 = pitches[2];

    resultShower.fPitch_0 = pitches[0];
    resultShower.fPitch_1 = pitches[1];
    resultShower.fPitch_2 = pitches[2];

    if (_verbose)
    {
      std::cout << "pitches : "
                << pitches[0] << ", "
                << pitches[1] << ", "
                << pitches[2] << std::endl;
    }

    // loop through planes
    for (size_t n = 0; n < clusters.size(); n++)
    {
      auto const& clus = clusters.at(n);

      // get the hits associated with this cluster
      auto const& hits = clus._hits;
      // get the plane associated with this cluster
      auto const& pl = clus._plane;
      // get start point on pllane
      auto& start2D = clus._start;
      auto& end2D = clus._end;

      if (pl == 0)
      {
        resultShower.fClStart_t_0 = start2D.t;
        resultShower.fClStart_w_0 = start2D.w;
        resultShower.fClEnd_t_0 = end2D.t;
        resultShower.fClEnd_w_0 = end2D.w;
      }
      else if (pl == 1)
      {
        resultShower.fClStart_t_1 = start2D.t;
        resultShower.fClStart_w_1 = start2D.w;
        resultShower.fClEnd_t_1 = end2D.t;
        resultShower.fClEnd_w_1 = end2D.w;
      }
      else if (pl == 2)
      {
        resultShower.fClStart_t_2 = start2D.t;
        resultShower.fClStart_w_2 = start2D.w;
        resultShower.fClEnd_t_2 = end2D.t;
        resultShower.fClEnd_w_2 = end2D.w;
      }

      float start_angle = clus._start_angle_2d;
      float end_angle = clus._end_angle_2d;

      if (_verbose)
      {
        std::cout << "PLANE : " << pl << std::endl;
        std::cout << "start angle : " << start_angle << std::endl;
        std::cout << "end angle : " << end_angle << std::endl;
        std::cout << "simple angle : " << atan2((end2D.t - start2D.t), (end2D.w - start2D.w)) << std::endl;
      }

      double pitch = pitches[pl];

      std::vector<double> cluster_axis;
      std::vector<double> cluster_start;
      std::vector<double> cluster_end;

      if (pitch >= 0)
      {
        cluster_axis = {cos(start_angle),
                        sin(start_angle)};
        cluster_start = {start2D.w - _rectangle_tolerance * cos(start_angle),
                         start2D.t - _rectangle_tolerance * sin(start_angle)};
        cluster_end = {end2D.w, end2D.t};
      }
      else
      {
        cluster_axis = {-1. * cos(start_angle),
                        -1. * sin(start_angle)};
        cluster_start = {end2D.w + _rectangle_tolerance * cos(start_angle),
                         end2D.t + _rectangle_tolerance * sin(start_angle)};
        cluster_end = {start2D.w, start2D.t};
      }

      // Build rectangle 4 x 1 cm around the cluster axis
      std::vector<std::vector<double>> points;
      buildRectangle(_rectangle_length, _rectangle_width, cluster_start, cluster_axis, points);
      std::vector<double> dedx_v;

      if (_verbose)
      {
        std::cout << "pitch : " << pitch << std::endl;
        std::cout << "cluster_axis : " << cluster_axis[0] << " , " << cluster_axis[1] << std::endl;
        std::cout << "cluster_start : " << cluster_start[0] << " , " << cluster_start[1] << std::endl;
        std::cout << "cluster_end : " << cluster_end[0] << " , " << cluster_end[1] << std::endl;
      }

      for (auto &hit : hits)
      {
        std::vector<double> hit_pos = {hit.w, hit.t};

        bool is_within = isInside(hit_pos, points);

        if (is_within)
        {
          dedx_v.push_back(hit.charge / fabs(pitch));
        }

        if (_verbose)
        {
          std::cout << "hit : " << hit.w << " , " << hit.t << std::endl;
          std::cout << "is_within : " << is_within << std::endl;
          std::cout << "dEdx : " << hit.charge / fabs(pitch) << std::endl;
        }
      }

      int nhits = dedx_v.size();
      double dedx = median(dedx_v);

      if (_verbose)
      {
        for (auto const& aa : dedx_v)
        {
          std::cout << "dedx Module : \t dedx = " << aa << std::endl;
        }
        std::cout << "dedx Module : Final dEdx = " << dedx << std::endl;
      }

      if (pl == 0)
      {
        _nhits0 = nhits;
        _dedx0_v = dedx_v;
        _dedx0 = dedx;
        _ntot0 = hits.size();
        resultShower.fdEdxBox_0 = dedx;
      }
      if (pl == 1)
      {
        _nhits1 = nhits;
        _dedx1_v = dedx_v;
        _dedx1 = dedx;
        _ntot1 = hits.size();
        resultShower.fdEdxBox_1 = dedx;
      }
      if (pl == 2)
      {
        _nhits2 = nhits;
        _dedx2_v = dedx_v;
        _dedx2 = dedx;
        _ntot2 = hits.size();
        resultShower.fdEdxBox_2 = dedx;
      }

      // resultShower.fdEdxBox_v.at(pl) = dedx;
      // if (pl == 2)
      // {
      //   resultShower.fBestdEdxBoxPlane = pl;
      //   resultShower.fBestdEdxBox   = dedx;
      // }
    }
    _dedx_tree->Fill();
    return;
  }

  double dEdxBoxModule::median(std::vector<double> values)
  {
    size_t size = values.size();

    if (size == 0)
    {
      return -1;
    }
    else
    {
      sort(values.begin(), values.end());
      if (size % 2 == 0)
      {
        return (values[size / 2 - 1] + values[size / 2]) / 2;
      }
      else
      {
        return values[size / 2];
      }
    }
  }

  bool dEdxBoxModule::isInside( std::vector<double> P, std::vector< std::vector<double> > V)
  {
    int nvert = (int)V.size();
    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++)
    {
      if ( ((V[i][1]>P[1]) != (V[j][1]>P[1])) &&
      (P[0] < (V[j][0]-V[i][0]) * (P[1]-V[i][1]) / (V[j][1]-V[i][1]) + V[i][0]) )
      c = !c;
    }
    return c;
  }

  double dEdxBoxModule::getPitch(const TVector3 &direction, const int &pl)
  {
    // prepare a direction vector for the plane
    TVector3 wireDir = {0., 0., 0.};
    // the direction of the plane is the vector uniting two consecutive wires
    // such that this vector is perpendicular to both wires
    // basically this is the vector perpendicular to the wire length direction,
    // and still in the wire-plane direction
    if (pl == 0)
      wireDir = {0., -sqrt(3) / 2., 1 / 2.};
    else if (pl == 1)
      wireDir = {0., sqrt(3) / 2., 1 / 2.};
    else if (pl == 2)
      wireDir = {0., 0., 1.};

    // cosine between shower direction and plane direction gives the factor
    // by which to divide 0.3, the minimum wire-spacing
    double minWireSpacing = 0.3;
    double cos = wireDir.Dot(direction);

    cos /= (wireDir.Mag() * direction.Mag());
    // if cosine is 0 the direction is perpendicular and the wire-spacing is
    // infinite
    if (cos == 0)
      return std::numeric_limits<double>::max();
    double pitch = minWireSpacing / cos;
    return pitch;
  }

  void dEdxBoxModule::buildRectangle(double length, double width,
                                      std::vector<double> &start,
                                      std::vector<double> &axis,
                                      std::vector<std::vector<double>> &points)
  {
    double perp_axis[2] = {-axis[1], axis[0]};

    std::vector<double> p1 = {start[0] + perp_axis[0] * width / 2,
                              start[1] + perp_axis[1] * width / 2};
    std::vector<double> p2 = {p1[0] + axis[0] * length, p1[1] + axis[1] * length};

    std::vector<double> p3 = {start[0] - perp_axis[0] * width / 2,
                              start[1] - perp_axis[1] * width / 2};
    std::vector<double> p4 = {p3[0] + axis[0] * length, p3[1] + axis[1] * length};

    points.insert(points.end(), {p1, p2, p4, p3});
  }

DEFINE_ART_CLASS_TOOL(dEdxBoxModule)
} //showerreco

#endif
