#ifndef DEDXMODULE_CXX
#define DEDXMODULE_CXX

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

namespace showerreco {

  class dEdxModule : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    dEdxModule(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~dEdxModule(){};

    void configure(const fhicl::ParameterSet& pset);
    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
    void initialize();
    
  protected:

    // distance along which to calculate dEdx
    double _dtrunk;

    /*
    float ChargeCorrection(const double& q, const double& w, const double& t, const TVector3& dir, const TVector3& vtx,
			   const int& pl, const lariov::TPCEnergyCalibProvider& energyCalibProvider);
    */

    // debugging tree
    TTree* _dedx_tree;
    double _dedx0, _dedx1, _dedx2;
    std::vector<double> _dedx0_v, _dedx1_v, _dedx2_v;
    std::vector<double> _dist0_v, _dist1_v, _dist2_v;
    int    _pl0, _pl1, _pl2;
    double _pitch0, _pitch1, _pitch2;
    int    _nhits0, _nhits1, _nhits2;
    int    _ntot0, _ntot1, _ntot2;
    double _px, _py, _pz;
    
  };
  
  dEdxModule::dEdxModule(const fhicl::ParameterSet& pset)
  {
    _name = "dEdxModule";
    configure(pset);

    art::ServiceHandle<art::TFileService> tfs;
    _dedx_tree = tfs->make<TTree>("_dedx_tree","dE/dx TTree");

    _dedx_tree->Branch("_px",&_px,"px/D");
    _dedx_tree->Branch("_py",&_py,"py/D");
    _dedx_tree->Branch("_pz",&_pz,"pz/D");

    _dedx_tree->Branch("_pl0",&_pl0,"pl0/I");
    _dedx_tree->Branch("_pitch0",&_pitch0,"pitch0/D");
    _dedx_tree->Branch("_ntot0",&_ntot0,"ntot0/I");
    _dedx_tree->Branch("_nhits0",&_nhits0,"nhits0/I");
    _dedx_tree->Branch("_dedx0",&_dedx0,"dedx0/D");
    _dedx_tree->Branch("_dedx0_v","std::vector<double>",&_dedx0_v);
    _dedx_tree->Branch("_dist0_v","std::vector<double>",&_dist0_v);

    _dedx_tree->Branch("_pl1",&_pl1,"pl1/I");
    _dedx_tree->Branch("_pitch1",&_pitch1,"pitch1/D");
    _dedx_tree->Branch("_ntot1",&_ntot1,"ntot1/I");
    _dedx_tree->Branch("_nhits1",&_nhits1,"nhits1/I");
    _dedx_tree->Branch("_dedx1",&_dedx1,"dedx1/D");
    _dedx_tree->Branch("_dedx1_v","std::vector<double>",&_dedx1_v);
    _dedx_tree->Branch("_dist1_v","std::vector<double>",&_dist1_v);

    _dedx_tree->Branch("_pl2",&_pl2,"pl2/I");
    _dedx_tree->Branch("_pitch2",&_pitch2,"pitch2/D");
    _dedx_tree->Branch("_ntot2",&_ntot2,"ntot2/I");
    _dedx_tree->Branch("_nhits2",&_nhits2,"nhits2/I");
    _dedx_tree->Branch("_dedx2",&_dedx2,"dedx2/D");
    _dedx_tree->Branch("_dedx2_v","std::vector<double>",&_dedx2_v);
    _dedx_tree->Branch("_dist2_v","std::vector<double>",&_dist2_v);

  }

  void dEdxModule::configure(const fhicl::ParameterSet& pset)
  {
    _dtrunk = pset.get<double>("dtrunk");
    _verbose   = pset.get<bool>("verbose",false);
  }
  
  void dEdxModule::initialize()
  {
    
    return;
  }
  
  void dEdxModule::do_reconstruction(const ::protoshower::ProtoShower & proto_shower, Shower_t & resultShower) {

    //handle to tpc energy calibration provider
    //const lariov::TPCEnergyCalibProvider& energyCalibProvider  = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
    
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()){
      std::stringstream ss;
      ss << "aaa";
      throw ShowerRecoException(ss.str());
    }
    
    auto& clusters = proto_shower.clusters();
    
    // grab shower direction
    auto const& dir3D = resultShower.fDCosStart;

    std::cout << "3D shower direction : " << dir3D[0] << ", " << dir3D[1] << ", " << dir3D[2] << std::endl;

    _px = dir3D[0];
    _py = dir3D[1];
    _pz = dir3D[2];

    // unused auto const& geomH = ::util::GeometryUtilities();

    // loop through planes
    for (size_t n = 0; n < clusters.size(); n++) {
      
      auto const& clus = clusters.at(n);
      
      // get the hits associated with this cluster
      auto const& hits = clus._hits;
      
      // get the plane associated with this cluster
      auto const& pl = clus._plane;

      // get start point on pllane
      auto& start2D = clus._start;
      
      std::cout << std::endl << "PLANE : " << pl << std::endl;

      auto const* geom = ::lar::providerFrom<geo::Geometry>();
      const geo::WireGeo& wire = geom->TPC().Plane(pl).MiddleWire();
      TVector3 wireunitperp = wire.Direction();//(wire.GetStart()-wire.GetEnd()).Unit();
      // rotate by 90 degrees around x
      TVector3 wireunit = {wireunitperp[0], -wireunitperp[2], wireunitperp[1]}; 
      std::cout << "wire unit on plane : " << pl << " is " << wireunit[0] << ", " << wireunit[1] << ", " << wireunit[2] << std::endl;
      double cosPlane = fabs(cos(wireunit.Angle(dir3D)));

      std::vector<double> dedx_v;
      std::cout << "dtrunk is " << _dtrunk << std::endl;
      dedx_v.resize(3 * _dtrunk);
      for (size_t jj=0; jj < dedx_v.size(); jj++) 
	dedx_v.at(jj) = 0.;
      std::cout << "dedx_v size is " << dedx_v.size() << std::endl;
      //double dedx;
      int nhits = 0;
      double pitch = 0.3 / cosPlane;
      std::cout << " dEdx Module : pitch = " << pitch << " from function" <<  std::endl;      

      // loop through hits and find those within some radial distance of the start point
      // loop over hits
      for (auto const &h : hits) {
	
	double d2D = sqrt( pow(h.w - start2D.w, 2) + pow(h.t - start2D.t, 2) );
	double d3D = d2D / cosPlane;
	size_t d3Delement = (size_t)(d3D * 3);
	double dE = h.charge;// * ChargeCorrection(h.charge, h.w, h.t, resultShower.fDCosStart, resultShower.fXYZStart, pl, energyCalibProvider);
	double dEdx = dE / pitch;

	if (d3Delement >= dedx_v.size()) continue;

	std::cout << "\t d2D : " << d2D << "\t d3D : " << d3D << " \t d3D int : " << d3Delement 
		  << "\t dEdx : " << dEdx
		  << std::endl;

	dedx_v.at( d3Delement ) += dEdx;
	nhits += 1;
	
      }// loop over all hits

      std::vector<double> dedx_empty_v;
      double dedx;

      for (size_t n=0; n < dedx_v.size(); n++) {
	if (dedx_v.at(n) != 0){
	  std::cout << "\t adding an element..." << std::endl;
	  dedx_empty_v.push_back(dedx_v.at(n));
	  resultShower.fdEdx_v_v.at(pl).push_back(dedx_v.at(n));
	  std::cout << "\t added..." << std::endl;
	}
      }// for all dedx values
      std::cout << "done erasing" << std::endl;
      
      std::cout << "number of dedx points :  " << dedx_empty_v.size() << std::endl;
      
      if (dedx_empty_v.size() == 0)
	dedx = 0.;

      else {
	std::sort( dedx_empty_v.begin(), dedx_empty_v.end() );
	//std::nth_element(dedx_empty_v.begin(), dedx_empty_v.end(), dedx_empty_v.end() );
	dedx = dedx_empty_v[dedx_empty_v.size()/2.];
      }

      for (auto const& aa : dedx_empty_v)
	std::cout << "dedx Module : \t dedx = " << aa << std::endl;
      std::cout << "dedx Module : Final dEdx = " << dedx << std::endl;

      if (pl == 0) {
	_pitch0 = pitch;
	_nhits0 = nhits;
	_dedx0_v = dedx_v;
	_dedx0 = dedx;
	_ntot0 = hits.size();
      }
      if (pl == 1) {
	_pitch1 = pitch;
	_nhits1 = nhits;
	_dedx1_v = dedx_v;
	_dedx1 = dedx;
	_ntot1 = hits.size();
      }
      if (pl == 2) {
	_pitch2 = pitch;
	_nhits2 = nhits;
	_dedx2_v = dedx_v;
	_dedx2 = dedx;
	_ntot2 = hits.size();
      }

      resultShower.fdEdx_v.at(pl) = dedx;
      if (pl == 2) {
	resultShower.fBestdEdxPlane = pl;
	resultShower.fBestdEdx   = dedx;
      }
      
    }// for all clusters (planes)

    _dedx_tree->Fill();
    
    return;
  }

  /*
  float dEdxModule::ChargeCorrection(const double& q, const double& w, const double& t, const TVector3& dir, const TVector3& vtx,
				     const int& pl, const lariov::TPCEnergyCalibProvider& energyCalibProvider){

    // find 3D position of hit                                                                                                        
    double z = w;
    double x = t;

    // get 2D distance of hit to vtx                                                                                                  
    double r2D = sqrt( ( (z-vtx.Z()) * (z-vtx.Z()) ) + ( (x-vtx.X()) * (x-vtx.X()) ) );
    double r3D = r2D/dir[1];

    auto xyz = vtx + dir * r3D;

    float yzcorrection = energyCalibProvider.YZdqdxCorrection(pl, xyz.Y(), xyz.Z());
    float xcorrection  = energyCalibProvider.XdqdxCorrection(pl,  xyz.X());

    if (!yzcorrection) yzcorrection = 1.;
    if (!xcorrection ) xcorrection  = 1.;

    return yzcorrection * xcorrection;
  }
  */

  DEFINE_ART_CLASS_TOOL(dEdxModule)
} //showerreco

#endif
