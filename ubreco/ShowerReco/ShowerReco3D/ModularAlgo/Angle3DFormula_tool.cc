#ifndef ANGLE3DFORMULA_CXX
#define ANGLE3DFORMULA_CXX

#include <iostream>
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShowerRecoModuleBase.h"
/**
   \class ShowerRecoModuleBase
   User defined class ShowerRecoModuleBase ... these comments are used to generate
   doxygen documentation!
 */

namespace showerreco {

  class Angle3DFormula : public ShowerRecoModuleBase {
    
  public:
    
    /// Default constructor
    Angle3DFormula(const fhicl::ParameterSet& pset);
    
    /// Default destructor
    ~Angle3DFormula(){};

    void configure(const fhicl::ParameterSet& pset);
    
    /// set the maximum angle that is allowed as an error (radians)
    void setMaxAngleError(double err) { _max_err = err; }
    
    /// set whether we should validate the direction using info from reco'd vertex
    void setValidateDirection(bool on) { _validate_dir = on; }
    
    /// set minimum dot product between shower direction and vtx -> strt direction
    void setMinDotProduct(double d) { _dot_min = d; }
    
    void do_reconstruction(const ::protoshower::ProtoShower &, Shower_t &);
    
  private:
    
    // maximum error in the angle that is allowed
    double _max_err;
    
    // validate direction by using reconstructed vertex
    bool _validate_dir;
    
    // minimum dot product value
    double _dot_min;
    
  };

  Angle3DFormula::Angle3DFormula(const fhicl::ParameterSet& pset)
  {
    configure(pset);
    _name = "Angle3DFormula";
  }

  void Angle3DFormula::configure(const fhicl::ParameterSet& pset)
  {
    _verbose   = pset.get<bool>("verbose",false);
    _max_err = pset.get<double>("max_err");
    _validate_dir = pset.get<bool>("validate_dir");
    _dot_min = pset.get<double>("dot_min");
  }
  
  void Angle3DFormula::do_reconstruction(const ::protoshower::ProtoShower & proto_shower,
					 Shower_t& resultShower) {
    
    
    //if the module does not have 2D cluster info -> fail the reconstruction
    if (!proto_shower.hasCluster2D()){
      std::stringstream ss;
      ss << "Fail @ algo " << this->name() << " due to missing 2D cluster";
      throw ShowerRecoException(ss.str());
    }
    
    auto & clusters = proto_shower.clusters();
    
    auto const& geomH = ::util::GeometryUtilities();
    
    // use geometry helper function to go from 2 2D angles to 1 3D angle
    // do this for all plane combinations available
    
    if (_verbose)
      std::cout << std::endl;

    // store theta & phi angles
    std::vector<double> theta_v;
    std::vector<double> phi_v;
    
    for (size_t i = 0; i < clusters.size(); i++) {
      for (size_t j = i + 1; j < clusters.size(); j++) {
	
	auto const& cl1 = clusters[i];
	auto const& cl2 = clusters[j];
	
	double theta, phi;
	
	if (_verbose)
	  std::cout << "expamining cluster pair (" << i << ", " << j << ")" << std::endl;
	
	
	if (_verbose) {
	  std::cout << "input angles : " << std::endl
		    << "\t " << cl1._angle_2d << " @ plane " << cl1._plane << std::endl
		    << "\t " << cl2._angle_2d << " @ plane " << cl2._plane << std::endl;
	}
	
	geomH.Get3DaxisN((int)cl1._plane, (int)cl2._plane,
			 cl1._angle_2d, cl2._angle_2d,
			 phi, theta);

	theta *= (3.1415/180.);
	phi   *= (3.1415/180.);
	
	if (_verbose)
	  std::cout << "return angle 3D : phi -> " << phi << ", theta -> " << theta << std::endl;
	
	// if we are returned the garbage values that the function returns...
	// very ugly...
	if ( (phi == 0) or (theta == -999) )
	  continue;
	
	phi_v.push_back(phi);
	theta_v.push_back(theta);
	
      }// for index j
    }// for index i

    
    // get the average theta & phi
    double theta = 0., phi = 0.;
    double theta_err = 0., phi_err = 0.;
    // if no values actually used -> return
    if (theta_v.size() == 0)
      return;
    // if only 1 pair of clusters was used -> use their return values
    else if (theta_v.size() == 1) {
      theta = theta_v[0];
      phi   = phi_v[0];
    }
    // if 2 pairs were used -> average
    else if (theta_v.size() == 2) {
      // average the 2 results
      theta = (theta_v[0] + theta_v[1]) / 2.;
      phi   = (phi_v[0] + phi_v[1]) / 2.;
      theta_err = sqrt ( (theta_v[0] - theta) * (theta_v[0] - theta) +
			 (theta_v[1] - theta) * (theta_v[1] - theta) );
      phi_err   = sqrt ( (phi_v[0] - phi) * (phi_v[0] - phi) +
			 (phi_v[1] - phi) * (phi_v[1] - phi) );
    }// if 2 pairs used
    else if (theta_v.size() == 3) {
      // create new vectors where to store the angles measured (if they are good)
      std::vector<double> theta_good;
      std::vector<double> phi_good;
      std::vector<double> theta_err_good;
      std::vector<double> phi_err_good;
      for (size_t i = 0; i < theta_v.size(); i++) {
	for (size_t j = i + 1; j < theta_v.size(); j++) {
	  // average the 2 results
	  theta = (theta_v[i] + theta_v[j]) / 2.;
	  phi   = (phi_v[i] + phi_v[j]) / 2.;
	  theta_err = sqrt ( (theta_v[i] - theta) * (theta_v[i] - theta) +
			     (theta_v[j] - theta) * (theta_v[j] - theta) );
	  phi_err   = sqrt ( (phi_v[i] - phi) * (phi_v[i] - phi) +
			     (phi_v[j] - phi) * (phi_v[j] - phi) );
	  // if the angle errors are within tolerance
	  if ( (theta_err < _max_err) and (phi_err < _max_err) ) {
	    theta_good.push_back(theta);
	    phi_good.push_back(phi);
	    theta_err_good.push_back(theta_err);
	    phi_err_good.push_back(phi_err);
	  }// if the angle errors are within tolerance
	}// for all j elements
      }// for all i elements
      // if the new vectors are empty -> return
      if (theta_good.size() == 0)
	return;
      // otherwise, average the results
      theta = phi = theta_err = phi_err = 0.;
      for (size_t n = 0; n < theta_good.size(); n++) {
	theta += theta_good[n];
	phi   += phi_good[n];
      }
      theta /= theta_good.size();
      phi   /= phi_good.size();
      for (size_t n = 0; n < theta_good.size(); n++) {
	theta_err += (theta_good[n] - theta) * (theta_good[n] - theta);
	phi_err += (phi_good[n] - phi) * (phi_good[n] - phi);
      }
      theta_err = sqrt(theta_err);
      phi_err   = sqrt(phi_err);
    }// if using 3 pairs
    
    if (_verbose)
      std::cout << "final angles : phi -> " << phi << ", theta -> " << theta << std::endl;
    
    resultShower.fDCosStart[1] = sin(theta);
    resultShower.fDCosStart[0] = cos(theta) * sin(phi);
    resultShower.fDCosStart[2] = cos(theta) * cos(phi);

    // calculate impact parameter with the vertex
    if (_validate_dir){
      
      // get the 3D start point of the shower
      auto const& start3D = resultShower.fXYZStart;
      
      if (proto_shower.hasVertex() == false){
	std::cout << "Number of vertices is not one!" << std::endl;
	return;
      }
      
      // get the proto-shower 3D vertex
      // take the vtx -> start point direction as the 3D direction
      auto const& vtx = proto_shower.vertex();
      
      std::vector<double> dir3D = {0,0,0};
      
      dir3D[0] = start3D[0] - vtx[0];
      dir3D[1] = start3D[1] - vtx[1];
      dir3D[2] = start3D[2] - vtx[2];
      
      // normalize
      double mag = sqrt( dir3D[0]*dir3D[0] + dir3D[1]*dir3D[1] + dir3D[2]*dir3D[2] );
      dir3D[0] /= mag;
      dir3D[1] /= mag;
      dir3D[2] /= mag;
      
      // ask that dot-product between vertex -> start and 3D direction is within some tolerance
      double dot = ( (dir3D[0] * resultShower.fDCosStart[0]) +
		     (dir3D[1] * resultShower.fDCosStart[1]) +
		     (dir3D[2] * resultShower.fDCosStart[2]) );

      if (_verbose){
	std::cout << "3D vtx   @   [" << vtx[0] << ", " << vtx[1] << ", " << vtx[2] << "]" << std::endl;
	std::cout << "3D start @   [" << start3D[0] << ", " << start3D[1] << ", " << start3D[2] << "]" << std::endl;
	std::cout << "Direction is [" << dir3D[0] << ", " << dir3D[1] << ", " << dir3D[2] << "]" << std::endl;
	std::cout << "Shr Dir   is [" << resultShower.fDCosStart[0] << ", " << resultShower.fDCosStart[1] << ", " << resultShower.fDCosStart[2] << "]" << std::endl;
	std::cout << "Dot product is : " << dot << std::endl;
      }

      if ( dot < _dot_min ){
	std::stringstream ss;
	ss << "Fail @ algo " << this->name() << " due to 3D dir mis-aligned with vertex";
	throw ShowerRecoException(ss.str());
      }
      
    }// if we should validate the direction
    
    return;
    
  }

  DEFINE_ART_CLASS_TOOL(Angle3DFormula)    
} //showerreco

#endif
