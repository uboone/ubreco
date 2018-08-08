#ifndef SHOWERRECO_SHOWERRECOMANAGER_CXX
#define SHOWERRECO_SHOWERRECOMANAGER_CXX

#include "ShowerRecoManager.h"
#include <iomanip>

namespace showerreco {
  
  ShowerRecoManager::ShowerRecoManager()
  { }
  
  void ShowerRecoManager::Initialize()
  {
    for (auto & alg : _alg_v) {
      alg->initialize();
      _alg_time_v.push_back(0.);
      _alg_ctr_v.push_back(0);
    }
    
    return;
  }
  
  void ShowerRecoManager::Reset()
  {
    _proto_showers.clear();
    
    return;
  }
  
  void ShowerRecoManager::Reconstruct(std::vector<showerreco::Shower_t>& showers)
  {
    
    showers.clear();
    showers.reserve(_proto_showers.size());
    
    // for all pfparticle proto-showers
    for (auto const& proto_shower : _proto_showers) 
      showers.push_back(RecoOneShower(proto_shower));
    
    // Check that the showers reconstructed are the same length as the proto_showers vector
    if (showers.size() != _proto_showers.size()) {
      throw ShowerRecoException("ERROR: number of reconstructed showers doesn't match input list!!");
    }
    
    return;
  }

  ::showerreco::Shower_t ShowerRecoManager::RecoOneShower(const ::protoshower::ProtoShower& proto_shower)
  {
    
    // reset product shoer
    Shower_t result;
    Reset(result);
    
    // make a local copy of the shower to track differences
    Shower_t localCopy = result;


    // loop through reconstruction modules
    for (size_t n = 0; n < _alg_v.size(); n++) {
      
      _watch.Start();
      
      try {
	_alg_v[n] -> do_reconstruction(proto_shower, result);
      }// if reco succeeds
      catch (ShowerRecoException e) {
	//catch (std::exception e) {
	result.fPassedReconstruction = false;
	std::cout << e.what() << std::endl;
	return result;
      }// if reco fails
      _alg_time_v[n] += _watch.RealTime();
      _alg_ctr_v[n] += 1;
      if (_debug && _verbose) {
	printChanges(localCopy, result, _alg_v[n]->name());
	localCopy = result;
      }// if verbose
    }// for all reconstruction modules
  
    // if we made it this far, the shower is good!
    result.fPassedReconstruction = true;    

    return result;
  }
  
  void ShowerRecoManager::Reset(Shower_t& result) {
    
    size_t nPlanes = 3;
    
    result.fTotalEnergy_v.resize(nPlanes);
    result.fSigmaTotalEnergy_v.resize(nPlanes);
    result.fdEdx_v.resize(nPlanes);
    result.fdQdx_v.resize(nPlanes);
    result.fSigmadEdx_v.resize(nPlanes);
    result.fSigmadQdx_v.resize(nPlanes);
    result.fHitdQdx_v.resize(nPlanes);
    
    result.fShoweringLength.resize(nPlanes); // resizing Showering Length Vector
    result.fTotalMIPEnergy_v.resize(nPlanes);
    result.fSigmaTotalMIPEnergy_v.resize(nPlanes);
    
    result.fPlaneIsBad.resize(nPlanes);

    result.fPassedReconstruction = true;

    return;
  }

void ShowerRecoManager::PrintModuleList() {
  
  std::cout << "Print the list of modules to run in Shower Reco Alg Modular:\n";
  int i = 0;
  for (auto & alg : _alg_v) {
    std::cout << "\t" << i << ") " << alg -> name() << "\n";
    i++;
  }
  
}

  void ShowerRecoManager::printChanges(const Shower_t & localCopy,
				     const Shower_t result,
				     std::string moduleName) {
  

  bool changed;

  // Look at each value of Shower_t and if it has changed, print out that change
  std::cout << "\nPrinting the list of changes made by module " << moduleName << std::endl;

  // Compare vectors by x/y/z/ values
  // Doing comparisons in the order of values in Shower_t
  // Cos Start:
  if (localCopy.fDCosStart.X() != result.fDCosStart.X() ||
      localCopy.fDCosStart.Y() != result.fDCosStart.Y() ||
      localCopy.fDCosStart.X() != result.fDCosStart.X() ) {
    std::cout << "\tfDCosStart has changed from ("
              << localCopy.fDCosStart.X() << ", "
              << localCopy.fDCosStart.Y() << ", "
              << localCopy.fDCosStart.Z() << ") to ("
              << result.fDCosStart.X() << ", "
              << result.fDCosStart.Y() << ", "
              << result.fDCosStart.Z() << ") "
              << std::endl;
  }

  // Sigma Cos Start:
  if (localCopy.fSigmaDCosStart.X() != result.fSigmaDCosStart.X() ||
      localCopy.fSigmaDCosStart.Y() != result.fSigmaDCosStart.Y() ||
      localCopy.fSigmaDCosStart.X() != result.fSigmaDCosStart.X() ) {
    std::cout << "\tfSigmaDCosStart has changed from ("
              << localCopy.fSigmaDCosStart.X() << ", "
              << localCopy.fSigmaDCosStart.Y() << ", "
              << localCopy.fSigmaDCosStart.Z() << ") to ("
              << result.fSigmaDCosStart.X() << ", "
              << result.fSigmaDCosStart.Y() << ", "
              << result.fSigmaDCosStart.Z() << ") "
              << std::endl;
  }

  // XYZ Start:
  if (localCopy.fXYZStart.X() != result.fXYZStart.X() ||
      localCopy.fXYZStart.Y() != result.fXYZStart.Y() ||
      localCopy.fXYZStart.X() != result.fXYZStart.X() ) {
    std::cout << "\tfXYZStart has changed from ("
              << localCopy.fXYZStart.X() << ", "
              << localCopy.fXYZStart.Y() << ", "
              << localCopy.fXYZStart.Z() << ") to ("
              << result.fXYZStart.X() << ", "
              << result.fXYZStart.Y() << ", "
              << result.fXYZStart.Z() << ") "
              << std::endl;
  }

  // Sigma XYZ Start
  if (localCopy.fSigmaXYZStart.X() != result.fSigmaXYZStart.X() ||
      localCopy.fSigmaXYZStart.Y() != result.fSigmaXYZStart.Y() ||
      localCopy.fSigmaXYZStart.X() != result.fSigmaXYZStart.X() ) {
    std::cout << "\tfSigmaXYZStart has changed from ("
              << localCopy.fSigmaXYZStart.X() << ", "
              << localCopy.fSigmaXYZStart.Y() << ", "
              << localCopy.fSigmaXYZStart.Z() << ") to ("
              << result.fSigmaXYZStart.X() << ", "
              << result.fSigmaXYZStart.Y() << ", "
              << result.fSigmaXYZStart.Z() << ") "
              << std::endl;
  }

  // Length
  if (localCopy.fLength != result.fLength) {
    std::cout << "\tfDCosStart has changed from " << localCopy.fLength
              << " to "  << result.fLength << std::endl;
  }

  // BestdEdx
  if (localCopy.fBestdEdx != result.fBestdEdx) {
    std::cout << "\tfBestdEdx has changed from " << localCopy.fBestdEdx
              << " to "  << result.fBestdEdx << std::endl;
  }

  // Opening Angle
  if (localCopy.fOpeningAngle != result.fOpeningAngle) {
    std::cout << "\tfDCosStart has changed from " << localCopy.fOpeningAngle
              << " to "  << result.fOpeningAngle << std::endl;
  }

  // Since these should be length = nplanes, checking each plane for change
  // If changed, print out

  // Total Energy
  changed = false;
  if (localCopy.fTotalEnergy != result.fTotalEnergy)
    changed = true;
  if (changed) {
    std::cout << "\tfTotalEnergy has changed from "
              << localCopy.fTotalEnergy << " to "
              << result.fTotalEnergy << std::endl;
  }

  //Total Energy vector
  changed = false;
  for (unsigned int i = 0; i < localCopy.fTotalEnergy_v.size(); i++) {
    if (localCopy.fTotalEnergy_v[i] != result.fTotalEnergy_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfTotalEnergy_v has changed from (";
    for (auto & val : localCopy.fTotalEnergy_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fTotalEnergy_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }



  // Sigma Total Energy
  changed = false;
  if (localCopy.fSigmaTotalEnergy != result.fSigmaTotalEnergy)
    changed = true;
  if (changed) {
    std::cout << "\tfSigmaTotalEnergy has changed from "
              << localCopy.fSigmaTotalEnergy << " to "
              << result.fSigmaTotalEnergy << std::endl;
  }

  // Sigma Total Energy vector
  changed = false;
  for (unsigned int i = 0; i < localCopy.fSigmaTotalEnergy_v.size(); i++) {
    if (localCopy.fSigmaTotalEnergy_v[i] != result.fSigmaTotalEnergy_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfSigmaTotalEnergy_v has changed from (";
    for (auto & val : localCopy.fSigmaTotalEnergy_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fSigmaTotalEnergy_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }

  // dEdx
  changed = false;
  if (localCopy.fdEdx != result.fdEdx)
    changed = true;
  if (changed) {
    std::cout << "\tfdEdx has changed from "
              << localCopy.fdEdx << " to "
              << result.fdEdx << std::endl;
  }

  // dEdx_v
  changed = false;
  for (unsigned int i = 0; i < localCopy.fdEdx_v.size(); i++) {
    if (localCopy.fdEdx_v[i] != result.fdEdx_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfdEdx_v has changed from (";
    for (auto & val : localCopy.fdEdx_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fdEdx_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }

  // dQdx
  changed = false;
  if (localCopy.fdQdx_v != result.fdQdx_v)
    changed = true;
  if (changed) {
    std::cout << "\tfdQdx has changed from"
              << localCopy.fdQdx << " to "
              << result.fdQdx << std::endl;
  }

  // dQdx_v
  changed = false;
  for (unsigned int i = 0; i < localCopy.fdQdx_v.size(); i++) {
    if (localCopy.fdQdx_v[i] != result.fdQdx_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfdQdx_v has changed from (";
    for (auto & val : localCopy.fdQdx_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fdQdx_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }

  // HitdQdx_v
  changed = false;
  for (unsigned int i = 0; i < localCopy.fHitdQdx_v.size(); i++) {
    for (unsigned int j = 0; j < localCopy.fHitdQdx_v[i].size(); j++) {
      if (localCopy.fHitdQdx_v[i][j] != result.fHitdQdx_v[i][j]) {
        changed = true;
        break;
      }
    }
  }
  if (changed) {
    std::cout << "\tfHitdQdx_v has changed from (";
    for (auto & val : localCopy.fHitdQdx_v )
      for (auto & v : val) std::cout << v << " ";
    std::cout << ") to (";
    for (auto & val : result.fHitdQdx_v )
      for (auto & v : val )std::cout << v << " ";
    std::cout << ")" << std::endl;
  }


  // sigma dEdx
  changed = false;
  if (localCopy.fSigmadEdx != result.fSigmadEdx)
    changed = true;
  if (changed) {
    std::cout << "\tfSigmadEdx has changed from "
              << localCopy.fSigmadEdx << " to "
              << result.fSigmadEdx << std::endl;
  }

  // sigma dEdx_v
  changed = false;
  for (unsigned int i = 0; i < localCopy.fSigmadEdx_v.size(); i++) {
    if (localCopy.fSigmadEdx_v[i] != result.fSigmadEdx_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfSigmadEdx_v has changed from (";
    for (auto & val : localCopy.fSigmadEdx_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fSigmadEdx_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }


  // sigma dQdx
  changed = false;
  if (localCopy.fSigmadQdx != result.fSigmadQdx)
    changed = true;
  if (changed) {
    std::cout << "\tfSigmadQdx has changed from "
              << localCopy.fSigmadQdx << " to "
              << result.fSigmadQdx << std::endl;
  }

  // sigma dQdx_v
  changed = false;
  for (unsigned int i = 0; i < localCopy.fSigmadQdx_v.size(); i++) {
    if (localCopy.fSigmadQdx_v[i] != result.fSigmadQdx_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfSigmadQdx_v has changed from (";
    for (auto & val : localCopy.fSigmadQdx_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fSigmadQdx_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }

  // Total MIP Energy
  changed = false;
  for (unsigned int i = 0; i < localCopy.fTotalMIPEnergy_v.size(); i++) {
    if (localCopy.fTotalMIPEnergy_v[i] != result.fTotalMIPEnergy_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfTotalMIPEnergy has changed from (";
    for (auto & val : localCopy.fTotalMIPEnergy_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fTotalMIPEnergy_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }

  // Sigma Total MIP Energy
  changed = false;
  for (unsigned int i = 0; i < localCopy.fSigmaTotalMIPEnergy_v.size(); i++) {
    if (localCopy.fSigmaTotalMIPEnergy_v[i] != result.fSigmaTotalMIPEnergy_v[i]) {
      changed = true;
      break;
    }
  }
  if (changed) {
    std::cout << "\tfSigmaTotalMIPEnergy has changed from (";
    for (auto & val : localCopy.fSigmaTotalMIPEnergy_v ) std::cout << val << " ";
    std::cout << ") to (";
    for (auto & val : result.fSigmaTotalMIPEnergy_v ) std::cout << val << " ";
    std::cout << ")" << std::endl;
  }

  if (localCopy.fBestPlane != result.fBestPlane) {
    std::cout << "\tfBestPlane has changed from " << localCopy.fBestPlane.Plane
              << " to "  << result.fBestPlane.Plane << std::endl;
  }


  std::cout << std::endl;

}


// finalize function
void ShowerRecoManager::Finalize(TFile* fout)
{

  // loop through algos and evaluate time-performance
  std::cout << std::endl
            << "=================== Time Report =====================" << std::endl;
  for (size_t n = 0; n < _alg_v.size(); n++) {
    double alg_time = _alg_time_v[n] / ((double)_alg_ctr_v[n]);
    std::cout <<  std::setw(25) << _alg_v[n]->name() << "\t Algo Time: "
              << std::setw(10) << alg_time * 1.e6     << " [us/proto_shower]"
              << "\t Proto-Showers Scanned: " << _alg_ctr_v[n] << std::endl;
  }

  std::cout << "=====================================================" << std::endl
            << std::endl;

  return;
}



}

#endif
