#include "NuSelectionTruth.h"
#include "vector"

using namespace nsm;

  NuSelectionTruth::NuSelectionTruth(){
    _isCC=false;
    _isEligible=false;
    _isFC=false;
    _vtxInside=false;
    _nuPdg=-1;
    _vtxX=-1.;
    _vtxY=-1.;
    _vtxZ=-1.;
    _time=-1.;
    _nuEnergy=-1.;
    _energyInside=-1.;
    _electronInside=-1.;
  }

  NuSelectionTruth::NuSelectionTruth( bool isCC,
				      bool isEligible,
				      bool isFC,
				      bool vtxInside,
				      int nuPdg,
				      float vtxX,
				      float vtxY,
				      float vtxZ,
				      float time,
				      float nuEnergy,
				      float energyInside,
				      float electronInside ){
    _isCC=isCC;
    _isEligible=isEligible;
    _isFC=isFC;
    _vtxInside=vtxInside;
    _nuPdg=nuPdg;
    _vtxX=vtxX;
    _vtxY=vtxY;
    _vtxZ=vtxZ;
    _time=time;
    _nuEnergy=nuEnergy;
    _energyInside=energyInside;
    _electronInside=electronInside;
  }

  void NuSelectionTruth::SetIsCC(bool isCC){ this->_isCC = isCC; }
  void NuSelectionTruth::SetIsEligible(bool isEligible){ this->_isEligible = isEligible; }
  void NuSelectionTruth::SetIsFC(bool isFC){ this->_isFC = isFC; }
  void NuSelectionTruth::SetVtxInside(bool vtxInside){ this->_vtxInside = vtxInside; }
  void NuSelectionTruth::SetNuPdg(int nuPdg){ this->_nuPdg = nuPdg; }
  void NuSelectionTruth::SetVtxX(float vtxX){ this->_vtxX = vtxX; }
  void NuSelectionTruth::SetVtxY(float vtxY){ this->_vtxY = vtxY; }
  void NuSelectionTruth::SetVtxZ(float vtxZ){ this->_vtxZ = vtxZ; }
  void NuSelectionTruth::SetTime(float time){ this->_time = time; }
  void NuSelectionTruth::SetNuEnergy(float nuEnergy){ this->_nuEnergy = nuEnergy; }
  void NuSelectionTruth::SetEnergyInside(float energyInside){ this->_energyInside = energyInside; }
  void NuSelectionTruth::SetElectronInside(float electronInside){ this->_electronInside = electronInside; }

  const bool & NuSelectionTruth::GetIsCC() const { return this->_isCC; }
  const bool & NuSelectionTruth::GetIsEligible() const { return this->_isEligible; }
  const bool & NuSelectionTruth::GetIsFC() const { return this->_isFC; }
  const bool & NuSelectionTruth::GetIsVtxInside() const { return this->_vtxInside; }
  const int & NuSelectionTruth::GetNuPdg() const { return this->_nuPdg; }
  const float & NuSelectionTruth::GetVtxX() const { return this->_vtxX; }
  const float & NuSelectionTruth::GetVtxY() const { return this->_vtxY; }
  const float & NuSelectionTruth::GetVtxZ() const { return this->_vtxZ; }
  const float & NuSelectionTruth::GetTime() const { return this->_time; }
  const float & NuSelectionTruth::GetNuEnergy() const { return this->_nuEnergy; }
  const float & NuSelectionTruth::GetEnergyInside() const { return this->_energyInside; }
  const float & NuSelectionTruth::GetElectronInside() const { return this->_electronInside; }


