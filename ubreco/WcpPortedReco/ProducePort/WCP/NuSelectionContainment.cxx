#include "NuSelectionContainment.h"
#include <vector>

#include <iostream>

  nsm::NuSelectionContainment::NuSelectionContainment(){
    _flash_found=false;
    _flash_time=-1.;
    _flash_measPe=-1.;
    _flash_predPe=-1.;
    _found=false;
    _match_type=0;
    _isFC=false;
    _isTGM=false;
    _notFC_FV=false;
    _notFC_SP=false;
    _notFC_DC=false;
  }

  nsm::NuSelectionContainment::NuSelectionContainment(bool flash_found,
						 float flash_time,
						 float flash_measPe,
						 float flash_predPe,
						 bool found,
						 unsigned int match_type,
						 bool isFC,
						 bool isTGM,
						 bool notFC_FV,
						 bool notFC_SP,
						 bool notFC_DC){
    _flash_found=flash_found;
    _flash_time=flash_time;
    _flash_measPe=flash_measPe;
    _flash_predPe=flash_predPe;
    _found=found;
    _match_type=match_type;
    _isFC=isFC;
    _isTGM=isTGM;
    _notFC_FV=notFC_FV;
    _notFC_SP=notFC_SP;
    _notFC_DC=notFC_DC;
  }

  void nsm::NuSelectionContainment::SetFlashFound(bool flash_found){ this->_flash_found = flash_found; }
  void nsm::NuSelectionContainment::SetFlashTime(float flash_time){ this->_flash_time = flash_time; }
  void nsm::NuSelectionContainment::SetFlashMeasPe(float flash_measPe){ this->_flash_measPe = flash_measPe; }
  void nsm::NuSelectionContainment::SetFlashPredPe(float flash_predPe){ this->_flash_predPe = flash_predPe; }
  void nsm::NuSelectionContainment::SetMatchType(unsigned int match_type){ this->_match_type = match_type; }
  void nsm::NuSelectionContainment::SetMatchFound(bool found){ this->_found = found; }
  void nsm::NuSelectionContainment::SetIsFC(bool isFC){ this->_isFC = isFC; }
  void nsm::NuSelectionContainment::SetIsTGM(bool isTGM){ this->_isTGM = isTGM; }
  void nsm::NuSelectionContainment::SetNotFCFV(bool notFC_FV){ this->_notFC_FV = notFC_FV; }
  void nsm::NuSelectionContainment::SetNotFCSP(bool notFC_SP){ this->_notFC_SP = notFC_SP; }
  void nsm::NuSelectionContainment::SetNotFCDC(bool notFC_DC){ this->_notFC_DC = notFC_DC; }
  
  const bool  & nsm::NuSelectionContainment::GetFlashFound() const { return this->_flash_found; }
  const float & nsm::NuSelectionContainment::GetFlashTime() const { return this->_flash_time; }
  const float & nsm::NuSelectionContainment::GetFlashMeasPe() const { return this->_flash_measPe; }
  const float & nsm::NuSelectionContainment::GetFlashPredPe() const { return this->_flash_predPe; }
  const unsigned int & nsm::NuSelectionContainment::GetMatchType() const { return this->_match_type; }
  const bool & nsm::NuSelectionContainment::GetMatchFound() const { return this->_found; }
  const bool & nsm::NuSelectionContainment::GetIsFC() const { return this->_isFC; }
  const bool & nsm::NuSelectionContainment::GetIsTGM() const { return this->_isTGM; }
  const bool & nsm::NuSelectionContainment::GetNotFCFV() const { return this->_notFC_FV; }
  const bool & nsm::NuSelectionContainment::GetNotFCSP() const { return this->_notFC_SP; }
  const bool & nsm::NuSelectionContainment::GetNotFCDC() const { return this->_notFC_DC; }


