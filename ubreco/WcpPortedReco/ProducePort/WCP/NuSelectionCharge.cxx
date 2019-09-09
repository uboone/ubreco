#include "NuSelectionCharge.h"
#include "vector"

using namespace nsm;

  NuSelectionCharge::NuSelectionCharge() {
    _charge_measU=-1.;
    _charge_measV=-1.;
    _charge_measY=-1.;
  }

  NuSelectionCharge::NuSelectionCharge(float charge_measU,
				       float charge_measV,
				       float charge_measY){
    _charge_measU=charge_measU;
    _charge_measV=charge_measV;
    _charge_measY=charge_measY;
  }

  void NuSelectionCharge::SetChargeU(float charge_measU){ this->_charge_measU = charge_measU; }
  void NuSelectionCharge::SetChargeV(float charge_measV){ this->_charge_measV = charge_measV; }
  void NuSelectionCharge::SetChargeY(float charge_measY){ this->_charge_measY = charge_measY; }
  
  const float & NuSelectionCharge::GetChargeU() const { return this->_charge_measU; }
  const float & NuSelectionCharge::GetChargeV() const { return this->_charge_measV; }
  const float & NuSelectionCharge::GetChargeY() const { return this->_charge_measY; }
