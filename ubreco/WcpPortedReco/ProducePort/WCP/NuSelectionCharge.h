#ifndef NuSelectionCharge_h
#define NuSelectionCharge_h

namespace nsm{

  class NuSelectionCharge{
  public:
    NuSelectionCharge();
    NuSelectionCharge(float charge_measU,
		      float charge_measV,
		      float charge_measY);

    void SetChargeU(float);
    void SetChargeV(float);
    void SetChargeY(float);

    const float & GetChargeU() const;
    const float & GetChargeV() const;
    const float & GetChargeY() const;

  private:
    float _charge_measU;
    float _charge_measV;
    float _charge_measY;
  };
}

#endif
