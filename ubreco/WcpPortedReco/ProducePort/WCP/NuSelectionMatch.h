#ifndef NuSelectionMatch_h
#define NuSelectionMatch_h

#include <vector>

namespace nsm{

  class NuSelectionMatch{
  public:
    NuSelectionMatch();
    NuSelectionMatch(float completeness,
		     float completeness_energy,
		     float purity,
		     float purity_xy,
		     float purity_xz,
		     float charge,
		     float energy);
    //virtual ~NuSelectionMatch();

    void SetCompleteness(float);
    void SetCompletenessEnergy(float);
    void SetPurity(float);
    void SetPurityXY(float);
    void SetPurityXZ(float);
    void SetCharge(float);
    void SetEnergy(float);

    const float & GetCompleteness() const;
    const float & GetCompletenessEnergy() const;
    const float & GetPurity() const;
    const float & GetPurityXY() const;
    const float & GetPurityXZ() const;
    const float & GetCharge() const;
    const float & GetEnergy() const;

  private:
    float _completeness;
    float _completeness_energy;
    float _purity;
    float _purity_xy;
    float _purity_xz;
    float _charge;
    float _energy;
  };
}

#endif

