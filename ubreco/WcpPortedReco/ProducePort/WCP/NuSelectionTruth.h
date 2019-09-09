#ifndef NuSelectionTruth_h
#define NuSelectionTruth_h

namespace nsm{

  class NuSelectionTruth{
  public:
    NuSelectionTruth();
    NuSelectionTruth( bool isCC,
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
		      float electronInside );

    void SetIsCC(bool);
    void SetIsEligible(bool);
    void SetIsFC(bool);
    void SetVtxInside(bool);
    void SetNuPdg(int);
    void SetVtxX(float);
    void SetVtxY(float);
    void SetVtxZ(float);
    void SetTime(float);
    void SetNuEnergy(float);
    void SetEnergyInside(float);
    void SetElectronInside(float);

    const bool & GetIsCC() const;
    const bool & GetIsEligible() const;
    const bool & GetIsFC() const;
    const bool & GetIsVtxInside() const;
    const int & GetNuPdg() const;
    const float & GetVtxX() const;
    const float & GetVtxY() const;
    const float & GetVtxZ() const;
    const float & GetTime() const;
    const float & GetNuEnergy() const;
    const float & GetEnergyInside() const;
    const float & GetElectronInside() const;

  private:
    bool _isCC;
    bool _isEligible;
    bool _isFC;
    bool _vtxInside;
    int _nuPdg;
    float _vtxX;
    float _vtxY;
    float _vtxZ;
    float _time;
    float _nuEnergy;
    float _energyInside;
    float _electronInside;
  };
}

#endif

