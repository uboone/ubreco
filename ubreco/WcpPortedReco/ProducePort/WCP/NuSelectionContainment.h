#ifndef NuSelectionContainment_h
#define NuSelectionContainment_h

#include <vector>

namespace nsm{

  class NuSelectionContainment{
  public:
    NuSelectionContainment();
    NuSelectionContainment(
			   bool flash_found,
			   float flash_time,
			   float flash_measPe,
			   float flash_predPe,
			   bool found,
			   unsigned int match_type,
			   bool isFC,
			   bool isTGM,
			   bool notFC_FV,
			   bool notFC_SP,
			   bool notFC_DC);

    void SetFlashFound(bool flash_found);
    void SetFlashTime(float flash_time);
    void SetFlashMeasPe(float flash_measPe);
    void SetFlashPredPe(float flash_predPe);
    void SetMatchFound(bool found);
    void SetMatchType(unsigned int match_type);
    void SetIsFC(bool isFC);
    void SetIsTGM(bool isTGM);
    void SetNotFCFV(bool notFC_FV);
    void SetNotFCSP(bool notFC_SP);
    void SetNotFCDC(bool notFC_DC);

    const bool & GetFlashFound() const;
    const float & GetFlashTime() const;
    const float & GetFlashMeasPe() const;
    const float & GetFlashPredPe() const;
    const bool & GetMatchFound() const;
    const unsigned int & GetMatchType() const;
    const bool & GetIsFC() const;
    const bool & GetIsTGM() const;
    const bool & GetNotFCFV() const;
    const bool & GetNotFCSP() const;
    const bool & GetNotFCDC() const;

  private:
    bool _flash_found;
    float _flash_time;
    float _flash_measPe;
    float _flash_predPe;
    bool _found;
    unsigned int _match_type;
    bool _isFC;
    bool _isTGM;
    bool _notFC_FV;
    bool _notFC_SP;
    bool _notFC_DC;
  };
}

#endif
