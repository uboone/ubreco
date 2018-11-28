#include "UBEventWaveform.h"

namespace wcopreco {

   UBEventWaveform::UBEventWaveform()  {


  }

  void UBEventWaveform::addWaveform( UBOpWaveformForm_t type, const OpWaveform& wfm ) {
    push_back_wfm( (int)type, wfm );
  }

}
