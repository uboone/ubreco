#ifndef UBEVENTWAVEFORM_H
#define UBEVENTWAVEFORM_H

//WCOpReco includes
#include "EventOpWaveforms.h"
#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "EventOpWaveforms.h"

//root includes
#include "TObject.h"
#include "TH1S.h"
#include "TF1.h"
#include "TVirtualFFT.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TTimeStamp.h"

//c++ includes
#include <string>
#include <iostream>
#include <sstream>


namespace wcopreco {

  typedef enum {kbeam_hg = 0, kbeam_lg, kcosmic_hg, kcosmic_lg, kbeam_merged, kcosmic_merged, kNumTypes } UBOpWaveformForm_t;

  class UBEventWaveform : public EventOpWaveforms {
  public:


    UBEventWaveform();
    virtual ~UBEventWaveform() {};

    void addWaveform( UBOpWaveformForm_t type, const OpWaveform& wfm );
    std::vector<float> get_op_gain() {return op_gain;}
    void set_op_gain(std::vector<float> gains_v) {op_gain = gains_v;}

    std::vector<float> get_op_gainerror() {return op_gainerror;}
    void set_op_gainerror(std::vector<float> gainserror_v) {op_gainerror = gainserror_v;}

  protected:
    std::vector<float> op_gain;
    std::vector<float> op_gainerror;


  };

}

#endif
