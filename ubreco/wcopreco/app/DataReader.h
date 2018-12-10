#ifndef DATAREADER_H
#define DATAREADER_H

//WCOpReco includes
//data
#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "EventOpWaveforms.h"
#include "UBEventWaveform.h"


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
#include "TCanvas.h"

//c++ includes
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

namespace wcopreco  {

  // This class is designed to perform a read-in of microboone data from a
  // root file and organize it in datastructures designed by WCOpReco.
  class DataReader {
  public:
    DataReader(std::string *filepath);
    ~DataReader() ;

    UBEventWaveform Reader(int event_num);

    void LoopThroughWfms(std::vector<short> ch,
      std::vector<double> timestamp,
      TClonesArray Eventwaveform,
      int type,
      OpWaveformCollection &wfm_collection);

    //Make a bunch of root reading datamembers:
    Int_t nevents;
    TFile *file;
    TTree *tree;
    //These will be branches we'll need to read
    std::vector<short> * cosmic_hg_opch;
    std::vector<short> * cosmic_lg_opch;
    std::vector<short> * beam_hg_opch;
    std::vector<short> * beam_lg_opch;

    std::vector<double> * cosmic_hg_timestamp;
    std::vector<double> * cosmic_lg_timestamp;
    std::vector<double> * beam_hg_timestamp;
    std::vector<double> * beam_lg_timestamp;

    std::vector<float> * op_gain;
    std::vector<float> * op_gainerror;
    double triggerTime;
    double test;

    TClonesArray * cosmic_hg_wf;
    TClonesArray * cosmic_lg_wf;
    TClonesArray * beam_hg_wf;
    TClonesArray * beam_lg_wf;
    int eventNo;


    //General Datamembers
    UBEventWaveform _UB_Ev_wfm;
    int type;

  protected:

  };

}
#endif
