#ifndef COphit_h
#define COphit_h

#include "OpWaveform.h"
#include "Config_Params.h"
#include <vector>
#include <math.h>

#include "Config_COpHit.h"


namespace wcopreco{

  class COphit{
  public:
    COphit(int ch_no, OpWaveform *wfm, double time, double gain, double gain_err, const Config_COpHit &config);
    ~COphit();

    double get_time() const {return time;};
    double get_baseline() const {return baseline;};
    double get_peak() const {return peak;};
    double get_integral() const {return integral;};
    double get_gain() const {return gain;}
    int get_ch_no() const {return channel_no;};

    // derived
    double get_PE() const {return PE;};
    double get_PE_err() const {return PE_err;};
    bool get_type()const {return good_baseline;};


  protected:
    double cal_integral(double peak);
    bool good_baseline;

    int channel_no; // FEM channel number
    double time; // start time in us
    double gain; // to be set
    double gain_err;

    double baseline; // baseline
    double peak; // maximum PE
    double integral; // integral

    double PE;
    double PE_err;

    Config_COpHit _cfg;

  };

  typedef std::vector<COphit*> COphitSelection;

}

#endif
