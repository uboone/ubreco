#ifndef CONFIG_UB_RC_H
#define CONFIG_UB_RC_H

namespace wcopreco{

  class Config_UB_rc {
  public:
    Config_UB_rc();
    ~Config_UB_rc() {};

    //UB_rc - parameter for ch 28
    double  _rc_tau_badch;
    //UB_rc - parameter for "good channels"
    double  _rc_tau_goodch;

    void _set_rc_tau_badch(double value) {_rc_tau_badch = value;}
    double _get_rc_tau_badch() {return _rc_tau_badch;}

    void _set_rc_tau_goodch(double value) {_rc_tau_goodch = value;}
    double _get_rc_tau_goodch() {return _rc_tau_goodch;}


  protected:

  };

}

#endif
