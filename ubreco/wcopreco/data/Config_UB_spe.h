#ifndef CONFIG_UB_SPE_H
#define CONFIG_UB_SPE_H

namespace wcopreco{

  class Config_UB_spe {
  public:
    Config_UB_spe();
    ~Config_UB_spe() {};

    //UB_spe paramters
    double  _spe_p0;
    double  _spe_p1;

    void _set_spe_p0(double value) {_spe_p0 = value;}
    double _get_spe_p0() {return _spe_p0;}

    void _set_spe_p1(double value) {_spe_p1 = value;}
    double _get_spe_p1() {return _spe_p1;}


  protected:

  };

}

#endif
