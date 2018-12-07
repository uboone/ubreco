#ifndef CONFIG_HITFINDER_COSMIC_H
#define CONFIG_HITFINDER_COSMIC_H

namespace wcopreco{

  class Config_Hitfinder_Cosmic {
  public:
    Config_Hitfinder_Cosmic();
    ~Config_Hitfinder_Cosmic() {};

    //min time diff to be considered a separate hit (in microseconds)
    double  _ophit_group_t_diff_max;

    void _set_ophit_group_t_diff_max(double value) {_ophit_group_t_diff_max = value;}
    double _get_ophit_group_t_diff_max() {return _ophit_group_t_diff_max;}



  protected:

  };

}

#endif
