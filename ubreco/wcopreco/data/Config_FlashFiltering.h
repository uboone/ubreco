#ifndef CONFIG_FLASHFILTERING_H
#define CONFIG_FLASHFILTERING_H

namespace wcopreco{

  class Config_FlashFiltering {
  public:
    Config_FlashFiltering();
    ~Config_FlashFiltering() {};

    //time threshold for filtering
    double  _flash_filter_time_thresh;
    //pe threshold for filtering
    double  _flash_filter_pe_thresh;
    //true if need to swap channels
    bool    _do_swap_channels;

    void _set_flash_filter_time_thresh(double value) {_flash_filter_time_thresh = value;}
    double _get_flash_filter_time_thresh() {return _flash_filter_time_thresh;}

    void _set_flash_filter_pe_thresh(double value) {_flash_filter_pe_thresh = value;}
    double _get_flash_filter_pe_thresh() {return _flash_filter_pe_thresh;}
    
    void _set_do_swap_channels(bool do_swap) {_do_swap_channels = do_swap;}
    bool _get_do_swap_channels() {return _do_swap_channels;}


  protected:

  };

}

#endif
