#ifndef Opflash_h
#define Opflash_h

#include "COphit.h"
#include "Config_Opflash.h"
#include <math.h>
#include <set>

namespace wcopreco{
  class Opflash{
  public:
    Opflash(COphitSelection &ophits, const Config_Opflash &configOpF);
    Opflash(const std::vector<std::vector<double>> &vec_v, double start_time, int start_bin, int end_bin, const Config_Opflash &configOpF);
    ~Opflash();

    void Add_l1info(std::vector<double>* vec1, std::vector<double> *vec2, double start_time , int start_bin, int end_bin, const Config_Opflash &configOpF);

    void set_flash_id(int value){flash_id = value;};
    int get_flash_id() const {return flash_id;};

    double get_time() const {return time;};
    double get_total_PE() const {return total_PE;};
    double get_PE(int ch) const {return PE[ch];};
    std::vector<double> get_pe_v() const {return PE;};
    std::vector<double> get_pe_v_nocor() const {return PEnc;};
    double get_PE_err(int ch) const {return PE_err[ch];};
    bool get_fired(int ch);
    int get_num_fired() const {return fired_channels.size();};

    int get_type() const {return type;}
    double get_low_time() const {return low_time;};
    double get_high_time() const {return high_time;};

    std::vector<double>& get_l1_fired_time() {return l1_fired_time;};
    std::vector<double>& get_l1_fired_pe() {return l1_fired_pe;};

    void swap_channels();



  protected:

    int type;
    int flash_id;
    Config_Opflash _cfgOpF;
    double low_time;
    double high_time;

    double time;
    double total_PE;
    std::vector<int> fired_channels;
    std::vector<double> PE;
    std::vector<double> PEnc;
    std::vector<double> PE_err;

    std::vector<double> l1_fired_time;
    std::vector<double> l1_fired_pe;
  };

  struct OpFlashCompare{
    bool operator() (Opflash *a, Opflash *b) const{
      if (a->get_time() < b->get_time()){
	return true;
      }else{
	return false;
      }
    }
  };

  typedef std::vector<Opflash*> OpflashSelection;
  typedef std::set<Opflash*, OpFlashCompare> OpFlashSet;

}

#endif
