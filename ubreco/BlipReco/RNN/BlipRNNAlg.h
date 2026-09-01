#ifndef BLIPRNNALG_H
#define BLIPRNNALG_H

// Framework, libtorch, etc
#include "fhiclcpp/ParameterSet.h"
#include <torch/script.h>
#include <algorithm>

// Blip data types
#include "ubreco/BlipReco/Utils/BlipCore.h"


typedef std::vector<float>      vfloat_t;
typedef std::vector<vfloat_t>   vvfloat_t;

namespace blip {

  //--------------------------------------------
  class BlipRNNAlg {
    
    public:
    
      //Constructor/destructor
      BlipRNNAlg( fhicl::ParameterSet const& pset );
      BlipRNNAlg();
      ~BlipRNNAlg();
      
      bool loadModel();
      vfloat_t predict(vvfloat_t const&);
      vfloat_t predict(blipobj::Blip const&, std::vector<blipobj::HitInfo> const&);

      bool model_loaded;
      std::string loaded_model_path;

    private:
     
      // FCL configs
      std::string fModelFilename;
     
      // PyTorch things
      torch::jit::script::Module model;
      torch::Device device;
  
  };

}

#endif
