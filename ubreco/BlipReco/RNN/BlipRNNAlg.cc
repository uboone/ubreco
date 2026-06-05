#include "ubreco/BlipReco/RNN/BlipRNNAlg.h"

namespace blip {

  //###########################################################
  // Constructor
  //###########################################################
  BlipRNNAlg::BlipRNNAlg( fhicl::ParameterSet const& pset ):
    device(torch::kCPU)
  {
    fModelFilename   = pset.get<std::string>   ("ModelFilename", "model_blipRNN.pt");
    model_loaded=false;
    loadModel();
  }
  
  //====================================================
  // Destructor
  BlipRNNAlg::~BlipRNNAlg()
  {
  }
  
  //====================================================
  // Load model
  bool BlipRNNAlg::loadModel() {
    // Skip if already loaded with same path
    if (model_loaded && loaded_model_path == fModelFilename ) {
      std::cout << "Model already loaded, skipping...\n";
      return true;
    }
    try {
      cet::search_path sp("FW_SEARCH_PATH");
      std::cout<<"Searching for file...\n";
      std::cout<<sp.find_file(fModelFilename)<<"\n";
      model = torch::jit::load(sp.find_file(fModelFilename));
      model.to(device);
      model.eval();
      model_loaded = true;
      loaded_model_path = fModelFilename;
      std::cout << "Model loaded: " << loaded_model_path << std::endl;
      return true;
    } catch (...) {
        std::cout << "Failed to load model "<<fModelFilename<<"\n";
        return false;
      }
  }
  
  //====================================================
  // Predict blip 3D direction
  std::vector<float> BlipRNNAlg::predict(vvfloat_t const& blip_hit_seq){
    if( !model_loaded ) return {};
    try{
      
      int nFeatures=5;
      int nHits=blip_hit_seq.size();
      auto bliphits = torch::zeros({1, nHits, nFeatures});
      for(int i=0; i<nHits; i++){
        for(int j=0; j<nFeatures; j++) 
          bliphits[0][i][j] = blip_hit_seq[i][j];
      }
     
      // create input
      auto length = torch::tensor({nHits}, torch::kLong);
      torch::NoGradGuard no_grad;
      std::vector<torch::jit::IValue> input;
      input.push_back(bliphits);
      input.push_back(length);
      
      // feed it to the model
      at::Tensor output = model.forward(input).toTensor();

      // get the result and normalize magnitude
      auto accessor = output.accessor<float, 2>();
      TVector3 vdir(accessor[0][0],accessor[0][1],accessor[0][2]);
      vdir.SetMag(1.);
      
      return { (float)vdir.X(), (float)vdir.Y(), (float)vdir.Z() };
    
    } catch (const std::exception& e) {
      std::cerr << "RNN Prediction failed: "<<e.what()<<std::endl;
      return {};
    }
  }

  //============================================================
  // This function takes in a blipobj::Blip object and pre-processed the
  // hit information to feed into the function above
  std::vector<float> BlipRNNAlg::predict(blipobj::Blip const& blp, std::vector<blipobj::HitInfo> const& hitinfo){
    
    // 3D direction (initialized to dummy)
    std::vector<float> dir{-9,-9,-9};

    // Basic cuts: require blip be valid
    if( !blp.isValid ) return dir;

    // Conversion factor for casting ticks into units of
    // wire spacings (TODO: make this not hard-coded?)
    float tick_to_us  = 0.5; // microsecond
    float drift_speed = 1.1; // mm per microsecond
    float wire_sep    = 3.0; // mm
    float convfact    = tick_to_us * drift_speed / wire_sep;
    
    // Collections of hit info info per plane
    std::array<std::vector<float>, kNplanes > blip_pl_wire;   // [plane] -> vector of wires
    std::array<std::vector<float>, kNplanes > blip_pl_time;
    std::array<std::vector<float>, kNplanes > blip_pl_rms;
    std::array<std::vector<float>, kNplanes > blip_pl_amp;

    // Loop over the planes
    int nhits = 0;
    float sum_driftT = 0;
    float min_driftT = 999999;
    for(size_t i=0; i<kNplanes; i++){
      auto& clust = blp.clusters[i];
      
      // loop this blip's hits on this plane
      for(auto hid : clust.HitIDs ){
        float t = hitinfo[hid].driftTime * convfact;
        float rms = hitinfo[hid].rms * convfact;
        float amp = hitinfo[hid].amp;
        blip_pl_wire[i].push_back( hitinfo[hid].wire );
        blip_pl_time[i].push_back( t );
        blip_pl_rms[i] .push_back( rms );
        blip_pl_amp[i] .push_back( amp );
        nhits++;
        sum_driftT += t;
        if( t < min_driftT ) min_driftT = t;
      }
      
      // basic sanity check
      if( blip_pl_wire[i].empty() ) continue;

      // make corrections to wires on each plane
      float max_wire = *std::max_element(blip_pl_wire[i].begin(), blip_pl_wire[i].end());
      for(auto &wn : blip_pl_wire[i]) wn = max_wire - wn;
    }

    // If there were 2 hits or less, forget about it
    if( nhits <= 2 ) return dir;
    
    // Ok we got all the hits. Shift the drift times
    for(size_t i=0; i<kNplanes; i++){
      for(auto& t : blip_pl_time[i] ) t -= min_driftT;
    }
    
    // create vector of hit info needed for RNN inference
    // order of hit features: rms, amp, wire, plane, time
    vvfloat_t blip_hit_seq;
    for(size_t i=0; i<kNplanes; i++){
      for(size_t j=0; j<blip_pl_wire[i].size(); j++){
        std::vector<float> hit_features;
        hit_features.push_back( blip_pl_rms[i][j] );
        hit_features.push_back( blip_pl_amp[i][j] );
        hit_features.push_back( blip_pl_wire[i][j] );
        hit_features.push_back( i );
        hit_features.push_back( blip_pl_time[i][j] );
        blip_hit_seq.push_back( hit_features );
      }
    }
    
    return predict(blip_hit_seq);
  
  }

}


