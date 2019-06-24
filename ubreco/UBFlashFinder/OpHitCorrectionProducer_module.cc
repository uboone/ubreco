#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RecoBase/OpHit.h"

#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"

#include "ubevt/Utilities/PMTRemapProvider.h"
#include "ubevt/Utilities/PMTRemapService.h"

#include "ubevt/Database/LightYieldService.h"
#include "ubevt/Database/LightYieldProvider.h"
#include "ubevt/Database/UbooneLightYieldProvider.h"

#include <memory>

class OpHitCorrectionProducer;


class OpHitCorrectionProducer : public art::EDProducer {
public:
  explicit OpHitCorrectionProducer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpHitCorrectionProducer(OpHitCorrectionProducer const &) = delete;
  OpHitCorrectionProducer(OpHitCorrectionProducer &&) = delete;
  OpHitCorrectionProducer & operator = (OpHitCorrectionProducer const &) = delete;
  OpHitCorrectionProducer & operator = (OpHitCorrectionProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:
  art::InputTag inputTag_;
  std::string outputProducer;
};


OpHitCorrectionProducer::OpHitCorrectionProducer(fhicl::ParameterSet const & p)
{
  inputTag_ = p.get<art::InputTag>("OpHitsInputTag");
  outputProducer = p.get<std::string>("OutputProducer","ophitcorrectionBeam");
  //
  produces<std::vector<recob::OpHit> >(outputProducer);

}

void OpHitCorrectionProducer::produce(art::Event & e)
{
  auto output = std::make_unique<std::vector<recob::OpHit> >();
  const auto& input = e.getValidHandle<std::vector<recob::OpHit> >(inputTag_);
  
  const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
  const lariov::LightYieldProvider& ly_provider = art::ServiceHandle<lariov::LightYieldService>()->GetProvider();
  const ::util::PMTRemapProvider &pmtremap_provider = art::ServiceHandle<util::PMTRemapService>()->GetProvider();

  for (const auto& oph : *input){

    auto oldch = pmtremap_provider.OriginalOpChannel(oph.OpChannel());
    float gaincor = gain_provider.ExtraInfo(oldch%100).GetFloatData("amplitude_gain")/20.;
    float lycor = ly_provider.LYScaling(oldch%100);
    output->emplace_back( recob::OpHit(oph.OpChannel(), 
				       oph.PeakTime(), 
				       oph.PeakTimeAbs(), 
				       oph.Frame(), 
				       oph.Width(), 
				       oph.Area(), 
				       oph.Amplitude(), 
				       oph.PE()*gaincor*lycor, 
				       oph.FastToTotal() ) );

  }
  e.put(std::move(output));
}

DEFINE_ART_MODULE(OpHitCorrectionProducer)
