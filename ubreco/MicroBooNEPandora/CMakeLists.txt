cet_make_library(
  SOURCE
  DummyVertexParticleCreationAlgorithm.cc
  MicroBooNEContent.cc
  MicroBooNEExampleAlgorithm.cc
  MicroBooNEExternalVertexCreationAlgorithm.cc
  MicroBooNEMasterAlgorithm.cc
  MicroBooNEPreProcessingAlgorithm.cc
  LIBRARIES
  PUBLIC
  larpandoracontent::LArPandoraContent
  PandoraPFA::PandoraSDK
)

cet_build_plugin(
  MicroBooNEPandora art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::MicroBooNEPandora
  larpandora::LArPandoraInterface
)

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)

