cet_make_library(
  SOURCE
  DummyVertexParticleCreationAlgorithm.cc
  MicroBooNEContent.cc
  MicroBooNEDLContent.cc
  MicroBooNEExampleAlgorithm.cc
  MicroBooNEExternalVertexCreationAlgorithm.cc
  MicroBooNEMasterAlgorithm.cc
  MicroBooNEDLMasterAlgorithm.cc
  MicroBooNEPreProcessingAlgorithm.cc
  LIBRARIES
  PUBLIC
  larpandoracontent::LArPandoraContent
  larpandoracontent::LArPandoraDLContent
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

