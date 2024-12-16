
# Define all the base Pandora files
set(SOURCE_FILES
  DummyVertexParticleCreationAlgorithm.cc
  MicroBooNEContent.cc
  MicroBooNEExampleAlgorithm.cc
  MicroBooNEExternalVertexCreationAlgorithm.cc
  MicroBooNEMasterAlgorithm.cc
  MicroBooNEPreProcessingAlgorithm.cc
)
set(LIBRARIES larpandoracontent::LArPandoraContent)

# See if Pandora was built with DL support, and include extra files if so.
find_package(Torch QUIET EXPORT)

if (Torch_FOUND)
  message("Building MicroBooNE Pandora Wrapper with DL Support!")
  list(APPEND SOURCE_FILES MicroBooNEDLContent.cc MicroBooNEDLMasterAlgorithm.cc)
  list(APPEND LIBRARIES larpandoracontent::LArPandoraDLContent)
else()
  message("Building MicroBooNE Pandora Wrapper without DL Support!")
endif()

cet_make_library(
  SOURCE
  ${SOURCE_FILES}
  LIBRARIES
  PUBLIC
  ${LIBRARIES}
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

