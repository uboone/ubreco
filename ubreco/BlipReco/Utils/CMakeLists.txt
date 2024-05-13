cet_make_library(
  SOURCE
  BlipUtils.cc
  LIBRARIES
  PUBLIC
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  lardata::Utilities
  lardata::DetectorPropertiesService
  nusimdata::SimulationBase
  canvas::canvas
)

art_dictionary(
  DICTIONARY_LIBRARIES
  lardataobj::RecoBase
  nusimdata::SimulationBase
)

install_headers()
install_fhicl()
install_source()
