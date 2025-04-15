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

install_headers()
install_fhicl()
install_source()
