cet_make_library(
  SOURCE
  BlipRecoAlg.cc
  LIBRARIES
  PUBLIC
  ubobj::DataOverlay
  ubevt::Database
  ubreco::BlipReco_Utils
  larevt::CalibrationDBI_IOVData
  larreco::Calorimetry
  lardata::DetectorClocksService
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

install_headers()
install_fhicl()
install_source()
