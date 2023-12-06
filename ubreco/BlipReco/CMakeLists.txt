cet_build_plugin(
  BlipAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::BlipReco_Alg
  larreco::Calorimetry
  larsim::MCCheater_BackTrackerService_service
  larevt::SpaceChargeService
  lardata::DetectorClocksService  
  larcore::ServiceUtil
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  BlipRecoProducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::BlipReco_Alg
  larreco::Calorimetry
  larsim::MCCheater_BackTrackerService_service
  larevt::CalibrationDBI_IOVData
  lardata::Utilities
  lardata::LArPropertiesService
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  ParticleDump art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::BlipReco_Alg
  larsim::MCCheater_BackTrackerService_service
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  TrackMasker art::EDProducer
  LIBRARIES
  PRIVATE
  larsim::MCCheater_BackTrackerService_service
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)



add_subdirectory(Utils)
add_subdirectory(Alg)
add_subdirectory(job)
install_fhicl()
install_source()
install_headers()

file(GLOB txt_file *.txt)
install_fw(LIST ${txt_file})
