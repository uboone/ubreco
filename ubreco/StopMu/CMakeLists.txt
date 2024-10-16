cet_build_plugin(
  StopMuMCFilter art::EDFilter
  LIBRARIES
  PRIVATE
  nusimdata::SimulationBase
)

cet_build_plugin(
  StopMuMCGen art::EDProducer
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  nutools::EventGeneratorBase_CRY
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  nugen::EventGeneratorBase_GENIE
  CRY::CRY
  ROOT::EG
)

cet_build_plugin(
  StopMu art::EDProducer
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  ThroughMu art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_fhicl()
install_headers()
install_source()
