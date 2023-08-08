cet_build_plugin(
  RecoEfficiency art::EDAnalyzer
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

install_fhicl()
install_headers()
install_source()

