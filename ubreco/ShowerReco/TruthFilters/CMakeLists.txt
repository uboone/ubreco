cet_build_plugin(
  NeutrinoFilter art::EDFilter
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_headers()
install_source()
