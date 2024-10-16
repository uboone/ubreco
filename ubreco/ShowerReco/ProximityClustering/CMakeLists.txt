cet_build_plugin(
  CCVertex art::EDFilter
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  CosmicFilter art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  ProximityClustering art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_ProximityClustering_Algorithms
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  TruthVertex art::EDFilter
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
)

add_subdirectory(Algorithms)

install_fhicl()
install_headers()
install_source()
