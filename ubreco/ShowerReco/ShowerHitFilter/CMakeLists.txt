cet_build_plugin(
  NuGraphShowerHits art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  NuGraphShowerHitsByPfp art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  VertexSaver art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
)

install_fhicl()
install_headers()
install_source()

