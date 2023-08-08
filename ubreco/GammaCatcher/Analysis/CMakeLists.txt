cet_build_plugin(
  ExampleAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  LoadCCVertex art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
)

install_fhicl()
install_headers()
install_source()

