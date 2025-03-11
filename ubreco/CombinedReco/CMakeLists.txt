cet_build_plugin(
  EmptyFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_fhicl()
install_headers()
install_source()
