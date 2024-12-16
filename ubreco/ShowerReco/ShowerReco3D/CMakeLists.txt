cet_build_plugin(
  ShowerTruthStudies art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::MCBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  ShrReco3D art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_ShowerReco3D_Base
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorClocksService
  lardataobj::RecoBase
  lardataobj::MCBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
  larreco::RecoAlg
  Eigen3::Eigen
)

add_subdirectory(ProtoShower)
add_subdirectory(Base)
add_subdirectory(ModularAlgo)

install_fhicl()
install_headers()
install_source()

