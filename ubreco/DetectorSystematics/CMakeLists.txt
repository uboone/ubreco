cet_build_plugin(
  DummyRawDigitCreator art::EDProducer
  LIBRARIES
  PRIVATE
  lardataobj::RawData
)

cet_build_plugin(
  HitVariationsAnalyzer art::EDAnalyzer
  LIBRARIES
  PRIVATE
  larevt::SpaceChargeService
  larcore::ServiceUtil
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  WireModifier art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()

