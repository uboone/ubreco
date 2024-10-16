cet_build_plugin(
  OpNoiseCreateMask art::EDProducer
  LIBRARIES
  PRIVATE
  larreco::HitFinder
  larevt::SpaceChargeService
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

install_headers()
install_fhicl()
install_source()


