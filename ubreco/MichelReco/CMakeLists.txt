cet_build_plugin(
  MichelRecoDriver art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::MichelReco_Algo
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
)

add_subdirectory(Fmwk)
add_subdirectory(Algo)

install_headers()
install_fhicl()
install_source()

