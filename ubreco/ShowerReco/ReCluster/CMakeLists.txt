cet_build_plugin(
  PhotonMerge art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_ClusterMerging_CMToolBase
  lardata::Utilities
  lardataobj::RecoBase
)

install_fhicl()
install_headers()
install_source()
