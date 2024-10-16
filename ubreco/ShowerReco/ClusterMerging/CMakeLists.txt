cet_build_plugin(
  ClusterMatcher art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_ClusterMerging_CMToolBase
)

cet_build_plugin(
  ClusterMerger art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_ClusterMerging_CMToolApp
  lardata::DetectorClocksService
)

art_dictionary(
  DICTIONARY_LIBRARIES
  lardataobj::RecoBase
 )

add_subdirectory(CMToolBase)
add_subdirectory(CMTAlgMerge)
add_subdirectory(CMTAlgMatch)
add_subdirectory(CMTAlgPriority)
add_subdirectory(CMToolApp)

install_fhicl()
install_headers()
install_source()
