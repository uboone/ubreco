cet_make_library(
  SOURCE
  CMergeHelper.cxx
  LIBRARIES
  PUBLIC
  ubreco::ShowerReco_ClusterMerging_CMToolBase
)

install_fhicl()
install_headers()
install_source()
