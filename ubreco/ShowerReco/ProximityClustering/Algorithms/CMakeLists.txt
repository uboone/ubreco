cet_make_library(
  SOURCE
  ProximityClusterer.cxx
  LIBRARIES
  PUBLIC
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art::Framework_Principal
  messagefacility::MF_MessageLogger
)

install_fhicl()
install_headers()
install_source()
