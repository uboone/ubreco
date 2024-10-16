cet_make_library(
  LOCAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE
  FlashFiltering.cxx
  Flashes_beam.cxx
  Flashes_cosmic.cxx
  HitFinder_beam.cxx
  HitFinder_cosmic.cxx
  LIBRARIES
  PUBLIC
  ubreco::wcopreco_deconv
  ROOT::Hist
)

install_headers()
#install_fhicl()
install_source()
