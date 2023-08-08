include_directories(../data)
include_directories(../deconv)

cet_make_library(
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
