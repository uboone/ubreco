include_directories(../data)

cet_make_library(
  SOURCE
  Deconvolver.cxx
  kernel_fourier.cxx
  kernel_fourier_container.cxx
  LIBRARIES
  PUBLIC
  ubreco::wcopreco_data
  ROOT::Hist
)

install_headers()
#install_fhicl()
install_source()
