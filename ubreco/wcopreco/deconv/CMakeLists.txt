cet_make_library(
  LOCAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
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
