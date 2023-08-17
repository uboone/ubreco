cet_make_library(
  LOCAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE
  DataReader.cxx
  Saturation_Merger.cxx
  UBAlgo.cxx
  UBEventWaveform.cxx
  UB_rc.cxx
  UB_spe.cxx
  LIBRARIES
  PUBLIC
  ubreco::wcopreco_algo
  ROOT::Hist
)

install_headers()
#install_fhicl()
install_source()
