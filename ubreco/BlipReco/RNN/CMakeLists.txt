cet_make_library(
  SOURCE BlipRNNAlg.cc
  LIBRARIES
    PUBLIC
    art_root_io::TFileService_service
    ubobj::DataOverlay
    lardata::Utilities
    nusimdata::SimulationBase
    TorchScatter::TorchScatter
)

install_headers()
install_fhicl()
install_source()

file(GLOB model_file *.ps)
install_fw(LIST ${model_file})
