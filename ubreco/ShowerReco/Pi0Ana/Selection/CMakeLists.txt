cet_make_library(
  SOURCE
  SelectionAlg.cxx
  TruncMean.cxx
  LIBRARIES
  PUBLIC
  lardata::Utilities
  lardataobj::RecoBase
  art::Framework_Principal
)

install_headers()
install_source()

