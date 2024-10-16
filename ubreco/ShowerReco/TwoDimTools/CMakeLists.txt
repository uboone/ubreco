cet_make_library(
  SOURCE
  Linearity.cxx
  Poly2D.cxx
  LIBRARIES
  PUBLIC
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  canvas::canvas
  ROOT::MathCore
)

install_fhicl()
install_headers()
install_source()
