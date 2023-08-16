cet_make_library(
  SOURCE
  ShowerAnaBase.cxx
  ShowerRecoAlgBase.cxx
  ShowerRecoException.cxx
  ShowerRecoUtils.cxx_dep
  ShrRecoManager.cxx
  LIBRARIES
  PUBLIC
  lardata::Utilities
  larcore::Geometry_Geometry_service
  ROOT::Physics
)

cet_make_library(
  LIBRARY_NAME ShowerRecoModuleBase INTERFACE
  SOURCE ShowerRecoModuleBase.h
  LIBRARIES INTERFACE
  lardata::Utilities
  larcore::Geometry_Geometry_service
)

cet_write_plugin_builder(
  uboone_tool::ShowerRecoModuleBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubreco::ShowerRecoModuleBase
)

install_fhicl()
install_headers()
install_source()

