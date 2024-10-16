cet_make_library(LIBRARY_NAME ProtoShowerAlgBase INTERFACE
  SOURCE ProtoShowerAlgBase.h
  LIBRARIES INTERFACE
  lardata::Utilities
  lardataobj::RecoBase
  canvas::canvas
)

cet_write_plugin_builder(
  uboone_tool::ProtoShowerAlgBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubreco::ProtoShowerAlgBase
)

include(uboone_tool::ProtoShowerAlgBase)

cet_build_plugin(
  ProtoShowerCMTool uboone_tool::ProtoShowerAlgBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  ProtoShowerPandora uboone_tool::ProtoShowerAlgBase
  LIBRARIES
  PRIVATE
)

install_fhicl()
install_headers()
install_source()

