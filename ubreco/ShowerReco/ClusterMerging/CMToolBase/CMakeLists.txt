cet_make_library(
  SOURCE
  CMManagerBase.cxx
  CMTException.cxx
  CMatchBookKeeper.cxx
  CMatchManager.cxx
  CMergeBookKeeper.cxx
  CMergeManager.cxx
  Cluster.cxx
  ClusterMaker.cxx
  LIBRARIES
  PUBLIC
  ubreco::ShowerReco_TwoDimTools
  lardata::Utilities
  larcore::Geometry_Geometry_service
  art::Utilities
  art_root_io::TFileService_service
  ROOT::Core
)

cet_make_library(
  LIBRARY_NAME CBoolAlgoBase INTERFACE
  SOURCE CBoolAlgoBase.h
  LIBRARIES
  INTERFACE
  art::Utilities
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_write_plugin_builder(
  uboone_tool::CBoolAlgoBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubreco::CBoolAlgoBase
)

cet_make_library(
  LIBRARY_NAME CFloatAlgoBase INTERFACE
  SOURCE CFloatAlgoBase.h
  LIBRARIES
  INTERFACE
  art::Utilities
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_write_plugin_builder(
  uboone_tool::CFloatAlgoBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubreco::CFloatAlgoBase
)

cet_make_library(
  LIBRARY_NAME CPriorityAlgoBase INTERFACE
  SOURCE CPriorityAlgoBase.h
  LIBRARIES
  INTERFACE
  art::Utilities
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_write_plugin_builder(
  uboone_tool::CPriorityAlgoBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubreco::CPriorityAlgoBase
)

install_fhicl()
install_headers()
install_source()
