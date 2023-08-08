include(uboone_tool::ShowerRecoModuleBase)

cet_build_plugin(
  Angle3DFormula uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
  ROOT::Physics
)

cet_build_plugin(
  Angle3DFromVtxQweighted uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  Angle3DFromVtx uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  EmptyModule uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  FillLength uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  FilterPFPart uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  FilterShowers uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  LinearEnergy uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  StartPoint3DModule uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  StartPoint3DfromVtx uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  StartPointfromY2D uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  YPlaneStartPoint3D uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  dEdxModule uboone_tool::ShowerRecoModuleBase
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
  ROOT::Tree
)

install_fhicl()
install_headers()
install_source()
