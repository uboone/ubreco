include(uboone_tool::CFloatAlgoBase)

cet_build_plugin(
  CFAlgoDevel uboone_tool::CFloatAlgoBase
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
)

cet_build_plugin(
  CFAlgoIoU uboone_tool::CFloatAlgoBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  CFAlgoTimeOverlap uboone_tool::CFloatAlgoBase
  LIBRARIES
  PRIVATE
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
)

install_fhicl()
install_headers()
install_source()
