include(uboone_tool::CBoolAlgoBase)

cet_build_plugin(
  CBAlgoAngleCompat uboone_tool::CBoolAlgoBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  CBAlgoFake uboone_tool::CBoolAlgoBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  CBAlgoPolar uboone_tool::CBoolAlgoBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  CBAlgoVtxAlign uboone_tool::CBoolAlgoBase
  LIBRARIES
  PRIVATE
  Eigen3::Eigen
)

install_fhicl()
install_headers()
install_source()

