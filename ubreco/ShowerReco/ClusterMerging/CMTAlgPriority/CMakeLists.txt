include(uboone_tool::CPriorityAlgoBase)

cet_build_plugin(
  CPAlgoNHits uboone_tool::CPriorityAlgoBase
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  CPAlgoQSum uboone_tool::CPriorityAlgoBase
  LIBRARIES
  PRIVATE
)

install_fhicl()
install_headers()
install_source()
