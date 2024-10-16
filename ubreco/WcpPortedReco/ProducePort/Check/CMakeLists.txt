cet_build_plugin(
  CheckFlash art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  CheckHitsAgainst art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  CheckMetricsPlus art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  art_root_io::TFileService_service
)

cet_build_plugin(
  CheckMetrics art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  art_root_io::TFileService_service
)

cet_build_plugin(
  CheckWires art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
)

cet_build_plugin(
  WCPcheckout art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  larcore::ServiceUtil
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  WCPselection art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  larsim::EventWeight_Base
  larevt::SpaceChargeService
  larcore::ServiceUtil
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

install_fhicl()
install_headers()
install_source()
