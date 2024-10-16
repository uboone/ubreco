cet_build_plugin(
  NuSelectionMetrics art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  ROOT::Tree
)

cet_build_plugin(
  NuSelectionPlus art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
)

cet_build_plugin(
  PortedFlash art::EDProducer
  LIBRARIES
  PRIVATE
  ubcore::Geometry
  lardata::Utilities
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
)

cet_build_plugin(
  PortedHits art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
)

cet_build_plugin(
  PortedSpacePoints art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
)

cet_build_plugin(
  SelectWires art::EDProducer
  LIBRARIES
  PRIVATE
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
)

cet_build_plugin(
  WCPHybrid art::EDProducer
  LIBRARIES
  PRIVATE
  larpandora::LArPandoraInterface
  lardata::Utilities
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
)

cet_build_plugin(
  WireCellPF art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::WcpPort
  nusimdata::SimulationBase
)

add_subdirectory(Check)
add_subdirectory(job)

install_fhicl()
install_headers()
install_source()
