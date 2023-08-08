cet_build_plugin(
  MuCSDT art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::MuCS
  lardata::LArPropertiesService
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  ifdh::ifdh
)

cet_build_plugin(
  MuCSExtrapolate art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::MuCS
  lardataobj::RawData
  ROOT::Graf
)

cet_build_plugin(
  MuCSMCGen art::EDProducer
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  nutools::EventGeneratorBase_CRY
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  CRY::CRY
  ROOT::EG
)

cet_build_plugin(
  MuCSMerger art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::MuCS
  lardataobj::RawData
  ifdh::ifdh
)

cet_build_plugin(
  MuCSReader art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubobj::MuCS
  larreco::RecoAlg
  larsim::Simulation
  lardata::LArPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  art_root_io::TFileService_service
  nusimdata::SimulationBase
)

cet_build_plugin(
  MuCSRecoToMC art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::MuCS
  art_root_io::TFileService_service
  nutools::EventGeneratorBase_CRY
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  ROOT::EG
)

cet_build_plugin(
  MuCSReco art::EDProducer
  LIBRARIES
  PRIVATE
  ubobj::MuCS
  lardataobj::RawData
  SQLite::SQLite3
)

cet_build_plugin(
  MuCSTrackTagger art::EDProducer
  LIBRARIES
  PRIVATE
  larevt::PmtGainService
  lardata::LArPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_fhicl()
install_source()
install_headers()
