cet_build_plugin(
  AnalyzeFlashMatchChi2 art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  FilteredHitsProducer art::EDProducer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  FilteredHitsProducerByPfp art::EDProducer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  FilteredHitsTruthProducer art::EDProducer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  NuSliceReBuilderProducer art::EDProducer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  nusimdata::SimulationBase
)

cet_build_plugin(
  PandoraEventAnalyser art::EDAnalyzer
  LIBRARIES
  PRIVATE
  art_root_io::TFileService_service
)

cet_build_plugin(
  StoreFlashMatchChi2 art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubevt::Database
  ubobj::DataOverlay
  larpandora::LArPandoraInterface
  larevt::CalibrationDBI_IOVData
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

include(larpandora::SliceIdTool)

cet_build_plugin(
  FlashNeutrinoId larpandora::SliceIdTool
  LIBRARIES
  PRIVATE
  ubreco::LLSelectionTool_OpT0Finder_Algorithms
  ubreco::PandoraEventBuildingFlashID_HitCosmicTag_Algorithms
  ubevt::Database
  ubobj::DataOverlay
  larpandora::LArPandoraEventBuilding
  larreco::RecoAlg
  larevt::CalibrationDBI_IOVData
  art_root_io::TFileService_service
)

add_subdirectory(job)
add_subdirectory(HitCosmicTag)

install_headers()
install_fhicl()
install_source()
