cet_build_plugin(
  BackTracking art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  CosmicBackgrounds art::EDFilter
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  PerfectClustering art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  Pi0AnalyzerDATA art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_Pi0Ana_Selection
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  Pi0Analyzer art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_Pi0Ana_Selection
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::MCBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  Pi0EventFilter art::EDFilter
  LIBRARIES
  PRIVATE
  ROOT::Tree
)

cet_build_plugin(
  Pi0Filter art::EDFilter
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_Pi0Ana_Selection
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  Pi0PhysicsDATA art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_Pi0Ana_Selection
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  Pi0Physics art::EDAnalyzer
  LIBRARIES
  PRIVATE
  ubreco::ShowerReco_Pi0Ana_Selection
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

add_subdirectory(Selection)

install_fhicl()
install_headers()
install_source()

