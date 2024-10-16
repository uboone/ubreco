cet_build_plugin(
  CheckNeutrino art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  ClusterTrackDistance art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  Gamma3D art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  SelectNeutrino art::EDFilter
  LIBRARIES
  PRIVATE
)

cet_build_plugin(
  TruthStudies art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  nusimdata::SimulationBase
  ROOT::Tree
)

cet_build_plugin(
  gammacorrelation art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardata::DetectorPropertiesService
  lardataobj::RecoBase
  art_root_io::TFileService_service
  ROOT::Tree
)

add_subdirectory(Analysis)
add_subdirectory(HitFinding)
add_subdirectory(CosmicBackground)
add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

