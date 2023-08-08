cet_build_plugin(
  ACPTtrigMCFilter art::EDFilter
  LIBRARIES
  PRIVATE
  nusimdata::SimulationBase
)

cet_build_plugin(
  ACPTtrigMCGen art::EDProducer
  LIBRARIES
  PRIVATE
  larcore::Geometry_Geometry_service
  art_root_io::TFileService_service
  nutools::EventGeneratorBase_CRY
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  nugen::EventGeneratorBase_GENIE
  CRY::CRY
  ROOT::EG
)

cet_build_plugin(
  ACPTtrig art::EDProducer
  LIBRARIES
  PRIVATE
  ubevt::Utilities
  larevt::SpaceChargeService
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  CosmicTaggingAnodeCathodePiercing art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
)

cet_build_plugin(
  T0RecoAnodeCathodePiercingAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  T0RecoAnodeCathodePiercingValidationAna art::EDAnalyzer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

cet_build_plugin(
  T0RecoAnodeCathodePiercing art::EDProducer
  LIBRARIES
  PRIVATE
  lardata::Utilities
  larcore::Geometry_Geometry_service
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  art_root_io::TFileService_service
  ROOT::Tree
)

install_headers()
install_fhicl()
install_source()

