cet_make_library(
  SOURCE
  FlashAlgoBase.cxx
  FlashAlgoFactory.cxx
  FlashFinderFMWKInterface.cxx
  FlashFinderManager.cxx
  PECalib.cxx
  SimpleFlashAlgo.cxx
  LIBRARIES
  PUBLIC
  ubcore::Geometry
  larcore::Geometry_Geometry_service
  art::Framework_Core
  art_root_io::TFileService_service
)

cet_build_plugin(
  OpHitCorrectionProducer art::EDProducer
  LIBRARIES
  PRIVATE
  ubevt::Utilities
  larevt::CalibrationDBI_IOVData
  larevt::PmtGainService
  lardataobj::RecoBase
)

cet_build_plugin(
  UBFlashFinder art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::UBFlashFinder
  ubevt::Utilities
  larevt::PmtGainService
  lardata::Utilities
  lardataobj::RecoBase
)

cet_build_plugin(
  UBWCFlashFinder art::EDProducer
  LIBRARIES
  PRIVATE
  ubreco::UBFlashFinder
  ubreco::wcopreco_app
  larevt::PmtGainService
  lardata::DetectorClocksService
  larcore::Geometry_Geometry_service
  lardataobj::RawData
  lardataobj::RecoBase
  art_root_io::TFileService_service
  Eigen3::Eigen
  ROOT::Tree
)

install_headers()
install_source()
install_fhicl()
add_subdirectory(job)
