include_directories(../wcopreco/data)
include_directories(../wcopreco/deconv)
include_directories(../wcopreco/algo)
include_directories(../wcopreco/app)
include_directories($ENV{EIGEN_INC})

art_make( BASENAME_ONLY
          LIBRARY_NAME  uboone_UBFlashFinder
          LIB_LIBRARIES
                   larcorealg_Geometry
                   larcore_Geometry_Geometry_service
                   ubcore::UBOpReadoutMap_service
                   lardataobj_RecoBase
                   ${UBOONECODE_LIB}
                   art::Framework_Core
                   art::Framework_Principal
                   art::Framework_Services_Registry
                   art_root_io::tfile_support ROOT::Core
                   art_root_io::TFileService_service
                   art::Persistency_Common
                   art::Persistency_Provenance
                   art::Utilities canvas::canvas
                   messagefacility::MF_MessageLogger
                   fhiclcpp::fhiclcpp
                   cetlib::cetlib
                   cetlib_except::cetlib_except
                   ROOT::Geom
                   ROOT::XMLIO
                   ROOT::Gdml
          MODULE_LIBRARIES
                   uboone_UBFlashFinder
                   larcorealg_Geometry
                   larcore_Geometry_Geometry_service
                   larevt_CalibrationDBI_IOVData
                   larevt_CalibrationDBI_Providers
                   ubreco_wcopreco_data
                   ubreco_wcopreco_deconv
                   ubreco_wcopreco_algo
                   ubreco_wcopreco_app
                   ubevt_Database
                   ubevt_Utilities
                   lardataobj_RecoBase
                   ${UBOONECODE_LIB}
                   art::Framework_Core
                   art::Framework_Principal
                   art::Framework_Services_Registry
                   art_root_io::TFileService_service
                   art::Persistency_Common
                   art::Persistency_Provenance
                   art::Utilities canvas::canvas
                   messagefacility::MF_MessageLogger
                   fhiclcpp::fhiclcpp
                   cetlib::cetlib
                   cetlib_except::cetlib_except
                   ROOT::Geom
                   ROOT::XMLIO
                   ROOT::Gdml
                   ROOT::Tree
        )

install_headers()
install_source()
install_fhicl()
add_subdirectory(job)
