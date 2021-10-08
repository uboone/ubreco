include_directories( $ENV{EIGEN_INC} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
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
                   ${LARDATA_LIB}
                   ${ART_FRAMEWORK_CORE}
                   ${ART_FRAMEWORK_PRINCIPAL}

                   ${ART_FRAMEWORK_SERVICES_REGISTRY}
                   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                   art_Persistency_Common
                   art_Persistency_Provenance
                   art_Utilities canvas
                   ${MF_MESSAGELOGGER}
                   ${FHICLCPP}
                   cetlib cetlib_except
                   ${ROOT_GEOM}
                   ${ROOT_XMLIO}
                   ${ROOT_GDML}
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
		   #${UBOONECODE_LIB}
		   ${LARDATA_LIB}
		   ${ART_FRAMEWORK_CORE}
		   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_BASIC}
		   ${ART_FRAMEWORK_SERVICES_REGISTRY}
		   ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		   ${ART_ROOT_IO_TFILESERVICE_SERVICE}
		   art_Persistency_Common canvas
		   art_Persistency_Provenance canvas
		   art_Utilities canvas
		   ${MF_MESSAGELOGGER}
		   ${FHICLCPP}
		   cetlib cetlib_except
		   ${ROOT_GEOM}
		   ${ROOT_XMLIO}
		   ${ROOT_GDML}
		   ${ROOT_BASIC_LIB_LIST}

        )							   

install_headers()
install_source()
install_fhicl()
add_subdirectory(job)
