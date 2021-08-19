art_make(
          LIB_LIBRARIES
			ubreco_ShowerReco_TwoDimTools
			art_Framework_Core
			art_Framework_Services_Registry
			art_Framework_Principal
			lardataobj_RecoBase
			${ART_FRAMEWORK_CORE}
			larcorealg_Geometry
			lardata_Utilities # added now
		   	uboone_Geometry
			larcore_Geometry_Geometry_service
			art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB} cetlib_except
			${Boost_FILESYSTEM_LIBRARY}

		   	${ART_FRAMEWORK_SERVICES_REGISTRY}
		   	${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
          TOOL_LIBRARIES
	  		  ubreco_ShowerReco_ClusterMerging_CMToolBase
			  ubreco_ShowerReco_TwoDimTools
			  ubevt::SignalShapingServiceMicroBooNE_service
                          lardata_Utilities_LArFFT_service
                          larcorealg_Geometry
			  larevt_CalibrationDBI_IOVData
			  larevt_CalibrationDBI_Providers
			  lardataobj_RecoBase
			  larcore_Geometry_Geometry_service
			  lardata_Utilities
			  nurandom_RandomUtils_NuRandomService_service
  			  ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			  art_Persistency_Common
			  art_Persistency_Provenance
			  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
			  ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
			  art_Utilities
			  canvas
			  ${MF_MESSAGELOGGER}
			  ${MF_UTILITIES}
			  ${FHICLCPP}
			  ${CETLIB} cetlib_except
			  ${Boost_FILESYSTEM_LIBRARY}
                          ${CLHEP}


        )

install_fhicl()
install_headers()
install_source()

