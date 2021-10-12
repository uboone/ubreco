art_make(
          LIB_LIBRARIES
                        ubreco_ShowerReco_ShowerReco3D_Base
                        #ubreco_ShowerReco_ShowerReco3D_ProtoShower
                        ${ART_FRAMEWORK_CORE}
                        art_Framework_Services_Registry
                        larcorealg_Geometry
                        uboone_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larreco_Calorimetry
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB} cetlib_except

                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
          TOOL_LIBRARIES
                          ubreco_ShowerReco_ClusterMerging_CMToolBase
                          ubreco_ShowerReco_TwoDimTools
                          ubevt::SignalShapingServiceMicroBooNE_service
                          lardata_Utilities_LArFFT_service
                          larcorealg_Geometry
                          larreco_Calorimetry
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
                          ${ROOT_GEOM}
                          ${Boost_FILESYSTEM_LIBRARY}
                          CLHEP::CLHEP
                          ROOT::Tree
        )

install_fhicl()
install_headers()
install_source()
