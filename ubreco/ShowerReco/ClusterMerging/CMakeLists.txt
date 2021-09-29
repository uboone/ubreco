art_make(
          LIB_LIBRARIES
                        ${ART_FRAMEWORK_CORE}
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        cetlib::cetlib
                        larcorealg_Geometry
                        lardata_Utilities
                        lardataobj_RecoBase
                        lardataobj_MCBase
                        larcore_Geometry_Geometry_service
          MODULE_LIBRARIES
                        ubreco_ShowerReco_ClusterMerging_CMToolBase
                        ubreco_ShowerReco_ClusterMerging_CMToolApp
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardataobj_RecoBase
                        lardataobj_MCBase
                        lardataobj_AnalysisBase
                        lardata_RecoObjects
                        larreco_RecoAlg
                        lardata_Utilities
                        larevt_Filters
                        nusimdata_SimulationBase
                        ${LARRECO_LIB}
                        ${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Persistency_Common
                        art_Persistency_Provenance
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        cetlib::cetlib
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
           DICT_LIBRARIES nusimdata::SimulationBase
        )

add_subdirectory(CMToolBase)
add_subdirectory(CMTAlgMerge)
add_subdirectory(CMTAlgMatch)
add_subdirectory(CMTAlgPriority)
add_subdirectory(CMToolApp)

install_fhicl()
install_headers()
install_source()
