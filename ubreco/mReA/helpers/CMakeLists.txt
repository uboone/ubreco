
art_make(
  # Needed for standalone c++ helper files
  LIBRARY_NAME
    ubreco_mReA_helpers
  LIB_LIBRARIES
    larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation
    nug4_ParticleNavigation
    lardataobj_Simulation
    lardata_Utilities
    larevt_Filters
    lardataobj_RawData
    lardataobj_RecoBase 
    lardata_RecoObjects 
    larreco_RecoAlg
    larpandora_LArPandoraInterface
    SignalShapingServiceMicroBooNE_service
    nusimdata_SimulationBase
    ubobj_Trigger
    ubobj_Optical
    ${LARRECO_LIB}
    ${LARDATA_LIB}
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    art_Persistency_Common canvas
    art_Persistency_Provenance canvas
    art_Utilities canvas
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib 
    cetlib_except
    ${ROOT_GEOM}
    ${ROOT_XMLIO}
    ${ROOT_GDML}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
  )

install_headers()
#install_fhicl()
install_source()
