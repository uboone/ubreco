
art_make(
  MODULE_LIBRARIES 
    larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation
    nug4::ParticleNavigation
    lardataobj_Simulation
    lardata_Utilities
    larevt_Filters
    lardataobj_RawData
    lardataobj_RecoBase 
    lardata_RecoObjects 
    larreco_RecoAlg
    larpandora_LArPandoraInterface
    larpandora_LArPandoraEventBuilding
    SignalShapingServiceMicroBooNE_service
    nusimdata::SimulationBase
    ubobj_Trigger
    ubobj_Optical
    ubreco_mReA_helpers
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib 
    cetlib_except::cetlib_except
    ROOT::Geom
    ROOT::XMLIO
    ROOT::Gdml
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.21.01) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

add_subdirectory(helpers)
add_subdirectory(job)

install_source()
