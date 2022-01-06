art_make(
          LIB_LIBRARIES
                        art::Framework_Core
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
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
                        nusimdata::SimulationBase
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art_root_io::TFileService_service
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
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
