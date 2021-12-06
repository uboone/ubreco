art_make(
          LIB_LIBRARIES
                        ubreco_ShowerReco_ShowerReco3D_Base
                        #ubreco_ShowerReco_ShowerReco3D_ProtoShower
                        art::Framework_Core
                        art::Framework_Services_Registry
                        larcorealg_Geometry
                        uboone_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larreco_Calorimetry
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except::cetlib_except

                        art::Framework_Services_Registry
                        art_root_io::tfile_support ROOT::Core
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
                          nurandom::RandomUtils_NuRandomService_service
                          art::Framework_Core
                          art::Framework_Principal
                          art::Framework_Services_Registry
                          art::Persistency_Common
                          art::Persistency_Provenance
                          art_root_io::tfile_support ROOT::Core
                          art::Framework_Services_Optional_RandomNumberGenerator_service
                          art_root_io::TFileService_service
                          ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                          art::Utilities
                          canvas::canvas
                          messagefacility::MF_MessageLogger
                          
                          fhiclcpp::fhiclcpp
                          cetlib::cetlib cetlib_except::cetlib_except
                          ROOT::Geom
                          Boost::filesystem
                          CLHEP::CLHEP
                          ROOT::Tree
        )

install_fhicl()
install_headers()
install_source()
