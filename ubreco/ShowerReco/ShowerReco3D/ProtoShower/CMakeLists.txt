
cet_find_library( PANDORASDK NAMES PandoraSDK PATHS ENV PANDORA_LIB NO_DEFAULT_PATH)

art_make(
          LIB_LIBRARIES
			art_Framework_Core
			art::Framework_Services_Registry
			art::Framework_Principal
			art::Persistency_Provenance
			art::Framework_Core
			larcorealg_Geometry
			lardataobj_RecoBase
		   	larcorealg_Geometry
			#larpandora_LArPandoraObjects
                        #LArPandoraContent
                        #larpandoraobj
            		PANDORASDK
			larcore_Geometry_Geometry_service
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
			#larpandora_LArPandoraObjects
                        #LArPandoraContent
                        #larpandoraobj
            		PANDORASDK
			  nurandom::RandomUtils_NuRandomService_service
  			  art::Framework_Core
			  art::Framework_Principal
			  art::Framework_Services_Registry
			  art::Persistency_Common
			  art::Persistency_Provenance
			  art_root_io::tfile_support ROOT::Core
                          art::Framework_Services_Optional_RandomNumberGenerator_service
			  art_root_io::TFileService_service
			  art::Utilities
			  canvas::canvas
			  messagefacility::MF_MessageLogger
			  fhiclcpp::fhiclcpp
			  cetlib::cetlib cetlib_except::cetlib_except
			  ROOT::Geom
			  Boost::filesystem
                          CLHEP::CLHEP
	MODULE_LIBRARIES larpandora_LArPandoraInterface
        )

install_fhicl()
install_headers()
install_source()

