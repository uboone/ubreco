art_make(
          LIB_LIBRARIES
			ubreco_ShowerReco_TwoDimTools
			art::Framework_Core
			art::Framework_Services_Registry
			art::Framework_Principal
			lardataobj_RecoBase
			art::Framework_Core
			larcorealg_Geometry
			lardata_Utilities # added now
		   	uboone_Geometry
			larcore_Geometry_Geometry_service
			art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except::cetlib_except
			Boost::filesystem

		   	art::Framework_Services_Registry
		   	art_root_io::tfile_support ROOT::Core
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
			  Boost::filesystem
                          CLHEP::CLHEP


        )

install_fhicl()
install_headers()
install_source()

