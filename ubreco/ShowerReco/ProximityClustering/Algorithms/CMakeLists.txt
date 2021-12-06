art_make(
          LIB_LIBRARIES
		   	larcorealg_Geometry
			lardata_Utilities
			larcore_Geometry_Geometry_service
                        lardataalg_DetectorInfo
			lardataobj_RecoBase
			art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except::cetlib_except

		   	art::Framework_Services_Registry
		   	art_root_io::tfile_support ROOT::Core

        )

install_fhicl()
install_headers()
install_source()

