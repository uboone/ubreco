art_make(
	  LIB_LIBRARIES
			art_Framework_Core
			art::Framework_Services_Registry
			art::Framework_Principal art::Framework_Core
			larcorealg_Geometry lardataobj_RecoBase
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
                        lardataalg_DetectorInfo
			art::Utilities canvas messagefacility::MF_MessageLogger
			 cetlib::cetlib

			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			messagefacility::MF_MessageLogger
			cetlib::cetlib cetlib_except

        )

install_fhicl()
install_headers()
install_source()
