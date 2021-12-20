art_make(
          LIB_LIBRARIES
			ubreco_ShowerReco_ClusterMerging_CMToolBase
			art_Framework_Core
			art::Framework_Services_Registry
			art::Framework_Principal
			lardataobj_RecoBase
			art::Framework_Core
			larcorealg_Geometry
			lardata_Utilities # added now
		   	#uboone_Geometry
			larcore_Geometry_Geometry_service
			art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        
                        cetlib::cetlib cetlib_except::cetlib_except

		   	art::Framework_Services_Registry
		   	art_root_io::tfile_support ROOT::Core

        )

install_fhicl()
install_headers()
install_source()

