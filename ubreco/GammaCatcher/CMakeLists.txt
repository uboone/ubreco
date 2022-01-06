
art_make(	
	  LIB_LIBRARIES
			art::Framework_Core
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib
                        cetlib_except::cetlib_except
          MODULE_LIBRARIES 
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_AnalysisBase	
			lardata_RecoObjects
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata::SimulationBase
                        art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art::Framework_Services_Optional_RandomNumberGenerator_service
			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
        )

add_subdirectory(Analysis)
add_subdirectory(HitFinding)
add_subdirectory(CosmicBackground)
add_subdirectory(job)

install_headers()
install_fhicl()
install_source()

