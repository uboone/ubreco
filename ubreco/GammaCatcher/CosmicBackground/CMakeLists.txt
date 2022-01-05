
art_make(	
	  LIB_LIBRARIES
			art::Framework_Core
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			
			cetlib::cetlib cetlib_except
                   	#larpandora_LArPandoraInterface	
                        
			${ROOT_BASIC_LIB_LIST}
			larcorealg_Geometry
			lardata_Utilities
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_Simulation	
			larcore_Geometry_Geometry_service
          MODULE_LIBRARIES 
	  		ubreco_ShowerReco_Pi0Ana_Selection
		   	#larpandora_LArPandoraInterface	
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_Simulation	
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
			cetlib::cetlib cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
			${ROOT_BASIC_LIB_LIST}
        )
install_fhicl()
install_headers()
install_source()

