
art_make(	
	  LIB_LIBRARIES
			art::Framework_Core
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib
                        cetlib_except::cetlib_except
                   	#larpandora_LArPandoraInterface	
			${ROOT_BASIC_LIB_LIST}
			larcorealg_Geometry
			lardata_Utilities
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_Simulation	
			lardataobj_MCBase
			larcore_Geometry_Geometry_service
          MODULE_LIBRARIES 
		   	#larpandora_LArPandoraInterface	
			larevt_CalibrationDBI_IOVData
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_Simulation	
			lardataobj_AnalysisBase	
			lardata_RecoObjects
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata::SimulationBase
			nurandom::RandomUtils_NuRandomService_service
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
			ROOT::Physics
      			ROOT::Geom
      			ROOT::XMLIO
      			ROOT::Gdml
			${ROOT_BASIC_LIB_LIST}
			ROOT::Core
			ROOT::EG
			CRY::CRY
        )


install_fhicl()
install_headers()
install_source()
