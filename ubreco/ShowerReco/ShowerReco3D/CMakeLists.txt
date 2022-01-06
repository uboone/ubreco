
art_make(
	  LIB_LIBRARIES
			art::Framework_Core
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib
                        cetlib_except::cetlib_except
                   	#larpandora_LArPandoraInterface
			larcorealg_Geometry
			lardata_Utilities
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_Simulation
			larcore_Geometry_Geometry_service
          MODULE_LIBRARIES
	  		ubreco_ShowerReco_ShowerReco3D_Base
	  		#ubreco_ShowerReco_ShowerReco3D_ProtoShower
		   	#larpandora_LArPandoraInterface
                        lardataalg_DetectorInfo
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_AnalysisBase
			lardataobj_Simulation
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
			art::Utilities
                        canvas::canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml

        )

add_subdirectory(ProtoShower)
add_subdirectory(Base)
add_subdirectory(ModularAlgo)

install_fhicl()
install_headers()
install_source()

