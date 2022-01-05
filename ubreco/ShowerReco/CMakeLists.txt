
art_make(
	  LIB_LIBRARIES
			art::Framework_Core
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
			
			cetlib::cetlib cetlib_except
                   	#larpandora_LArPandoraInterface

          MODULE_LIBRARIES
		   	#larpandora_LArPandoraInterface
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_MCBase
			lardataobj_AnalysisBase
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
			art::Persistency_Common canvas::canvas
			art::Persistency_Provenance canvas::canvas
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
			
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml

        )

add_subdirectory(TwoDimTools)
add_subdirectory(ShowerReco3D)
add_subdirectory(ClusterMerging)
add_subdirectory(ProximityClustering)
add_subdirectory(ReCluster)
add_subdirectory(Pi0Ana)
add_subdirectory(TruthFilters)

install_headers()
install_source()

