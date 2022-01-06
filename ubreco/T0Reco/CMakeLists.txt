
art_make(
	  LIB_LIBRARIES
			art::Framework_Core
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
	                cetlib::cetlib
                        cetlib_except::cetlib_except

          MODULE_LIBRARIES
                        ubevt_Utilities
                        lardataalg_DetectorInfo
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_AnalysisBase
			lardata_RecoObjects
			larreco_RecoAlg
			lardata_Utilities
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			larcoreobj_SummaryData
			larevt_Filters
                        nurandom::RandomUtils_NuRandomService_service
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
	                cetlib::cetlib
                        cetlib_except::cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ROOT::EG
                        CRY::CRY

        )

install_headers()
install_fhicl()
install_source()

