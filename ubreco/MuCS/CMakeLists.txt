
# for classes_def.xml!!

art_make( BASENAME_ONLY
      MODULE_LIBRARIES
			larevt_CalibrationDBI_IOVData
      			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larcoreobj_SummaryData
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			larsim_MCCheater_ParticleInventoryService_service
			larsim_MCCheater_BackTrackerService_service
			lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			lardataobj_RecoBase
			lardata_RecoObjects
			larreco_RecoAlg
			ubobj_MuCS
			ubevt::SignalShapingServiceMicroBooNE_service
			lardataobj_AnalysisBase
      			nusimdata::SimulationBase
      			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
			art::Framework_Services_Optional_RandomNumberGenerator_service
			art::Persistency_Common canvas::canvas
			art::Persistency_Provenance canvas::canvas
			art::Utilities canvas::canvas
			nurandom::RandomUtils_NuRandomService_service
			messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
			cetlib::cetlib
                        cetlib_except::cetlib_except
			CLHEP::CLHEP
			ROOT::Physics
      			ROOT::Geom
      			ROOT::XMLIO
      			ROOT::Gdml
			ROOT::Core
			ROOT::EG
			SQLite::SQLite3
			CRY::CRY
			ifdh::ifdh
)

install_fhicl()
install_source()
install_headers()
