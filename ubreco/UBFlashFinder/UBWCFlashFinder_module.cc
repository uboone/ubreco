////////////////////////////////////////////////////////////////////////
// Class:       UBWCFlashFinder
// Module Type: producer
// File:        UBWCFlashFinder_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft 
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//DB
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Interface/PmtGainProvider.h"
#include "larevt/CalibrationDBI/IOVData/CalibrationExtraInfo.h"
//RawDigits
#include "lardataobj/RawData/raw.h"
//#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RawData/OpDetWaveform.h"
//OpFlash
#include "lardataobj/RecoBase/OpFlash.h"
#include <memory>
#include <string>
#include "FlashFinderFMWKInterface.h" //pmtana
//WCOpReco
#include "Config_Params.h"
#include "OpWaveform.h"
#include "OpWaveformCollection.h"
#include "UBEventWaveform.h"
#include "UB_rc.h"
#include "UB_spe.h"
#include "UBAlgo.h"
#include "kernel_fourier_container.h"

class UBWCFlashFinder;

class UBWCFlashFinder : public art::EDProducer {
public:
  explicit UBWCFlashFinder(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  UBWCFlashFinder(UBWCFlashFinder const &) = delete;
  UBWCFlashFinder(UBWCFlashFinder &&) = delete;
  UBWCFlashFinder & operator = (UBWCFlashFinder const &) = delete;
  UBWCFlashFinder & operator = (UBWCFlashFinder &&) = delete;

  // Required functions.
  void produce(art::Event & evt) override;


private:

  // Declare member data here.
  typedef enum {kBeamHighGain=0,kBeamLowGain,kCosmicHighGain,kCosmicLowGain} OpDiscrTypes;
  std::string _OpDataProducer;
  std::string _OpSatDataProducer;
  std::vector<std::string> _OpDataTypes;
  std::string _TriggerProducer;
  std::vector<float> pmt_gain;
  std::vector<float> pmt_gainerr;
  bool _usePmtGainDB;
  bool _remap_ch;
  bool _useExtSat;
  float _OpDetFreq;
  std::vector<std::string> _flashProducts;
  std::vector<std::string> _saturationProducts;
  ::wcopreco::Config_Params flash_pset;
  ::wcopreco::UBAlgo flash_algo;
  bool _saveAnaTree;
  //for ana output
  TTree* _outtree;
  Int_t   fEventID;
  Int_t   fFlashID;
  Float_t fFlashTime; 
  Float_t fFlashWidth;
  Float_t fAbsTime;
  bool    fInBeamFrame;
  int     fOnBeamTime;
  Float_t fTotalPE;
  Int_t   fFlashFrame;
  Int_t   fFlashType; //0=undefined 1=cosmic 2=beam
  std::vector< double > fPEPerCh; //per flash    
  Float_t fYCenter;
  Float_t fYWidth;
  Float_t fZCenter;
  Float_t fZWidth;


  void reco_default(art::Event &evt, double &triggerTime);
  void reco_external_sat(art::Event &evt, double &triggerTime);
  ::wcopreco::UBEventWaveform fill_evt_wf(double triggerTime,
					  std::vector<raw::OpDetWaveform> bgh, 
  					  std::vector<raw::OpDetWaveform> blg,
  					  std::vector<raw::OpDetWaveform> chg,
  					  std::vector<raw::OpDetWaveform> clg,
  					  std::vector<float> pmt_gain, std::vector<float> pmt_gainerr );
  void fill_wfmcollection(double triggerTime,
			  std::vector<raw::OpDetWaveform> opwfms,
			  int type,
			  ::wcopreco::OpWaveformCollection &wfm_collection);
  std::vector<::wcopreco::kernel_fourier_container> fill_kernel_container(std::vector<float> pmt_gain);
  void GetFlashLocation(std::vector<double>, double&, double&, double&, double&);
  void fill_ana_tree(recob::OpFlash flash,int idx, int type);
};


UBWCFlashFinder::UBWCFlashFinder(fhicl::ParameterSet const & p)
: EDProducer(p)
// Initialize member data here.
{
  _OpDataProducer   = p.get<std::string>("OpDataProducer", "pmtreadout" );   // Waveform Module name, to get waveforms
  _OpSatDataProducer= p.get<std::string>("OpSatDataProducer", "saturation" );   // Saturation corrected waveforms
  _OpDataTypes      = p.get<std::vector<std::string> >("OpDataTypes");
  _flashProducts    = p.get<std::vector<std::string> >("FlashProducts");
  _saturationProducts = p.get<std::vector<std::string> >("SaturationProducts");
  _TriggerProducer  = p.get<std::string>("TriggerProducer","daq");
  pmt_gain          = p.get<std::vector<float> >("PMTGains");
  pmt_gainerr       = p.get<std::vector<float> >("PMTGainErrors");
  _usePmtGainDB     = p.get<bool>("usePmtGainDB");
  _remap_ch         = p.get<bool>("RemapCh");
  _useExtSat        = p.get<bool>("ExtSaturation",false);
  _OpDetFreq        = p.get<float>("OpDetFreq");
  _saveAnaTree      = p.get<bool>("SaveAnaTree");

  // configure
  flash_pset.set_do_swap_channels(_remap_ch);
  flash_pset.set_tick_width_us(1./_OpDetFreq*1.e6);
  flash_pset.Check_common_parameters();
  flash_algo.Configure(flash_pset);

  //make ana tree
  if(_saveAnaTree){
    art::ServiceHandle< art::TFileService > tfs;
    _outtree = tfs->make<TTree>("outtree","per flash tree");
    _outtree->Branch("EventID",    &fEventID,     "EventID/I");
    _outtree->Branch("FlashID",    &fFlashID,     "FlashID/I");
    _outtree->Branch("FlashType",  &fFlashType,   "FlashType/I");
    _outtree->Branch("YCenter",    &fYCenter,     "YCenter/F");
    _outtree->Branch("ZCenter",    &fZCenter,     "ZCenter/F");
    _outtree->Branch("YWidth",     &fYWidth,      "YWidth/F");
    _outtree->Branch("ZWidth",     &fZWidth,      "ZWidth/F");
    _outtree->Branch("FlashTime",  &fFlashTime,   "FlashTime/F");
    _outtree->Branch("FlashWidth", &fFlashWidth,  "FlashWidth/F");
    _outtree->Branch("AbsTime",    &fAbsTime,     "AbsTime/F");
    _outtree->Branch("FlashFrame", &fFlashFrame,  "FlashFrame/I");
    _outtree->Branch("InBeamFrame",&fInBeamFrame, "InBeamFrame/B");
    _outtree->Branch("OnBeamTime", &fOnBeamTime,  "OnBeamTime/I");
    _outtree->Branch("TotalPE",    &fTotalPE,     "TotalPE/F");
    _outtree->Branch("PEPerCh", &fPEPerCh);
  }


  for ( unsigned int cat=0; cat<2; cat++ ) {
    produces< std::vector<recob::OpFlash> >( _flashProducts[cat] );
    produces< std::vector<raw::OpDetWaveform> >( _saturationProducts[cat] );
  }

}

void UBWCFlashFinder::produce(art::Event & evt)
{
  std::unique_ptr< std::vector<recob::OpFlash> > opflashes_beam(new std::vector<recob::OpFlash>);
  std::unique_ptr< std::vector<recob::OpFlash> > opflashes_cosmic(new std::vector<recob::OpFlash>);
  std::unique_ptr< std::vector<raw::OpDetWaveform> > saturation_beam(new std::vector<raw::OpDetWaveform>);
  std::unique_ptr< std::vector<raw::OpDetWaveform> > saturation_cosmic(new std::vector<raw::OpDetWaveform>);

  // initialize data handles and services
  art::ServiceHandle<geo::Geometry> geo;
  auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
  double triggerTime = ts->TriggerTime();

  //get gains from database if requested
  if(_usePmtGainDB){
    pmt_gain.clear();
    pmt_gainerr.clear();
    const lariov::PmtGainProvider& gain_provider = art::ServiceHandle<lariov::PmtGainService>()->GetProvider();
    for (unsigned int i=0; i!= geo->NOpDets(); ++i) {
      if (geo->IsValidOpChannel(i) && i<32) {
	pmt_gain.push_back(gain_provider.Gain(i));
	pmt_gain.push_back(gain_provider.GainErr(i));
	//pmt_gain.push_back(gain_provider.ExtraInfo(i).GetFloatData("amplitude_gain"));
	//pmt_gainerr.push_back(gain_provider.ExtraInfo(i).GetFloatData("amplitude_gain_err"));
      }
    }
  }

  //reconstruct
  flash_algo.Configure(flash_pset);
  if(!_useExtSat) reco_default(evt, triggerTime);
  else reco_external_sat(evt, triggerTime);

  //get saturation corrected wf
  auto const satb_v = flash_algo.get_merged_beam();
  auto const satc_v = flash_algo.get_merged_cosmic();
  for(const auto& isat : satb_v){
    std::vector<unsigned short> wf;
    for(unsigned int c=0; c<isat.size(); c++){
      wf.push_back(isat[c]);
    }
    raw::TimeStamp_t t = isat.get_time_from_trigger()+triggerTime;
    raw::Channel_t ch = isat.get_ChannelNum();
    raw::OpDetWaveform wfm(t,ch,wf);
    saturation_beam->emplace_back(std::move(wfm));
    wf.clear();
  }
  for(const auto& isat : satc_v){
    std::vector<unsigned short> wf;
    for(unsigned int c=0; c<isat.size(); c++){
      wf.push_back(isat[c]);
    }
    raw::TimeStamp_t t = isat.get_time_from_trigger()+triggerTime;
    raw::Channel_t ch = isat.get_ChannelNum();
    raw::OpDetWaveform wfm(t,ch,wf);
    saturation_cosmic->emplace_back(std::move(wfm));
    wf.clear();

  }

  //get flashes
  auto const flash_v = flash_algo.get_flashes();

  int idx=0;
  for(const auto& lflash :  flash_v) {
    double Ycenter, Zcenter, Ywidth, Zwidth;
    GetFlashLocation(lflash->get_pe_v(), Ycenter, Zcenter, Ywidth, Zwidth);

    recob::OpFlash flash(lflash->get_time(), lflash->get_high_time()-lflash->get_low_time(),
			 triggerTime + lflash->get_time(),
			 (triggerTime + lflash->get_time()) / 1600.,
			 lflash->get_pe_v(),
                         0, 0, 1, // this are just default values
                         Ycenter, Ywidth, Zcenter, Zwidth);
    //fill ana tree if requested
    if(_saveAnaTree) fill_ana_tree(flash,idx,lflash->get_type());
    idx++;

    if(lflash->get_type()==2) opflashes_beam->emplace_back(std::move(flash));
    else if(lflash->get_type()==1) opflashes_cosmic->emplace_back(std::move(flash));
    else{
    }
  }

  flash_algo.clear_flashes();

  evt.put(std::move(saturation_beam),_saturationProducts[0]);
  evt.put(std::move(saturation_cosmic),_saturationProducts[1]);
  evt.put(std::move(opflashes_beam),_flashProducts[0]);
  evt.put(std::move(opflashes_cosmic),_flashProducts[1]);
}

//-----------------------------------------//
// run algo with default configuration
// (using internal saturation correction)
//----------------------------------------//
void UBWCFlashFinder::reco_default(art::Event &evt, double &triggerTime){

  art::Handle< std::vector< raw::OpDetWaveform > > wfCHGHandle;
  art::Handle< std::vector< raw::OpDetWaveform > > wfBHGHandle;
  art::Handle< std::vector< raw::OpDetWaveform > > wfCLGHandle;
  art::Handle< std::vector< raw::OpDetWaveform > > wfBLGHandle;

  evt.getByLabel( _OpDataProducer, _OpDataTypes[kBeamHighGain], wfBHGHandle);
  std::vector<raw::OpDetWaveform> const& opwfms_bhg(*wfBHGHandle);
  evt.getByLabel( _OpDataProducer, _OpDataTypes[kBeamLowGain], wfBLGHandle);
  std::vector<raw::OpDetWaveform> const& opwfms_blg(*wfBLGHandle);
  evt.getByLabel( _OpDataProducer, _OpDataTypes[kCosmicHighGain], wfCHGHandle);
  std::vector<raw::OpDetWaveform> const& opwfms_chg(*wfCHGHandle);
  evt.getByLabel( _OpDataProducer, _OpDataTypes[kCosmicLowGain], wfCLGHandle);
  std::vector<raw::OpDetWaveform> const& opwfms_clg(*wfCLGHandle);

  std::vector<raw::OpDetWaveform> sort_blg;
  std::vector<raw::OpDetWaveform> sort_clg;
  if( /*evt.run() <= 3984*/ opwfms_blg.size()>opwfms_bhg.size() ){
    for (unsigned i=0; i<opwfms_blg.size(); i++){
      if(opwfms_blg.at(i).size()>=1500) sort_blg.push_back(opwfms_blg.at(i));
      if(opwfms_blg.at(i).size()<1500) sort_clg.push_back(opwfms_blg.at(i));
    }
  }
  ::wcopreco::UBEventWaveform UB_evt_wf;
  if(sort_blg.size()>0) UB_evt_wf = fill_evt_wf(triggerTime, opwfms_bhg, sort_blg, opwfms_chg, sort_clg, pmt_gain,pmt_gainerr);
  else if(sort_blg.size()<=0) UB_evt_wf = fill_evt_wf(triggerTime, opwfms_bhg, opwfms_blg, opwfms_chg, opwfms_clg, pmt_gain,pmt_gainerr);

  std::vector<wcopreco::kernel_fourier_container> kernel_container_v = fill_kernel_container( pmt_gain);
  flash_algo.SaturationCorrection(&UB_evt_wf);
  flash_algo.Run(&pmt_gain, &pmt_gainerr, &kernel_container_v);
  
}

//----------------------------------------//
// run algo with external saturation corr.
//----------------------------------------//
void UBWCFlashFinder::reco_external_sat(art::Event &evt, double &triggerTime){

  art::Handle< std::vector< raw::OpDetWaveform > > wfCHGHandle;
  art::Handle< std::vector< raw::OpDetWaveform > > wfBHGHandle;

  evt.getByLabel( _OpSatDataProducer, _OpDataTypes[kBeamHighGain], wfBHGHandle);
  std::vector<raw::OpDetWaveform> const& opwfms_bhg(*wfBHGHandle);
  evt.getByLabel( _OpSatDataProducer, _OpDataTypes[kCosmicHighGain], wfCHGHandle);
  std::vector<raw::OpDetWaveform> const& opwfms_chg(*wfCHGHandle);

  ::wcopreco::OpWaveformCollection BHG_wfm_collection;
  BHG_wfm_collection.set_op_gain(pmt_gain);
  BHG_wfm_collection.set_op_gainerror(pmt_gainerr);

  ::wcopreco::OpWaveformCollection CHG_wfm_collection;
  CHG_wfm_collection.set_op_gain(pmt_gain);
  CHG_wfm_collection.set_op_gainerror(pmt_gainerr);

  //Fill up wfm collections
  fill_wfmcollection(triggerTime, opwfms_bhg, ::wcopreco::kbeam_merged, BHG_wfm_collection);
  fill_wfmcollection(triggerTime, opwfms_chg, ::wcopreco::kcosmic_merged, CHG_wfm_collection);

  std::vector<wcopreco::kernel_fourier_container> kernel_container_v = fill_kernel_container( pmt_gain);
  flash_algo.set_merged_beam(BHG_wfm_collection);
  flash_algo.set_merged_cosmic(CHG_wfm_collection);
  flash_algo.Run(&pmt_gain, &pmt_gainerr, &kernel_container_v);

}

//--------------------------------------//
// fill ana tree
//--------------------------------------//
void UBWCFlashFinder::fill_ana_tree(recob::OpFlash flash, int idx, int type){
  
  fFlashTime   = flash.Time();
  fFlashWidth  = flash.TimeWidth(); 
  fFlashID     = idx;
  fFlashType   = type;
  fYCenter     = flash.YCenter();
  fZCenter     = flash.ZCenter();
  fYWidth      = flash.YWidth();
  fZWidth      = flash.ZWidth();
  fInBeamFrame = flash.InBeamFrame();
  fOnBeamTime  = flash.OnBeamTime();
  fAbsTime     = flash.AbsTime();
  fFlashFrame  = flash.Frame();
  fTotalPE     = flash.TotalPE();
  fPEPerCh     = flash.PEs();

  _outtree->Fill();
}

//------------------------------------------------//
// fill UBEventWaveform 
// fills beam high/low gain, cosmic high/low gain
//------------------------------------------------//
::wcopreco::UBEventWaveform UBWCFlashFinder::fill_evt_wf(double triggerTime,
							 std::vector<raw::OpDetWaveform> bhg,
							 std::vector<raw::OpDetWaveform> blg,
							 std::vector<raw::OpDetWaveform> chg,
							 std::vector<raw::OpDetWaveform> clg, 
							 std::vector<float> pmt_gain, std::vector<float> pmt_gainerr){
  
  ::wcopreco::UBEventWaveform _UB_evt_wfm;
  std::vector<wcopreco::OpWaveformCollection> empty_vec;
  _UB_evt_wfm.set_wfm_v( empty_vec );
  _UB_evt_wfm.set_op_gain(pmt_gain);
  _UB_evt_wfm.set_op_gainerror(pmt_gainerr);
  
  ::wcopreco::OpWaveformCollection CHG_wfm_collection;
  CHG_wfm_collection.set_op_gain(pmt_gain);
  CHG_wfm_collection.set_op_gainerror(pmt_gainerr);
  
  ::wcopreco::OpWaveformCollection CLG_wfm_collection;
  CLG_wfm_collection.set_op_gain(pmt_gain);
  CLG_wfm_collection.set_op_gainerror(pmt_gainerr);

  ::wcopreco::OpWaveformCollection BHG_wfm_collection;
  BHG_wfm_collection.set_op_gain(pmt_gain);
  BHG_wfm_collection.set_op_gainerror(pmt_gainerr);

  ::wcopreco::OpWaveformCollection BLG_wfm_collection;
  BLG_wfm_collection.set_op_gain(pmt_gain);
  BLG_wfm_collection.set_op_gainerror(pmt_gainerr);

  //Fill up wfm collections
  fill_wfmcollection(triggerTime, bhg, ::wcopreco::kbeam_hg, BHG_wfm_collection);
  fill_wfmcollection(triggerTime, blg, ::wcopreco::kbeam_lg, BLG_wfm_collection);
  fill_wfmcollection(triggerTime, chg, ::wcopreco::kcosmic_hg, CHG_wfm_collection);
  fill_wfmcollection(triggerTime, clg, ::wcopreco::kcosmic_lg, CLG_wfm_collection);

  _UB_evt_wfm.add_entry(BHG_wfm_collection, ::wcopreco::kbeam_hg );
  _UB_evt_wfm.add_entry(BLG_wfm_collection, ::wcopreco::kbeam_lg );
  _UB_evt_wfm.add_entry(CHG_wfm_collection, ::wcopreco::kcosmic_hg );
  _UB_evt_wfm.add_entry(CLG_wfm_collection, ::wcopreco::kcosmic_lg );

  return _UB_evt_wfm;

}

//--------------------------------------//
// fill the 4 unmerged waveform collections
// merged waveforms filled by algo.Run()
//-------------------------------------//
void UBWCFlashFinder::fill_wfmcollection(double triggerTime,
					 std::vector<raw::OpDetWaveform> opwfms,
					 int type,
					 ::wcopreco::OpWaveformCollection &wfm_collection) {
  
  for(auto &opwfm : opwfms)  {
    int ch = opwfm.ChannelNumber();
    double timestamp = opwfm.TimeStamp();

    if ( ch%100>=32 )
      continue;

    if ( (type == ::wcopreco::kbeam_hg)||(type==::wcopreco::kbeam_lg) || (type==::wcopreco::kbeam_merged) ){
      ::wcopreco::OpWaveform wfm(ch%100, timestamp-triggerTime, type, flash_pset._get_cfg_deconvolver()._get_nbins_beam());
      for (int bin=0; bin<flash_pset._get_cfg_deconvolver()._get_nbins_beam(); bin++) {
	wfm[bin]=(double)opwfm[bin];
      }
      wfm_collection.add_waveform(wfm);
    }

    if ( (type == ::wcopreco::kcosmic_hg)||(type==::wcopreco::kcosmic_lg) || (type==::wcopreco::kcosmic_merged) ){
      ::wcopreco::OpWaveform wfm(ch%100, timestamp-triggerTime, type, flash_pset._get_cfg_cophit()._get_nbins_cosmic());
      for (int bin=0; bin<flash_pset._get_cfg_cophit()._get_nbins_cosmic(); bin++) {
	wfm[bin]=(double)opwfm[bin];
      }
      wfm_collection.add_waveform(wfm);
    }

  }
  return;
}
//--------------------------------//
// make kernels for deconvolution
//--------------------------------//
std::vector<wcopreco::kernel_fourier_container> UBWCFlashFinder::fill_kernel_container(std::vector<float> pmt_gain){

  std::vector<wcopreco::kernel_fourier_container> kernel_container_v;
  kernel_container_v.resize(flash_pset._get_cfg_deconvolver()._get_num_channels());

  for (int i =0 ; i<flash_pset._get_cfg_deconvolver()._get_num_channels(); i++){ 
    ::wcopreco::UB_rc *rc; 
    if(i==28) rc=new wcopreco::UB_rc(true, true, flash_pset._get_cfg_ub_rc());
    else{ rc=new wcopreco::UB_rc(true, false, flash_pset._get_cfg_ub_rc());}
    ::wcopreco::UB_spe *spe = new wcopreco::UB_spe(true, pmt_gain.at(i), flash_pset._get_cfg_ub_spe()); //Place UB_spe on heap, so object not deleted
    kernel_container_v.at(i).add_kernel(spe);
    kernel_container_v.at(i).add_kernel(rc);
  }

  return kernel_container_v;
}

//----------------------------------------//
// calculate PE-weighted flash location
// this frunction is from UBFlashFinder
//---------------------------------------//
void UBWCFlashFinder::GetFlashLocation(std::vector<double> pePerOpChannel, 
				       double& Ycenter, 
				       double& Zcenter, 
				       double& Ywidth, 
				       double& Zwidth){

  // Reset variables
  Ycenter = Zcenter = 0.;
  Ywidth  = Zwidth  = -999.;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {
    if (opch > 31 ){
      continue;
    }

    // Get physical detector location for this opChannel
    double PMTxyz[3];
    ::pmtana::OpDetCenterFromOpChannel(opch, PMTxyz);
    // Add up the position, weighting with PEs
    sumy    += pePerOpChannel[opch]*PMTxyz[1];
    sumy2   += pePerOpChannel[opch]*PMTxyz[1]*PMTxyz[1];
    sumz    += pePerOpChannel[opch]*PMTxyz[2];
    sumz2   += pePerOpChannel[opch]*PMTxyz[2]*PMTxyz[2];

    totalPE += pePerOpChannel[opch];
  }

  Ycenter = sumy/totalPE;
  Zcenter = sumz/totalPE;

  // This is just sqrt(<x^2> - <x>^2)
  if ( (sumy2*totalPE - sumy*sumy) > 0. ) 
    Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;
  
  if ( (sumz2*totalPE - sumz*sumz) > 0. ) 
    Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
}

DEFINE_ART_MODULE(UBWCFlashFinder)
