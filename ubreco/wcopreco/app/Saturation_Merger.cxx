#include "Saturation_Merger.h"


namespace wcopreco {

  Saturation_Merger::Saturation_Merger(UBEventWaveform UB_Ev, const Config_Saturation_Merger &cfg)
  :_cfg(cfg) {

    OpWaveformCollection BHG_WFs;
    OpWaveformCollection BLG_WFs;
    OpWaveformCollection CHG_WFs;
    OpWaveformCollection CLG_WFs;

    BHG_WFs = (UB_Ev.get_wfm_v() [kbeam_hg]);
    BLG_WFs = (UB_Ev.get_wfm_v() [kbeam_lg]);
    CHG_WFs = (UB_Ev.get_wfm_v() [kcosmic_hg]);
    CLG_WFs = (UB_Ev.get_wfm_v() [kcosmic_lg]);

    scale_lowgains(&BLG_WFs,&CLG_WFs);

    //Set data member merged versions of waveforms
    merged_beam = *beam_merger(&BHG_WFs, &BLG_WFs);
    merged_cosmic = *cosmic_merger(&CHG_WFs, &CLG_WFs);

    //Make all the individual waveforms the new type (merged beam or merged cosmic (5 and 6))
    for (size_t n = 0; n< merged_beam.size(); n++){
      merged_beam.at(n).set_type(kbeam_merged);
    }
    for (size_t n = 0; n< merged_cosmic.size(); n++){
      merged_cosmic.at(n).set_type(kcosmic_merged);
    }


    UB_Ev_Merged.add_entry(merged_beam,   kbeam_merged );
    UB_Ev_Merged.add_entry(merged_cosmic, kcosmic_merged );
    UB_Ev_Merged.set_op_gain(   UB_Ev.get_op_gain()   );
    UB_Ev_Merged.set_op_gainerror( UB_Ev.get_op_gainerror()   );

  }//End of Class Constructor

  void Saturation_Merger::scale_lowgains(OpWaveformCollection *BLG_WFs, OpWaveformCollection *CLG_WFs){
      //First lets do the Beam Low Gain Waveform Rescaling
      double baseline;
      double temp_baseline;
      float scalefactor;
      int channel=-1;
      int nbins;
      //Start with a loop through the beam lg waveforms
      for (size_t n=0; n < BLG_WFs->size() ; n++){
        baseline =_cfg._baseline_default;
        nbins = BLG_WFs->at(n).size();
        temp_baseline = findBaselineLg(&BLG_WFs->at(n), nbins);


        if (fabs(temp_baseline-baseline) <= _cfg._baseline_difference_max ) {

          baseline = temp_baseline;

        }
        //Now let's do the actual rescaling by looping through the current waveform's bins:
        channel =  BLG_WFs->at(n).get_ChannelNum();
        scalefactor = findScaling(channel);

        for (int bin = 0; bin < nbins; bin++){

          BLG_WFs->at(n).at(bin) = floor(( (BLG_WFs->at(n).at(bin)-baseline)*scalefactor ) + baseline) ;

        }

      }

      //Now let's scale the Cosmic Lowgains,  easier
      for (size_t n=0; n < CLG_WFs->size() ; n++){
        nbins = CLG_WFs->at(n).size();
        //Baseline is first value at first tick
        baseline = CLG_WFs->at(n).at(0);
        //Now let's do the actual rescaling by looping through the current waveform's bins:
        channel =  CLG_WFs->at(n).get_ChannelNum();
        scalefactor = findScaling(channel);
        for (int bin = 0; bin < nbins; bin++){
          CLG_WFs->at(n).at(bin) = floor(( (CLG_WFs->at(n).at(bin)-baseline)*scalefactor ) + baseline) ;
        }
      }


  }//End of Function

  double wcopreco::Saturation_Merger::findBaselineLg(OpWaveform *wfm, int nbin){
    TH1F *h = new TH1F("h","",1000,_cfg._low_bound_baseline_search-0.5,_cfg._high_bound_baseline_search-0.5);
    double baseline=0;
    for(int i=0; i!=_cfg._nbins_baseline_search; i++){
      double content = wfm->at(i);
      //    baseline += content;
      if(content>_cfg._low_bound_baseline_search && content<_cfg._high_bound_baseline_search){ h->Fill(content); }
    }
    //  baseline /= 6.;
    baseline = h->GetBinCenter(h->GetMaximumBin()+1);
    delete h;
    return baseline;
  }//End of Function

OpWaveformCollection* Saturation_Merger::cosmic_merger(OpWaveformCollection* CHG, OpWaveformCollection* CLG){
  short saturation_threshold =_cfg._sat_threshold;
  int tick_window = _cfg._cosmic_tick_window;
  float tick = _cfg._tick_width_us;

  /*
  Loop through high gain discriminators, if a waveform is saturated:
  -- if isolated from all low gain discriminators add to output
  -- if not isolated add the corresponding low gain discriminator
  */

  //Create a map to ensure we only use each lowgain wfm at most once.
  std::map<int,bool> is_used;
  for (size_t i =0;i<CLG->size();i++){
    is_used[i] = false;
  }
  bool is_saturated;
  //bool is_paired;
  int count_bin_sat;

  for (size_t idx_chg = 0; idx_chg<CHG->size(); idx_chg++){
    is_saturated = false;
    count_bin_sat =0;
    //Loop through all bins in this waveform to see if saturated
    for (size_t bin =0; bin<CHG->at(idx_chg).size(); bin++){
      if (CHG->at(idx_chg).at(bin) > saturation_threshold) {
        count_bin_sat++ ;
        if(count_bin_sat >= 3){
          //Consider it saturated
          is_saturated = true ;
          break;
        }
      }
    }
    //If it's saturated then we need to look for a low gain pair
    if (is_saturated){
      //is_paired =0;
      //Get Channel and timestamp of CHG for matching
      float time_High = CHG->at(idx_chg).get_time_from_trigger();
      short ch_High = CHG->at(idx_chg).get_ChannelNum();

      //Loop through all LG looking for a friend
      for (size_t idx_clg = 0; idx_clg<CLG->size(); idx_clg++){
        //If this lowgain waveform has already been used, skip to the next one
        if (is_used[idx_clg]) {
          continue;
        }
        //Get potential friend's channel and timestamp
        float time_Low = CLG->at(idx_clg).get_time_from_trigger();
        short ch_Low = CLG->at(idx_clg).get_ChannelNum();

        if (ch_Low == ch_High && fabs(time_High-time_Low) < (float)tick_window*(float)tick){
          CHG->at(idx_chg) = CLG->at(idx_clg);
          is_used[idx_clg] = true; //set this LG wfm to true so we don't use it again.
          //is_paired = true;
          break; //Found pair, don't have to keep looking
        }

      }//End of loop through low gains
      // If not paired
      // Do nothing to change this high gain waveform, is saturated and unfixable
      // Keep it in the collection to be returned
    }
    //Since this is a boolean it should always be true or false, but else if for check
    else if (is_saturated ==false){
      //It was not saturated.
      //Keep it in the collection to be returned.
      //Go Through Low Gains, and mark any that match it as used so they
      //don't go in the end collection as well.

      for (size_t idx_clg=0; idx_clg<CLG->size(); idx_clg++){
        if (is_used[idx_clg] ==true) {continue;} //Skip index if already used, save computation
        float time_High = CHG->at(idx_chg).get_time_from_trigger();
        float time_Low= CLG->at(idx_clg).get_time_from_trigger();
        short ch_High = CHG->at(idx_chg).get_ChannelNum();
        short ch_Low = CLG->at(idx_clg).get_ChannelNum();
        if ((ch_Low == ch_High) && (fabs(time_High-time_Low) < tick_window*tick)){

          //std::cout << "Time Difference: " <<  fabs(time_High-time_Low) <<std::endl;
          is_used[idx_clg] =true;
          break; //Found a lowgain waveform that matched, stop searching
        }
      }


    }
    else {
      std::cout << "XXXXXXXXXXXXXXXXXXXXXXXxxxxxxxx PROBLEM HG meets no criteria! xxxxxxxxxxxxxXxXXXXXXXXXXXX";
    }
  }
  //Now we add in all the unpaired low gain waveforms.
  for (size_t idx_clg=0; idx_clg<CLG->size(); idx_clg++){
    if (is_used[idx_clg] ==false ){
      //If this LG isn't used yet then add it in to the merged collection (CHG)
      CHG->add_waveform( CLG->at(idx_clg) );
    }
  }
  return CHG;
}

OpWaveformCollection* Saturation_Merger::beam_merger(OpWaveformCollection* BHG, OpWaveformCollection* BLG) {
  OpWaveformCollection* merged_beam = NULL;
  if ((BHG->size() != BLG->size() ) && (BHG->size()>0)) {
    std::cout << "Beam High Gain Collection and Beam Low Gain Collection do not have the same number of entries. Returning Empy merge\n";
    return merged_beam;
  }
  short saturation_threshold =_cfg._sat_threshold;
  bool is_saturated;
  int count_bin_sat=0;

  int idx_ch_hg ;
  int idx_ch_lg ;

  //Loop through each waveform
  //(generally one per channel so 32 or 36 waveforms for each gain of beam)
  for(int i =0; i<_cfg._num_channels; i++){

    //These calls get the index of the wfm that occurs at the channel i
    //for each gain, in case the two Opwaveformcollections aren't ordered the same

    idx_ch_hg = BHG->get_channel2index(i).at(0);
    idx_ch_lg = BLG->get_channel2index(i).at(0);

    is_saturated = false;
    count_bin_sat=0;
    //Now loop through HG wfm to see if saturated (3+ ticks at value >4050)
    for (size_t bin =0; bin<BHG->at(idx_ch_hg).size(); bin++){
      if (BHG->at(idx_ch_hg).at(bin) > saturation_threshold) {
        count_bin_sat++ ;
        if(count_bin_sat >= _cfg._nbins_saturation_threshold){
          //Consider it saturated
          is_saturated = true ;
          break;
        }
      }
    }
    if(is_saturated) {
          //If high gain waveform is saturated, replace the saturated sections with part of the low gain wfm
          std::vector<std::pair<short,short> > tickVec = findSaturationTick(&BHG->at(idx_ch_hg), saturation_threshold);
          BHG->at(idx_ch_hg) = replaceSaturatedBin( (BHG->at(idx_ch_hg)), (BLG->at(idx_ch_lg)), tickVec);

    }
  }

  return BHG;
}


std::vector<std::pair<short,short> > Saturation_Merger::findSaturationTick(OpWaveform *wfm, short saturation_threshold){
  std::vector<std::pair<short,short> > result;
  bool saturatedStatus = false;
  std::pair<short,short> tempPair;

  for(int i=0; i<(int)wfm->size(); i++){
    if(wfm->at(i)>saturation_threshold){
      if(saturatedStatus == false){ tempPair.first = i; }
      saturatedStatus = true;
    }
    if(wfm->at(i)<saturation_threshold && saturatedStatus == true){
      saturatedStatus = false;
      tempPair.second = i;
      result.push_back(tempPair);
    }
  }
  return result;
}

OpWaveform wcopreco::Saturation_Merger::replaceSaturatedBin(OpWaveform &high, OpWaveform &low, std::vector<std::pair<short,short>> saturation_ranges){

  for(int i=0; i<(int)saturation_ranges.size(); i++){
    for(short j=saturation_ranges.at(i).first; j<saturation_ranges.at(i).second; j++){
      high.at(j) = low.at(j);
    }
  }
  return high;
}

  float Saturation_Merger::findScaling(size_t channel){

    if(channel > (_cfg._scaling_by_channel.size() - 1)) {std::cout << "    Error in saturation_merger::findScaling, you you input channel bigger than 31:  " << channel << std::endl;
    return 0;}
    else {return _cfg._scaling_by_channel[channel];}
  }//End of Function

}//end of namespace
