#include "Deconvolver.h"

namespace wcopreco {


  wcopreco::Deconvolver::Deconvolver(OpWaveformCollection &merged_beam, bool with_filters, std::vector<kernel_fourier_container> &input_k_container_v, const Config_Deconvolver & cfg)
  :_cfg(cfg)
  {

    //int type = merged_beam.at(0).get_type();
    op_gain = merged_beam.get_op_gain();
    kernel_container_v = &input_k_container_v;
    filter_status = with_filters;

  }


  OpWaveformCollection wcopreco::Deconvolver::Deconvolve_Collection(OpWaveformCollection & merged_beam)

    {
      //Process the Beam:
      //Note that the following code is supposed to only deal with beam waveforms, 32 channels and 1500 bin wfms.
      OpWaveform wfm(0,0,0,0);

      for (int ch=0; ch<_cfg._num_channels; ch++){
        wfm = merged_beam.at(ch);
        nbins = wfm.size();
	
	//get the gain for this ch: need to decide whether to do deconvolution
	float chgain = merged_beam.get_op_gain().at(ch);
	
	//remove baselines (baseline here are determined by the start of the waveform)
	Remove_Baseline_Leading_Edge(wfm);

	Remove_Baseline_Secondary(wfm);
	
	if(chgain>0){	  
	  //Do deconvolution (need to add a way to incorporate kernels)
	  deconvolved_collection.add_waveform(Deconvolve_One_Wfm(wfm, kernel_container_v->at(wfm.get_ChannelNum())));
	}
	else{
	  deconvolved_collection.add_waveform(wfm);
	}
      }
      
      return deconvolved_collection;
    }//End of Deconvolve_Collection


    // void Deconvolver::add_kernel_container_entry(kernel_fourier *kernel, int channel) {
    //
    //   //If channel is ommitted from the function arguments then we add the kernel to all
    //   if (channel <0) {
    //     for (int i =0 ; i<_cfg._num_channels; i++){
    //       kernel_container_v->at(i).add_kernel(kernel);
    //     }
    //   }
    //   else if (channel < _cfg._num_channels) {
    //     kernel_container_v->at(channel).add_kernel(kernel);
    //   }
    //   else
    //   {std::cout << "You're asking to add a kernel to a nonexistant channel. Options are 0-"<<_cfg._num_channels-1<<" for individual beam channels or input channel < 0 to apply to all channels (by default params )\n";}
    //
    // }

    double Deconvolver::HighFreqFilter(double frequency)
    {
      double freq_filter = exp(-pow(frequency/_cfg._high_freq_p0,_cfg._high_freq_p1));
      return freq_filter;
    }

    double Deconvolver::LateLightFilter(double frequency2)
    {

      double latelight_filter = (1-exp(-pow(frequency2/_cfg._latelight_filter_p0,2)))*exp(-pow(frequency2/_cfg._latelight_filter_p1 , _cfg._latelight_filter_p2));
      return latelight_filter;
     }

     void Deconvolver::Remove_Baseline_Leading_Edge(OpWaveform &wfm)
     {
       int size = wfm.size();
       double baseline = wfm.at(0);
       for (int i=0; i<size; i++) {
         wfm.at(i) = wfm.at(i) - baseline;
       }
       return;
     }

     void Deconvolver::Remove_Baseline_Secondary(OpWaveform &wfm)
     {
        TH1F h1("h1","h1",200,-100,100);

        h1.Reset();
        for (int j=0;j!=_cfg._nbins_baseline_search;j++){
            h1.Fill(wfm.at(j));
        }
        double baseline = h1.GetMaximumBin()-(_cfg._baseline_safety_subtraction);
        if (fabs(baseline)>=_cfg._baseline_difference_max) {baseline = 0;}
        for (size_t j=0;j!=wfm.size();j++){
            wfm.at(j) = wfm.at(j) - baseline;
        }

        return;
     }

     std::pair<double,double> Deconvolver::cal_mean_rms(std::vector<double> wfm, int nbin)
     {
        //calculate the mean and rms values
        TH1F *h4 = new TH1F("h4","h4",2000,-10,10);
        double mean, rms;
        for (int i=0;i!=nbin;i++){
          double content = wfm.at(i);
          if (fabs(content)<10)
          h4->Fill(content);
        }
        mean = h4->GetBinCenter(h4->GetMaximumBin()+1);

        double arg = _cfg._xq;
        double par[3];

	if(h4->Integral() !=0.){
	  h4->GetQuantiles(1,&par[1],&arg);
	  arg = _cfg._xq - _cfg._xq_diff;
	  
	  h4->GetQuantiles(1,&par[0],&arg);
	  arg = _cfg._xq + _cfg._xq_diff;
	  
	  h4->GetQuantiles(1,&par[2],&arg);
	  
	  rms = sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);
	}
	else{mean = rms = 0;}
        delete h4;
        return std::make_pair(mean,rms);
     }


     double Deconvolver::KS_maxdiff(int n, double *array1, double *array2){

       for (int index =1; index<n; index++){
         array1[index] += array1[index-1];
         array2[index] += array2[index-1];
       }
       double maxdiff=0;
       double thisdiff=0;
       for (int index=0; index<n; index++){
         array1[index] = array1[index]/array1[n-1];
         array2[index] = array2[index]/array2[n-1];
         thisdiff =fabs(array1[index]-array2[index]);
         if (thisdiff > maxdiff) {maxdiff = thisdiff;}
       }
       return maxdiff;
     }

     OpWaveform Deconvolver::Deconvolve_One_Wfm(OpWaveform & wfm, const kernel_fourier_container & kernel_container) {
       //BEGIN DECONVOLUTION MARKER
       std::vector<double> wfm_doubles(wfm.begin(), wfm.end());
       float bin_width = (_cfg._tick_width_us*1e-6 ); // e-6 to go from microseconds to seconds
       int nbins = wfm.size();

       //get power spectrum of data
       //Create Mag and Phase vectors
       std::vector<double> mag_raw;
       std::vector<double> phase_raw;
       mag_raw.resize(nbins);
       phase_raw.resize(nbins);

       //Get the input to the fourier transform ready
       std::vector<double> power_spec_d(nbins,0);
       power_spec_d = wfm_doubles;
       double *in = power_spec_d.data();
       //Start the fourier transform (real to complex)
       TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &nbins, "R2C");
       fftr2c->SetPoints(in);
       fftr2c->Transform();
       double *re = new double[nbins]; //Real -> Magnitude
       double *im = new double[nbins]; //Imaginary -> Phase

       fftr2c->GetPointsComplex(re, im); //Put the values in the arrays

       //Copy those array values into vectors passed in by reference.
       // This is inefficient, but makes for an easier user interface.
       double im_i =0;
       double re_i = 0;
       for (int i =0; i<nbins; i++){
         //Calculate the phase_raw
         im_i = im[i];
         re_i = re[i];
         double ph = 0;
         if (TMath::Abs(re_i) > 1e-13){
            ph = atan2(im_i , re_i);

         }
         phase_raw.at(i) = ph;

         //End of phase_raw calc

         //Calculate the mag_v
         double magnitude = TMath::Sqrt(re[i]*re[i]+im[i]*im[i]);
         mag_raw[i] = magnitude;
         //End of mag_v calc
       }


       delete fftr2c;
       //double max_freq_MHz = (1/_cfg._tick_width_us)*_cfg._nbins_beam*2*TMath::Pi();

       double *value_re = new double[nbins];
       double *value_im = new double[nbins];
       double *value_re1 = new double[nbins];
       double *value_im1 = new double[nbins];

       //Set up the kernels to be deconvolved out:
      int channel = wfm.get_ChannelNum();
      int num_kernels = kernel_container_v->at(channel).size();
      std::vector<std::vector<double>> mag_kernel;
      mag_kernel.resize(num_kernels);
      std::vector<std::vector<double>> phase_kernel;
      phase_kernel.resize(num_kernels);


      for (int n=0; n < num_kernels; n++ ) {
        kernel_container_v->at(channel).at(n)->Get_pow_spec(nbins, bin_width, &mag_kernel.at(n), &phase_kernel.at(n));
      }


       for (int i=0;i<nbins;i++){
         double freq;
         if (i<=nbins/2){
     	     freq = ((double)i/(double)nbins*2.)*1.0;
         }
         else{
     	     freq = (((double)nbins-(double)i)/(double)nbins*2.)*1.0;
         }

         double rho = mag_raw.at(i);

         double phi = phase_raw.at(i);
         for (int n=0;n<num_kernels;n++){
           rho = rho / mag_kernel.at(n).at(i);
           phi = phi - phase_kernel.at(n).at(i);
         }


         if (i==0) rho = 0;
         //Perform Deconv with Filters
         if (filter_status){
           value_re[i] = rho * (cos(phi)/nbins) * LateLightFilter(freq);
           value_im[i] = rho * (sin(phi)/nbins) * LateLightFilter(freq);
           value_re1[i] = rho * (cos(phi)/nbins) * HighFreqFilter(freq);
           value_im1[i] = rho * (sin(phi)/nbins) * HighFreqFilter(freq);
         }
         //Perform Deconv without Filters
         else{
           value_re[i] = rho * (cos(phi)/nbins) ;
           value_im[i] = rho * (sin(phi)/nbins) ;
           value_re1[i] = rho * (cos(phi)/nbins) ;
           value_im1[i] = rho * (sin(phi)/nbins) ;
         }

       }

       // ROI finding
       TVirtualFFT *ifft = TVirtualFFT::FFT(1,&nbins,"C2R");
       ifft->SetPointsComplex(value_re,value_im);
       ifft->Transform();

       double *re_inv = new double[nbins]; //Real -> Magnitude
       double *im_inv = new double[nbins]; //Imaginary -> Phase

       ifft->GetPointsComplex(re_inv, im_inv);

       std::vector<double> inverse_res;
       inverse_res.resize(nbins);
       for (int i = 0; i < nbins ; i++){
         inverse_res.at(i) = re_inv[i];

       };
       delete ifft;

       // calculate rms and mean
       std::pair<double,double> results = cal_mean_rms(inverse_res, nbins);
       std::vector<double> hflag;
       hflag.resize(nbins);
       for (int i=0;i<nbins;i++){
         double content = inverse_res.at(i);
         if (fabs(content-results.first)>5*results.second){
	   for (int j=-20;j!=20;j++){
	     double flag =1.0;
	     if((i+j) >= 0 && (i+j) < _cfg._nbins_beam) hflag.at(i+j) = flag;
	   }
	 }
       }

       // solve for baseline
       TVirtualFFT *ifft1 = TVirtualFFT::FFT(1,&nbins,"C2R");
       ifft1->SetPointsComplex(value_re1,value_im1);
       ifft1->Transform();

       double *re_inv1 = new double[nbins]; //Real -> Magnitude
       double *im_inv1 = new double[nbins]; //Imaginary -> Phase

       ifft1->GetPointsComplex(re_inv1, im_inv1);

       OpWaveform inverse_res1(channel,wfm.get_time_from_trigger(), wfm.get_type(), nbins);
       for (int i = 0; i < nbins ; i++){
         double value = re_inv1[i];
         inverse_res1.at(i) = value;
       }

       double A11 = 0, A12 = 0, A21=0, A22=0;
       double B1 = 0, B2 = 0;
       double a=0, b=0;
       for (int i=0;i!=_cfg._nbins_beam;i++){
         double bincenter = i+.5;
         if (hflag.at(i)==0){
         	B2 += inverse_res1.at(i);
         	B1 += inverse_res1.at(i) * bincenter;
         	A11 += pow(bincenter,2);
         	A12 += bincenter;
         	A21 += bincenter;
         	A22 += 1;
         }
       }

       if (A22>0){
         a = (B1*A22-B2*A12)/(A11*A22-A21*A12);
         b = (B1*A21-B2*A11)/(A22*A11-A12*A21);
       }
       for (int i=0;i!=_cfg._nbins_beam;i++){
         double bincenter = i+.5;
         inverse_res1.at(i) = inverse_res1.at(i) - a * bincenter -b;
       }
       results = cal_mean_rms(inverse_res1, nbins);
       for (int i=0;i!=_cfg._nbins_beam;i++){
         if (i<_cfg._nbins_beam-_cfg._n_bins_end_wfm){
            inverse_res1.at(i) = inverse_res1.at(i) -results.first+_cfg._small_content_bump;
         }else{
   	       inverse_res1.at(i) = 0;
         }
       }
       delete[] re;
       delete[] im;
       delete[] re_inv;
       delete[] im_inv;
       delete[] re_inv1;
       delete[] im_inv1;
       delete[] value_re;
       delete[] value_im;
       delete[] value_re1;
       delete[] value_im1;
       //END OF DECONVOLUTION
       return inverse_res1;
     }

     // void Deconvolver::clear_kernels() {
     //   for (int i=0; i< kernel_container_v->size() ; i++) {
     //     kernel_container_v->at(i).clear_kernels_v();
     //   }
     // }



}
