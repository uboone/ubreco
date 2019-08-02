#include "SEAviewer.h"
namespace seaview{
    //default
    SEAviewer::SEAviewer(std::string intag, geo::GeometryCore const * ingeom, detinfo::DetectorProperties const * intheDetector ): tag(intag), geom(ingeom), theDetector(intheDetector){
        chan_max = {-9999,-9999,-9999};
        chan_min = {9999,9999,9999};
        tick_max = -99999;
        tick_min = 99999;

        plot_true_vertex = false;
        vertex_tick.resize(3); 
        vertex_chan.resize(3); 
        vertex_graph.resize(3); 

        true_vertex_tick.resize(3); 
        true_vertex_chan.resize(3); 
        true_vertex_graph.resize(3); 

        tick_shift = 350;
        chan_shift = 100;

        n_showers=0;
        n_pfps = 0;
        has_been_clustered = false;
        hit_threshold = -10;

        rangen = new TRandom3(0);
    }


    int SEAviewer::setBadChannelList(std::vector<std::pair<int,int>> &in){
        m_bad_channel_list = in;
        return 0;
    }

    int SEAviewer::addSliceHits(std::vector<art::Ptr<recob::Hit>>& hits){
        slice_hits = hits;   
        for(auto &h: slice_hits){
            map_unassociated_hits[h] = true;
            map_slice_hits[h] = true;
        }
        return 0;
    }

    int SEAviewer::setHitThreshold(double h){
        hit_threshold = h;
        return 0;    
    }
    int SEAviewer::addAllHits(std::vector<art::Ptr<recob::Hit>>& hits){
        all_hits = hits;

        std::vector<std::vector<double>>  vec_tick(3);
        std::vector<std::vector<double>>  vec_chan(3);

        for(auto&h:all_hits){
            if(map_slice_hits.count(h)==0){
                double wire = (double)h->WireID().Wire;
                vec_chan[(int)h->View()].push_back(wire);
                vec_tick[(int)h->View()].push_back((double)h->PeakTime());
                //tick_max = std::max(tick_max, (double)h->PeakTime());
                //tick_min = std::min(tick_min, (double)h->PeakTime());
                //chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                //chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);
            }
        }

        for(int i=0; i<3; i++){
            vec_all_graphs.emplace_back(vec_tick[i].size(),&vec_chan[i][0],&vec_tick[i][0]); 
        }

        vec_all_ticks = vec_tick;
        vec_all_chans = vec_chan;

        return 0;
    }

    int SEAviewer::calcUnassociatedHits(){
        int n_assoc=0;

        std::vector<std::vector<double>>  vec_tick(3);
        std::vector<std::vector<double>>  vec_chan(3);
        std::vector<std::vector<std::vector<double>>> vec_pts(3);
        std::vector<std::vector<art::Ptr<recob::Hit>>> vec_hits(3);

        for(auto&h:slice_hits){
            if(map_unassociated_hits[h]){


                if(h->SummedADC() < hit_threshold) continue;

                n_assoc++;
                double wire = (double)h->WireID().Wire;
                vec_chan[(int)h->View()].push_back(wire);
                vec_tick[(int)h->View()].push_back((double)h->PeakTime());

                vec_pts[(int)h->View()].push_back({wire,(double)h->PeakTime()});
                vec_hits[(int)h->View()].push_back(h);
                //   tick_max = std::max(tick_max, (double)h->PeakTime());
                //   tick_min = std::min(tick_min, (double)h->PeakTime());
                //   chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                //   chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

            }

        }

        for(int i=0; i<3; i++){
            vec_unass_graphs.emplace_back(vec_tick[i].size(),&vec_chan[i][0],&vec_tick[i][0]); 
        }

        vec_unass_ticks = vec_tick;
        vec_unass_chans = vec_chan;
        vec_unass_pts = vec_pts;
        vec_unass_hits = vec_hits;

        return n_assoc;
    }

        
    int SEAviewer::addPFParticleHits(std::vector<art::Ptr<recob::Hit>>& hits, std::string legend){
        n_pfps++;

        vec_pfp_legend.push_back(legend);

        std::vector<std::vector<double>>  vec_tick(3);
        std::vector<std::vector<double>>  vec_chan(3);

        for(auto &h: hits){
            double wire = (double)h->WireID().Wire;
            vec_chan[(int)h->View()].push_back(wire);
            vec_tick[(int)h->View()].push_back((double)h->PeakTime());
            tick_max = std::max(tick_max, (double)h->PeakTime());
            tick_min = std::min(tick_min, (double)h->PeakTime());
            chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
            chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

            //remove from unassociated hits
            map_unassociated_hits[h] = false;
        }

        std::vector<TGraph> t_graphs;
        for(int i=0; i<3; i++){
            t_graphs.emplace_back(vec_tick[i].size(),&vec_chan[i][0],&vec_tick[i][0]); 
        }
        vec_graphs.push_back(t_graphs);
        vec_ticks.push_back(vec_tick);
        vec_chans.push_back(vec_chan);
        return 0;
    }

    int SEAviewer::addShower(art::Ptr<recob::Shower>&shr){
        n_showers++;

        vec_showers.push_back(shr); 
    
        return 0;
    }


    std::vector<std::vector<double>> SEAviewer::to2D(std::vector<double> & threeD){

        auto const TPC = (*geom).begin_TPC();
        auto ID = TPC.ID();
        int fCryostat = ID.Cryostat;
        int fTPC = ID.TPC;
        
        std::vector<std::vector<double>> ans(3);

        for(int i=0; i<3; i++){
            double wire = (double)calcWire(threeD[1], threeD[2], i, fTPC, fCryostat, *geom);
            double time = calcTime(threeD[0], i, fTPC,fCryostat, *theDetector);

            ans[i] = {wire,time};
        }

        return ans;
    }



    int SEAviewer::loadVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z){

        auto const TPC = (*geom).begin_TPC();
        auto ID = TPC.ID();
        int fCryostat = ID.Cryostat;
        int fTPC = ID.TPC;

        for(int i=0; i<3; i++){

            std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *geom)};
            std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, *theDetector)};

            vertex_tick[i] = time[0];
            vertex_chan[i] = wire[0];

            chan_max[i] = std::max( chan_max[i],wire[0]);
            chan_min[i] = std::min( chan_min[i],wire[0]);

            TGraph gtmp(1, &wire[0], &time[0]); 
            vertex_graph[i] = gtmp;
        }

        return 0;
    }


    int SEAviewer::addTrueVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z){

        plot_true_vertex = true;

        auto const TPC = (*geom).begin_TPC();
        auto ID = TPC.ID();
        int fCryostat = ID.Cryostat;
        int fTPC = ID.TPC;

        for(int i=0; i<3; i++){

            std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *geom)};
            std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, *theDetector)};

            true_vertex_tick[i] = time[0];
            true_vertex_chan[i] = wire[0];

            chan_max[i] = std::max( chan_max[i],wire[0]);
            chan_min[i] = std::min( chan_min[i],wire[0]);

            TGraph gtmp(1, &wire[0], &time[0]); 
            true_vertex_graph[i] = gtmp;
        }

        return 0;
    }


    int SEAviewer::Print(){


        std::string print_name = "SEAview_"+tag;
        TCanvas *can=new TCanvas(print_name.c_str(),print_name.c_str(),3000,800);
        can->Divide(4,1,0,0.1);

        double plot_point_size=0.4;        

        //******************************* First plot "Vertex" ***************************************

        for(int i=0; i<3; i++){
            TPad * pader = (TPad*)can->cd(i+1);

            if(i==0 || i ==4 || i == 8) pader->SetLeftMargin(0.1);


            vertex_graph[i].SetMarkerStyle(29);
            vertex_graph[i].SetMarkerSize(2);
            vertex_graph[i].SetMarkerColor(kMagenta-3);
            vertex_graph[i].GetYaxis()->SetRangeUser(tick_min-tick_shift,tick_max+tick_shift);
            vertex_graph[i].GetXaxis()->SetLimits(chan_min[i]-chan_shift,chan_max[i]+chan_shift);
            vertex_graph[i].SetTitle(("Plane " +std::to_string(i)).c_str());
            vertex_graph[i].GetYaxis()->SetTitle("Peak Hit Time Tick");
            vertex_graph[i].GetXaxis()->SetTitle( ("Wire Number Plane " +std::to_string(i)).c_str());
            vertex_graph[i].Draw("ap");

            if(i>0){
                vertex_graph[i].GetYaxis()->SetLabelOffset(999);
                vertex_graph[i].GetYaxis()->SetLabelSize(0);
            }
        }

        /********************************* Non Slice  Hits ****************************/


        for(int i=0; i<3; i++){
            can->cd(i+1);
            if(vec_all_graphs[i].GetN()>0){//need a check in case this track has no hits on this plane.
                vec_all_graphs[i].Draw("p same"); 
                vec_all_graphs[i].SetMarkerColor(kGray);
                vec_all_graphs[i].SetFillColor(kWhite);
                vec_all_graphs[i].SetMarkerStyle(20);
                vec_all_graphs[i].SetMarkerSize(plot_point_size*0.75);
            }
        }

        //******************************** DeadWireRegions********************************************
        for(size_t i=0; i< m_bad_channel_list.size(); i++){
            int badchan = m_bad_channel_list[i].first;                                       
            int ok = m_bad_channel_list[i].second;       

            if(ok>1)continue;
            auto hs = geom->ChannelToWire(badchan);

            int thisp = (int)hs[0].Plane;
            double bc = hs[0].Wire;

            if(chan_min[thisp]-chan_shift < bc && bc < chan_max[thisp]+chan_shift ){
                can->cd(thisp+1);
                TLine *l = new TLine(bc,tick_min-tick_shift,bc,tick_max+tick_shift);
                l->SetLineColor(kGray+1);
                l->Draw("same");
                can->cd(thisp+5);
                l->Draw("same");
                can->cd(thisp+9);
                l->Draw("same");
            }
        }


        ///******************************** Plotting all PFP's *********************************8

        std::vector<int> tcols = {kRed-7, kBlue-7, kGreen-3, kOrange-3, kCyan-3, kMagenta-3, kGreen+1 , kRed+1};
        int used_col=0;

        if(n_pfps > (int)tcols.size()){
            for(int i =0; i< (int)(n_pfps +2); i++){
                //tcols.push_back(tcols[(int)rangen->Uniform(0,7)]+(int)rangen->Uniform(-5,5));
                tcols.push_back(kRed);
            }
        }


        for(int p=0; p<n_pfps; p++){

            int tcol = tcols[used_col];
            used_col++;

            for(int i=0; i<3; i++){
                can->cd(i+1);
                if(vec_graphs[p][i].GetN()>0){//need a check in case this track has no hits on this plane.

                    vec_graphs[p][i].Draw("p same"); 
                    vec_graphs[p][i].SetMarkerColor(tcol);
                    vec_graphs[p][i].SetFillColor(tcol);
                    vec_graphs[p][i].SetMarkerStyle(20);
                    vec_graphs[p][i].SetMarkerSize(plot_point_size);
                }
            }
        }


        //Plot all Shower lines. Might need a bit of work here..
        std::vector<TLine*> lines;  

        for(size_t s=0; s<vec_showers.size(); ++s){
            std::vector<double> shr_start_3D= {vec_showers[s]->ShowerStart().X(), vec_showers[s]->ShowerStart().Y(),vec_showers[s]->ShowerStart().Z()};
            std::vector<double> shr_other_3D=  {vec_showers[s]->ShowerStart().X()+vec_showers[s]->Direction().X(),vec_showers[s]->ShowerStart().Y()+vec_showers[s]->Direction().Y(), vec_showers[s]->ShowerStart().Z()+vec_showers[s]->Direction().Z()};

            //std::cout<<" "<<shr_start_3D[0]<<" "<<shr_start_3D[1]<<" "<<shr_start_3D[2]<<std::endl;
            //std::cout<<" "<<shr_other_3D[0]<<" "<<shr_other_3D[1]<<" "<<shr_other_3D[2]<<std::endl;

            std::vector<std::vector<double>> start_pt =   this->to2D(shr_start_3D);
            std::vector<std::vector<double>> other_pt =   this->to2D(shr_other_3D);
          
            for(int i=0; i<3; i++){
                //std::cout<<start_pt[i][0]<<" "<<start_pt[i][1]<<" "<<other_pt[i][0]<<" "<<other_pt[i][1]<<std::endl;
                can->cd(i+1);
                double slope = (start_pt[i][1]-other_pt[i][1])/(start_pt[i][0]-other_pt[i][0]);
                double inter = start_pt[i][1]-slope*start_pt[i][0];

                double x1_plot = other_pt[i][0];//chan_min[i]-chan_shift;
                double y1_plot = slope*x1_plot+inter;

                double x2_plot;
                if(other_pt[i][0]<start_pt[i][0]){
                    x2_plot = chan_max[i]+chan_shift;
                }else{
                    x2_plot = chan_min[i]-chan_shift;    
                }
                double y2_plot = slope*x2_plot+inter;
            
                can->cd(i+1);
                TLine *l = new TLine(x1_plot, y1_plot, x2_plot, y2_plot);
                lines.push_back(l);
                l->SetLineColorAlpha(tcols[s],0.5);
                l->SetLineWidth(1);
                l->SetLineStyle(2);
                l->Draw();
            
            }

        }




        //If its be clusterized, plot clusters here. Lets try a color surrounding the black.

        if(has_been_clustered){ 

            std::vector<int> cluster_colors(vec_clusters.size()+1,0);
            std::vector<int> base_col = {632,416, 600, 400, 616,  432,  800, 820,  840, 860,  880, 900};

            for(size_t j=0; j< vec_clusters.size()+1; j++){
                int b = (int)rangen->Uniform(0,11);
                int mod = (int)rangen->Uniform(-10,+3);

                cluster_colors[j] = base_col[b]+mod;
            }
            int c_offset = 0;


            for(auto &c: vec_clusters){
                int pl = c.getPlane();
                can->cd(pl+1);    
                    if (c.getGraph()->GetN()>0){
                        c.getGraph()->Draw("p same");
                        c.getGraph()->SetMarkerColor(cluster_colors[c_offset]);
                        c.getGraph()->SetFillColor(cluster_colors[c_offset]);
                        c.getGraph()->SetMarkerStyle(20);
                        c.getGraph()->SetMarkerSize(plot_point_size*2.0);
                        c_offset++;
                    }
            }
        }//end clusters


        /********************************* Unassociated Hits ****************************/
        for(int i=0; i<3; i++){
            can->cd(i+1);
            if(vec_unass_graphs[i].GetN()>0){//need a check in case this track has no hits on this plane.

                vec_unass_graphs[i].Draw("p same"); 
                vec_unass_graphs[i].SetMarkerColor(kBlack);
                vec_unass_graphs[i].SetFillColor(kBlack);
                vec_unass_graphs[i].SetMarkerStyle(20);
                vec_unass_graphs[i].SetMarkerSize(plot_point_size);
            }
        }

        //****** just plto vertex again with elipse;
        for(int i=0; i<3; i++){
            can->cd(i+1);
            vertex_graph[i].Draw("p same");

            double rad_cm = 12.0;
            TEllipse ell_p(vertex_chan[i],vertex_tick[i],rad_cm/0.3,rad_cm*25);
            ell_p.SetLineColor(kRed);
            ell_p.SetFillStyle(0);
            ell_p.Draw("same");

        }



        //**************************** INFO ***************************/
        TPad *p_top_info = (TPad*)can->cd(4);
        p_top_info->cd();

        /*TLatex pottex;
          pottex.SetTextSize(0.045);
          pottex.SetTextAlign(13);  //align at top
          pottex.SetNDC();
          std::string pot_draw = "Run: "+std::to_string(m_run_number)+" SubRun: "+std::to_string(m_subrun_number)+" Event: "+std::to_string(m_event_number);
          pottex.DrawLatex(.1,.94, pot_draw.c_str());
          */
        TLegend l_top(0.1,0.1,0.9,0.9);

        for(int p=0; p<n_pfps; p++){


            if(vec_graphs[p][0].GetN()>0){
                l_top.AddEntry(&vec_graphs[p][0],vec_pfp_legend[p].c_str(),"f");
            }else if(vec_graphs[p][1].GetN()>0){
                l_top.AddEntry(&vec_graphs[p][1],vec_pfp_legend[p].c_str(),"f");
            }else if(vec_graphs[p][2].GetN()>0){
                l_top.AddEntry(&vec_graphs[p][2],vec_pfp_legend[p].c_str(),"f");
            }

        }
        l_top.SetHeader(print_name.c_str(),"C");
        l_top.SetLineWidth(0);
        l_top.SetLineColor(kWhite);
        l_top.Draw("same");





        can->Update();
        can->SaveAs((print_name+".pdf").c_str(),"pdf");




        return 0;
    }

    int SEAviewer::runDBSCAN(double min_pts, double eps){

        has_been_clustered = true;
        num_clusters = {0,0,0};
        cluster_labels.resize(3);

        for(int i=0; i<3; i++){

            //      std::cout<<"SinglePhoton::SSS\t||\tStarting to run DBSCAN for plane: "<<i<<" has "<<pts_to_recluster[i].size()<<" pts to do using eps: "<<eps<<" and min_pts: "<<min_pts<<std::endl; 
            DBSCAN ReCluster(eps,min_pts);
            cluster_labels[i] =  ReCluster.Scan2D(vec_unass_pts[i]);

            for(auto &c: cluster_labels[i]){
                num_clusters[i] = std::max(c,num_clusters[i]);

            }




            std::cout << "DBSCAN\t||\tOn this plane "<<i<<" DBSCAN found: "<<num_clusters[i]<<" clusters"<<std::endl;
        }

        //Now fill the clusters

        for(int i=0; i<3; i++){
           
            for(int c=0; c<num_clusters[i]; c++){
   
                std::vector<std::vector<double>> pts;
                std::vector<art::Ptr<recob::Hit>> hitz;
                for(size_t p=0; p< vec_unass_pts[i].size(); p++){
                    if(cluster_labels[i][p] == 0) continue;//noise 
                    if(cluster_labels[i][p] == c){
                        
                        pts.push_back(vec_unass_pts[i][p]);
                        hitz.push_back(vec_unass_hits[i][p]);
                    }
                
                }
                vec_clusters.emplace_back(c,i,pts,hitz);
            }
        }

        return 0;

    }




}
