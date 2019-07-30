#include "SEAviewer.h"
namespace seaview{
    //default
    SEAviewer::SEAviewer(std::string intag, geo::GeometryCore const * ingeom, detinfo::DetectorProperties const * intheDetector ): tag(intag), geom(ingeom), theDetector(intheDetector){
        std::cout<<tag<<std::endl;
        chan_max = {-9999,-9999,-9999};
        chan_min = {9999,9999,9999};
        tick_max = -99999;
        tick_min = 99999;

        vertex_tick.resize(3); 
        vertex_chan.resize(3); 
        vertex_graph.resize(3); 
        n_pfps = 0;
    }


    int SEAviewer::setBadChannelList(std::vector<std::pair<int,int>> &in){
        m_bad_channel_list = in;
        return 0;
    }

    int SEAviewer::addSliceHits(std::vector<art::Ptr<recob::Hit>>& hits){
        slice_hits = hits;   
        for(auto &h: slice_hits){
            map_unassociated_hits[h] = true;
        }
        return 0;
    }

    int SEAviewer::calcUnassociatedHits(){
        int n_assoc=0;

        std::vector<std::vector<double>>  vec_tick(3);
        std::vector<std::vector<double>>  vec_chan(3);


        for(auto&h:slice_hits){
            if(map_unassociated_hits[h]){
                n_assoc++;
                
                double wire = (double)h->WireID().Wire;
                vec_chan[(int)h->View()].push_back(wire);
                vec_tick[(int)h->View()].push_back((double)h->PeakTime());
                tick_max = std::max(tick_max, (double)h->PeakTime());
                tick_min = std::min(tick_min, (double)h->PeakTime());
                chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);
     
                }

        }
        
        for(int i=0; i<3; i++){
            vec_unass_graphs.emplace_back(vec_tick[i].size(),&vec_chan[i][0],&vec_tick[i][0]); 
        }

        vec_unass_ticks = vec_tick;
        vec_unass_chans = vec_chan;

        return n_assoc;
    }

    int SEAviewer::addPFParticleHits(std::vector<art::Ptr<recob::Hit>>& hits){
        n_pfps++;

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


    int SEAviewer::Print(){


        std::string print_name = "SEAview_"+tag;
        TCanvas *can=new TCanvas(print_name.c_str(),print_name.c_str(),3000,800);
        can->Divide(4,1,0,0.1);


        //******************************* First plot "Vertex" ***************************************

        for(int i=0; i<3; i++){
            TPad * pader = (TPad*)can->cd(i+1);

            if(i==0 || i ==4 || i == 8) pader->SetLeftMargin(0.1);


            vertex_graph[i].SetMarkerStyle(29);
            vertex_graph[i].SetMarkerSize(4);
            vertex_graph[i].SetMarkerColor(kMagenta-3);
            vertex_graph[i].GetYaxis()->SetRangeUser(tick_min*0.98,tick_max*1.02);
            vertex_graph[i].GetXaxis()->SetLimits(chan_min[i]*0.98,chan_max[i]*1.02);
            vertex_graph[i].SetTitle(("Plane " +std::to_string(i)).c_str());
            vertex_graph[i].GetYaxis()->SetTitle("Peak Hit Time Tick");
            vertex_graph[i].GetXaxis()->SetTitle( ("Wire Number Plane " +std::to_string(i)).c_str());
            vertex_graph[i].Draw("ap");

            if(i>0){
                vertex_graph[i].GetYaxis()->SetLabelOffset(999);
                vertex_graph[i].GetYaxis()->SetLabelSize(0);
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

            if(chan_min[thisp]*0.98 < bc && bc < chan_max[thisp]*1.02 ){
                can->cd(thisp+1);
                TLine *l = new TLine(bc,tick_min*0.98,bc,tick_max*1.02);
                l->SetLineColor(kGray+1);
                l->Draw("same");
                can->cd(thisp+5);
                l->Draw("same");
                can->cd(thisp+9);
                l->Draw("same");
            }
        }


        ///******************************** Plotting all PFP's *********************************8
        double plot_point_size=0.6;        

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




         can->Update();
         can->SaveAs((print_name+".pdf").c_str(),"pdf");




        return 0;
    }

}
