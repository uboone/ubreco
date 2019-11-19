//
// SelectNeutrino class
//
#include <fstream>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace filter {

  class SelectNeutrino : public art::EDFilter  {

  public:

    explicit SelectNeutrino(fhicl::ParameterSet const& );

    bool filter(art::Event& evt) override;

    std::vector < unsigned int > SetOfBadEvents()   const { return fBadEvents;}
    std::vector < unsigned int > SetOfBadRuns()     const { return fBadRuns;  }

  private:

    std::vector < unsigned int >            fBadEvents;
    std::vector < unsigned int >            fBadRuns;

    std::vector < unsigned int >            fSelEvents;
    std::vector < unsigned int >            fSelRuns;
    std::vector < unsigned int >            fSelSubRuns;
    std::string fEventList;
    int         fSelection; //0: reject events based on input
    //>0: accept events based on txt file
    //<0: reject events based on txt file


  }; //class SelectNeutrino
}


filter::SelectNeutrino::SelectNeutrino(fhicl::ParameterSet const& pset)
: EDFilter{pset}
{
  fBadEvents  = pset.get < std::vector <unsigned int> >("BadEvents");
  fBadRuns    = pset.get < std::vector <unsigned int> >("BadRuns");

  fSelection = pset.get< int >("Selection");
  fEventList = pset.get< std::string >("EventList");
  fSelEvents.clear();
  fSelRuns.clear();
  std::ifstream in;
  in.open(fEventList.c_str());
  char line[1024];
  while(1){
    in.getline(line,1024);
    if (!in.good()) break;
    unsigned int n0, n1, n2;
    sscanf(line,"%u %u %u",&n0,&n1,&n2);
    fSelRuns.push_back(n0);
    fSelSubRuns.push_back(n1);
    fSelEvents.push_back(n2);
  }
  in.close();
}

bool filter::SelectNeutrino::filter(art::Event &evt)
{
  unsigned int evtNo = (unsigned int) evt.id().event();
  unsigned int runNo = (unsigned int) evt.run();
  unsigned int subrunNo = (unsigned int) evt.subRun();

  std::cout << "SELNU run " << runNo << " evt " << evtNo << std::endl;

  if (fSelection==0){
    std::vector <unsigned int> sobe = SetOfBadEvents();
    std::vector <unsigned int> sobr = SetOfBadRuns();
    if (sobe.size() != sobr.size()) {
      throw cet::exception("SelectNeutrino.cxx: ") << " BadEvent and BadRun list must be same length. Line " <<__LINE__ << ", " << __FILE__ << "\n";
    }

    for (unsigned int ii=0; ii<sobe.size(); ++ii){
      if(sobe.at(ii)==evtNo && sobr.at(ii)==runNo)
      {
        mf::LogInfo("SelectNeutrino: ") << "\t\n Skipping run/event " << runNo <<"/"<< evtNo << " by request.\n";
        return false;
      }
    }
    return true;
  }
  else{
    for (unsigned int ii = 0; ii<fSelRuns.size(); ii++){
      if (fSelRuns[ii] == runNo && fSelSubRuns[ii] == subrunNo && fSelEvents[ii] == evtNo){
        //std::cout<<"true"<<std::endl;
        if (fSelection>0){
          return true;
        }
        else{
          return false;
        }
      }
    }
    if (fSelection>0){
      return false;
    }
    else {
      return true;
    }
  }
}

namespace filter {

  DEFINE_ART_MODULE(SelectNeutrino)

} //namespace filt
