//
// Original Author:  thomas.mccauley@cern.ch
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <iostream>
#include <string>
#include <fstream>

class ZmumuFilter : public edm::EDFilter {
public:
  explicit ZmumuFilter(const edm::ParameterSet&);
  ~ZmumuFilter();

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::InputTag muonInputTag_;

  std::ofstream csvOut_;
  std::string csvFileName_;
  double minMuonPt_;
  double maxMuonEta_;
  double maxRelIso_;

  double invariantMassMin_;
  double invariantMassMax_;
};

ZmumuFilter::ZmumuFilter(const edm::ParameterSet& iConfig)
  : muonInputTag_(iConfig.getParameter<edm::InputTag>("muonInputTag")),
    csvFileName_(iConfig.getParameter<std::string>("csvFileName")),
    minMuonPt_(iConfig.getParameter<double>("minMuonPt")),
    maxMuonEta_(iConfig.getParameter<double>("maxMuonEta")),
    maxRelIso_(iConfig.getParameter<double>("maxRelIso")),
    invariantMassMin_(iConfig.getParameter<double>("invariantMassMin")),
    invariantMassMax_(iConfig.getParameter<double>("invariantMassMax"))
{          
  csvOut_.open(csvFileName_.c_str());
}
  
ZmumuFilter::~ZmumuFilter()
{}

bool
ZmumuFilter::filter(edm::Event& event, const edm::EventSetup& setup)
{
  edm::Handle<reco::MuonCollection> muons;
  event.getByLabel(muonInputTag_, muons);

  if ( ! muons.isValid() )
  {  
    std::cerr<<"ZmumuFilter: invalid collection"<<std::endl;
    return false;
  }
  
  int charge;
  double pt, eta, phi;
  double last_pt, last_eta, last_phi;
  float iso, dxy;
  float last_iso, last_dxy;

  // Only examine if there are precisely 2 muons in the event (for simplicity)
  if ( muons->size() != 2 )
    return false;

  int last_charge = 0;
 
  for ( reco::MuonCollection::const_iterator it = muons->begin(), end = muons->end(); 
        it != end; ++it) 
  {
    // We are only looking for global muons
    if ( ! (*it).isGlobalMuon() )
      return false;

    pt = (*it).globalTrack()->pt();
    phi = (*it).globalTrack()->phi();
    eta = (*it).globalTrack()->eta();

    if ( eta > maxMuonEta_ )
      return false;

    charge = (*it).charge();
    iso = (*it).isolationR03().sumPt;
    dxy = (*it).globalTrack()->dxy();

    if ( last_charge == 0 ) // i.e. this is the first of the pair of muons
    {
      last_charge = charge;
      last_pt = pt;
      last_eta = eta;
      last_phi = phi;
      last_iso = iso;
      last_dxy = dxy;
    }
    
    else // we are on the second muon of the pair and can compare charge and calculate invariant mass
    { 
      double M = 2*last_pt*pt*(cosh(last_eta-eta)-cos(last_phi-phi));
      M = sqrt(M);

      if ( M < invariantMassMin_ || M > invariantMassMax_ )
        return false;

      csvOut_<< event.id().run() <<","<< event.id().event() <<","
             << last_pt <<","<< last_eta <<","<< last_phi <<","<< last_charge <<","
             << last_dxy <<","<< last_iso <<","
             << pt <<","<< eta <<","<< phi <<","<< charge <<","
             << dxy <<","<< iso <<std::endl;
    } 
  }
  
  return true;
}

void
ZmumuFilter::beginJob()
{
  csvOut_<<"Run,Event,pt1,eta1,phi1,Q1,dxy1,iso1,pt2,eta2,phi2,Q2,dxy2,iso2"<<std::endl;
}

void 
ZmumuFilter::endJob() {
}

DEFINE_FWK_MODULE(ZmumuFilter);
