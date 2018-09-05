#include <memory>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/EventBase.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoEgamma/EgammaTools/interface/ECALPositionCalculator.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "/mnt/t3nfs01/data01/shome/mratti/cmssw_workarea/EcalDPG/CMSSW_10_0_0/src/ECAL/ECALLowLevelComparison/interface/ECALNoiseStudy.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TVector3.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <map>

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

//
// constructors and destructor
ECALNoiseStudy::ECALNoiseStudy(const edm::ParameterSet& ps)
{
  // collections
  vertexToken_               = consumes<reco::VertexCollection>(ps.getParameter<edm::InputTag>("PVTag"));

  genParticleCollection_     = consumes<reco::GenParticleCollection>(ps.getParameter<edm::InputTag>("genParticleCollection"));

  recHitCollection_EB_       = consumes<EcalRecHitCollection> (ps.getParameter<edm::InputTag>("recHitCollection_EB"));
  recHitCollection_EE_       = consumes<EcalRecHitCollection> (ps.getParameter<edm::InputTag>("recHitCollection_EE"));

  PFrecHitCollection_        = consumes<reco::PFRecHitCollection> (ps.getParameter<edm::InputTag>("PFrecHitCollection"));

  PFclusterCollection_       = consumes<reco::PFClusterCollection> (ps.getParameter<edm::InputTag>("PFclusterCollection"));

  superClusterCollection_EB_ = consumes<reco::SuperClusterCollection> (ps.getParameter<edm::InputTag>("superClusterCollection_EB"));
  superClusterCollection_EE_ = consumes<reco::SuperClusterCollection> (ps.getParameter<edm::InputTag>("superClusterCollection_EE"));

  beamSpot_                  = consumes<reco::BeamSpot> (ps.getParameter<edm::InputTag>("beamSpot"));

  // thresholds
  ethrEB_                    = ps.getParameter<double>("ethrEB");
  ethrEE_                    = ps.getParameter<double>("ethrEE");
  scEtThrEB_                 = ps.getParameter<double>("scEtThrEB");
  scEtThrEE_                 = ps.getParameter<double>("scEtThrEE");

  naiveId_ = 0;

  // configurations for plots in delta eta bins
  std::map<TString, Double_t> start_eta;
  start_eta["EB"]=-1.5;
  start_eta["EEM"]=-3.0;
  start_eta["EEP"]=1.5;

  std::map<TString,Double_t> delta_eta;
  delta_eta["EB"]=0.1; //  can be overwritten by job option
  delta_eta["EEM"]=0.1;
  delta_eta["EEP"]=0.1;

  std::map<TString, Int_t> nBins_eta;
  nBins_eta["EB"]= (int) (3.0/delta_eta["EB"]);
  nBins_eta["EEM"]= (int) (1.5/delta_eta["EEM"]);
  nBins_eta["EEP"]= (int) (1.5/delta_eta["EEP"]);


  for (TString region : regions){
    // create the keys for the eta bins
    for (Int_t i=0; i<nBins_eta[region]; i++){
      Float_t low_edge = (start_eta[region] + i * delta_eta[region]);
      Float_t up_edge =  (start_eta[region] + (i+1) * delta_eta[region]);
      //std::cout << TString::Format("%.2f_%.2f", low_edge, up_edge) << std::endl;
      TString key = TString::Format("%.2f_%.2f", low_edge, up_edge);
      if(key.Contains(".")) key.ReplaceAll(".", "dot");
      if(key.Contains("-")) key.ReplaceAll("-", "n");
      eta_keys[region].push_back(key);
      eta_edges[region][key].first = low_edge;
      eta_edges[region][key].second = up_edge;
    }
  }

  // configurations for plots in gen energy bins
  for (TString region: regions){
    for (Int_t i=1; i<10; i++){
      TString key = TString::Format("%d_%d", i, i+1 );
      Et_keys[region].push_back(key);
      Et_edges[region][key].first = i;
      Et_edges[region][key].second = i+1;
    }
  }

  // histos
  edm::Service<TFileService> fs;
  TFileDirectory genDir = fs->mkdir( "general" );
  TFileDirectory recHitsDir = fs->mkdir( "recHits" );
  TFileDirectory PFrecHitsDir = fs->mkdir( "PFrecHits" );
  TFileDirectory PFclustersDir = fs->mkdir( "PFClusters" );
  TFileDirectory etaBinnedDir = fs->mkdir( "etaBinnedQuantities" );
  TFileDirectory EtBinnedDir = fs->mkdir( "EtBinnedQuantities" );
  TFileDirectory eventDir = fs->mkdir( "event" );


  // Vertices
  h_numberOfEvents = genDir.make<TH1D>("h_numberOfEvents","h_numberOfEvents",10,0,10);
  h_nPVs = genDir.make<TH1D>("h_nPVs","h_nPVs",80,0,80);

  // gen particles
  h_genP_pt                 = genDir.make<TH1D>("h_genP_pt", "h_genP_pt", 100, 0., 100. );
  h_genP_eta                = genDir.make<TH1D>("h_genP_eta", "h_genP_eta", 120, -3., 3. );
  h_genP_phi                = genDir.make<TH1D>("h_genP_phi", "h_genP_phi", 128, -3.2, 3.2);
  h_genP_status             = genDir.make<TH1D>("h_genP_status", "h_genP_status", 10, 0., 10.);
  h_genP_pdgid              = genDir.make<TH1D>("h_genP_pdgid", "h_genP_pdgid", 40, 0., 40.);
  h_genP_pt_EB                = genDir.make<TH1D>("h_genP_pt_EB", "h_genP_pt_EB", 100, 0., 100. );
  h_genP_pt_EEP               = genDir.make<TH1D>("h_genP_pt_EEP", "h_genP_pt_EEP", 100, 0., 100. );
  h_genP_pt_EEM               = genDir.make<TH1D>("h_genP_pt_EEM", "h_genP_pt_EEM", 100, 0., 100. );

  // Rechits, barrel
  h_recHits_EB_size          = recHitsDir.make<TH1D>("h_recHits_EB_size", "h_recHitsEB_size", 100, 500, 3500 );
  h_recHits_EB_eta           = recHitsDir.make<TH1D>("h_recHits_EB_eta","h_recHits_EB_eta",148,-1.48,1.48);
  h_recHits_EB_maxEneEta     = recHitsDir.make<TH1D>("h_recHits_EB_maxEneEta","h_recHits_EB_maxEneEta",148,-1.48,1.48);
  h_recHits_EB_energy        = recHitsDir.make<TH1D>("h_recHits_EB_energy","h_recHitsEB_energy",1000,0,20);
  h_recHits_EB_energyMax     = recHitsDir.make<TH1D>("h_recHits_EB_energyMax","h_recHitsEB_energyMax",100,0,20);
  h_recHits_EB_time          = recHitsDir.make<TH1D>("h_recHits_EB_time","h_recHits_EB_time",400,-100,100);
  h_recHits_EB_Chi2          = recHitsDir.make<TH1D>("h_recHits_EB_Chi2","h_recHits_EB_Chi2",500,0,50);
  h_recHits_EB_OutOfTimeChi2 = recHitsDir.make<TH1D>("h_recHits_EB_OutOfTimeChi2","h_recHits_EB_OutOfTimeChi2",1000,0,100);
  h_recHits_EB_E1oE4         = recHitsDir.make<TH1D>("h_recHits_EB_E1oE4","h_recHitsEB_E1oE4",148, 0, 1.48);
  h_recHits_EB_iPhiOccupancy = recHitsDir.make<TH1D>("h_recHits_EB_iPhiOccupancy","h_recHits_EB_iPhiOccupancy",360,1.,361. );
  h_recHits_EB_iEtaOccupancy = recHitsDir.make<TH1D>("h_recHits_EB_iEtaOccupancy","h_recHits_EB_iEtaOccupancy",172,-86.,86.);
  h_recHits_EB_occupancy     = recHitsDir.make<TH2D>("h_recHits_EB_occupancy","h_recHits_EB_occupancy",360,1.,361.,172,-86.,86. );
  h_recHits_EB_occupancy_gt10 = recHitsDir.make<TH2D>("h_recHits_EB_occupancy_gt10","h_recHits_EB_occupancy_gt10",360,1.,361.,172,-86.,86. );
  h_recHits_EB_occupancy_lt10 = recHitsDir.make<TH2D>("h_recHits_EB_occupancy_lt10","h_recHits_EB_occupancy_lt10",360,1.,361.,172,-86.,86. );
  h_recHits_EB_eneVSieta     = recHitsDir.make<TH2D>("h_recHits_EB_eneVSieta", "h_recHits_EB_eneVSieta", 100,0,20, 172,-86.,86.);
  h_recHits_EB_energy_spike  = recHitsDir.make<TH1D>("h_recHits_EB_energy_spike","h_recHitsEB_energy_spike",2000,0,500);

  // Rechits, barrel (with spike cleaning)
  h_recHits_EB_size_cleaned          = recHitsDir.make<TH1D>("h_recHits_EB_size_cleaned", "h_recHitsEB_size_cleaned", 1000, 0, 10000 );
  h_recHits_EB_energy_cleaned        = recHitsDir.make<TH1D>("h_recHits_EB_energy_cleaned","h_recHitsEB_energy_cleaned",11000,-50,500);
  h_recHits_EB_time_cleaned          = recHitsDir.make<TH1D>("h_recHits_EB_time_cleaned","h_recHits_EB_time_cleaned",400,-100,100);
  h_recHits_EB_Chi2_cleaned          = recHitsDir.make<TH1D>("h_recHits_EB_Chi2_cleaned","h_recHits_EB_Chi2_cleaned",1000,0,100);

  // Rechits, endcap
  h_recHits_EE_size           = recHitsDir.make<TH1D>("h_recHits_EE_size","h_recHits_EE_size",100,0,1000);
  h_recHits_EEP_size          = recHitsDir.make<TH1D>("h_recHits_EEP_size","h_recHits_EEP_size",100,0,1000);
  h_recHits_EEP_eta           = recHitsDir.make<TH1D>("h_recHits_EEP_eta","h_recHits_EEP_eta",74,1.48,3.);
  h_recHits_EEP_maxEneEta     = recHitsDir.make<TH1D>("h_recHits_EEP_maxEneEta","h_recHits_EEP_maxEneEta",74,1.48,3.);
  h_recHits_EEP_energy        = recHitsDir.make<TH1D>("h_recHits_EEP_energy","h_recHits_EEP_energy",50,0,100);
  h_recHits_EEP_energy_gt25   = recHitsDir.make<TH1D>("h_recHits_EEP_energy_gt25","h_recHits_EEP_energy_gt25",50,0,100);
  h_recHits_EEP_energyMax     = recHitsDir.make<TH1D>("h_recHits_EEP_energyMax","h_recHitsEEP_energyMax",50,0,100);
  h_recHits_EEP_time          = recHitsDir.make<TH1D>("h_recHits_EEP_time","h_recHits_EEP_time",400,-100,100);
  h_recHits_EEP_Chi2          = recHitsDir.make<TH1D>("h_recHits_EEP_Chi2","h_recHits_EEP_Chi2",500,0,50);
  h_recHits_EEP_OutOfTimeChi2 = recHitsDir.make<TH1D>("h_recHits_EEP_OutOfTimeChi2","h_recHits_EEP_OutOfTimeChi2",500,0,50);
  h_recHits_EEP_E1oE4         = recHitsDir.make<TH1D>("h_recHits_EEP_E1oE4","h_recHitsEEP_E1oE4",150, 0, 1.5);
  h_recHits_EEP_iXoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEP_iXoccupancy","h_recHits_EEP_iXoccupancy",100,0.,100.);
  h_recHits_EEP_iYoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEP_iYoccupancy","h_recHits_EEP_iYoccupancy",100,0.,100.);
  h_recHits_EEP_occupancy     = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy","h_recHits_EEP_occupancy",100,0.,100.,100,0.,100. );
  h_recHits_EEP_occupancy_etaphi = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy_etaphi","h_recHits_EEP_occupancy_etaphi",40,1.0,3.0,100,-3.2,3.2 );
  h_recHits_EEP_occupancy_gt10 = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy_gt10","h_recHits_EEP_occupancy_gt10",100,0.,100.,100,0.,100. );
  h_recHits_EEP_occupancy_lt10 = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy_lt10","h_recHits_EEP_occupancy_lt10",100,0.,100.,100,0.,100. );

  h_recHits_EEM_size          = recHitsDir.make<TH1D>("h_recHits_EEM_size","h_recHits_EEM_size",100,0,1000);
  h_recHits_EEM_eta           = recHitsDir.make<TH1D>("h_recHits_EEM_eta","h_recHits_EEM_eta",74,-3.,-1.48);
  h_recHits_EEM_maxEneEta     = recHitsDir.make<TH1D>("h_recHits_EEM_maxEneEta","h_recHits_EEM_maxEneEta",74,-3.,-1.48);
  h_recHits_EEM_energy        = recHitsDir.make<TH1D>("h_recHits_EEM_energy","h_recHits_EEM_energy",50,0,100);
  h_recHits_EEM_energy_gt25   = recHitsDir.make<TH1D>("h_recHits_EEM_energy_gt25","h_recHits_EEM_energy_gt25",50,0,100);
  h_recHits_EEM_energyMax     = recHitsDir.make<TH1D>("h_recHits_EEM_energyMax","h_recHitsEEM_energyMax",50,0,100);
  h_recHits_EEM_time          = recHitsDir.make<TH1D>("h_recHits_EEM_time","h_recHits_EEM_time",400,-100,100);
  h_recHits_EEM_Chi2          = recHitsDir.make<TH1D>("h_recHits_EEM_Chi2","h_recHits_EEM_Chi2",500,0,50);
  h_recHits_EEM_OutOfTimeChi2 = recHitsDir.make<TH1D>("h_recHits_EEM_OutOfTimeChi2","h_recHits_EEM_OutOfTimeChi2",500,0,50);
  h_recHits_EEM_E1oE4         = recHitsDir.make<TH1D>("h_recHits_EEM_E1oE4","h_recHitsEEM_E1oE4",150, 0, 1.5);
  h_recHits_EEM_iXoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEM_iXoccupancy","h_recHits_EEM_iXoccupancy",100,0.,100.);
  h_recHits_EEM_iYoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEM_iYoccupancy","h_recHits_EEM_iYoccupancy",100,0.,100.);
  h_recHits_EEM_occupancy     = recHitsDir.make<TH2D>("h_recHits_EEM_occupancy","h_recHits_EEM_occupancy",100,0.,100.,100,0.,100. );
  h_recHits_EEM_occupancy_gt10 = recHitsDir.make<TH2D>("h_recHits_EEM_occupancy_gt10","h_recHits_EEM_occupancy_gt10",100,0.,100.,100,0.,100. );
  h_recHits_EEM_occupancy_lt10 = recHitsDir.make<TH2D>("h_recHits_EEM_occupancy_lt10","h_recHits_EEM_occupancy_lt10",100,0.,100.,100,0.,100. );

  // full eta distributions
  h_recHits_eta = recHitsDir.make<TH1D>("h_recHits_eta","h_recHits_eta",300,-3.,3.);

  // energy of neighbours for maximal energy deposit in given eta bin
  h_recHits_EEP_neighbourEt_eta20 = recHitsDir.make<TH2D>("h_recHits_EEP_neighbourEt_eta20","h_recHits_EEP_neighbourEt_eta20",20,-10.,10.,20,-10.,10. );
  h_recHits_EEP_neighbourEt_eta24 = recHitsDir.make<TH2D>("h_recHits_EEP_neighbourEt_eta24","h_recHits_EEP_neighbourEt_eta24",20,-10.,10.,20,-10.,10. );
  h_recHits_EEP_sumneighbourEt_eta20 = recHitsDir.make<TH1D>("h_recHits_EEP_sumneighbourEt_eta20","h_recHits_EEP_sumneighbourEt_eta20",100,0.,10. );
  h_recHits_EEP_sumneighbourEt_eta24 = recHitsDir.make<TH1D>("h_recHits_EEP_sumneighbourEt_eta24","h_recHits_EEP_sumneighbourEt_eta24",100,0.,10. );

  // --------- PF rechits
  h_PFrecHits_EB_eta           = PFrecHitsDir.make<TH1D>("h_PFrecHits_EB_eta","h_PFrecHits_EB_eta",148,-1.48,1.48);
  h_PFrecHits_EB_energy        = PFrecHitsDir.make<TH1D>("h_PFrecHits_EB_energy","h_PFrecHitsEB_energy",100,0,20);
  h_PFrecHits_EB_time          = PFrecHitsDir.make<TH1D>("h_PFrecHits_EB_time","h_PFrecHits_EB_time",400,-100,100);
  h_PFrecHits_EB_occupancy     = PFrecHitsDir.make<TH2D>("h_PFrecHits_EB_occupancy","h_PFrecHits_EB_occupancy",360,1.,361.,172,-86.,86. );
  h_PFrecHits_EB_eneVSieta     = PFrecHitsDir.make<TH2D>("h_PFrecHits_EB_eneVSieta", "h_PFrecHits_EB_eneVSieta", 100,0,20, 172,-86.,86.);

  h_PFrecHits_EEP_eta           = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEP_eta","h_PFrecHits_EEP_eta",74,1.48,3);
  h_PFrecHits_EEP_energy        = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEP_energy","h_PFrecHits_EEP_energy",50,0,100);
  h_PFrecHits_EEP_time          = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEP_time","h_PFrecHits_EEP_time",400,-100,100);
  h_PFrecHits_EEP_occupancy     = PFrecHitsDir.make<TH2D>("h_PFrecHits_EEP_occupancy","h_PFrecHits_EEP_occupancy",100,0.,100.,100,0.,100. );

  h_PFrecHits_EEM_eta           = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEM_eta","h_PFrecHits_EEM_eta",74,-3.,-1.48);
  h_PFrecHits_EEM_energy        = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEM_energy","h_PFrecHits_EEM_energy",50,0,100);
  h_PFrecHits_EEM_time          = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEM_time","h_PFrecHits_EEM_time",400,-100,100);
  h_PFrecHits_EEM_occupancy     = PFrecHitsDir.make<TH2D>("h_PFrecHits_EEM_occupancy","h_PFrecHits_EEM_occupancy",100,0.,100.,100,0.,100. );

  // --------- PF clusters
  h_PFclusters_EB_size    = PFclustersDir.make<TH1D>("h_PFclusters_EB_size","h_PFclusters_EB_size",100,0.,100.);
  h_PFclusters_EB_nXtals  = PFclustersDir.make<TH1D>("h_PFclusters_EB_nXtals","h_PFclusters_EB_nXtals",50,0.,50.);
  h_PFclusters_EB_energy  = PFclustersDir.make<TH1D>("h_PFclusters_EB_energy","h_PFclusters_EB_energy",200,0.,10.);
  h_PFclusters_EB_et      = PFclustersDir.make<TH1D>("h_PFclusters_EB_et","h_PFclusters_EB_et",200,0.,10.);
  h_PFclusters_EB_eta     = PFclustersDir.make<TH1D>("h_PFclusters_EB_eta","h_PFclusters_EB_eta",148,-1.48,1.48);
  h_PFclusters_EB_phi     = PFclustersDir.make<TH1D>("h_PFclusters_EB_phi","h_PFclusters_EB_phi",128,-3.2,3.2);
  //h_PFclusters_EB_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_EB_eOverEtrue","h_PFclusters_EB_eOverEtrue",100,0.,2.);

  h_PFclusters_EEP_size   = PFclustersDir.make<TH1D>("h_PFclusters_EEP_size","h_PFclusters_EEP_size",100,0.,100.);
  h_PFclusters_EEP_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_EEP_nXtals","h_PFclusters_EEP_nXtals",50,0.,50.);
  h_PFclusters_EEP_energy = PFclustersDir.make<TH1D>("h_PFclusters_EEP_energy","h_PFclusters_EEP_energy",200,0.,10.);
  h_PFclusters_EEP_et     = PFclustersDir.make<TH1D>("h_PFclusters_EEP_et","h_PFclusters_EEP_et",200,0.,10.);
  h_PFclusters_EEP_eta    = PFclustersDir.make<TH1D>("h_PFclusters_EEP_eta","h_PFclusters_EEP_eta",300,-3.,3.);
  h_PFclusters_EEP_phi    = PFclustersDir.make<TH1D>("h_PFclusters_EEP_phi","h_PFclusters_EEP_phi",128,-3.2,3.2);
  //h_PFclusters_EEP_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_EEP_eOverEtrue","h_PFclusters_EEP_eOverEtrue",100,0.,2.);

  h_PFclusters_EEM_size   = PFclustersDir.make<TH1D>("h_PFclusters_EEM_size","h_PFclusters_EEM_size",100,0.,100.);
  h_PFclusters_EEM_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_EEM_nXtals","h_PFclusters_EEM_nXtals",50,0.,50.);
  h_PFclusters_EEM_energy = PFclustersDir.make<TH1D>("h_PFclusters_EEM_energy","h_PFclusters_EEM_energy",200,0.,10.);
  h_PFclusters_EEM_et     = PFclustersDir.make<TH1D>("h_PFclusters_EEM_et","h_PFclusters_EEM_et",200,0.,10.);
  h_PFclusters_EEM_eta    = PFclustersDir.make<TH1D>("h_PFclusters_EEM_eta","h_PFclusters_EEM_eta",300,-3.,3.);
  h_PFclusters_EEM_phi    = PFclustersDir.make<TH1D>("h_PFclusters_EEM_phi","h_PFclusters_EEM_phi",128,-3.2,3.2);
  //h_PFclusters_EEM_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_EEM_eOverEtrue","h_PFclusters_EEM_eOverEtrue",100,0.,2.);

  h_PFclusters_eta        = PFclustersDir.make<TH1D>("h_PFclusters_eta",   "h_PFclusters_eta",   300,-3.,3.);

  h_PFclusters_deltaR_gen = PFclustersDir.make<TH1D>("h_PFclusters_deltaR_gen", "h_PFclusters_deltaR_gen", 1000., 0., 10.);
  h_PFclusters_deltaR_gen_EB = PFclustersDir.make<TH1D>("h_PFclusters_deltaR_gen_EB", "h_PFclusters_deltaR_gen_EB", 1000., 0., 10.);
  h_PFclusters_deltaR_gen_EEP = PFclustersDir.make<TH1D>("h_PFclusters_deltaR_gen_EEP", "h_PFclusters_deltaR_gen_EEP", 1000., 0., 10.);
  h_PFclusters_deltaR_gen_EEM = PFclustersDir.make<TH1D>("h_PFclusters_deltaR_gen_EEM", "h_PFclusters_deltaR_gen_EEM", 1000., 0., 10.);

  h_PFclusters500_deltaR_gen = PFclustersDir.make<TH1D>("h_PFclusters500_deltaR_gen", "h_PFclusters500_deltaR_gen", 1000., 0., 10.);
  h_PFclusters500_deltaR_gen_EB = PFclustersDir.make<TH1D>("h_PFclusters500_deltaR_gen_EB", "h_PFclusters500_deltaR_gen_EB", 1000., 0., 10.);
  h_PFclusters500_deltaR_gen_EEP = PFclustersDir.make<TH1D>("h_PFclusters500_deltaR_gen_EEP", "h_PFclusters500_deltaR_gen_EEP", 1000., 0., 10.);
  h_PFclusters500_deltaR_gen_EEM = PFclustersDir.make<TH1D>("h_PFclusters500_deltaR_gen_EEM", "h_PFclusters500_deltaR_gen_EEM", 1000., 0., 10.);

  h_PFclusters1000_deltaR_gen = PFclustersDir.make<TH1D>("h_PFclusters1000_deltaR_gen", "h_PFclusters1000_deltaR_gen", 1000., 0., 10.);
  h_PFclusters1000_deltaR_gen_EB = PFclustersDir.make<TH1D>("h_PFclusters1000_deltaR_gen_EB", "h_PFclusters1000_deltaR_gen_EB", 1000., 0., 10.);
  h_PFclusters1000_deltaR_gen_EEP = PFclustersDir.make<TH1D>("h_PFclusters1000_deltaR_gen_EEP", "h_PFclusters1000_deltaR_gen_EEP", 1000., 0., 10.);
  h_PFclusters1000_deltaR_gen_EEM = PFclustersDir.make<TH1D>("h_PFclusters1000_deltaR_gen_EEM", "h_PFclusters1000_deltaR_gen_EEM", 1000., 0., 10.);

  // --------- PF clusters - gen matched
  h_PFclusters_genMatched_EB_size    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_size","h_PFclusters_genMatched_EB_size",100,0.,100.);
  h_PFclusters_genMatched_EB_nXtals  = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_nXtals","h_PFclusters_genMatched_EB_nXtals",50,0.,50.);
  h_PFclusters_genMatched_EB_energy  = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_energy","h_PFclusters_genMatched_EB_energy",200,0.,10.);
  h_PFclusters_genMatched_EB_et  = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_et","h_PFclusters_genMatched_EB_et",200,0.,10.);
  h_PFclusters_genMatched_EB_eta     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_eta","h_PFclusters_genMatched_EB_eta",148,-1.48,1.48);
  h_PFclusters_genMatched_EB_phi     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_phi","h_PFclusters_genMatched_EB_phi",128,-3.2,3.2);
  h_PFclusters_genMatched_EB_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_eOverEtrue","h_PFclusters_genMatched_EB_eOverEtrue",100,0.,2.);

  h_PFclusters_genMatched_EEP_size   = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_size","h_PFclusters_genMatched_EEP_size",100,0.,100.);
  h_PFclusters_genMatched_EEP_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_nXtals","h_PFclusters_genMatched_EEP_nXtals",50,0.,50.);
  h_PFclusters_genMatched_EEP_energy = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_energy","h_PFclusters_genMatched_EEP_energy",200,0.,10.);
  h_PFclusters_genMatched_EEP_et     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_et","h_PFclusters_genMatched_EEP_et",200,0.,10.);
  h_PFclusters_genMatched_EEP_eta    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_eta","h_PFclusters_genMatched_EEP_eta",300,-3.,3.);
  h_PFclusters_genMatched_EEP_phi    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_phi","h_PFclusters_genMatched_EEP_phi",128,-3.2,3.2);
  h_PFclusters_genMatched_EEP_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_eOverEtrue","h_PFclusters_genMatched_EEP_eOverEtrue",100,0.,2.);

  h_PFclusters_genMatched_EEM_size   = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_size","h_PFclusters_genMatched_EEM_size",100,0.,100.);
  h_PFclusters_genMatched_EEM_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_nXtals","h_PFclusters_genMatched_EEM_nXtals",50,0.,50.);
  h_PFclusters_genMatched_EEM_energy = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_energy","h_PFclusters_genMatched_EEM_energy",200,0.,10.);
  h_PFclusters_genMatched_EEM_et     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_et","h_PFclusters_genMatched_EEM_et",200,0.,10.);
  h_PFclusters_genMatched_EEM_eta    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_eta","h_PFclusters_genMatched_EEM_eta",300,-3.,3.);
  h_PFclusters_genMatched_EEM_phi    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_phi","h_PFclusters_genMatched_EEM_phi",128,-3.2,3.2);
  h_PFclusters_genMatched_EEM_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_eOverEtrue","h_PFclusters_genMatched_EEM_eOverEtrue",100,0.,2.);





  // Super Clusters ----------------------------------------------
  // ... barrel
  h_superClusters_EB_size      = fs->make<TH1D>("h_superClusters_EB_size","h_superClusters_EB_size",15,0.,15.);
  h_superClusters_EB_nXtals    = fs->make<TH1D>("h_superClusters_EB_nXtals","h_superClusters_EB_nXtals",30,0.,150.);
  h_superClusters_EB_nBC       = fs->make<TH1D>("h_superClusters_EB_nBC","h_superClusters_EB_nBC",20,0.,20.);
  h_superClusters_EB_energy    = fs->make<TH1D>("h_superClusters_EB_energy","h_superClusters_EB_energy",50,0.,200.);
  h_superClusters_EB_rawEnergy = fs->make<TH1D>("h_superClusters_EB_rawEnergy","h_superClusters_EB_rawEnergy",50,0.,200.);
  h_superClusters_EB_et        = fs->make<TH1D>("h_superClusters_EB_et","h_superClusters_EB_et",50,0.,200.);

  // ... endcap
  h_superClusters_EEP_size   = fs->make<TH1D>("h_superClusters_EEP_size","h_superClusters_EEP_size",15,0.,15.);
  h_superClusters_EEP_nXtals = fs->make<TH1D>("h_superClusters_EEP_nXtals","h_superClusters_EEP_nXtals",30,0.,60.);
  h_superClusters_EEP_nBC    = fs->make<TH1D>("h_superClusters_EEP_nBC","h_superClusters_EEP_nBC",20,0.,20.);
  h_superClusters_EEP_energy = fs->make<TH1D>("h_superClusters_EEP_energy","h_superClusters_EEP_energy",50,0.,200.);
  h_superClusters_EEP_rawEnergy = fs->make<TH1D>("h_superClusters_EEP_rawEnergy","h_superClusters_EEP_rawEnergy",50,0.,200.);
  h_superClusters_EEP_et     = fs->make<TH1D>("h_superClusters_EEP_et","h_superClusters_EEP_et",50,0.,200.);

  h_superClusters_EEM_size   = fs->make<TH1D>("h_superClusters_EEM_size","h_superClusters_EEM_size",15,0.,15.);
  h_superClusters_EEM_nXtals = fs->make<TH1D>("h_superClusters_EEM_nXtals","h_superClusters_EEM_nXtals",30,0.,60.);
  h_superClusters_EEM_nBC    = fs->make<TH1D>("h_superClusters_EEM_nBC","h_superClusters_EEM_nBC",20,0.,20.);
  h_superClusters_EEM_energy = fs->make<TH1D>("h_superClusters_EEM_energy","h_superClusters_EEM_energy",50,0.,200.);
  h_superClusters_EEM_rawEnergy = fs->make<TH1D>("h_superClusters_EEM_rawEnergy","h_superClusters_EEM_rawEnergy",50,0.,200.);
  h_superClusters_EEM_et     = fs->make<TH1D>("h_superClusters_EEM_et","h_superClusters_EEM_et",50,0.,200.);

  h_superClusters_eta        = fs->make<TH1D>("h_superClusters_eta","h_superClusters_eta",      300,-3.,3.);
  h_superClusters_EB_eta     = fs->make<TH1D>("h_superClusters_EB_eta","h_superClusters_EB_eta",148,-1.48,1.48);
  h_superClusters_EE_eta     = fs->make<TH1D>("h_superClusters_EE_eta","h_superClusters_EE_eta",300,-3.,3.);

  // check golden fraction
  h_superClusters_nBCvsEta = fs->make<TH2D>("h_superClusters_nBCvsEta","h_superClusters_nBCvsEta",300,-3.,3.,20,0.,20.);
  h_superClusters_nBC_0to1     = fs->make<TH1D>("h_superClusters_nBC_0to1",    "h_superClusters_nBC_0to1",    20,0.,20.);
  h_superClusters_nBC_1to1d5   = fs->make<TH1D>("h_superClusters_nBC_1to1d5",  "h_superClusters_nBC_1to1d5",  20,0.,20.);
  h_superClusters_nBC_1d5to1d8 = fs->make<TH1D>("h_superClusters_nBC_1d5to1d8","h_superClusters_nBC_1d5to1d8",20,0.,20.);
  h_superClusters_nBC_1d8to2d1 = fs->make<TH1D>("h_superClusters_nBC_1d8to2d1","h_superClusters_nBC_1d8to2d1",20,0.,20.);
  h_superClusters_nBC_2d1to2d5 = fs->make<TH1D>("h_superClusters_nBC_2d1to2d5","h_superClusters_nBC_2d1to2d5",20,0.,20.);
  h_superClusters_nBC_2d5to3   = fs->make<TH1D>("h_superClusters_nBC_2d5to3",  "h_superClusters_nBC_2d5to3",  20,0.,20.);


  // --------- Rechits vs eta
  for (TString region : regions){
    for (TString key : eta_keys[region]){
      TString histo_name = "h_RecHits_" + region + "_energy_" + key;
      h_recHits_energy_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,1000,0,10);
      histo_name = "h_RecHits_" + region + "_et_" + key;
      h_recHits_et_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,1000,0,10);
    }
  }

  // --------- PFRechits vs eta
  for (TString region : regions){
    for (TString key : eta_keys[region]){
      TString histo_name = "h_PfRecHits_" + region + "_energy_" + key;
      h_PFrecHits_energy_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,1000,0,10);
    }
  }

  // --------- E (PF clusters) over E true binned in eta
  for (TString region : regions){
    for (TString key : eta_keys[region]){
      TString histo_name = "h_PFclusters_genMatched_" + region + "_eOverEtrue_" + key;
      h_PFclusters_genMatched_eOverEtrue_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,100,0.,2.);
    }
  }

  // --------- E (PF clusters) over E true binned in et
  for (TString region : regions){
    for (TString key : Et_keys[region]){
      TString histo_name = "h_PFclusters_genMatched_" + region + "_eOverEtrue_" + key;
      h_PFclusters_genMatched_eOverEtrue_EtBinned[region][key] = EtBinnedDir.make<TH1F>(histo_name,histo_name,100,0.,2.);
      histo_name = "h_genP_" + region + "_nEvts_" + key;
      h_genP_nEvts_EtBinned[region][key] = EtBinnedDir.make<TH1F>(histo_name,histo_name,1,0.,1.);
     }
  }

  // --------- event by event diagnostic plots
  for (int i=0; i<100; i++){

    TString histo_name = "h_genP_etaVsPhi_"  + TString::Format("%d", i);
    h_genP_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_RecHits_etaVsPhi_" + TString::Format("%d", i);
    h_recHits_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_PFclusters_etaVsPhi_" + TString::Format("%d", i);
    h_PFclusters_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_PFclusters_genMatched_etaVsPhi_" + TString::Format("%d", i);
    h_PFclusters_genMatched_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
  }


}



ECALNoiseStudy::~ECALNoiseStudy() {}

// ------------ method called to for each event  ------------
void ECALNoiseStudy::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  TH1::StatOverflows(kTRUE);

  // Get vertices
  edm::Handle<reco::VertexCollection> vtx_h;
  ev.getByToken(vertexToken_, vtx_h);
  if ( ! vtx_h.isValid() ) std::cout << "ECALNoiseStudy: vtx collection not found" << std::endl;

  // Get the BS position
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  ev.getByToken(beamSpot_,recoBeamSpotHandle);
  if ( ! recoBeamSpotHandle.isValid() ) std::cout << "ECALNoiseStudy: BS collection not found" << std::endl;
  const reco::BeamSpot::Point& BSPosition = recoBeamSpotHandle->position();

  int nvertices(0);
  for(reco::VertexCollection::const_iterator it = vtx_h->begin(); it !=vtx_h->end() ; ++it){
    if(it->isValid() && it->ndof() > 4. && it->position().Rho() < 2. && fabs(it->position().Z() - BSPosition.z()) < 24){
      nvertices++;
    }
  }
  h_nPVs->Fill(nvertices);

  // Gen particles
  edm::Handle<reco::GenParticleCollection> genParticles_handle;
  ev.getByToken( genParticleCollection_, genParticles_handle );
  if ( ! genParticles_handle.isValid() ) std::cout << "ECALNoiseStudy::analyze --> genParticles not found" << std::endl;
  const reco::GenParticleCollection* genParticles = genParticles_handle.product ();

  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin (); genParticle != genParticles->end () ;++genParticle) {
    //std::cout << naiveId_ << " pdgid=" << genParticle->pdgId() << " status=" <<  genParticle->status() << " pt=" << genParticle->pt() << std::endl;
    h_genP_pt->Fill(genParticle->pt());
    h_genP_eta->Fill(genParticle->eta());
    h_genP_phi->Fill(genParticle->phi());
    h_genP_status->Fill(genParticle->status());
    h_genP_pdgid->Fill(genParticle->pdgId());
    if(naiveId_<100) h_genP_etaVsPhi.at(naiveId_)->Fill(genParticle->eta(), genParticle->phi(), genParticle->energy());

    // FIXME: this is a poor way of doing this, please automatise
    if (fabs(genParticle->eta()) < 1.48) {
      h_genP_pt_EB->Fill(genParticle->pt());
      for(TString key: Et_keys["EB"]){
        if (genParticle->pt() >= Et_edges["EB"][key].first && genParticle->pt() < Et_edges["EB"][key].second){
          h_genP_nEvts_EtBinned["EB"][key]->Fill(0.5);
          break; // when you found it, do not loop over the other keys
        }
      }
    }
    else if (genParticle->eta()  > 1.48) {
      h_genP_pt_EEP->Fill(genParticle->pt());
      for(TString key: Et_keys["EEP"]){
        // FIXME: in reality these edges are the same between EB and EE... but ok
        if (genParticle->pt() >= Et_edges["EEP"][key].first && genParticle->pt() < Et_edges["EEP"][key].second){
          h_genP_nEvts_EtBinned["EEP"][key]->Fill(0.5);
          break; // when you found it, do not loop over the other keys
        }
      }
    }
    else {
      h_genP_pt_EEM->Fill(genParticle->pt());
      for(TString key: Et_keys["EEM"]){
        if (genParticle->pt() >= Et_edges["EEM"][key].first && genParticle->pt() < Et_edges["EEM"][key].second){
          h_genP_nEvts_EtBinned["EEM"][key]->Fill(0.5);
          break; // when you found it, do not loop over the other keys
        }
      }
    }

  } // end loop gen particles

  // calo geometry
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  const CaloGeometry *geometry = pGeometry.product();


  // --- REC HITS, barrel -------------------------------------------------------------------------------------
  edm::Handle<EcalRecHitCollection> recHitsEB;
  ev.getByToken( recHitCollection_EB_, recHitsEB );
  if ( ! recHitsEB.isValid() ) std::cout << "ECALNoiseStudy::analyze --> recHitsEB not found" << std::endl;
  const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product ();

  float maxERecHitEB_ene = -999.;
  float maxERecHitEB_eta = -999.;
  int sizeEB_cleaned = 0;
  for ( EcalRecHitCollection::const_iterator itr = theBarrelEcalRecHits->begin (); itr != theBarrelEcalRecHits->end () ;++itr) {

    EBDetId ebid( itr -> id() );
    GlobalPoint mycell = geometry -> getPosition(DetId(itr->id()));

    // max energy rec hit
    if (itr -> energy() > maxERecHitEB_ene ) {
      maxERecHitEB_ene = itr -> energy() ;
      maxERecHitEB_eta = mycell.eta();
    }

    // distributions for all RH above energy thresholds
    if ( itr -> energy() > ethrEB_ ){
      h_recHits_EB_energy        -> Fill( itr -> energy() );
      //std::cout << "energy=" << itr->energy() << " et=" << itr->pt() << std::endl;
      h_recHits_EB_time          -> Fill( itr -> time() );
      h_recHits_EB_Chi2          -> Fill( itr -> chi2() );
      h_recHits_EB_eneVSieta     -> Fill( itr->energy() , ebid.ieta() );
      h_recHits_EB_occupancy     -> Fill( ebid.iphi() , ebid.ieta() );
      if (itr->energy() > 10) h_recHits_EB_occupancy_gt10 -> Fill( ebid.iphi() , ebid.ieta() );
      if (itr->energy() < 10) h_recHits_EB_occupancy_lt10 -> Fill( ebid.iphi() , ebid.ieta() );
      h_recHits_EB_iPhiOccupancy -> Fill( ebid.iphi() );
      h_recHits_EB_iEtaOccupancy -> Fill( ebid.ieta() );
      h_recHits_EB_eta           -> Fill( mycell.eta() );
      h_recHits_eta              -> Fill( mycell.eta() );
      if (naiveId_<100) h_recHits_etaVsPhi.at(naiveId_)-> Fill( mycell.eta(), mycell.phi(), itr->energy() );

      for(TString key : eta_keys["EB"]){
        //std::cout << key << "  " << itr->energy() << std::endl;
        if( mycell.eta() >= eta_edges["EB"][key].first && mycell.eta() < eta_edges["EB"][key].second){
          h_recHits_energy_etaBinned["EB"][key]->Fill(itr -> energy());
          Double_t et = itr -> energy() *  TMath::Sin(2*TMath::ATan(TMath::Exp(-mycell.eta())));
          h_recHits_et_etaBinned["EB"][key]->Fill(et);
          break; // when you found it, exit
        }
      }

    }

    // spikes
    double et = itr -> energy()*mycell.perp()/mycell.mag();
    float R4 = EcalTools::swissCross( ebid, *theBarrelEcalRecHits, 0. );
    if ( itr -> energy() > 3. && abs(ebid.ieta())!=85 )  h_recHits_EB_E1oE4-> Fill( R4 );
    if ( R4 > 0.95 ) h_recHits_EB_energy_spike -> Fill( itr -> energy() );

    // spike cleaning
    if ( R4 > 0.95 && et > 3.) continue;
    sizeEB_cleaned++;

    if ( itr -> energy() > ethrEB_ ){
      h_recHits_EB_energy_cleaned -> Fill( itr -> energy() );
      h_recHits_EB_time_cleaned   -> Fill( itr -> time() );
      h_recHits_EB_Chi2_cleaned   -> Fill( itr -> chi2() );
    }
  }

  h_recHits_EB_size          -> Fill( recHitsEB->size() );
  h_recHits_EB_size_cleaned  -> Fill( sizeEB_cleaned );
  h_recHits_EB_energyMax     -> Fill( maxERecHitEB_ene );
  h_recHits_EB_maxEneEta     -> Fill( maxERecHitEB_eta );




  // ----------------
  // ... endcap
  edm::Handle<EcalRecHitCollection> recHitsEE;
  ev.getByToken( recHitCollection_EE_, recHitsEE );
  if ( ! recHitsEE.isValid() ) std::cout << "ECALNoiseStudy::analyze --> recHitsEE not found" << std::endl;
  const EcalRecHitCollection* theEndcapEcalRecHits = recHitsEE.product () ;

  int nHitsEEP = 0;
  int nHitsEEM = 0;
  float maxERecHitEEP_ene = -999.;
  float maxERecHitEEM_ene = -999.;
  float maxERecHitEEP_eta = -999.;
  float maxERecHitEEM_eta = -999.;

  float maxERecHit20_ene = -999;
  float maxERecHit20_ix = -999;
  float maxERecHit20_iy = -999;
  float maxERecHit24_ene = -999;
  float maxERecHit24_ix = -999;
  float maxERecHit24_iy = -999;

  for ( EcalRecHitCollection::const_iterator itr = theEndcapEcalRecHits->begin (); itr != theEndcapEcalRecHits->end () ; ++itr) {

    EEDetId eeid( itr -> id() );
    GlobalPoint mycell = geometry->getPosition(itr->detid());

    // EE+
    if ( eeid.zside() > 0 ){

	    nHitsEEP++;

	    // max energy rec hit
	    if (itr -> energy() > maxERecHitEEP_ene){
	      maxERecHitEEP_ene = itr -> energy() ;
	      maxERecHitEEP_eta = mycell.eta();
	    }

      // max energy rec hit in eta 2.0-2.1
      if (itr->energy()>maxERecHit20_ene and mycell.eta()>=2.0 and  mycell.eta()<2.1 ){
        maxERecHit20_ene = itr->energy();
        maxERecHit20_ix = eeid.ix();
        maxERecHit20_iy = eeid.iy();
      } // max energy rec hit in eta 2.4-2.5
      if (itr->energy()>maxERecHit24_ene and mycell.eta()>=2.4 and  mycell.eta()<2.5 ){
        maxERecHit24_ene = itr->energy();
        maxERecHit24_ix = eeid.ix();
        maxERecHit24_iy = eeid.iy();
      }

	    // distributions for all RH above energy thresholds
    	if ( itr -> energy() > ethrEE_ ){
    	  h_recHits_EEP_energy        -> Fill( itr -> energy() );
    	  if(fabs(mycell.eta())>2.5) h_recHits_EEP_energy_gt25   -> Fill( itr -> energy() );
    	  h_recHits_EEP_time          -> Fill( itr -> time() );
    	  h_recHits_EEP_Chi2          -> Fill( itr -> chi2() );
    	  h_recHits_EEP_occupancy     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5 );
        h_recHits_EEP_occupancy_etaphi -> Fill (mycell.eta(), mycell.phi());
    	  if (itr->energy() >10) h_recHits_EEP_occupancy_gt10     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5 );
    	  if (itr->energy() <10) h_recHits_EEP_occupancy_lt10     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5 );
    	  h_recHits_EEP_iXoccupancy   -> Fill( eeid.ix() - 0.5 );
    	  h_recHits_EEP_iYoccupancy   -> Fill( eeid.iy() - 0.5 );
    	  h_recHits_EEP_eta           -> Fill( mycell.eta() );
    	  h_recHits_eta               -> Fill( mycell.eta() );
        if (naiveId_<100) h_recHits_etaVsPhi.at(naiveId_)-> Fill( mycell.eta(), mycell.phi(), itr->energy() );

        for(TString key : eta_keys["EEP"]){
          //std::cout << key << "  " << itr->energy() << std::endl;
          if( mycell.eta() >= eta_edges["EEP"][key].first && mycell.eta() < eta_edges["EEP"][key].second){
             //std::cout << "Filling " << std::endl;
             h_recHits_energy_etaBinned["EEP"][key]->Fill(itr -> energy());
             Double_t et = itr -> energy() *  TMath::Sin(2*TMath::ATan(TMath::Exp(-mycell.eta())));
             h_recHits_et_etaBinned["EEP"][key]->Fill(et);
             break; // when you found it, exit
          }
        }
    	}
    }

    // EE-
    if ( eeid.zside() < 0 ){

	    nHitsEEM++;

	    // max energy rec hit
	    if (itr -> energy() > maxERecHitEEM_ene){
	      maxERecHitEEM_ene = itr -> energy() ;
	      maxERecHitEEM_eta = mycell.eta();
	    }

          // distributions for all RH above energy threshold
      if (  itr -> energy() > ethrEE_ ) {
        h_recHits_EEM_energy        -> Fill( itr -> energy() );
        if(fabs(mycell.eta())>2.5) h_recHits_EEM_energy_gt25   -> Fill( itr -> energy() );
        h_recHits_EEM_time          -> Fill( itr -> time() );
        h_recHits_EEM_Chi2          -> Fill( itr -> chi2() );
        h_recHits_EEM_occupancy     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
        if(itr->energy() >10) h_recHits_EEM_occupancy_gt10     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
        if(itr->energy() <10) h_recHits_EEM_occupancy_lt10     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
        h_recHits_EEM_iXoccupancy   -> Fill( eeid.ix() - 0.5 );
        h_recHits_EEM_iYoccupancy   -> Fill( eeid.iy() - 0.5 );
        h_recHits_EEM_eta           -> Fill( mycell.eta() );
        h_recHits_eta               -> Fill( mycell.eta() );
        if (naiveId_<100) h_recHits_etaVsPhi.at(naiveId_)-> Fill( mycell.eta(), mycell.phi(), itr->energy() );

        for(TString key : eta_keys["EEM"]){
          //std::cout << key << "  " << itr->energy() << std::endl;
          if( mycell.eta() >= eta_edges["EEM"][key].first && mycell.eta() < eta_edges["EEM"][key].second){
             //std::cout << "Filling " << std::endl;
             h_recHits_energy_etaBinned["EEM"][key]->Fill(itr -> energy());
             Double_t et = itr -> energy() *  TMath::Sin(2*TMath::ATan(TMath::Exp(-mycell.eta())));
             h_recHits_et_etaBinned["EEM"][key]->Fill(et);
             break; // when you found it, exit
           }
         }
       }
     }
   } // end loop over EE rec hits


  // size
  h_recHits_EE_size    -> Fill( recHitsEE->size() );
  h_recHits_EEP_size   -> Fill( nHitsEEP );
  h_recHits_EEM_size   -> Fill( nHitsEEM );

  h_recHits_EEP_energyMax -> Fill(maxERecHitEEP_ene);
  h_recHits_EEM_energyMax -> Fill(maxERecHitEEM_ene);
  h_recHits_EEP_maxEneEta -> Fill(maxERecHitEEP_eta);
  h_recHits_EEM_maxEneEta -> Fill(maxERecHitEEM_eta);



  // make another loop over rechits and have a look at the ones in the vicinity of the rechit with maximal energy
  float sumEnergy_20=0;
  float sumEnergy_24=0;

  for ( EcalRecHitCollection::const_iterator itr = theEndcapEcalRecHits->begin (); itr != theEndcapEcalRecHits->end () ; ++itr) {
    EEDetId eeid( itr -> id() );
    GlobalPoint mycell = geometry->getPosition(itr->detid());

    float delta20_ix = maxERecHit20_ix - eeid.ix();
    float delta20_iy = maxERecHit20_iy - eeid.iy();

    float delta24_ix = maxERecHit24_ix - eeid.ix();
    float delta24_iy = maxERecHit24_iy - eeid.iy();

    Double_t et = itr -> energy() *  TMath::Sin(2*TMath::ATan(TMath::Exp(-mycell.eta())));


    if(fabs(delta20_ix)<=10 && fabs(delta20_iy)<=10){
      if(mycell.eta()>=2.0 and mycell.eta()<2.1){
        h_recHits_EEP_neighbourEt_eta20->Fill(delta20_ix, delta20_iy, et);
        sumEnergy_20 += et;
      }
    }
    if(fabs(delta24_ix)<=10 && fabs(delta24_iy)<=10){
      if(mycell.eta()>=2.4 and mycell.eta()<2.5){
        h_recHits_EEP_neighbourEt_eta24->Fill(delta24_ix, delta24_iy, et);
        sumEnergy_24 += et;
      }
    } // end condition of vicinity






  } // end loop over rechits

  h_recHits_EEP_sumneighbourEt_eta20->Fill(sumEnergy_20);
  h_recHits_EEP_sumneighbourEt_eta24->Fill(sumEnergy_24);

  // -------------------------------------------------------------------------
  // --- PF rechits ---
  edm::Handle<reco::PFRecHitCollection> PFrecHits_handle;
  ev.getByToken( PFrecHitCollection_, PFrecHits_handle );
  if ( ! PFrecHits_handle.isValid() ) std::cout << "ECALNoiseStudy::analyze --> PFrecHits not found" << std::endl;
  const reco::PFRecHitCollection* PFrecHits = PFrecHits_handle.product ();

  for (reco::PFRecHitCollection::const_iterator itr = PFrecHits->begin(); itr != PFrecHits->end(); itr++ ) {

    GlobalPoint mycell = geometry -> getPosition(DetId(itr->detId()));

    //std::cout << "id=" << itr->detId() << " eta="  << mycell.eta() << " energy=" << itr->energy() <<  " ishigher=" << (itr->detId()> 872420480) << " isbarrel=" << (fabs(mycell.eta())<1.45)<< std::endl;
    EBDetId ebid( itr -> detId() );

    if ( itr -> energy() > ethrEB_ ){

      // BARREL
      if ( fabs(mycell.eta())<=1.48 ) {
        //std::cout << "found pfrechit in barrel" << std::endl;

        h_PFrecHits_EB_time          -> Fill( itr -> time() );
        h_PFrecHits_EB_eneVSieta     -> Fill( itr->energy() , ebid.ieta() );
        h_PFrecHits_EB_occupancy     -> Fill( ebid.iphi() , ebid.ieta() );
        h_PFrecHits_EB_eta           -> Fill( mycell.eta() );
        h_PFrecHits_EB_energy        -> Fill( itr->energy() );

        for(TString key : eta_keys["EB"]){
          if( mycell.eta() >= eta_edges["EB"][key].first && mycell.eta() < eta_edges["EB"][key].second){
            h_PFrecHits_energy_etaBinned["EB"][key]->Fill(itr -> energy());
            break; // when you found it, exit
          }
        }
      }

      // End-caps
      else{
        //std::cout << "found pfrechit in endcap" << std::endl;
        EEDetId eeid( itr -> detId() );
        // TODO make sure that this is correct !
        // EEP
        if (mycell.eta() > 0){

          h_PFrecHits_EEP_time          -> Fill( itr -> time() );
          h_PFrecHits_EEP_occupancy     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
          h_PFrecHits_EEP_eta           -> Fill( mycell.eta() );
          h_PFrecHits_EEP_energy        -> Fill( itr->energy() );

          for(TString key : eta_keys["EEP"]){
            if( mycell.eta() >= eta_edges["EEP"][key].first && mycell.eta() < eta_edges["EEP"][key].second){
              h_PFrecHits_energy_etaBinned["EEP"][key]->Fill(itr -> energy());
              break; // when you found it, exit
            }
          }
        }
        // EEM
        else {

          h_PFrecHits_EEM_time          -> Fill( itr -> time() );
          h_PFrecHits_EEM_occupancy     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
          h_PFrecHits_EEM_eta           -> Fill( mycell.eta() );
          h_PFrecHits_EEM_energy        -> Fill( itr->energy() );

          for(TString key : eta_keys["EEM"]){
            //std::cout << key << "  " << itr->energy() << std::endl;
            if( mycell.eta() >= eta_edges["EEM"][key].first && mycell.eta() < eta_edges["EEM"][key].second){
              h_PFrecHits_energy_etaBinned["EEM"][key]->Fill(itr -> energy());
              break; // when you found it, exit
            }
          }
        }
      } // end end-caps
    } // end if threshold
  } // end loop over pfrechits

  // -------------------------------------------------------------------------
  // --- PF CLUSTERS ---
  edm::Handle<reco::PFClusterCollection> PFclusters_handle;
  ev.getByToken( PFclusterCollection_, PFclusters_handle );
  if ( ! PFclusters_handle.isValid() ) std::cout << "ECALNoiseStudy::analyze --> PFclusters not found" << std::endl;
  const reco::PFClusterCollection* PFclusters = PFclusters_handle.product ();

  int size_PFclusters = 0; // how many pf clusters in EB in the event
  int size_PFclusters_EB = 0; // how many pf clusters in EB in the event
  int size_PFclusters_EEP = 0; // ""
  int size_PFclusters_EEM = 0; // ""

  int size_PFclusters_genMatched = 0; // how many gen matched pf clusters in EB in the event
  int size_PFclusters_genMatched_EB = 0; // how many gen matched pf clusters in EB in the event
  int size_PFclusters_genMatched_EEP = 0; // ""
  int size_PFclusters_genMatched_EEM = 0; // ""

  for (reco::PFClusterCollection::const_iterator itr = PFclusters->begin(); itr != PFclusters->end(); itr++ ) {
    size_PFclusters++;
    //std::cout << size_PFclusters << "  pt=" << itr->pt() << std::endl;
    h_PFclusters_eta -> Fill( itr->eta() );

    // barrel
    if (fabs(itr->eta()) < 1.48){
      size_PFclusters_EB++;
      h_PFclusters_EB_nXtals -> Fill( (*itr).hitsAndFractions().size() );
      h_PFclusters_EB_energy -> Fill( itr->energy() );
      h_PFclusters_EB_et     -> Fill( itr->pt() );
      h_PFclusters_EB_eta    -> Fill( itr->eta() );
      h_PFclusters_EB_phi    -> Fill( itr->phi() );
    }
    // end-caps
    else{
      // EEP
      if (itr->eta() > 0){
        size_PFclusters_EEP++;
        h_PFclusters_EEP_nXtals -> Fill( (*itr).hitsAndFractions().size() );
        h_PFclusters_EEP_energy -> Fill( itr->energy() );
        h_PFclusters_EEP_et     -> Fill( itr->pt() );
        h_PFclusters_EEP_eta    -> Fill( itr->eta() );
        h_PFclusters_EEP_phi    -> Fill( itr->phi() );
      }
      // EEM
      else {
        size_PFclusters_EEM++;
        h_PFclusters_EEM_nXtals -> Fill( (*itr).hitsAndFractions().size() );
        h_PFclusters_EEM_energy -> Fill( itr->energy() );
        h_PFclusters_EEM_et     -> Fill( itr->pt() );
        h_PFclusters_EEM_eta    -> Fill( itr->eta() );
        h_PFclusters_EEM_phi    -> Fill( itr->phi() );
      }
    } // end eta conditions

    // cluster matching with gen particle
    // delta R
    for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin (); genParticle != genParticles->end () ;++genParticle) {

      // matchin only with photons of status 1
      if(genParticle->pdgId()!=22 or genParticle->status()!= 1) continue;
      //if(genParticle->pt()<9) continue; //  only consider clusters matched to 9-10 GeV photons

      double deltaPhi = TVector2::Phi_mpi_pi( genParticle->phi() - itr->phi());
      double deltaEta = genParticle->eta() - itr->eta();
      double deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

      h_PFclusters_deltaR_gen->Fill(deltaR);
      if (fabs(genParticle->eta()) < 1.48) h_PFclusters_deltaR_gen_EB->Fill(deltaR);
      else if (genParticle->eta()  > 1.48) h_PFclusters_deltaR_gen_EEP->Fill(deltaR);
      else                                 h_PFclusters_deltaR_gen_EEM->Fill(deltaR);

      if (itr->pt()>0.5) { // 500 MeV
        h_PFclusters500_deltaR_gen->Fill(deltaR);
        if (fabs(genParticle->eta()) < 1.48) h_PFclusters500_deltaR_gen_EB->Fill(deltaR);
        else if (genParticle->eta()  > 1.48) h_PFclusters500_deltaR_gen_EEP->Fill(deltaR);
        else                                 h_PFclusters500_deltaR_gen_EEM->Fill(deltaR);
      }
      if (itr->pt()>1) { // 1 GeV
        h_PFclusters1000_deltaR_gen->Fill(deltaR);
        if (fabs(genParticle->eta()) < 1.48) h_PFclusters1000_deltaR_gen_EB->Fill(deltaR);
        else if (genParticle->eta()  > 1.48) h_PFclusters1000_deltaR_gen_EEP->Fill(deltaR);
        else                                 h_PFclusters1000_deltaR_gen_EEM->Fill(deltaR);
      }

      if(deltaR < 1.41*2*0.0174 && itr->pt() > 0.4 ) { // FIXME:  THRESHOLD and itr->pt() > 1.
        size_PFclusters_genMatched++;
        //if(size_PFclusters_genMatched>1) std::cout << "More than one cluster matched to gen particle" << std::endl;

        if(fabs(itr->eta()) < 1.48) {
          size_PFclusters_genMatched_EB++;
          h_PFclusters_genMatched_EB_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_genMatched_EB_energy -> Fill( itr->energy() );
          h_PFclusters_genMatched_EB_et     -> Fill( itr->pt() );
          h_PFclusters_genMatched_EB_eta    -> Fill( itr->eta() );
          h_PFclusters_genMatched_EB_phi    -> Fill( itr->phi() );
          h_PFclusters_genMatched_EB_eOverEtrue->Fill(itr->energy()/genParticle->energy());

          for(TString key : eta_keys["EB"]){
            if( itr->eta() >= eta_edges["EB"][key].first && itr->eta() < eta_edges["EB"][key].second){
              h_PFclusters_genMatched_eOverEtrue_etaBinned["EB"][key]->Fill(itr->energy()/genParticle->energy());
              break; // when you found it, do not loop over the other keys
            }
          }
          for(TString key: Et_keys["EB"]){
            // please note that here we bin based on the gen particle NOT reco particle
            if (genParticle->pt() >= Et_edges["EB"][key].first && genParticle->pt() < Et_edges["EB"][key].second){
              h_PFclusters_genMatched_eOverEtrue_EtBinned["EB"][key]->Fill(itr->energy()/genParticle->energy());
              break; // when you found it, do not loop over the other keys
            }
          }
        }
        else if(itr->eta()  > 1.48) {
          size_PFclusters_genMatched_EEP++;
          h_PFclusters_genMatched_EEP_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_genMatched_EEP_energy -> Fill( itr->energy() );
          h_PFclusters_genMatched_EEP_et     -> Fill( itr->pt() );
          h_PFclusters_genMatched_EEP_eta    -> Fill( itr->eta() );
          h_PFclusters_genMatched_EEP_phi    -> Fill( itr->phi() );
          h_PFclusters_genMatched_EEP_eOverEtrue->Fill(itr->energy()/genParticle->energy());
          for(TString key : eta_keys["EEP"]){
            if( itr->eta() >= eta_edges["EEP"][key].first && itr->eta() < eta_edges["EEP"][key].second){
              h_PFclusters_genMatched_eOverEtrue_etaBinned["EEP"][key]->Fill(itr->energy()/genParticle->energy());
              break; // when you found it, exit
            }
          }
          for(TString key: Et_keys["EEP"]){
            // please note that here we bin based on the gen particle NOT reco particle
            if (genParticle->pt() >= Et_edges["EEP"][key].first && genParticle->pt() < Et_edges["EEP"][key].second){
              h_PFclusters_genMatched_eOverEtrue_EtBinned["EEP"][key]->Fill(itr->energy()/genParticle->energy());
              break; // when you found it, do not loop over the other keys
            }
          }
        } //
        else if(itr->eta()  < -1.48) {
          size_PFclusters_genMatched_EEM++;
          h_PFclusters_genMatched_EEM_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_genMatched_EEM_energy -> Fill( itr->energy() );
          h_PFclusters_genMatched_EEM_et     -> Fill( itr->pt() );
          h_PFclusters_genMatched_EEM_eta    -> Fill( itr->eta() );
          h_PFclusters_genMatched_EEM_phi    -> Fill( itr->phi() );
          h_PFclusters_genMatched_EEM_eOverEtrue->Fill(itr->energy()/genParticle->energy());
          for(TString key : eta_keys["EEM"]){
            if( itr->eta() >= eta_edges["EEM"][key].first && itr->eta() < eta_edges["EEM"][key].second){
              h_PFclusters_genMatched_eOverEtrue_etaBinned["EEM"][key]->Fill(itr->energy()/genParticle->energy());
              break; // when you found it, exit
            }
          }
          for(TString key: Et_keys["EEM"]){
            // please note that here we bin based on the gen particle NOT reco particle
            if (genParticle->pt() >= Et_edges["EEM"][key].first && genParticle->pt() < Et_edges["EEP"][key].second){
              h_PFclusters_genMatched_eOverEtrue_EtBinned["EEM"][key]->Fill(itr->energy()/genParticle->energy());
              break; // when you found it, do not loop over the other keys
            }
          }
        } // end if EEM
        //std::cout << "I am here for event " << naiveId_ << std::endl;
        if(naiveId_<100) h_PFclusters_genMatched_etaVsPhi.at(naiveId_)->Fill(itr->eta(),itr->phi(), itr->energy() );
      } // end if matching
    } // end loop over gen particles

  if(naiveId_<100) h_PFclusters_etaVsPhi.at(naiveId_)->Fill(itr->eta(),itr->phi(), itr->energy() );


  } // end loop on PFclusters

  h_PFclusters_EB_size->Fill(size_PFclusters_EB);
  h_PFclusters_EEP_size->Fill(size_PFclusters_EEP);
  h_PFclusters_EEM_size->Fill(size_PFclusters_EEM);
  h_PFclusters_genMatched_EB_size->Fill(size_PFclusters_genMatched_EB);
  h_PFclusters_genMatched_EEP_size->Fill(size_PFclusters_genMatched_EEP);
  h_PFclusters_genMatched_EEM_size->Fill(size_PFclusters_genMatched_EEM);

  // ---- Super Clusters ---------
  // ... barrel
  edm::Handle<reco::SuperClusterCollection> superClusters_EB_h;
  ev.getByToken( superClusterCollection_EB_, superClusters_EB_h );
  if ( ! superClusters_EB_h.isValid() ) std::cout << "ECALNoiseStudy::analyze --> superClusters_EB_h not found" << std::endl;
  const reco::SuperClusterCollection* theBarrelSuperClusters = superClusters_EB_h.product () ;

  for (reco::SuperClusterCollection::const_iterator itSC = theBarrelSuperClusters->begin(); itSC != theBarrelSuperClusters->end(); ++itSC ) {

    double scEt    = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
    if (scEt < scEtThrEB_ ) continue;
    h_superClusters_EB_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
    h_superClusters_EB_nBC    -> Fill( itSC -> clustersSize());
    h_superClusters_EB_energy -> Fill( itSC -> energy() );
    h_superClusters_EB_rawEnergy -> Fill( itSC -> rawEnergy() );
    h_superClusters_EB_et     -> Fill( scEt );
    h_superClusters_eta       -> Fill( itSC -> eta() );
    h_superClusters_EB_eta    -> Fill( itSC -> eta() );
    h_superClusters_nBCvsEta  -> Fill( itSC -> eta(), itSC -> clustersSize() );

    if ( fabs(itSC->eta())<1 ) h_superClusters_nBC_0to1 -> Fill(itSC->clustersSize());
    if ( fabs(itSC->eta())<1.5 && fabs(itSC->eta())>=1) h_superClusters_nBC_1to1d5 -> Fill(itSC->clustersSize());
  }
  h_superClusters_EB_size         -> Fill( superClusters_EB_h->size() );

  // ... endcap
  edm::Handle<reco::SuperClusterCollection> superClusters_EE_h;
  ev.getByToken( superClusterCollection_EE_, superClusters_EE_h );
  if ( ! superClusters_EE_h.isValid() ) std::cout << "ECALNoiseStudy::analyze --> superClusters_EE_h not found" << std::endl;
  const reco::SuperClusterCollection* theEndcapSuperClusters = superClusters_EE_h.product () ;

  int nSuperClustersEEP = 0;
  int nSuperClustersEEM = 0;
  for (reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin(); itSC != theEndcapSuperClusters->end(); ++itSC ) {

    double scEt = itSC -> energy() * sin(2.*atan( exp(- itSC->position().eta() )));
    if (scEt < scEtThrEE_ ) continue;
    h_superClusters_eta       -> Fill( itSC -> eta() );
    h_superClusters_EE_eta    -> Fill( itSC -> eta() );
    h_superClusters_nBCvsEta  -> Fill( itSC -> eta(), itSC -> clustersSize() );
    if ( fabs(itSC->eta())<1.8 && fabs(itSC->eta())>=1.5) h_superClusters_nBC_1d5to1d8 -> Fill(itSC->clustersSize());
    if ( fabs(itSC->eta())<2.1 && fabs(itSC->eta())>=1.8) h_superClusters_nBC_1d8to2d1 -> Fill(itSC->clustersSize());
    if ( fabs(itSC->eta())<2.5 && fabs(itSC->eta())>=2.1) h_superClusters_nBC_2d1to2d5 -> Fill(itSC->clustersSize());
    if ( fabs(itSC->eta())<=3  && fabs(itSC->eta())>=2.5) h_superClusters_nBC_2d5to3   -> Fill(itSC->clustersSize());

    if  ( itSC -> z() > 0 ){
      h_superClusters_EEP_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
      h_superClusters_EEP_nBC    -> Fill( itSC -> clustersSize() );
      h_superClusters_EEP_energy -> Fill( itSC -> energy() );
      h_superClusters_EEP_rawEnergy -> Fill( itSC -> rawEnergy() );
      h_superClusters_EEP_et -> Fill( scEt );
      nSuperClustersEEP++;
    }
    if  ( itSC -> z() < 0 ){
      h_superClusters_EEM_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
      h_superClusters_EEM_nBC    -> Fill( itSC -> clustersSize() );
      h_superClusters_EEM_energy -> Fill( itSC -> energy() );
      h_superClusters_EEM_rawEnergy -> Fill( itSC -> rawEnergy() );
      h_superClusters_EEM_et -> Fill( scEt );
      nSuperClustersEEM++;
    }
  }

  h_superClusters_EEP_size->Fill( nSuperClustersEEP );
  h_superClusters_EEM_size->Fill( nSuperClustersEEM );

  naiveId_++;

}


// ------------ method called once each job just before starting event loop  ------------
void  ECALNoiseStudy::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void ECALNoiseStudy::endJob() {

  h_numberOfEvents ->Fill(0.,naiveId_);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ECALNoiseStudy);
