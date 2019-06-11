// TODO: fix the bad tabs created by ATOM
// TODO: change the binning for energy plots: for EE energy can be as large as et*10
// TODO: exclude dead channels from genParticle counting

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
///#include "RecoEcal/EgammaCoreTools/interface/EcalRecHitLess.h"
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

#include "ECAL/ECALLowLevelComparison/interface/ECALNoiseStudy.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TVector3.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <map>
#include <tuple>

using namespace cms ;
using namespace edm ;
using namespace std ;
using namespace reco;

//
// constructors and destructor
ECALNoiseStudy::ECALNoiseStudy(const edm::ParameterSet& ps)
{
  // analysis configurations
  anaName_                    = ps.getParameter<std::string>("anaName");

  std::cout << "Setting up analysis for channel=" << anaName_ << std::endl;

  // collections to access
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
  scEtThrEB_                   = ps.getParameter<double>("scEtThrEB");
  scEtThrEE_                 = ps.getParameter<double>("scEtThrEE");

  naiveId_ = 0;

  // binning of rechit and pf rechit energy plots
  // the binning can vary from region to the other ( EB vs EE ), but also
  // within EE two regions are defined:
  // * low: |eta| < 2.4 , or roughly equivalently |ring| > 10
  // * high: the complementary region
  std::map<TString, std::map<TString,std::tuple<int, float, float>>> rH_e_bins;
  std::get<0>(rH_e_bins["EB"]["low"]) = 200;
  std::get<1>(rH_e_bins["EB"]["low"]) = 0.;
  std::get<2>(rH_e_bins["EB"]["low"]) = 2.;
  std::get<0>(rH_e_bins["EEP"]["low"]) = 1000;
  std::get<1>(rH_e_bins["EEP"]["low"]) = 0.;
  std::get<2>(rH_e_bins["EEP"]["low"]) = 10.;
  std::get<0>(rH_e_bins["EEM"]["low"]) = 1000;
  std::get<1>(rH_e_bins["EEM"]["low"]) = 0.;
  std::get<2>(rH_e_bins["EEM"]["low"]) = 10.;
  std::get<0>(rH_e_bins["EEP"]["high"]) = 1000;
  std::get<1>(rH_e_bins["EEP"]["high"]) = 0.;
  std::get<2>(rH_e_bins["EEP"]["high"]) = 50.;
  std::get<0>(rH_e_bins["EEM"]["high"]) = 1000;
  std::get<1>(rH_e_bins["EEM"]["high"]) = 0.;
  std::get<2>(rH_e_bins["EEM"]["high"]) = 50.;
  const float low_eta = 2.3;
  const float low_ring = 22.;

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

  // configurations for plots in ring bins
  std::map<TString, Double_t> start_ring;
  start_ring["EB"]=  -90.;//-86.;  // this is in ieta coordinates
  start_ring["EEP"] = 0.;  //(~11.5 -> ~52.) // this is in ir=sqrt((ix-50)**2+(iy-50)**2) coordinates
  start_ring["EEM"] = 0.;

  std::map<TString,Double_t> delta_ring;
  delta_ring["EB"]=  1; // in EB size size of a trigger tower is 5x5
  delta_ring["EEP"]= 1; // in end-caps fixed delta ring means variable number of crystal inside
  delta_ring["EEM"]= 1; // in end-caps fixed delta ring means variable number of crystal inside

  std::map<TString, Int_t> nBins_ring;
  nBins_ring["EB"]  = (int) (-start_ring["EB"]-start_ring["EB"])/delta_ring["EB"];
  nBins_ring["EEP"] = (int) (40-start_ring["EEP"])/delta_ring["EEP"];
  nBins_ring["EEM"] = nBins_ring["EEP"];

  for (TString region : regions){
    // create the keys for the ring bins
    for (Int_t i=0; i<nBins_ring[region]; i++){
      Float_t low_edge = (start_ring[region] + i * delta_ring[region]);
      Float_t up_edge =  (start_ring[region] + (i+1) * delta_ring[region]);
      //std::cout << TString::Format("%.2f_%.2f", low_edge, up_edge) << std::endl;
      TString key = TString::Format("%d_%d", int(low_edge), int(up_edge));
      //if(key.Contains(".")) key.ReplaceAll(".", "dot");
      if(key.Contains("-")) key.ReplaceAll("-", "n");
      ring_keys[region].push_back(key);
      ring_edges[region][key].first = low_edge;
      ring_edges[region][key].second = up_edge;
    }
  }

  // configurations for plots in gen energy bins
  /*// old configurations
  Et_keys.push_back("1_4");
  Et_keys.push_back("4_7");
  Et_keys.push_back("7_10");
  Et_edges["1_4"].first = 1.;
  Et_edges["1_4"].second = 4.;
  Et_edges["4_7"].first = 4.;
  Et_edges["4_7"].second = 7.;
  Et_edges["7_10"].first = 7.;
  Et_edges["7_10"].second = 10.;*/

  // new configurations
  Et_keys.push_back("1_20");
  Et_keys.push_back("20_40");
  Et_keys.push_back("40_60");
  Et_keys.push_back("60_80");
  Et_keys.push_back("80_100");
  Et_edges["1_20"].first = 1.;
  Et_edges["1_20"].second = 20.;
  Et_edges["20_40"].first = 20.;
  Et_edges["20_40"].second = 40.;
  Et_edges["40_60"].first = 40.;
  Et_edges["40_60"].second = 60.;
  Et_edges["60_80"].first = 60.;
  Et_edges["60_80"].second = 80.;
  Et_edges["80_100"].first = 80.;
  Et_edges["80_100"].second = 100.;

  Eta_keys.push_back("0p00_0p50");
  Eta_edges["0p00_0p50"].first = 0.;
  Eta_edges["0p00_0p50"].second = 0.5;
  Eta_keys.push_back("0p50_1p00");
  Eta_edges["0p50_1p00"].first = 0.5;
  Eta_edges["0p50_1p00"].second = 1.0;
  Eta_keys.push_back("1p00_1p48");
  Eta_edges["1p00_1p48"].first = 1.0;
  Eta_edges["1p00_1p48"].second = 1.479;
  Eta_keys.push_back("1p48_2p00");
  Eta_edges["1p48_2p00"].first = 1.479;
  Eta_edges["1p48_2p00"].second = 2.0;
  Eta_keys.push_back("2p00_2p50");
  Eta_edges["2p00_2p50"].first = 2.0;
  Eta_edges["2p00_2p50"].second = 2.5;
  Eta_keys.push_back("2p50_3p00");
  Eta_edges["2p50_3p00"].first = 2.5;
  Eta_edges["2p50_3p00"].second = 3.0;


    // histos (and a few graphs)
  edm::Service<TFileService> fs;
  TFileDirectory genDir = fs->mkdir( "general" );
  TFileDirectory recHitsDir = fs->mkdir( "recHits" );
  TFileDirectory PFrecHitsDir = fs->mkdir( "PFrecHits" );
  TFileDirectory PFclustersDir = fs->mkdir( "PFClusters" );
  TFileDirectory superClustersDir = fs->mkdir( "superClusters" );
  TFileDirectory etaBinnedDir = fs->mkdir( "etaBinnedQuantities" );
  TFileDirectory ringBinnedDir = fs->mkdir( "ringBinnedQuantities" );
  TFileDirectory EtaEtBinnedDir = fs->mkdir( "EtaEtBinnedQuantities" );
  TFileDirectory eventDir = fs->mkdir( "event" );


  // General
  h_numberOfEvents = genDir.make<TH1D>("h_numberOfEvents","h_numberOfEvents",10,0,10);
  h_nPVs = genDir.make<TH1D>("h_nPVs","h_nPVs",80,0,80);

  // gen particles
  h_genP_pt                 = genDir.make<TH1D>("h_genP_pt", "h_genP_pt", 100, 0., 100. );
  h_genP_energy             = genDir.make<TH1D>("h_genP_energy", "h_genP_energy", 100, 0., 100. );
  h_genP_eta                = genDir.make<TH1D>("h_genP_eta", "h_genP_eta", 120, -3., 3. );
  h_genP_phi                = genDir.make<TH1D>("h_genP_phi", "h_genP_phi", 128, -3.2, 3.2);
  h_genP_status             = genDir.make<TH1D>("h_genP_status", "h_genP_status", 10, 0., 10.);
  h_genP_pdgid              = genDir.make<TH1D>("h_genP_pdgid", "h_genP_pdgid", 40, 0., 40.);
  h_genP_pt_EB                = genDir.make<TH1D>("h_genP_pt_EB", "h_genP_pt_EB", 100, 0., 100. );
  h_genP_pt_EEP               = genDir.make<TH1D>("h_genP_pt_EEP", "h_genP_pt_EEP", 100, 0., 100. );
  h_genP_pt_EEM               = genDir.make<TH1D>("h_genP_pt_EEM", "h_genP_pt_EEM", 100, 0., 100. );

  // Rechits, barrel
  g_coord_EB_ieta_eta = recHitsDir.make<TGraph>(); //"h_coord_EB_ieta_eta", "h_coord_EB_ieta_eta", 172, -86., 86., 172, -1.479, 1.479); // here use exact values to get one-to-one correspondance
  g_coord_EB_iphi_phi = recHitsDir.make<TGraph>(); // "h_coord_EB_iphi_phi", "h_coord_EB_iphi_phi", 360, 1., 361., 360,-3.2,3.2); // here use exact values to get one-to-one correspondance
  g_coord_EE_ir_eta = recHitsDir.make<TGraph>();
  g_coord_EE_iphi_phi = recHitsDir.make<TGraph>();
  g_coord_EB_ieta_eta->SetTitle("g_coord_EB_ieta_eta");
  g_coord_EB_ieta_eta->SetName("g_coord_EB_ieta_eta");
  g_coord_EB_iphi_phi->SetTitle("g_coord_EB_iphi_phi");
  g_coord_EB_iphi_phi->SetName("g_coord_EB_iphi_phi");
  g_coord_EE_ir_eta->SetTitle("g_coord_EE_ir_eta");
  g_coord_EE_ir_eta->SetName("g_coord_EE_ir_eta");
  g_coord_EE_iphi_phi->SetTitle("g_coord_EE_iphi_phi");
  g_coord_EE_iphi_phi->SetName("g_coord_EE_iphi_phi");

  h_recHits_EB_size          = recHitsDir.make<TH1D>("h_recHits_EB_size", "h_recHitsEB_size", 100, 500, 3500 );
  h_recHits_EB_eta           = recHitsDir.make<TH1D>("h_recHits_EB_eta","h_recHits_EB_eta",148,-1.48,1.48);
  h_recHits_EB_energy        = recHitsDir.make<TH1D>("h_recHits_EB_energy","h_recHitsEB_energy",1000,0,20);
  h_recHits_EB_energyMax     = recHitsDir.make<TH1D>("h_recHits_EB_energyMax","h_recHitsEB_energyMax",100,0,20);
  h_recHits_EB_energyMaxEta  = recHitsDir.make<TH1D>("h_recHits_EB_energyMaxEta","h_recHits_EB_energyMaxEta",148,-1.48,1.48);
  h_recHits_EB_time          = recHitsDir.make<TH1D>("h_recHits_EB_time","h_recHits_EB_time",400,-100,100);
  h_recHits_EB_Chi2          = recHitsDir.make<TH1D>("h_recHits_EB_Chi2","h_recHits_EB_Chi2",500,0,50);
  //h_recHits_EB_OutOfTimeChi2 = recHitsDir.make<TH1D>("h_recHits_EB_OutOfTimeChi2","h_recHits_EB_OutOfTimeChi2",1000,0,100);
  h_recHits_EB_E1oE4         = recHitsDir.make<TH1D>("h_recHits_EB_E1oE4","h_recHitsEB_E1oE4",148, 0, 1.48);
  h_recHits_EB_iPhiOccupancy = recHitsDir.make<TH1D>("h_recHits_EB_iPhiOccupancy","h_recHits_EB_iPhiOccupancy",360,1.,361. );
  h_recHits_EB_iEtaOccupancy = recHitsDir.make<TH1D>("h_recHits_EB_iEtaOccupancy","h_recHits_EB_iEtaOccupancy",172,-86.,86.);
  h_recHits_EB_occupancy     = recHitsDir.make<TH2D>("h_recHits_EB_occupancy","h_recHits_EB_occupancy",360,1.,361.,172,-86.,86. );
  //h_recHits_EB_energy_etaphi     = recHitsDir.make<TH2D>("h_recHits_EB_energy_etaphi","h_recHits_EB_energy_etaphi",32,-3.2,3.2,30,-1.5,1.5 );
  //h_recHits_EB_energy_ietaiphi     = recHitsDir.make<TH2D>("h_recHits_EB_energy_ietaiphi","h_recHits_EB_energy_ietaiphi",360, 1.,361,172,-86.,86.);
  h_recHits_EB_occupancy_gt10 = recHitsDir.make<TH2D>("h_recHits_EB_occupancy_gt10","h_recHits_EB_occupancy_gt10",360,1.,361.,172,-86.,86. );
  h_recHits_EB_occupancy_lt10 = recHitsDir.make<TH2D>("h_recHits_EB_occupancy_lt10","h_recHits_EB_occupancy_lt10",360,1.,361.,172,-86.,86. );
  h_recHits_EB_energy_3D =  recHitsDir.make<TH3D>("h_recHits_EB_energy_3D","h_recHits_EB_energy_3D",360,1.,361.,172,-86.,86.,1000,0.,50.);

  //h_recHits_EB_eneVSieta     = recHitsDir.make<TH2D>("h_recHits_EB_eneVSieta", "h_recHits_EB_eneVSieta", 100,0,20, 172,-86.,86.);
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
  h_recHits_EEP_energy        = recHitsDir.make<TH1D>("h_recHits_EEP_energy","h_recHits_EEP_energy",50,0,100);
  h_recHits_EEP_energy_gt25   = recHitsDir.make<TH1D>("h_recHits_EEP_energy_gt25","h_recHits_EEP_energy_gt25",50,0,100);
  h_recHits_EEP_energyMax     = recHitsDir.make<TH1D>("h_recHits_EEP_energyMax","h_recHitsEEP_energyMax",50,0,100);
  h_recHits_EEP_energyMaxEta  = recHitsDir.make<TH1D>("h_recHits_EEP_energyMaxEta","h_recHits_EEP_energyMaxEta",74,1.48,3.);
  h_recHits_EEP_time          = recHitsDir.make<TH1D>("h_recHits_EEP_time","h_recHits_EEP_time",400,-100,100);
  h_recHits_EEP_Chi2          = recHitsDir.make<TH1D>("h_recHits_EEP_Chi2","h_recHits_EEP_Chi2",500,0,50);
  h_recHits_EEP_OutOfTimeChi2 = recHitsDir.make<TH1D>("h_recHits_EEP_OutOfTimeChi2","h_recHits_EEP_OutOfTimeChi2",500,0,50);
  h_recHits_EEP_E1oE4         = recHitsDir.make<TH1D>("h_recHits_EEP_E1oE4","h_recHitsEEP_E1oE4",150, 0, 1.5);
  h_recHits_EEP_iXoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEP_iXoccupancy","h_recHits_EEP_iXoccupancy",100,0.,100.);
  h_recHits_EEP_iYoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEP_iYoccupancy","h_recHits_EEP_iYoccupancy",100,0.,100.);
  h_recHits_EEP_occupancy     = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy","h_recHits_EEP_occupancy",100,0.,100.,100,0.,100. );
  //h_recHits_EEP_occupancy_etaphi = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy_etaphi","h_recHits_EEP_occupancy_etaphi",40,1.0,3.0,100,-3.2,3.2 );
  h_recHits_EEP_occupancy_gt10 = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy_gt10","h_recHits_EEP_occupancy_gt10",100,0.,100.,100,0.,100. );
  h_recHits_EEP_occupancy_lt10 = recHitsDir.make<TH2D>("h_recHits_EEP_occupancy_lt10","h_recHits_EEP_occupancy_lt10",100,0.,100.,100,0.,100. );
  h_recHits_EEP_energy_3D      = recHitsDir.make<TH3D>("h_recHits_EEP_energy_3D","h_recHits_EEP_energy_3D",100,0.,100.,100,0.,100.,1000,0.,50.);

  //h_recHits_EEP_energy_etaphi = recHitsDir.make<TH2D>("h_recHits_EEP_energy_etaphi","h_recHits_EEP_energy_etaphi",100,-3.2,3.2,40,1.0,3.0);
  //h_recHits_EEP_energy_ixiy   = recHitsDir.make<TH2D>("h_recHits_EEP_energy_ixiy","h_recHits_EEP_energy_ixiy",100,0.,100.,100,0.,100. );

  h_recHits_EEM_size          = recHitsDir.make<TH1D>("h_recHits_EEM_size","h_recHits_EEM_size",100,0,1000);
  h_recHits_EEM_eta           = recHitsDir.make<TH1D>("h_recHits_EEM_eta","h_recHits_EEM_eta",74,-3.,-1.48);
  h_recHits_EEM_energy        = recHitsDir.make<TH1D>("h_recHits_EEM_energy","h_recHits_EEM_energy",50,0,100);
  h_recHits_EEM_energy_gt25   = recHitsDir.make<TH1D>("h_recHits_EEM_energy_gt25","h_recHits_EEM_energy_gt25",50,0,100);
  h_recHits_EEM_energyMax     = recHitsDir.make<TH1D>("h_recHits_EEM_energyMax","h_recHitsEEM_energyMax",50,0,100);
  h_recHits_EEM_energyMaxEta  = recHitsDir.make<TH1D>("h_recHits_EEM_energyMaxEta","h_recHits_EEM_energyMaxEta",74,-3.,-1.48);
  h_recHits_EEM_time          = recHitsDir.make<TH1D>("h_recHits_EEM_time","h_recHits_EEM_time",400,-100,100);
  h_recHits_EEM_Chi2          = recHitsDir.make<TH1D>("h_recHits_EEM_Chi2","h_recHits_EEM_Chi2",500,0,50);
  h_recHits_EEM_OutOfTimeChi2 = recHitsDir.make<TH1D>("h_recHits_EEM_OutOfTimeChi2","h_recHits_EEM_OutOfTimeChi2",500,0,50);
  h_recHits_EEM_E1oE4         = recHitsDir.make<TH1D>("h_recHits_EEM_E1oE4","h_recHitsEEM_E1oE4",150, 0, 1.5);
  h_recHits_EEM_iXoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEM_iXoccupancy","h_recHits_EEM_iXoccupancy",100,0.,100.);
  h_recHits_EEM_iYoccupancy   = recHitsDir.make<TH1D>("h_recHits_EEM_iYoccupancy","h_recHits_EEM_iYoccupancy",100,0.,100.);
  h_recHits_EEM_occupancy     = recHitsDir.make<TH2D>("h_recHits_EEM_occupancy","h_recHits_EEM_occupancy",100,0.,100.,100,0.,100. );
  h_recHits_EEM_occupancy_gt10 = recHitsDir.make<TH2D>("h_recHits_EEM_occupancy_gt10","h_recHits_EEM_occupancy_gt10",100,0.,100.,100,0.,100. );
  h_recHits_EEM_occupancy_lt10 = recHitsDir.make<TH2D>("h_recHits_EEM_occupancy_lt10","h_recHits_EEM_occupancy_lt10",100,0.,100.,100,0.,100. );
  h_recHits_EEM_energy_3D      = recHitsDir.make<TH3D>("h_recHits_EEM_energy_3D","h_recHits_EEM_energy_3D",100,0.,100.,100,0.,100.,1000,0.,50.);

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
  h_PFrecHits_EB_energy_3D     = PFrecHitsDir.make<TH3D>("h_PFrecHits_EB_energy_3D","h_PFrecHits_EB_energy_3D",360,1.,361.,172,-86.,86.,1000,0.,50.);
  //h_PFrecHits_EB_eneVSieta     = PFrecHitsDir.make<TH2D>("h_PFrecHits_EB_eneVSieta", "h_PFrecHits_EB_eneVSieta", 100,0,20, 172,-86.,86.);

  h_PFrecHits_EEP_eta           = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEP_eta","h_PFrecHits_EEP_eta",74,1.48,3);
  h_PFrecHits_EEP_energy        = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEP_energy","h_PFrecHits_EEP_energy",50,0,100);
  h_PFrecHits_EEP_time          = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEP_time","h_PFrecHits_EEP_time",400,-100,100);
  h_PFrecHits_EEP_occupancy     = PFrecHitsDir.make<TH2D>("h_PFrecHits_EEP_occupancy","h_PFrecHits_EEP_occupancy",100,0.,100.,100,0.,100. );
  h_PFrecHits_EEP_energy_3D = PFrecHitsDir.make<TH3D>("h_PFrecHits_EEP_energy_3D","h_PFrecHits_EEP_energy_3D",100,0.,100.,100,0.,100.,1000,0.,50.);

  h_PFrecHits_EEM_eta           = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEM_eta","h_PFrecHits_EEM_eta",74,-3.,-1.48);
  h_PFrecHits_EEM_energy        = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEM_energy","h_PFrecHits_EEM_energy",50,0,100);
  h_PFrecHits_EEM_time          = PFrecHitsDir.make<TH1D>("h_PFrecHits_EEM_time","h_PFrecHits_EEM_time",400,-100,100);
  h_PFrecHits_EEM_occupancy     = PFrecHitsDir.make<TH2D>("h_PFrecHits_EEM_occupancy","h_PFrecHits_EEM_occupancy",100,0.,100.,100,0.,100. );
  h_PFrecHits_EEM_energy_3D     = PFrecHitsDir.make<TH3D>("h_PFrecHits_EEM_energy_3D","h_PFrecHits_EEM_energy_3D",100,0.,100.,100,0.,100.,1000,0.,50.);

  // --------- PF clusters
  h_PFclusters_EvsEta  = PFclustersDir.make<TH2D>("h_PFclusters_EvsEta", "h_PFclusters_EvsEta",    100, 0., 10., 300,-3.,3.);
  h_PFclusters_EtvsEta = PFclustersDir.make<TH2D>("h_PFclusters_EtvsEta", "h_PFclusters_EtvsEta",    100, 0., 10., 300,-3.,3.);

  h_PFclusters_EB_size    = PFclustersDir.make<TH1D>("h_PFclusters_EB_size","h_PFclusters_EB_size",100,0.,100.);
  h_PFclusters_EB_nXtals  = PFclustersDir.make<TH1D>("h_PFclusters_EB_nXtals","h_PFclusters_EB_nXtals",50,0.,50.);
  h_PFclusters_EB_energy  = PFclustersDir.make<TH1D>("h_PFclusters_EB_energy","h_PFclusters_EB_energy",500,0.,100.);
  h_PFclusters_EB_et      = PFclustersDir.make<TH1D>("h_PFclusters_EB_et","h_PFclusters_EB_et",500,0.,100.);
  h_PFclusters_EB_eta     = PFclustersDir.make<TH1D>("h_PFclusters_EB_eta","h_PFclusters_EB_eta",148,-1.48,1.48);
  h_PFclusters_EB_phi     = PFclustersDir.make<TH1D>("h_PFclusters_EB_phi","h_PFclusters_EB_phi",128,-3.2,3.2);
  //h_PFclusters_EB_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_EB_eOverEtrue","h_PFclusters_EB_eOverEtrue",100,0.,2.);

  h_PFclusters_EEP_size   = PFclustersDir.make<TH1D>("h_PFclusters_EEP_size","h_PFclusters_EEP_size",100,0.,100.);
  h_PFclusters_EEP_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_EEP_nXtals","h_PFclusters_EEP_nXtals",50,0.,50.);
  h_PFclusters_EEP_energy = PFclustersDir.make<TH1D>("h_PFclusters_EEP_energy","h_PFclusters_EEP_energy",500,0.,100.);
  h_PFclusters_EEP_et     = PFclustersDir.make<TH1D>("h_PFclusters_EEP_et","h_PFclusters_EEP_et",500,0.,100.);
  h_PFclusters_EEP_eta    = PFclustersDir.make<TH1D>("h_PFclusters_EEP_eta","h_PFclusters_EEP_eta",300,-3.,3.);
  h_PFclusters_EEP_phi    = PFclustersDir.make<TH1D>("h_PFclusters_EEP_phi","h_PFclusters_EEP_phi",128,-3.2,3.2);
  //h_PFclusters_EEP_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_EEP_eOverEtrue","h_PFclusters_EEP_eOverEtrue",100,0.,2.);

  h_PFclusters_EEM_size   = PFclustersDir.make<TH1D>("h_PFclusters_EEM_size","h_PFclusters_EEM_size",100,0.,100.);
  h_PFclusters_EEM_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_EEM_nXtals","h_PFclusters_EEM_nXtals",50,0.,50.);
  h_PFclusters_EEM_energy = PFclustersDir.make<TH1D>("h_PFclusters_EEM_energy","h_PFclusters_EEM_energy",500,0.,100.);
  h_PFclusters_EEM_et     = PFclustersDir.make<TH1D>("h_PFclusters_EEM_et","h_PFclusters_EEM_et",500,0.,100.);
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
  h_PFclusters_genMatched_EB_energy  = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_energy","h_PFclusters_genMatched_EB_energy",500,0.,100.);
  h_PFclusters_genMatched_EB_et  = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_et","h_PFclusters_genMatched_EB_et",500,0.,100.);
  h_PFclusters_genMatched_EB_eta     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_eta","h_PFclusters_genMatched_EB_eta",148,-1.48,1.48);
  h_PFclusters_genMatched_EB_phi     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_phi","h_PFclusters_genMatched_EB_phi",128,-3.2,3.2);
  h_PFclusters_genMatched_EB_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EB_eOverEtrue","h_PFclusters_genMatched_EB_eOverEtrue",100,0.,2.);

  h_PFclusters_genMatched_EEP_size   = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_size","h_PFclusters_genMatched_EEP_size",100,0.,100.);
  h_PFclusters_genMatched_EEP_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_nXtals","h_PFclusters_genMatched_EEP_nXtals",50,0.,50.);
  h_PFclusters_genMatched_EEP_energy = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_energy","h_PFclusters_genMatched_EEP_energy",500,0.,100.);
  h_PFclusters_genMatched_EEP_et     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_et","h_PFclusters_genMatched_EEP_et",500,0.,100.);
  h_PFclusters_genMatched_EEP_eta    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_eta","h_PFclusters_genMatched_EEP_eta",300,-3.,3.);
  h_PFclusters_genMatched_EEP_phi    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_phi","h_PFclusters_genMatched_EEP_phi",128,-3.2,3.2);
  h_PFclusters_genMatched_EEP_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEP_eOverEtrue","h_PFclusters_genMatched_EEP_eOverEtrue",100,0.,2.);

  h_PFclusters_genMatched_EEM_size   = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_size","h_PFclusters_genMatched_EEM_size",100,0.,100.);
  h_PFclusters_genMatched_EEM_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_nXtals","h_PFclusters_genMatched_EEM_nXtals",50,0.,50.);
  h_PFclusters_genMatched_EEM_energy = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_energy","h_PFclusters_genMatched_EEM_energy",500,0.,100.);
  h_PFclusters_genMatched_EEM_et     = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_et","h_PFclusters_genMatched_EEM_et",500,0.,100.);
  h_PFclusters_genMatched_EEM_eta    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_eta","h_PFclusters_genMatched_EEM_eta",300,-3.,3.);
  h_PFclusters_genMatched_EEM_phi    = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_phi","h_PFclusters_genMatched_EEM_phi",128,-3.2,3.2);
  h_PFclusters_genMatched_EEM_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_genMatched_EEM_eOverEtrue","h_PFclusters_genMatched_EEM_eOverEtrue",100,0.,2.);

  // Distributions for non-matched clusters
  // global plots for EB and EE
  h_PFclusters_noise_EvsEta  = PFclustersDir.make<TH2D>("h_PFclusters_noise_EvsEta", "h_PFclusters_noise_EvsEta",    100, 0., 10., 300,-3.,3.);
  h_PFclusters_noise_EtvsEta = PFclustersDir.make<TH2D>("h_PFclusters_noise_EtvsEta", "h_PFclusters_noise_EtvsEta", 100, 0., 10., 300,-3.,3.);

  h_PFclusters_noise_EB_size   = PFclustersDir.make<TH1D>("h_PFclusters_noise_EB_size","h_PFclusters_noise_EB_size",             100,0.,100.);
  h_PFclusters_noise_EB_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_noise_EB_nXtals","h_PFclusters_noise_EB_nXtals",         50,0.,50.);
  h_PFclusters_noise_EB_energy = PFclustersDir.make<TH1D>("h_PFclusters_noise_EB_energy", "h_PFclusters_noise_EB_energy",        500,0.,100.);
  h_PFclusters_noise_EB_et     = PFclustersDir.make<TH1D>("h_PFclusters_noise_EB_et", "h_PFclusters_noise_EB_et",                100, 0., 10.);
  h_PFclusters_noise_EB_eta    = PFclustersDir.make<TH1D>("h_PFclusters_noise_EB_eta", "h_PFclusters_noise_EB_eta",              300,-3.,3.);
  h_PFclusters_noise_EB_phi    = PFclustersDir.make<TH1D>("h_PFclusters_noise_EB_phi","h_PFclusters_noise_EB_phi",               128,-3.2,3.2);
  h_PFclusters_noise_EB_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_noise_EB_eOverEtrue","h_PFclusters_noise_EB_eOverEtrue",100,0.,2.);

  h_PFclusters_noise_EEP_size   = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEP_size","h_PFclusters_noise_EEP_size",100,0.,100.);
  h_PFclusters_noise_EEP_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEP_nXtals","h_PFclusters_noise_EEP_nXtals",50,0.,50.);
  h_PFclusters_noise_EEP_energy = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEP_energy","h_PFclusters_noise_EEP_energy",500,0.,100.);
  h_PFclusters_noise_EEP_et     = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEP_et","h_PFclusters_noise_EEP_et",500,0.,100.);
  h_PFclusters_noise_EEP_eta    = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEP_eta","h_PFclusters_noise_EEP_eta",300,-3.,3.);
  h_PFclusters_noise_EEP_phi    = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEP_phi","h_PFclusters_noise_EEP_phi",128,-3.2,3.2);
  h_PFclusters_noise_EEP_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEP_eOverEtrue","h_PFclusters_noise_EEP_eOverEtrue",100,0.,2.);

  h_PFclusters_noise_EEM_size   = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEM_size","h_PFclusters_noise_EEM_size",100,0.,100.);
  h_PFclusters_noise_EEM_nXtals = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEM_nXtals","h_PFclusters_noise_EEM_nXtals",50,0.,50.);
  h_PFclusters_noise_EEM_energy = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEM_energy","h_PFclusters_noise_EEM_energy",500,0.,100.);
  h_PFclusters_noise_EEM_et     = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEM_et","h_PFclusters_noise_EEM_et",500,0.,100.);
  h_PFclusters_noise_EEM_eta    = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEM_eta","h_PFclusters_noise_EEM_eta",300,-3.,3.);
  h_PFclusters_noise_EEM_phi    = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEM_phi","h_PFclusters_noise_EEM_phi",128,-3.2,3.2);
  h_PFclusters_noise_EEM_eOverEtrue = PFclustersDir.make<TH1D>("h_PFclusters_noise_EEM_eOverEtrue","h_PFclusters_noise_EEM_eOverEtrue",100,0.,2.);

  // Super Clusters ----------------------------------------------
  // ... barrel
  h_superClusters_EB_size      = superClustersDir.make<TH1D>("h_superClusters_EB_size","h_superClusters_EB_size",15,0.,15.);
  h_superClusters_EB_nXtals    = superClustersDir.make<TH1D>("h_superClusters_EB_nXtals","h_superClusters_EB_nXtals",30,0.,150.);
  h_superClusters_EB_nBC       = superClustersDir.make<TH1D>("h_superClusters_EB_nBC","h_superClusters_EB_nBC",20,0.,20.);
  h_superClusters_EB_energy    = superClustersDir.make<TH1D>("h_superClusters_EB_energy","h_superClusters_EB_energy",500,0.,100.);
  h_superClusters_EB_rawEnergy = superClustersDir.make<TH1D>("h_superClusters_EB_rawEnergy","h_superClusters_EB_rawEnergy",500,0.,100.);
  h_superClusters_EB_rawEt     = superClustersDir.make<TH1D>("h_superClusters_EB_rawEt","h_superClusters_EB_rawEt",500,0.,100.);

  // ... endcap
  h_superClusters_EEP_size   = superClustersDir.make<TH1D>("h_superClusters_EEP_size","h_superClusters_EEP_size",15,0.,15.);
  h_superClusters_EEP_nXtals = superClustersDir.make<TH1D>("h_superClusters_EEP_nXtals","h_superClusters_EEP_nXtals",30,0.,60.);
  h_superClusters_EEP_nBC    = superClustersDir.make<TH1D>("h_superClusters_EEP_nBC","h_superClusters_EEP_nBC",20,0.,20.);
  h_superClusters_EEP_energy = superClustersDir.make<TH1D>("h_superClusters_EEP_energy","h_superClusters_EEP_energy",500,0.,100.);
  h_superClusters_EEP_rawEnergy = superClustersDir.make<TH1D>("h_superClusters_EEP_rawEnergy","h_superClusters_EEP_rawEnergy",500,0.,100.);
  h_superClusters_EEP_rawEt     = superClustersDir.make<TH1D>("h_superClusters_EEP_rawEt","h_superClusters_EEP_rawEt",500,0.,100.);

  h_superClusters_EEM_size   = superClustersDir.make<TH1D>("h_superClusters_EEM_size","h_superClusters_EEM_size",15,0.,15.);
  h_superClusters_EEM_nXtals = superClustersDir.make<TH1D>("h_superClusters_EEM_nXtals","h_superClusters_EEM_nXtals",30,0.,60.);
  h_superClusters_EEM_nBC    = superClustersDir.make<TH1D>("h_superClusters_EEM_nBC","h_superClusters_EEM_nBC",20,0.,20.);
  h_superClusters_EEM_energy = superClustersDir.make<TH1D>("h_superClusters_EEM_energy","h_superClusters_EEM_energy",500,0.,100.);
  h_superClusters_EEM_rawEnergy = superClustersDir.make<TH1D>("h_superClusters_EEM_rawEnergy","h_superClusters_EEM_rawEnergy",500,0.,100.);
  h_superClusters_EEM_rawEt     = superClustersDir.make<TH1D>("h_superClusters_EEM_rawEt","h_superClusters_EEM_rawEt",500,0.,100.);

  h_superClusters_eta        = superClustersDir.make<TH1D>("h_superClusters_eta","h_superClusters_eta",      300,-3.,3.);
  h_superClusters_EB_eta     = superClustersDir.make<TH1D>("h_superClusters_EB_eta","h_superClusters_EB_eta",148,-1.48,1.48);
  h_superClusters_EE_eta     = superClustersDir.make<TH1D>("h_superClusters_EE_eta","h_superClusters_EE_eta",300,-3.,3.);

  // check golden fraction
  h_superClusters_nBCvsEta = superClustersDir.make<TH2D>("h_superClusters_nBCvsEta","h_superClusters_nBCvsEta",300,-3.,3.,20,0.,20.);
  h_superClusters_nBC_0to1     = superClustersDir.make<TH1D>("h_superClusters_nBC_0to1",    "h_superClusters_nBC_0to1",    20,0.,20.);
  h_superClusters_nBC_1to1d5   = superClustersDir.make<TH1D>("h_superClusters_nBC_1to1d5",  "h_superClusters_nBC_1to1d5",  20,0.,20.);
  h_superClusters_nBC_1d5to1d8 = superClustersDir.make<TH1D>("h_superClusters_nBC_1d5to1d8","h_superClusters_nBC_1d5to1d8",20,0.,20.);
  h_superClusters_nBC_1d8to2d1 = superClustersDir.make<TH1D>("h_superClusters_nBC_1d8to2d1","h_superClusters_nBC_1d8to2d1",20,0.,20.);
  h_superClusters_nBC_2d1to2d5 = superClustersDir.make<TH1D>("h_superClusters_nBC_2d1to2d5","h_superClusters_nBC_2d1to2d5",20,0.,20.);
  h_superClusters_nBC_2d5to3   = superClustersDir.make<TH1D>("h_superClusters_nBC_2d5to3",  "h_superClusters_nBC_2d5to3",  20,0.,20.);

  // delta R from truth particle
  h_superClusters_deltaR_gen = superClustersDir.make<TH1D>("h_superClusters_deltaR_gen", "h_superClusters_deltaR_gen", 1000., 0., 10.);
  h_superClusters_deltaR_gen_EB = superClustersDir.make<TH1D>("h_superClusters_deltaR_gen_EB", "h_superClusters_deltaR_gen_EB", 1000., 0., 10.);
  h_superClusters_deltaR_gen_EEP = superClustersDir.make<TH1D>("h_superClusters_deltaR_gen_EEP", "h_superClusters_deltaR_gen_EEP", 1000., 0., 10.);
  h_superClusters_deltaR_gen_EEM = superClustersDir.make<TH1D>("h_superClusters_deltaR_gen_EEM", "h_superClusters_deltaR_gen_EEM", 1000., 0., 10.);

  h_superClusters500_deltaR_gen = superClustersDir.make<TH1D>("h_superClusters500_deltaR_gen", "h_superClusters500_deltaR_gen", 1000., 0., 10.);
  h_superClusters500_deltaR_gen_EB = superClustersDir.make<TH1D>("h_superClusters500_deltaR_gen_EB", "h_superClusters500_deltaR_gen_EB", 1000., 0., 10.);
  h_superClusters500_deltaR_gen_EEP = superClustersDir.make<TH1D>("h_superClusters500_deltaR_gen_EEP", "h_superClusters500_deltaR_gen_EEP", 1000., 0., 10.);
  h_superClusters500_deltaR_gen_EEM = superClustersDir.make<TH1D>("h_superClusters500_deltaR_gen_EEM", "h_superClusters500_deltaR_gen_EEM", 1000., 0., 10.);

  h_superClusters1000_deltaR_gen = superClustersDir.make<TH1D>("h_superClusters1000_deltaR_gen", "h_superClusters1000_deltaR_gen", 1000., 0., 10.);
  h_superClusters1000_deltaR_gen_EB = superClustersDir.make<TH1D>("h_superClusters1000_deltaR_gen_EB", "h_superClusters1000_deltaR_gen_EB", 1000., 0., 10.);
  h_superClusters1000_deltaR_gen_EEP = superClustersDir.make<TH1D>("h_superClusters1000_deltaR_gen_EEP", "h_superClusters1000_deltaR_gen_EEP", 1000., 0., 10.);
  h_superClusters1000_deltaR_gen_EEM = superClustersDir.make<TH1D>("h_superClusters1000_deltaR_gen_EEM", "h_superClusters1000_deltaR_gen_EEM", 1000., 0., 10.);
  // gen matched super clusters
  h_superClusters_genMatched_EB_size    = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_size","h_superClusters_genMatched_EB_size",100,0.,100.);
  h_superClusters_genMatched_EB_nXtals  = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_nXtals","h_superClusters_genMatched_EB_nXtals",50,0.,50.);
  h_superClusters_genMatched_EB_energy  = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_energy","h_superClusters_genMatched_EB_energy",200,0.,10.);
  h_superClusters_genMatched_EB_rawEnergy  = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_rawEnergy","h_superClusters_genMatched_EB_rawEnergy",200,0.,10.);
  h_superClusters_genMatched_EB_rawEt  = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_rawEt","h_superClusters_genMatched_EB_rawEt",200,0.,10.);
  h_superClusters_genMatched_EB_eta     = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_eta","h_superClusters_genMatched_EB_eta",148,-1.48,1.48);
  h_superClusters_genMatched_EB_phi     = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_phi","h_superClusters_genMatched_EB_phi",128,-3.2,3.2);
  h_superClusters_genMatched_EB_eOverEtrue = superClustersDir.make<TH1D>("h_superClusters_genMatched_EB_eOverEtrue","h_superClusters_genMatched_EB_eOverEtrue",100,0.,2.);

  h_superClusters_genMatched_EEP_size   = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_size","h_superClusters_genMatched_EEP_size",100,0.,100.);
  h_superClusters_genMatched_EEP_nXtals = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_nXtals","h_superClusters_genMatched_EEP_nXtals",50,0.,50.);
  h_superClusters_genMatched_EEP_energy = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_energy","h_superClusters_genMatched_EEP_energy",200,0.,10.);
  h_superClusters_genMatched_EEP_rawEnergy = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_rawEnergy","h_superClusters_genMatched_EEP_rawEnergy",200,0.,10.);
  h_superClusters_genMatched_EEP_rawEt     = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_rawEt","h_superClusters_genMatched_EEP_rawEt",200,0.,10.);
  h_superClusters_genMatched_EEP_eta    = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_eta","h_superClusters_genMatched_EEP_eta",300,-3.,3.);
  h_superClusters_genMatched_EEP_phi    = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_phi","h_superClusters_genMatched_EEP_phi",128,-3.2,3.2);
  h_superClusters_genMatched_EEP_eOverEtrue = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEP_eOverEtrue","h_superClusters_genMatched_EEP_eOverEtrue",100,0.,2.);

  h_superClusters_genMatched_EEM_size   = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_size","h_superClusters_genMatched_EEM_size",100,0.,100.);
  h_superClusters_genMatched_EEM_nXtals = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_nXtals","h_superClusters_genMatched_EEM_nXtals",50,0.,50.);
  h_superClusters_genMatched_EEM_energy = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_energy","h_superClusters_genMatched_EEM_energy",200,0.,10.);
  h_superClusters_genMatched_EEM_rawEnergy = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_rawEnergy","h_superClusters_genMatched_EEM_rawEnergy",200,0.,10.);
  h_superClusters_genMatched_EEM_rawEt     = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_rawEt","h_superClusters_genMatched_EEM_rawEt",200,0.,10.);
  h_superClusters_genMatched_EEM_eta    = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_eta","h_superClusters_genMatched_EEM_eta",300,-3.,3.);
  h_superClusters_genMatched_EEM_phi    = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_phi","h_superClusters_genMatched_EEM_phi",128,-3.2,3.2);
  h_superClusters_genMatched_EEM_eOverEtrue = superClustersDir.make<TH1D>("h_superClusters_genMatched_EEM_eOverEtrue","h_superClusters_genMatched_EEM_eOverEtrue",100,0.,2.);

  // Distributions for non-matched clusters
  // FIXME: add distributions as for pf clusters
  h_superClusters_noise_rawEvsEta = superClustersDir.make<TH2D>("h_superClusters_noise_rawEvsEta", "h_superClusters_noise_rawEvsEta",    100, 0., 10., 300,-3.,3.);
  h_superClusters_noise_rawEtvsEta = superClustersDir.make<TH2D>("h_superClusters_noise_rawEtvsEta", "h_superClusters_noise_rawEtvsEta",    100, 0., 10., 300,-3.,3.);
  h_superClusters_noise_rawEnergy = superClustersDir.make<TH1D>("h_superClusters_noise_rawEnergy", "h_superClusters_noise_rawEnergy",    100, 0., 10.);
  h_superClusters_noise_rawEt = superClustersDir.make<TH1D>("h_superClusters_noise_rawEt", "h_superClusters_noise_rawEt",    100, 0., 10.);
  h_superClusters_noise_eta = superClustersDir.make<TH1D>("h_superClusters_noise_eta", "h_superClusters_noise_eta",             300,-3.,3.);


  // --------- Rechits vs eta
  for (TString region : regions){
    for (TString key : eta_keys[region]){
      TString histo_name = "h_recHits_" + region + "_energy_eta_" + key;
      TString subReg = "";
      if (fabs(eta_edges[region][key].first) < low_eta) subReg = "low";
      else subReg = "high";
      h_recHits_energy_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,std::get<0>(rH_e_bins[region][subReg]),std::get<1>(rH_e_bins[region][subReg]),std::get<2>(rH_e_bins[region][subReg]));
      histo_name = "h_recHits_" + region + "_et_" + key;
      h_recHits_et_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,std::get<0>(rH_e_bins[region][subReg]),std::get<1>(rH_e_bins[region][subReg]),std::get<2>(rH_e_bins[region][subReg]));
    }
  }

  // --------- Rechits vs ring
  for (TString region : regions){
    for (TString key : ring_keys[region]){
      TString histo_name = "h_recHits_" + region + "_energy_ring_" + key;
      TString subReg = "";
      if (fabs(ring_edges[region][key].first) > low_ring || region=="EB") subReg = "low";
      else subReg = "high";
      h_recHits_energy_ringBinned[region][key] = ringBinnedDir.make<TH1F>(histo_name,histo_name,std::get<0>(rH_e_bins[region][subReg]),std::get<1>(rH_e_bins[region][subReg]),std::get<2>(rH_e_bins[region][subReg]));
    }
  }

  // --------- PFRechits vs eta
  for (TString region : regions){
    for (TString key : eta_keys[region]){
      TString histo_name = "h_PFrecHits_" + region + "_energy_eta_" + key;
      TString subReg = "";
      if (fabs(eta_edges[region][key].first) < low_eta) subReg = "low";
      else subReg = "high";
      h_PFrecHits_energy_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,std::get<0>(rH_e_bins[region][subReg]),std::get<1>(rH_e_bins[region][subReg]),std::get<2>(rH_e_bins[region][subReg]));
    }
  }

  // --------- PFRechits vs ring
  for (TString region : regions){
    for (TString key : ring_keys[region]){
      TString histo_name = "h_PFrecHits_" + region + "_energy_ring_" + key;
      TString subReg = "";
      if (fabs(ring_edges[region][key].first) > low_ring || region=="EB") subReg = "low";
      else subReg = "high";
      h_PFrecHits_energy_ringBinned[region][key] = ringBinnedDir.make<TH1F>(histo_name,histo_name,std::get<0>(rH_e_bins[region][subReg]),std::get<1>(rH_e_bins[region][subReg]),std::get<2>(rH_e_bins[region][subReg]));
    }
  }

  // --------- E (PF / Super clusters) over E true binned in et and eta
  for (TString Et_key : Et_keys){
    for (TString Eta_key: Eta_keys){
      TString histo_name = "h_PFclusters_genMatched_eOverEtrue_Eta" + Eta_key + "_Et" + Et_key;
      h_PFclusters_genMatched_eOverEtrue_EtaEtBinned[Eta_key][Et_key] = EtaEtBinnedDir.make<TH1F>(histo_name,histo_name,100,0.,2.);
      histo_name = "h_superClusters_genMatched_eOverEtrue_Eta" + Eta_key + "_Et" + Et_key;
      h_superClusters_genMatched_eOverEtrue_EtaEtBinned[Eta_key][Et_key] = EtaEtBinnedDir.make<TH1F>(histo_name,histo_name,100,0.,2.);
      histo_name = "h_genP_nEvts_Eta" + Eta_key + "_Et" + Et_key;
      h_genP_nEvts_EtaEtBinned[Eta_key][Et_key] = EtaEtBinnedDir.make<TH1F>(histo_name,histo_name,1,0.,1.);

    }
  }
  /*
  // --------- E (PF clusters) over E true binned in eta
  for (TString region : regions){
    for (TString key : eta_keys[region]){
      TString histo_name = "h_PFclusters_genMatched_" + region + "_eOverEtrue_" + key;
      h_PFclusters_genMatched_eOverEtrue_etaBinned[region][key] = etaBinnedDir.make<TH1F>(histo_name,histo_name,100,0.,2.);
    }
  }
  // --------- E (PF clusters) over E true binned in et and eta
  for (TString region : regions){
    for (TString key : Et_keys[region]){
      TString histo_name = "h_PFclusters_genMatched_" + region + "_eOverEtrue_" + key;
      h_PFclusters_genMatched_eOverEtrue_EtBinned[region][key] = EtBinnedDir.make<TH1F>(histo_name,histo_name,100,0.,2.);
      histo_name = "h_genP_" + region + "_nEvts_" + key;
      h_genP_nEvts_EtBinned[region][key] = EtBinnedDir.make<TH1F>(histo_name,histo_name,1,0.,1.);
     }
  }*/

  // --------- event by event diagnostic plots
  for (int i=0; i<100; i++){

    // eta vs phi coordinates
    TString histo_name = "h_genP_etaVsPhi_"  + TString::Format("%d", i);
    h_genP_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_recHits_etaVsPhi_" + TString::Format("%d", i);
    h_recHits_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_PFclusters_etaVsPhi_" + TString::Format("%d", i);
    h_PFclusters_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_PFclusters_genMatched_etaVsPhi_" + TString::Format("%d", i);
    h_PFclusters_genMatched_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_superClusters_etaVsPhi_" + TString::Format("%d", i);
    h_superClusters_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));
    histo_name = "h_superClusters_genMatched_etaVsPhi_" + TString::Format("%d", i);
    h_superClusters_genMatched_etaVsPhi.push_back(eventDir.make<TH2D>(histo_name,  histo_name, 344,-3.0,3.0, 360,-3.14,3.14 ));

    // ix iy coordinates - PF rechits and PFclusters
    histo_name = "h_PFrecHits_EB_ietaiphi_" + TString::Format("%d", i);
    h_PFrecHits_EB_ietaiphi.push_back(eventDir.make<TH2D>(histo_name, histo_name, 172,-86.,86.,360,1.,361. ));
    histo_name = "h_PFclusters_EB_ietaiphi_" + TString::Format("%d", i);
    h_PFclusters_EB_ietaiphi.push_back(eventDir.make<TH2D>(histo_name, histo_name, 172,-86.,86.,360,1.,361. ));
    histo_name = "h_PFclusters_genMatched_EB_ietaiphi_" + TString::Format("%d", i);
    h_PFclusters_genMatched_EB_ietaiphi.push_back(eventDir.make<TH2D>(histo_name, histo_name, 172,-86.,86.,360,1.,361. ));

    histo_name = "h_PFrecHits_EEP_ixiy_" + TString::Format("%d", i);
    h_PFrecHits_EEP_ixiy.push_back(eventDir.make<TH2D>(histo_name, histo_name, 100,0.,100.,100,0.,100.));
    histo_name = "h_PFclusters_EEP_ixiy_" + TString::Format("%d", i);
    h_PFclusters_EEP_ixiy.push_back(eventDir.make<TH2D>(histo_name, histo_name, 100,0.,100.,100,0.,100.));
    histo_name = "h_PFclusters_genMatched_EEP_ixiy_" + TString::Format("%d", i);
    h_PFclusters_genMatched_EEP_ixiy.push_back(eventDir.make<TH2D>(histo_name, histo_name, 100,0.,100.,100,0.,100.));

    histo_name = "h_PFrecHits_EEM_ixiy_" + TString::Format("%d", i);
    h_PFrecHits_EEM_ixiy.push_back(eventDir.make<TH2D>(histo_name, histo_name, 100,0.,100.,100,0.,100.));
    histo_name = "h_PFclusters_EEM_ixiy_" + TString::Format("%d", i);
    h_PFclusters_EEM_ixiy.push_back(eventDir.make<TH2D>(histo_name, histo_name, 100,0.,100.,100,0.,100.));
    histo_name = "h_PFclusters_genMatched_EEM_ixiy_" + TString::Format("%d", i);
    h_PFclusters_genMatched_EEM_ixiy.push_back(eventDir.make<TH2D>(histo_name, histo_name, 100,0.,100.,100,0.,100.));
  }

}



ECALNoiseStudy::~ECALNoiseStudy() {}

// ------------ method called to for each event  ------------
void ECALNoiseStudy::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  TH1::StatOverflows(kTRUE);

/*  // Get vertices
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
*/
// MG: commented lines in order to run on 94X tags

  // ******************************
  // ******************************
  // ---- Gen particles ----
  // ******************************
  // ******************************
  edm::Handle<reco::GenParticleCollection> genParticles_handle;
  ev.getByToken( genParticleCollection_, genParticles_handle );
  if ( ! genParticles_handle.isValid() ) std::cout << "ECALNoiseStudy::analyze --> genParticles not found" << std::endl;
  const reco::GenParticleCollection* genParticles = genParticles_handle.product ();

  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin (); genParticle != genParticles->end () ;++genParticle) {

    //if (fabs(genParticle->pdgId())!=11 || genParticle->status()!=1 || genParticle->pt()<1.) continue; // FIXME: added only for quick tests on electrons
    //std::cout << naiveId_ << " pdgid=" << genParticle->pdgId() << " status=" <<  genParticle->status() << " pt=" << genParticle->pt() << std::endl;
    h_genP_pt->Fill(genParticle->pt());
    h_genP_energy->Fill(genParticle->energy());
    h_genP_eta->Fill(genParticle->eta());
    h_genP_phi->Fill(genParticle->phi());
    h_genP_status->Fill(genParticle->status());
    h_genP_pdgid->Fill(fabs(genParticle->pdgId()));
    if(naiveId_<100) h_genP_etaVsPhi.at(naiveId_)->Fill(genParticle->eta(), genParticle->phi(), genParticle->energy());

    for(TString Eta_key: Eta_keys){
      for(TString Et_key: Et_keys){
        if (     genParticle->pt() >= Et_edges[Et_key].first
              && genParticle->pt() < Et_edges[Et_key].second
              && fabs(genParticle->eta()) >= Eta_edges[Eta_key].first
              && fabs(genParticle->eta()) < Eta_edges[Eta_key].second
            ){
          h_genP_nEvts_EtaEtBinned[Eta_key][Et_key]->Fill(0.5);
        }
      }
    }

  } // end loop gen particles

  // calo geometry
  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  const CaloGeometry *geometry = pGeometry.product();


  // ******************************
  // ******************************
  // ---- Rec hits ----
  // ******************************
  // ******************************

  // --- REC HITS, barrel -------------------------------------------------------------------------------------
  edm::Handle<EcalRecHitCollection> recHitsEB;
  ev.getByToken( recHitCollection_EB_, recHitsEB );
  if ( ! recHitsEB.isValid() ) std::cout << "ECALNoiseStudy::analyze --> recHitsEB not found" << std::endl;
  const EcalRecHitCollection* theBarrelEcalRecHits = recHitsEB.product ();

  float maxERecHitEB_ene = -999.;
  float maxERecHitEB_eta = -999.;
  int sizeEB_cleaned = 0;
  int iEBrh = 0;
  for ( EcalRecHitCollection::const_iterator itr = theBarrelEcalRecHits->begin (); itr != theBarrelEcalRecHits->end () ;++itr) {

    EBDetId ebid( itr -> id() );
    GlobalPoint mycell = geometry -> getPosition(DetId(itr->id()));

    // fill geometry graphs
    iEBrh++;
    g_coord_EB_ieta_eta->SetPoint(iEBrh, ebid.ieta(), mycell.eta()); // even if points are filled more than once, it should not matter
    g_coord_EB_iphi_phi->SetPoint(iEBrh, ebid.iphi(), mycell.phi()); // even if points are filled more than once, it should not matter

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
      //h_recHits_EB_eneVSieta     -> Fill( itr->energy() , ebid.ieta() );
      h_recHits_EB_occupancy     -> Fill( ebid.iphi() , ebid.ieta() );
      h_recHits_EB_energy_3D     -> Fill (ebid.iphi() , ebid.ieta(), itr->energy());
      //h_recHits_EB_energy_etaphi     -> Fill( mycell.phi() , mycell.eta(), itr->energy() );
      //h_recHits_EB_energy_ietaiphi     -> Fill( ebid.iphi() ,  ebid.ieta(), itr->energy() );
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

      for(TString key : ring_keys["EB"]){
        if( ebid.ieta() >= ring_edges["EB"][key].first && ebid.ieta() < ring_edges["EB"][key].second){
          h_recHits_energy_ringBinned["EB"][key]->Fill(itr -> energy());
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
  h_recHits_EB_energyMaxEta  -> Fill( maxERecHitEB_eta );


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

      // fill geometry graphs
      g_coord_EE_ir_eta->SetPoint(nHitsEEP, TMath::Sqrt((eeid.ix()-50)*(eeid.ix()-50) + (eeid.iy()-50)*(eeid.iy()-50))-11., mycell.eta()); // even if points are filled more than once, it should not matter
      g_coord_EE_iphi_phi->SetPoint(nHitsEEP, TMath::ATan2((eeid.iy()-50),(eeid.ix()-50)), mycell.phi()); // even if points are filled more than once, it should not matter

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
          h_recHits_EEP_energy_3D     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5, itr->energy());
        //h_recHits_EEP_occupancy_etaphi -> Fill (mycell.eta(), mycell.phi());
        //h_recHits_EEP_energy_etaphi -> Fill (mycell.phi(), mycell.eta(), itr -> energy());
        //h_recHits_EEP_energy_ixiy -> Fill (eeid.ix() - 0.5, eeid.iy() - 0.5, itr -> energy());
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

        for(TString key : ring_keys["EEP"]){
          double r = TMath::Sqrt((eeid.ix()-50)*(eeid.ix()-50) + (eeid.iy()-50)*(eeid.iy()-50)) - 11.;
          if( r >= ring_edges["EEP"][key].first && r < ring_edges["EEP"][key].second){
             h_recHits_energy_ringBinned["EEP"][key]->Fill(itr -> energy());
             break; // when you found it, exit
          }
        }
      }
    } // end EE+

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
        h_recHits_EEM_occupancy     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5 );
        h_recHits_EEM_energy_3D     -> Fill( eeid.ix() - 0.5, eeid.iy() - 0.5, itr->energy());
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
        for(TString key : ring_keys["EEM"]){
          double r = TMath::Sqrt((eeid.ix()-50)*(eeid.ix()-50) + (eeid.iy()-50)*(eeid.iy()-50)) - 11.;
          if( r >= ring_edges["EEM"][key].first && r < ring_edges["EEM"][key].second){
             h_recHits_energy_ringBinned["EEM"][key]->Fill(itr -> energy());
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
  h_recHits_EEP_energyMaxEta -> Fill(maxERecHitEEP_eta);
  h_recHits_EEM_energyMaxEta -> Fill(maxERecHitEEM_eta);



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

  // ******************************
  // ******************************
  // ---- PF Rec hits ----
  // ******************************
  // ******************************

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
        //h_PFrecHits_EB_eneVSieta     -> Fill( itr->energy() , ebid.ieta() );
        h_PFrecHits_EB_occupancy     -> Fill( ebid.iphi() , ebid.ieta() );
        h_PFrecHits_EB_energy_3D     -> Fill( ebid.iphi() , ebid.ieta(), itr -> energy() > ethrEB_ );
        h_PFrecHits_EB_eta           -> Fill( mycell.eta() );
        h_PFrecHits_EB_energy        -> Fill( itr->energy() );

        if(naiveId_<100) h_PFrecHits_EB_ietaiphi.at(naiveId_)->Fill(ebid.ieta(), ebid.iphi(), itr->energy());

        for(TString key : eta_keys["EB"]){
          if( mycell.eta() >= eta_edges["EB"][key].first && mycell.eta() < eta_edges["EB"][key].second){
            h_PFrecHits_energy_etaBinned["EB"][key]->Fill(itr -> energy());
            break; // when you found it, exit
          }
        }
        for(TString key : ring_keys["EB"]){
          if( ebid.ieta() >= ring_edges["EB"][key].first && ebid.ieta() < ring_edges["EB"][key].second){
            h_PFrecHits_energy_ringBinned["EB"][key]->Fill(itr -> energy());
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
          h_PFrecHits_EEP_energy_3D     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5, itr->energy() );
          h_PFrecHits_EEP_eta           -> Fill( mycell.eta() );
          h_PFrecHits_EEP_energy        -> Fill( itr->energy() );

	  if(naiveId_<100) h_PFrecHits_EEP_ixiy.at(naiveId_)->Fill(eeid.ix()- 0.5, eeid.iy() - 0.5, itr->energy());

          for(TString key : eta_keys["EEP"]){
            if( mycell.eta() >= eta_edges["EEP"][key].first && mycell.eta() < eta_edges["EEP"][key].second){
              h_PFrecHits_energy_etaBinned["EEP"][key]->Fill(itr -> energy());
              break; // when you found it, exit
            }
          }
          for(TString key : ring_keys["EEP"]){
            double r = TMath::Sqrt((eeid.ix()-50)*(eeid.ix()-50) + (eeid.iy()-50)*(eeid.iy()-50)) - 11.;
            if( r >= ring_edges["EEP"][key].first && r < ring_edges["EEP"][key].second){
               h_PFrecHits_energy_ringBinned["EEP"][key]->Fill(itr -> energy());
               break; // when you found it, exit
            }
          }
        }
        // EEM
        else {

          h_PFrecHits_EEM_time          -> Fill( itr -> time() );
          h_PFrecHits_EEM_occupancy     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5 );
          h_PFrecHits_EEM_energy_3D     -> Fill( eeid.ix()- 0.5, eeid.iy() - 0.5, itr->energy() );
          h_PFrecHits_EEM_eta           -> Fill( mycell.eta() );
          h_PFrecHits_EEM_energy        -> Fill( itr->energy() );

	  if(naiveId_<100) h_PFrecHits_EEM_ixiy.at(naiveId_)->Fill(eeid.ix()- 0.5, eeid.iy() - 0.5, itr->energy());

          for(TString key : eta_keys["EEM"]){
            //std::cout << key << "  " << itr->energy() << std::endl;
            if( mycell.eta() >= eta_edges["EEM"][key].first && mycell.eta() < eta_edges["EEM"][key].second){
              h_PFrecHits_energy_etaBinned["EEM"][key]->Fill(itr -> energy());
              break; // when you found it, exit
            }
          }
          for(TString key : ring_keys["EEM"]){
            double r = TMath::Sqrt((eeid.ix()-50)*(eeid.ix()-50) + (eeid.iy()-50)*(eeid.iy()-50)) - 11.;
            if( r >= ring_edges["EEM"][key].first && r < ring_edges["EEM"][key].second){
               h_PFrecHits_energy_ringBinned["EEM"][key]->Fill(itr -> energy());
               break; // when you found it, exit
            }
          }
        } // end EE-
      } // end end-caps
    } // end if threshold
  } // end loop over pfrechits

  // ******************************
  // ******************************
  // ---- PF CLUSTERS ----
  // ******************************
  // ******************************
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

    //std::cout << "DEBUG components of cluster" << std::endl;
    //for (auto rhf : (*itr).recHitFractions()){
    //std::cout << " fraction of rechit owned by pfcluster: "<< rhf.fraction()  << std::endl;
    //std::cout << " pointer to rechit: detid=" << rhf.recHitRef()->detId() << " energy=" << rhf.recHitRef()->energy() << std::endl;
    //std::cout << " out " << rhf << std::endl;
    //}

    size_PFclusters++;
    //std::cout << size_PFclusters << "  pt=" << itr->pt() << std::endl;
    h_PFclusters_eta -> Fill( itr->eta() );
    h_PFclusters_EvsEta -> Fill(itr->energy(), itr->eta());
    h_PFclusters_EtvsEta -> Fill(itr->pt(), itr->eta());

    // barrel
    if (fabs(itr->eta()) < 1.48){
      size_PFclusters_EB++;
      h_PFclusters_EB_nXtals -> Fill( (*itr).hitsAndFractions().size() );
      h_PFclusters_EB_energy -> Fill( itr->energy() );
      h_PFclusters_EB_et     -> Fill( itr->pt() );
      h_PFclusters_EB_eta    -> Fill( itr->eta() );
      h_PFclusters_EB_phi    -> Fill( itr->phi() );
      if(naiveId_<100) {
        for (auto pfrhf : (*itr).recHitFractions()){ // object of type PFRecHitFraction
          EBDetId ebid( pfrhf.recHitRef()->detId() );
          h_PFclusters_EB_ietaiphi.at(naiveId_)->Fill(ebid.ieta(),ebid.iphi(), pfrhf.recHitRef()->energy() );
        }
      }
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

    if(naiveId_<100) h_PFclusters_etaVsPhi.at(naiveId_)->Fill(itr->eta(),itr->phi(), itr->energy() );

    // cluster matching with gen particle
    // delta R
    for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin (); genParticle != genParticles->end () ;++genParticle) {

      // matchin only with photons / electrons of status 1
      if(anaName_ ==      "DoublePhoton") {
        if (genParticle->pdgId()!=22 or genParticle->status()!= 1) continue;
      }
      else if(anaName_ == "DoubleElectron") {
        if (fabs(genParticle->pdgId())!=11 or genParticle->status()!= 1) continue;
      }
      else {
        std::cout << "**********************" << std::endl;
        std::cout << "WARNING: not going to select any particle status or pdg id: are you sure this is what you want ?" << std::endl;
        std::cout << "**********************" << std::endl;

      }
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

      if(deltaR < 1.41*2*0.0174){//1.41*2*0.0174) { // Delta R chosen to be ~ size of a crystal ~ = 1.41*2*0.0174 =~ 0.05 //  THRESHOLD and itr->pt() > 1.
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
        }
        else if(itr->eta()  > 1.48) {
          size_PFclusters_genMatched_EEP++;
          h_PFclusters_genMatched_EEP_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_genMatched_EEP_energy -> Fill( itr->energy() );
          h_PFclusters_genMatched_EEP_et     -> Fill( itr->pt() );
          h_PFclusters_genMatched_EEP_eta    -> Fill( itr->eta() );
          h_PFclusters_genMatched_EEP_phi    -> Fill( itr->phi() );
          h_PFclusters_genMatched_EEP_eOverEtrue->Fill(itr->energy()/genParticle->energy());
        } //
        else if(itr->eta()  < -1.48) {
          size_PFclusters_genMatched_EEM++;
          h_PFclusters_genMatched_EEM_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_genMatched_EEM_energy -> Fill( itr->energy() );
          h_PFclusters_genMatched_EEM_et     -> Fill( itr->pt() );
          h_PFclusters_genMatched_EEM_eta    -> Fill( itr->eta() );
          h_PFclusters_genMatched_EEM_phi    -> Fill( itr->phi() );
          h_PFclusters_genMatched_EEM_eOverEtrue->Fill(itr->energy()/genParticle->energy());
        } // end if EEM
        //std::cout << "I am here for event " << naiveId_ << std::endl;
        if(naiveId_<100) h_PFclusters_genMatched_etaVsPhi.at(naiveId_)->Fill(itr->eta(),itr->phi(), itr->energy() );


        for(TString Eta_key: Eta_keys){
          for(TString Et_key: Et_keys){
            if (     genParticle->pt() >= Et_edges[Et_key].first
                  && genParticle->pt() < Et_edges[Et_key].second
                  && fabs(genParticle->eta()) >= Eta_edges[Eta_key].first
                  && fabs(genParticle->eta()) < Eta_edges[Eta_key].second
                ){
              h_PFclusters_genMatched_eOverEtrue_EtaEtBinned[Eta_key][Et_key]->Fill(itr->energy()/genParticle->energy());
            }
          }
        }
        /*for(TString key : eta_keys["EEM"]){
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
        }*/


      } // end if matching
      else { // non matched clusters

        h_PFclusters_noise_EvsEta->Fill(itr->energy(), itr->eta());
        h_PFclusters_noise_EtvsEta->Fill(itr->pt(), itr->eta());

        if(fabs(itr->eta()) < 1.48) {
          h_PFclusters_noise_EB_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_noise_EB_energy -> Fill( itr->energy() );
          h_PFclusters_noise_EB_et     -> Fill( itr->pt() );
          h_PFclusters_noise_EB_eta    -> Fill( itr->eta() );
          h_PFclusters_noise_EB_phi    -> Fill( itr->phi() );
          h_PFclusters_noise_EB_eOverEtrue->Fill(itr->energy()/genParticle->energy());
        }
        else if(itr->eta()  > 1.48) {
          h_PFclusters_noise_EEP_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_noise_EEP_energy -> Fill( itr->energy() );
          h_PFclusters_noise_EEP_et     -> Fill( itr->pt() );
          h_PFclusters_noise_EEP_eta    -> Fill( itr->eta() );
          h_PFclusters_noise_EEP_phi    -> Fill( itr->phi() );
          h_PFclusters_noise_EEP_eOverEtrue->Fill(itr->energy()/genParticle->energy());
        } //
        else if(itr->eta()  < -1.48) {
          h_PFclusters_noise_EEM_nXtals -> Fill( (*itr).hitsAndFractions().size() );
          h_PFclusters_noise_EEM_energy -> Fill( itr->energy() );
          h_PFclusters_noise_EEM_et     -> Fill( itr->pt() );
          h_PFclusters_noise_EEM_eta    -> Fill( itr->eta() );
          h_PFclusters_noise_EEM_phi    -> Fill( itr->phi() );
          h_PFclusters_noise_EEM_eOverEtrue->Fill(itr->energy()/genParticle->energy());
        } // end if EEM
      } // end if not matching
    } // end loop over gen particles

  } // end loop on PFclusters

  h_PFclusters_EB_size->Fill(size_PFclusters_EB);
  h_PFclusters_EEP_size->Fill(size_PFclusters_EEP);
  h_PFclusters_EEM_size->Fill(size_PFclusters_EEM);
  h_PFclusters_genMatched_EB_size->Fill(size_PFclusters_genMatched_EB);
  h_PFclusters_genMatched_EEP_size->Fill(size_PFclusters_genMatched_EEP);
  h_PFclusters_genMatched_EEM_size->Fill(size_PFclusters_genMatched_EEM);
  h_PFclusters_noise_EB_size->Fill(size_PFclusters_EB-size_PFclusters_genMatched_EB);
  h_PFclusters_noise_EEP_size->Fill(size_PFclusters_EEP-size_PFclusters_genMatched_EEP);
  h_PFclusters_noise_EEM_size->Fill(size_PFclusters_EEM-size_PFclusters_genMatched_EEM);

  // ******************************
  // ******************************
  // ---- Super Clusters ---------
  // ******************************
  // ******************************
  // ... barrel
  int size_superClusters_genMatched = 0; // how many gen matched super clusters in EB in the event
  int size_superClusters_genMatched_EB = 0; // how many gen matched super clusters in EB in the event
  int size_superClusters_genMatched_EEP = 0; // ""
  int size_superClusters_genMatched_EEM = 0; // ""

  edm::Handle<reco::SuperClusterCollection> superClusters_EB_h;
  ev.getByToken( superClusterCollection_EB_, superClusters_EB_h );
  if ( ! superClusters_EB_h.isValid() ) std::cout << "ECALNoiseStudy::analyze --> superClusters_EB_h not found" << std::endl;
  const reco::SuperClusterCollection* theBarrelSuperClusters = superClusters_EB_h.product () ;

  for (reco::SuperClusterCollection::const_iterator itSC = theBarrelSuperClusters->begin(); itSC != theBarrelSuperClusters->end(); ++itSC ) {

    double rawSCet    = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));
    if (rawSCet < scEtThrEB_ ) continue;
    h_superClusters_EB_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
    h_superClusters_EB_nBC    -> Fill( itSC -> clustersSize());
    h_superClusters_EB_energy -> Fill( itSC -> energy() );
    h_superClusters_EB_rawEnergy -> Fill( itSC -> rawEnergy() );
    h_superClusters_EB_rawEt     -> Fill( rawSCet );
    h_superClusters_eta       -> Fill( itSC -> eta() ); // why using normal eta here and not geometric eta ?
    h_superClusters_EB_eta    -> Fill( itSC -> eta() );
    h_superClusters_nBCvsEta  -> Fill( itSC -> eta(), itSC -> clustersSize() );

    if ( fabs(itSC->eta())<1 ) h_superClusters_nBC_0to1 -> Fill(itSC->clustersSize());
    if ( fabs(itSC->eta())<1.5 && fabs(itSC->eta())>=1) h_superClusters_nBC_1to1d5 -> Fill(itSC->clustersSize());

    // superluster matching with gen particle
    // delta R
    for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin (); genParticle != genParticles->end () ;++genParticle) {

      // matchin only with photons / electrons of status 1
      if(anaName_ ==      "DoublePhoton") {
        if (genParticle->pdgId()!=22 or genParticle->status()!= 1) continue;
      }
      else if(anaName_ == "DoubleElectron") {
        if (fabs(genParticle->pdgId())!=11 or genParticle->status()!= 1) continue;
      }
      else {
        std::cout << "**********************" << std::endl;
        std::cout << "WARNING: not going to select any particle status or pdg id: are you sure this is what you want ?" << std::endl;
        std::cout << "**********************" << std::endl;

      }
      //if(genParticle->pt()<9) continue; //  only consider clusters matched to 9-10 GeV photons

      double deltaPhi = TVector2::Phi_mpi_pi( genParticle->phi() - itSC->phi());
      double deltaEta = genParticle->eta() - itSC->eta();
      double deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

      h_superClusters_deltaR_gen->Fill(deltaR);
      h_superClusters_deltaR_gen_EB->Fill(deltaR);

      if (rawSCet>0.5) { // 500 MeV
        h_superClusters500_deltaR_gen->Fill(deltaR);
        h_superClusters500_deltaR_gen_EB->Fill(deltaR);
      }
      if (rawSCet>1.0){
        h_superClusters1000_deltaR_gen->Fill(deltaR);
        h_superClusters1000_deltaR_gen_EB->Fill(deltaR);
      }

      // gen matched super clusters, DeltaR
      if(deltaR <0.25 ) { // Delta R chosen from DeltaR plot
        size_superClusters_genMatched++;
        //if(size_superClusters_genMatched>1) std::cout << "More than one super cluster matched to gen particle" << std::endl;

        size_superClusters_genMatched_EB++;
        h_superClusters_genMatched_EB_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
        h_superClusters_genMatched_EB_energy -> Fill( itSC->energy() );
        h_superClusters_genMatched_EB_rawEnergy -> Fill( itSC->rawEnergy() );
        h_superClusters_genMatched_EB_rawEt     -> Fill( rawSCet ); // calculated value
        h_superClusters_genMatched_EB_eta    -> Fill( itSC->eta() );
        h_superClusters_genMatched_EB_phi    -> Fill( itSC->phi() );
        h_superClusters_genMatched_EB_eOverEtrue->Fill(itSC->rawEnergy()/genParticle->energy());
        //std::cout << "I am here for event " << naiveId_ << std::endl;
        if(naiveId_<100) h_superClusters_genMatched_etaVsPhi.at(naiveId_)->Fill(itSC->eta(),itSC->phi(), itSC->rawEnergy() );


        for(TString Eta_key: Eta_keys){
          for(TString Et_key: Et_keys){
            if (     genParticle->pt() >= Et_edges[Et_key].first
                  && genParticle->pt() < Et_edges[Et_key].second
                  && fabs(genParticle->eta()) >= Eta_edges[Eta_key].first
                  && fabs(genParticle->eta()) < Eta_edges[Eta_key].second
                ){
              h_superClusters_genMatched_eOverEtrue_EtaEtBinned[Eta_key][Et_key]->Fill(itSC->energy()/genParticle->energy());
            }
          }
        }

      } // end if matching


    } // end loop over gen particles

    if(naiveId_<100) h_superClusters_etaVsPhi.at(naiveId_)->Fill(itSC->eta(),itSC->phi(), itSC->rawEnergy() );

  } // end loop on EB superclusters
  h_superClusters_EB_size                    -> Fill( superClusters_EB_h->size() );
  h_superClusters_genMatched_EB_size         -> Fill( size_superClusters_genMatched_EB );

  // ... endcap
  edm::Handle<reco::SuperClusterCollection> superClusters_EE_h;
  ev.getByToken( superClusterCollection_EE_, superClusters_EE_h );
  if ( ! superClusters_EE_h.isValid() ) std::cout << "ECALNoiseStudy::analyze --> superClusters_EE_h not found" << std::endl;
  const reco::SuperClusterCollection* theEndcapSuperClusters = superClusters_EE_h.product () ;

  int nSuperClustersEEP = 0;
  int nSuperClustersEEM = 0;
  for (reco::SuperClusterCollection::const_iterator itSC = theEndcapSuperClusters->begin(); itSC != theEndcapSuperClusters->end(); ++itSC ) {

    double rawSCet = itSC -> rawEnergy() * sin(2.*atan( exp(- itSC->position().eta() )));
    if (rawSCet < scEtThrEE_ ) continue;
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
      h_superClusters_EEP_rawEt -> Fill( rawSCet );
      nSuperClustersEEP++;
    }
    if  ( itSC -> z() < 0 ){
      h_superClusters_EEM_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
      h_superClusters_EEM_nBC    -> Fill( itSC -> clustersSize() );
      h_superClusters_EEM_energy -> Fill( itSC -> energy() );
      h_superClusters_EEM_rawEnergy -> Fill( itSC -> rawEnergy() );
      h_superClusters_EEM_rawEt -> Fill( rawSCet );
      nSuperClustersEEM++;
    }

    // trying to match with gen particles
    for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin (); genParticle != genParticles->end () ;++genParticle) {

      // matchin only with photons / electrons of status 1
      if(anaName_ ==      "DoublePhoton") {
        if (genParticle->pdgId()!=22 or genParticle->status()!= 1) continue;
      }
      else if(anaName_ == "DoubleElectron") {
        if (fabs(genParticle->pdgId())!=11 or genParticle->status()!= 1) continue;
      }
      else {
        std::cout << "**********************" << std::endl;
        std::cout << "WARNING: not going to select any particle status or pdg id: are you sure this is what you want ?" << std::endl;
        std::cout << "**********************" << std::endl;

      }

      double deltaPhi = TVector2::Phi_mpi_pi( genParticle->phi() - itSC->phi());
      double deltaEta = genParticle->eta() - itSC->eta();
      double deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

      h_superClusters_deltaR_gen->Fill(deltaR);
      if  ( itSC -> z() > 0 ) h_superClusters_deltaR_gen_EEP->Fill(deltaR);
      else                    h_superClusters_deltaR_gen_EEM->Fill(deltaR);

      if (rawSCet>0.5) { // 500 MeV
        h_superClusters500_deltaR_gen->Fill(deltaR);
        if  ( itSC -> z() > 0 ) h_superClusters500_deltaR_gen_EEP->Fill(deltaR);
        else                    h_superClusters500_deltaR_gen_EEM->Fill(deltaR);
      }
      if (rawSCet>1.0) { // 1 GeV
        h_superClusters1000_deltaR_gen->Fill(deltaR);
        if  ( itSC -> z() > 0 ) h_superClusters1000_deltaR_gen_EEP->Fill(deltaR);
        else                    h_superClusters1000_deltaR_gen_EEM->Fill(deltaR);
      }
      // gen matched super clusters, DeltaR
      if(deltaR <0.25 ) { // Delta R chosen from DeltaR plot

        size_superClusters_genMatched++;

        if  ( itSC -> z() > 0 ){

          size_superClusters_genMatched_EEP++;
          h_superClusters_genMatched_EEP_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
          h_superClusters_genMatched_EEP_energy -> Fill( itSC->energy() );
          h_superClusters_genMatched_EEP_rawEnergy -> Fill( itSC->rawEnergy() );
          h_superClusters_genMatched_EEP_rawEt     -> Fill( rawSCet ); // calculated value
          h_superClusters_genMatched_EEP_eta    -> Fill( itSC->eta() );
          h_superClusters_genMatched_EEP_phi    -> Fill( itSC->phi() );
          h_superClusters_genMatched_EEP_eOverEtrue->Fill(itSC->rawEnergy()/genParticle->energy());

        } else if ( itSC -> z() < 0 ){

          size_superClusters_genMatched_EEM++;
          h_superClusters_genMatched_EEM_nXtals -> Fill( (*itSC).hitsAndFractions().size() );
          h_superClusters_genMatched_EEM_energy -> Fill( itSC->energy() );
          h_superClusters_genMatched_EEM_rawEnergy -> Fill( itSC->rawEnergy() );
          h_superClusters_genMatched_EEM_rawEt     -> Fill( rawSCet); // calculated value
          h_superClusters_genMatched_EEM_eta    -> Fill( itSC->eta() );
          h_superClusters_genMatched_EEM_phi    -> Fill( itSC->phi() );
          h_superClusters_genMatched_EEM_eOverEtrue->Fill(itSC->rawEnergy()/genParticle->energy());

        }

        if(naiveId_<100) h_superClusters_genMatched_etaVsPhi.at(naiveId_)->Fill(itSC->eta(),itSC->phi(), itSC->rawEnergy() );

        for(TString Eta_key: Eta_keys){
          for(TString Et_key: Et_keys){
            if (     genParticle->pt() >= Et_edges[Et_key].first
                  && genParticle->pt() < Et_edges[Et_key].second
                  && fabs(genParticle->eta()) >= Eta_edges[Eta_key].first
                  && fabs(genParticle->eta()) < Eta_edges[Eta_key].second
                ){
              h_superClusters_genMatched_eOverEtrue_EtaEtBinned[Eta_key][Et_key]->Fill(itSC->energy()/genParticle->energy());
            }
          }
        }

      } // end if matching
      else{
        h_superClusters_noise_rawEvsEta->Fill(itSC->rawEnergy(), itSC->eta());
        h_superClusters_noise_rawEtvsEta->Fill(rawSCet, itSC->eta());
        h_superClusters_noise_rawEnergy->Fill(itSC->rawEnergy());
        h_superClusters_noise_rawEt->Fill(rawSCet);
        h_superClusters_noise_eta->Fill(itSC->eta());
      } // end if not matching
    } // end loop over gen particles

  if(naiveId_<100) h_superClusters_etaVsPhi.at(naiveId_)->Fill(itSC->eta(),itSC->phi(), itSC->rawEnergy() );

  } // end for super clusters end-caps

  h_superClusters_EEP_size->Fill( nSuperClustersEEP );
  h_superClusters_EEM_size->Fill( nSuperClustersEEM );
  h_superClusters_genMatched_EEP_size         -> Fill( size_superClusters_genMatched_EEP );
  h_superClusters_genMatched_EEM_size         -> Fill( size_superClusters_genMatched_EEM );

  naiveId_++;





}


// ------------ method called once each job just before starting event loop  ------------
void  ECALNoiseStudy::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void ECALNoiseStudy::endJob() {

  h_numberOfEvents ->Fill(0.,naiveId_);

  // Write graphs to file
  g_coord_EB_ieta_eta->Write();
  g_coord_EB_iphi_phi->Write();
  g_coord_EE_ir_eta->Write();
  g_coord_EE_iphi_phi->Write();
}


//define this as a plug-in
DEFINE_FWK_MODULE(ECALNoiseStudy);
