#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <iomanip> 

void DrawSlimValidationPlots(const Char_t* infile1, const Char_t* infile2, 
			     const Char_t* fileType = "png", 
			     const Char_t* dirName = ".")
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1); 
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleBorderSize(0);

  if (!infile1 || !infile2) {
    cout << " No input file specified !" << endl;
    return;
  }
  
  cout << "Producing validation plots for: 1) " << infile1 << " and 2) " << infile2 << endl;
  
  TFile* _f[2]; 
  cout << "opening the files in the read only mode" << endl;
  _f[0] = TFile::Open(infile2,"READ"); 
  _f[1] = TFile::Open(infile1,"READ");
  
  if(_f[0]) cout << "file0 is loaded" << endl;
  else cout << "not able to get the file0 " << endl;
  if(_f[1]) cout << "file1 is loaded" << endl;
  else cout << "not able to get the file1 " << endl;
  cout << "successfully opened the files" << endl;
  TH1D *h_numberOfEvents[2];
  for (int i=0;i<2;i++){
    h_numberOfEvents[i] = (TH1D*)_f[i]->Get("ecalslimvalidation/h_recHits_EEM_size") ; 
    cout<< "file[" << i << "] = " << h_numberOfEvents[i]->Integral() << " events " << endl;
  }
  if(h_numberOfEvents[0]) cout << "got histo" << endl;
  else{ 
    cout << "I am sorry histo not uploaded hence exiting ................" << endl;
    exit(1);
    
  }
  double nEvents0 = h_numberOfEvents[0]->Integral();
  double nEvents1 = h_numberOfEvents[1]->Integral();
  cout << std::setprecision(9) << "File blue is " << infile2 << endl;  
  cout << std::setprecision(9) << "File red is "  << infile1 << endl;
  cout << std::setprecision(9) << "Blue has " << nEvents0 << " events, Red has " << nEvents1 << " events" << endl;
  double s        = nEvents1/nEvents0;
  cout << "Ratio of number of events 0/1 "<<s << endl;

  const int nObj=57;
  const char *objName[nObj]={

   "ecalslimvalidation/h_recHits_EB_size",
   "ecalslimvalidation/h_recHits_EEP_size",
   "ecalslimvalidation/h_recHits_EEM_size",
   "ecalslimvalidation/h_recHits_EB_eta",  
   "ecalslimvalidation/h_recHits_EEP_eta",  
   "ecalslimvalidation/h_recHits_EEM_eta",  
   "ecalslimvalidation/h_recHits_EB_maxEneEta",  
   "ecalslimvalidation/h_recHits_EEP_maxEneEta",  
   "ecalslimvalidation/h_recHits_EEM_maxEneEta",  
   "ecalslimvalidation/h_recHits_EB_energy",
   "ecalslimvalidation/h_recHits_EEP_energy",
   "ecalslimvalidation/h_recHits_EEM_energy",
   "ecalslimvalidation/h_recHits_EB_energyMax",
   "ecalslimvalidation/h_recHits_EEP_energyMax",
   "ecalslimvalidation/h_recHits_EEM_energyMax",
   "ecalslimvalidation/h_recHits_eta",
   "ecalslimvalidation/h_recHits_EEP_energy_gt25",
   "ecalslimvalidation/h_recHits_EEM_energy_gt25",

   "ecalslimvalidation/h_basicClusters_EB_size",
   "ecalslimvalidation/h_basicClusters_EEP_size",
   "ecalslimvalidation/h_basicClusters_EEM_size",
   "ecalslimvalidation/h_basicClusters_EB_nXtals",
   "ecalslimvalidation/h_basicClusters_EEP_nXtals",
   "ecalslimvalidation/h_basicClusters_EEM_nXtals",
   "ecalslimvalidation/h_basicClusters_EB_energy",
   "ecalslimvalidation/h_basicClusters_EEP_energy",
   "ecalslimvalidation/h_basicClusters_EEM_energy",
   "ecalslimvalidation/h_basicClusters_eta",
   "ecalslimvalidation/h_basicClusters_EB_eta",
   "ecalslimvalidation/h_basicClusters_EE_eta",

   "ecalslimvalidation/h_superClusters_EB_size",
   "ecalslimvalidation/h_superClusters_EEP_size",
   "ecalslimvalidation/h_superClusters_EEM_size",
   "ecalslimvalidation/h_superClusters_EB_nXtals",
   "ecalslimvalidation/h_superClusters_EEP_nXtals",
   "ecalslimvalidation/h_superClusters_EEM_nXtals",
   "ecalslimvalidation/h_superClusters_EB_nBC",
   "ecalslimvalidation/h_superClusters_EEP_nBC",
   "ecalslimvalidation/h_superClusters_EEM_nBC",
   "ecalslimvalidation/h_superClusters_EB_energy",
   "ecalslimvalidation/h_superClusters_EEP_energy",
   "ecalslimvalidation/h_superClusters_EEM_energy",
   "ecalslimvalidation/h_superClusters_EB_rawEnergy",
   "ecalslimvalidation/h_superClusters_EEP_rawEnergy",
   "ecalslimvalidation/h_superClusters_EEM_rawEnergy",
   "ecalslimvalidation/h_superClusters_EB_et",
   "ecalslimvalidation/h_superClusters_EEP_et",
   "ecalslimvalidation/h_superClusters_EEM_et",
   "ecalslimvalidation/h_superClusters_eta",
   "ecalslimvalidation/h_superClusters_EB_eta",
   "ecalslimvalidation/h_superClusters_EE_eta",

   "ecalslimvalidation/h_superClusters_nBC_0to1",
   "ecalslimvalidation/h_superClusters_nBC_1to1d5",
   "ecalslimvalidation/h_superClusters_nBC_1d5to1d8",
   "ecalslimvalidation/h_superClusters_nBC_1d8to2d1",
   "ecalslimvalidation/h_superClusters_nBC_2d1to2d5",
   "ecalslimvalidation/h_superClusters_nBC_2d5to3"
};
  
const char *objTitle[nObj]={"Number of RecHits (EB)",
			    "Number of RecHits (EE+)",
			    "Number of RecHits (EE-)",
			    "RecHits Eta (EB)",  
			    "RecHits Eta (EEP)",  
			    "RecHits Eta (EEM)",  
			    "MaxEne RecHit Eta (EB)",  
			    "MaxEne RecHit Eta (EEP)",  
			    "MaxEne RecHit Eta (EEM)",  
			    "RecHits Energy (EB)",
			    "RecHits Energy (EE+)",
			    "RecHits Energy (EE-)",
			    "RecHits Max Energy (EB)",
			    "RecHits Max Energy (EE+)",
			    "RecHits Max Energy (EE-)",
			    "RecHits Eta",
			    "RecHits Energy (EE+, high eta)",
			    "RecHits Energy (EE-, high eta)",

			    "Number of Basic Clusters (EB)",
			    "Number of Basic Clusters (EE+)",
			    "Number of Basic Clusters (EE-)",
			    "Number of Crystals per Basic Cluster (EB)",
			    "Number of Crystals per Basic Cluster (EE+)",
			    "Number of Crystals per Basic Cluster (EE-)",
			    "Basic cluster Energy (EB)",
			    "Basic cluster Energy (EE+)",
			    "Basic cluster Energy (EE-)",
			    "Basic cluster Eta",  
			    "Basic cluster Eta (EB)",  
			    "Basic cluster Eta (EE)",  

			    "Number of Superclusters (EB)",
			    "Number of Superclusters (EE+)",
			    "Number of Superclusters (EE-)",
			    "Number of Crystals per Supercluster (EB)",
			    "Number of Crystals per Supercluster (EE+)",
			    "Number of Crystals per Supercluster (EE-)",
			    "Number of Basic Clusters  per Supercluster (EB)",
			    "Number of Basic Clusters  per Supercluster (EE-)",
			    "Number of Basic Clusters  per Supercluster (EE+)",
			    "Supercluster Energy (EB)",
			    "Supercluster Energy (EE+)",
			    "Supercluster Energy (EE-)",
			    "Supercluster RawEnergy (EB)",
			    "Supercluster RawEnergy (EE+)",
			    "Supercluster RawEnergy (EE-)",
			    "Supercluster Et (EB)",
			    "Supercluster Et (EE+)",
			    "Supercluster Et (EE-)",
			    "Superclusters Eta",
			    "Superclusters Eta (EB)",
			    "Superclusters Eta (EE)",

			    "Number of Basic Clusters  per Supercluster [0-1]",
			    "Number of Basic Clusters  per Supercluster [1-1.5]",
			    "Number of Basic Clusters  per Supercluster [1.5-1.8]",
			    "Number of Basic Clusters  per Supercluster [1.8-2.1]",
			    "Number of Basic Clusters  per Supercluster [2.1-2.5]",
			    "Number of Basic Clusters  per Supercluster [2.5-3.]"
};
 
const  char *labelX[nObj]={"Number of RecHits/Event",
			   "Number of RecHits/Event",
			   "Number of RecHits/Event",
			   "eta",
			   "eta",
			   "eta",
			   "max energy rechit eta",
			   "max energy rechit eta",
			   "max energy rechit eta",
			   "Energy (GeV)",
			   "Energy (GeV)",
			   "Energy (GeV)",
			   "Max Energy (GeV)",
			   "Max Energy (GeV)",
			   "Max Energy (GeV)",
			   "eta",
			   "BasicClusters/Event",
			   "BasicClusters/Event",
			   "BasicClusters/Event",
			   "Crystals/BasicCluster",
			   "Crystals/BasicCluster",
			   "Crystals/Basicluster",
			   "Energy (GeV)",
			   "Energy (GeV)",
			   "Energy (GeV)",
			   "eta",
			   "eta",
			   "eta",
			   "Energy (GeV)",
			   "Energy (GeV)",

			   "Superclusters/Event",
			   "Superclusters/Event",
			   "Superclusters/Event",
			   "Crystals/Supercluster",
			   "Crystals/Supercluster",
			   "Crystals/Supercluster",
			   "BasicClusters/Supercluster",
			   "BasicClusters/Supercluster",
			   "BasicClusters/Supercluster",
			   "Energy (GeV)",
			   "Energy (GeV)",
			   "Energy (GeV)",
			   "Raw Energy (GeV)",
			   "Raw Energy (GeV)",
			   "Raw Energy (GeV)",
			   "ET (GeV)",
			   "ET (GeV)",
			   "ET (GeV)",
			   "eta",
			   "eta",
			   "eta",
			   "BasicClusters/Supercluster",
			   "BasicClusters/Supercluster",
			   "BasicClusters/Supercluster",
			   "BasicClusters/Supercluster",
			   "BasicClusters/Supercluster",
			   "BasicClusters/Supercluster"};
                   
 const  char *labelY[1]={"Counts"};
 
 int reBin[nObj]  =  {1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      2,
		      1,1,   

		      1,1,1,
		      1,1,1,
		      4,4,4,
		      2,2,2,

		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      2,2,2,
		      1,1,1,
		      1,1,1};

 int optLogY[nObj] = {0,0,0,
		      1,1,1,
		      1,1,1,
		      0,0,0,
		      0,0,0,
		      1,
		      0,0,   

		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,

		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      1,1,1,
		      0,0,0,
		      0,0,0};

 TH1D* h[2][500]; 
 TCanvas *c[500];
 TPaveStats *st[500];
 
 bool isMissing = false;
 cout << "Finished with all declarations and now starting to loop over all the histograms" << endl;
 int iHisto = 0;
 while(iHisto<nObj){
   for (int ifile=0;ifile<2;ifile++){ 
     h[ifile][iHisto] = (TH1D*)_f[ifile]->Get(objName[iHisto]);
     if(h[ifile][iHisto] == 0)isMissing = true;
     if(h[ifile][iHisto] == 0) continue;
     h[ifile][iHisto]->Rebin(reBin[iHisto]);
     if (ifile == 0) {
       // open a new canvas
       c[iHisto] = new TCanvas(objName[iHisto],objName[iHisto],50+iHisto*20,50+iHisto*5,500,400);
       c[iHisto]->cd();
       // customize and plot
       h[ifile][iHisto]->GetYaxis()->SetTitle(labelY[0]);
       std::string histo_name = std::string(h[ifile][iHisto]->GetName());
       if(histo_name.find("Eta") != std::string::npos) h[ifile][iHisto]->GetYaxis()->SetTitle("ADC");
       h[ifile][iHisto]->GetXaxis()->SetTitle(labelX[iHisto]);
       h[ifile][iHisto]->SetMarkerColor(4);
       h[ifile][iHisto]->SetLineColor(4);
       h[ifile][iHisto]->SetFillStyle(3002);
       h[ifile][iHisto]->SetTitle(objTitle[iHisto]);
       //h[ifile][iHisto]->Scale(float(s));                          // normalize to the same number of events
       h[ifile][iHisto]->Scale(1./h[ifile][iHisto]->Integral());     // normalize to the same area       
       h[ifile][iHisto]->Draw("histE");
     }
     if (ifile == 1) {
       if(isMissing == true) continue;
       h[ifile][iHisto]->Scale(1./h[ifile][iHisto]->Integral());     // normalize to the same area       
       h[ifile][iHisto]->SetMarkerStyle(20);
       h[ifile][iHisto]->SetMarkerSize(0.7);
       h[ifile][iHisto]->SetMarkerColor(kRed);
       h[ifile][iHisto]->SetLineColor(kRed);
       std::string histo_name = std::string(h[ifile][iHisto]->GetName());
       if(histo_name.find("Eta") != std::string::npos) h[ifile][iHisto]->Draw("psameshist");
       else h[ifile][iHisto]->Draw("ephistsames");
       float maxy = max (h[ifile][iHisto]->GetMaximum(),h[ifile-1][iHisto]->GetMaximum() );
       float miny = h[ifile][iHisto]->GetMinimum();
       if (optLogY[iHisto]) {
	 c[iHisto]->SetLogy();
	 h[0][iHisto]->SetMaximum(maxy*2);
	 h[0][iHisto]->SetMinimum(0.00001);
       }
       else  h[0][iHisto]->SetMaximum(maxy*1.1);
       
       c[iHisto]->Update();
     }
     
   }//matches ifile.
 
   //cout << "finished drwaing 1D histograms " << endl;
   string temp = objName[iHisto];
   char *str = std::strtok ((char*)temp.c_str(),"/");
   str = std::strtok (NULL, "/");
   
   char myname[500];
   sprintf (myname,"%s/",dirName);
   strcat(myname,str);
   strcat(myname,".");
   strcat(myname,fileType);
   
   if(isMissing == false) c[iHisto]->Print(myname,fileType);
   isMissing = false;
   
   iHisto++;
   
   
 } //matches while loop. This is making sure that all histos are looped over.

}

