// script to print to pdf file all histograms of one ROOT directory
// example usage: root -l -b "loopdir.C(\"outputfiles/test_nuGun_v7_numEvent1000.root\", \"ecalnoisestudy\")"

void loopdir(TString filename="test_nuGun_v7_numEvent1000.root", TString dirname="ecalnoisestudy") {
   
   TFile *f = TFile::Open(filename);
   TDirectory *d = (TDirectory*)f->Get(dirname);

   TString outfilename = filename + ".pdf";

   TIter keyList(d->GetListOfKeys());
   TKey *key;
   TCanvas c1;
  
   //TString outfile="test_nuGun_v7_numEvent1000.pdf"
   c1.Print(outfilename + "[");
   
   while ((key = (TKey*)keyList())) {
      std::cout << "in file " << std::endl;
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      TH1 *h = (TH1*)key->ReadObj();
      std::cout << h->GetName() << std::endl;
      h->Draw("colz");
      c1.Print(outfilename);
   }
   c1.Print(outfilename + "]");
}

