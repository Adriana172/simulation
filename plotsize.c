{
#include <vector>

   using namespace std;
   TFile *myfile;
   vector<Double_t> *mypair = new std::vector<Double_t>();
   Double_t backrates[26], sizes[26];
   Int_t i = 0;
   TString bkg;

   //    sizes[0] = 100;
   /* new file */
   /* new file */

   for (Int_t k = 0; k < 1000; k += 100) {

      bkg = "";
      bkg += k;
      bkg += ".root";
      myfile = new TFile(bkg);

      gingko->SetBranchAddress("sizeandbkg", &mypair);
      gingko->GetEntry(0);

      backrates[i] = (*mypair)[0];
      sizes[i] = k;
      i++;
   }

   for (Int_t k = 1000; k < 10000; k += 1000) {

      bkg = "";
      bkg += k;
      bkg += ".root";
      myfile = new TFile(bkg);

      gingko->SetBranchAddress("sizeandbkg", &mypair);
      gingko->GetEntry(0);

      backrates[i] = (*mypair)[0];
      sizes[i] = k;
      i++;
   }

   for (Int_t k = 10000; k < 80000; k += 10000) {

      bkg = "";
      bkg += k;
      bkg += ".root";
      myfile = new TFile(bkg);

      gingko->SetBranchAddress("sizeandbkg", &mypair);
      gingko->GetEntry(0);

      backrates[i] = (*mypair)[0];
      sizes[i] = k;
      i++;
   }
   TGraph *gr = new TGraph(24, sizes, backrates);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("Efficiency=1.0; Bkgrate VS Memory in Bytes");
   gr->GetXaxis()->SetTitle("Background Rates");
   gr->GetYaxis()->SetTitle("Bytes of Memory");
   gr->Draw("ACP");
}
