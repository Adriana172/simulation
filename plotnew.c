{
#include <vector>
    
    TFile *ratiofile806;
    TFile *ratiofile1606;
    TFile *ratiofile3206;

    std::vector<Double_t> *hi = new std::vector<Double_t>();
    std::vector<Double_t> *hi16 = new std::vector<Double_t>();
    std::vector<Double_t> *hi32 = new std::vector<Double_t>();
    Int_t n = 26;
    Double_t x[26],y[26],y16[26],y32[26],ex[26];
    Int_t i = 0;
    /* new file */
    /* new file */
    // vector<double> N806 = {252, 257, 246, 244, 268, 232, 219, 237,245,249,259,252,238,239,237,235,250,251,232,244,248,241,243,251,242,253};
    // vector<double> N1606 = {253,241,244,229,259,246,256,258,250,236,248,231,248,247,247,245,260,260,246,270,229,247,251,261,255,259};
    // vector<double> N3206 = {267,258,246,258,231,261,249,239,252,261,234,245,242,233,230,264,248,234,231,235,237,243,269,251,259,234};
    vector<double> N808 =  {500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500};
    vector<double> N1608 = {500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500};
    vector<double> N3208 = {500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500};
    vector<double> Nnew1606= {991,996,1010,978,1000,1039,1038,1034,986,973,976,1006,1028,1011,957,965,992,999,995,1008,976,1017,1000,1010,956,1003}
    vector<double> Nnew3206= {1014,1000,951,983,994,1016,975,1000,993,989,1025,964,971,960,987,985,991,1006,980,953,1002,1014,973,978,991,945}
    vector<double> Nnew806= {1039,992,1002,1006,1001,996,1006,1032,973,973,989,990,988,987,1030,1002,1029,1008,997,1002,969}
    ratiofile808 = new TFile("ratio808.root");
    tree808->SetBranchAddress("vec",&hi);
    
    for (Int_t j = 0; j < 10; j++) {
        tree808->GetEntry(0);
        x[i] = j*100;
        y[i] = (*hi)[i]/N808[i];
        i++;
    }
    for (Int_t j = 1; j < 10; j++) {
        tree808->GetEntry(0);
        x[i] = j*1000;
        y[i] = (*hi)[i]/N808[i];
        i++;
    }
    for (Int_t j = 1; j < 8; j++) {
        tree808->GetEntry(0);
        x[i] = j*10000;
        y[i] = (*hi)[i]/N808[i];
        i++;
    }
    i=0;
    ratiofile1608 = new TFile("ratio1608.root");
    tree1608->SetBranchAddress("vec",&hi16);
    for (Int_t j = 0; j < 10; j++) {
        tree1608->GetEntry(0);
        y16[i] = (*hi16)[i]/N1608[i];
        i++;
    }
    for (Int_t j = 1; j < 10; j++) {
        tree1608->GetEntry(0);
        y16[i] = (*hi16)[i]/N1608[i];
        i++;
    }
    for (Int_t j = 1; j < 8; j++) {
        tree1608->GetEntry(0);
        y16[i] = (*hi16)[i]/N1608[i];
        i++;
    }

    i=0;
    ratiofile3208 = new TFile("ratio3208.root");
    tree3208->SetBranchAddress("vec",&hi32);
    
    for (Int_t j = 0; j < 10; j++) {
        tree3208->GetEntry(0);
        y32[i] = (*hi32)[i]/N3208[i];
        i++;
    }
    for (Int_t j = 1; j < 10; j++) {
        tree3208->GetEntry(0);
        y32[i] = (*hi32)[i]/N3208[i];
        i++;
    }
    for (Int_t j = 1; j < 8; j++) {
        tree3208->GetEntry(0);
        y32[i] = (*hi32)[i]/N3208[i];
        i++;
    }

    TGraph *gr1 = new TGraphErrors (n,x,y,ex);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(4);
    gr1->SetLineColor(4);
    gr1->Draw("ALP");
    gr1->GetYaxis()->SetLimits(0,1);
    gr1->SetMinimum(0.);
    gr1->SetMaximum(1.);

    TGraph *gr2 = new TGraphErrors (n,x,y16,ex);
    gr2->SetLineWidth(3);
    gr2->SetMarkerStyle(21);
    gr2->SetLineColor(2);
    gr2->Draw("CP");

    TGraph *gr3 = new TGraphErrors (n,x,y32,ex);
    gr3->SetLineWidth(2);
    gr3->SetMarkerStyle(20);
    gr3->SetLineColor(3);
    gr3->Draw("CP");
}
