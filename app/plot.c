{
#include <vector>

    TFile *file;
    std::vector<Double_t> *hi = new std::vector<Double_t>();
    Int_t n = 26;
    Double_t x[26], y[26];
    Int_t i = 0;
    Double_t N = 1000;
    file = new TFile("94.root");
    tree94->SetBranchAddress("vec94", &hi);

    for (Int_t j = 0; j < 10; j++)
    {
        tree94->GetEntry(0);
        x[i] = j * 100;
        y[i] = (*hi)[i] / N;
        i++;
    }
    for (Int_t j = 1; j < 10; j++)
    {
        tree94->GetEntry(0);
        x[i] = j * 1000;
        y[i] = (*hi)[i] / N;
        i++;
    }
    for (Int_t j = 1; j < 8; j++)
    {
        tree94->GetEntry(0);
        x[i] = j * 10000;
        y[i] = (*hi)[i] / N;
        i++;
    }
    TGraph *gr1 = new TGraph(n, x, y);
    gr1->SetTitle("x_road 8");
    gr1->GetXaxis()->SetTitle("Background Rate");
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->SetTitle("Trigger Efficiency");
    gr1->GetYaxis()->CenterTitle();
    gr1->GetXaxis()->SetRangeUser(0, 80000);
    gr1->GetYaxis()->SetRangeUser(0.7, 1);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(4);
    gr1->SetLineColor(4);

    gr1->Draw("ALP");
    auto legend = new TLegend(0.1, 0.3, 0.48, 0.2);
    legend->SetHeader("Detector efficiency", "C"); // option "C" allows to center the header
    legend->AddEntry(gr1, "E_D = 0.94", "l");
    legend->Draw();
}
