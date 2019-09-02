{
#include <vector>

    TFile *myfile;
    Double_t xReal;
    Double_t yReal;
    TString bkg;
    Int_t counter;
    Int_t nentries;
    std::vector<Double_t> *r94 = new std::vector<Double_t>();
    std::vector<Double_t> *xTrig = new std::vector<Double_t>();
    std::vector<Double_t> *yTrig = new std::vector<Double_t>();
    xlim = 5;
    ylim = 10;

    for (Int_t k = 0; k < 1000; k += 100)
    {

        bkg = "";
        bkg += k;
        bkg += ".root";
        myfile = new TFile(bkg);

        gingko->SetBranchAddress("real_x_muon", &xReal);
        gingko->SetBranchAddress("real_y_muon", &yReal);
        gingko->SetBranchAddress("trig_x", &xTrig);
        gingko->SetBranchAddress("trig_y", &yTrig);

        nentries = gingko->GetEntries();
        counter = 0;
        for (Int_t i = 0; i < nentries; i++)
        {
            gingko->GetEntry(i);
            for (Int_t j = 0; j < xTrig->size(); j++)
            {
                if ((abs(xReal - (*xTrig)[j]) < xlim) && (abs(yReal - (*yTrig)[j]) < ylim))
                {
                    counter++;
                    break;
                }
            }
        }
        r94->push_back(counter);
    }
    for (Int_t k = 1000; k < 10000; k += 1000)
    {

        bkg = "";
        bkg += k;
        bkg += ".root";
        myfile = new TFile(bkg);

        gingko->SetBranchAddress("real_x_muon", &xReal);
        gingko->SetBranchAddress("real_y_muon", &yReal);
        gingko->SetBranchAddress("trig_x", &xTrig);
        gingko->SetBranchAddress("trig_y", &yTrig);

        counter = 0;
        nentries = gingko->GetEntries();
        for (Int_t i = 0; i < nentries; i++)
        {
            gingko->GetEntry(i);
            for (Int_t j = 0; j < xTrig->size(); j++)
            {
                if ((abs(xReal - (*xTrig)[j]) < xlim) && (abs(yReal - (*yTrig)[j]) < ylim))
                {
                    counter++;
                    break;
                }
            }
        }
        r94->push_back(counter);
    }

    for (Int_t k = 10000; k < 80000; k += 10000)
    {

        bkg = "";
        bkg += k;
        bkg += ".root";
        myfile = new TFile(bkg);

        gingko->SetBranchAddress("real_x_muon", &xReal);
        gingko->SetBranchAddress("real_y_muon", &yReal);
        gingko->SetBranchAddress("trig_x", &xTrig);
        gingko->SetBranchAddress("trig_y", &yTrig);

        counter = 0;
        nentries = gingko->GetEntries();
        for (Int_t i = 0; i < nentries; i++)
        {
            gingko->GetEntry(i);
            for (Int_t j = 0; j < xTrig->size(); j++)
            {
                if ((abs(xReal - (*xTrig)[j]) < xlim) && (abs(yReal - (*yTrig)[j]) < ylim))
                {
                    counter++;
                    break;
                }
            }
        }
        r94->push_back(counter);
    }
    TFile *fout = new TFile("94.root", "RECREATE");
    TTree *tree94 = new TTree("tree94", "tree");
    tree94->Branch("vec94", &ratio96);
    tree94->Fill();
    fout->Write();
}

