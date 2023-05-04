void pid_check()
{
    // create pT bins
    const int nbin = 16;
    const double pt_min = 2.0;
    const double pt_max = 10.0;
    const double pt_interval = 0.5;
    double pt_bin[nbin + 1];
    for (int ibin = 0; ibin < nbin + 1; ibin++)
    {
        pt_bin[ibin] = pt_min + ibin * pt_interval;
    }

    // set up histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(11);

    TObjArray partpid;

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        partpid.Add(new TH1F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";particle PID;fraction", 2500, 0, 2500));
    }

    TFile *f = new TFile("Results/pythia/pythia.root");
    TTree *t = (TTree *)f->Get("ResultTree");

    Double_t weight;
    Int_t njets;
    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("njets", &njets);

    Double_t pt[100], pt1[100], pt2[100];
    Int_t b[100], pid1[100], pid2[100];
    t->SetBranchAddress("pt", pt);
    t->SetBranchAddress("b", b);
    t->SetBranchAddress("pt1", pt1);
    t->SetBranchAddress("pt2", pt2);
    t->SetBranchAddress("pid1", pid1);
    t->SetBranchAddress("pid2", pid2);

    int nev = t->GetEntries();
    cout << "number of events: " << nev << endl;

    for (int iev = 0; iev < nev; iev++)
    {
        t->GetEntry(iev);
        for (int ijet = 0; ijet < njets; ijet++)
        {
            if (b[ijet] == 0 || b[ijet] == -9)
                continue;
            for (int ibin = 0; ibin < nbin; ibin++)
            {
                if (pt_bin[ibin] < pt1[ijet] && pt1[ijet] < pt_bin[ibin + 1])
                {
                    ((TH1F *)partpid.At(ibin))->Fill(abs(pid1[ijet]), weight);
                }
                if (pt_bin[ibin] < pt2[ijet] && pt2[ijet] < pt_bin[ibin + 1])
                {
                    ((TH1F *)partpid.At(ibin))->Fill(abs(pid2[ijet]), weight);
                }
            }
        }
    }

    TCanvas *c = new TCanvas("c", "c");
    for (int ibin = 0; ibin < nbin; ibin++)
    {
        auto dist = ((TH1F *)partpid.At(ibin));
        Double_t integral = dist->Integral();
        dist->Scale(1. / integral);
        dist->Draw();
        dist->GetYaxis()->SetRangeUser(0., 1.);
        c->SaveAs(Form("plots/pid_%0.1f_pT_%0.1f.png", pt_bin[ibin], pt_bin[ibin + 1]), "png");
    }
}