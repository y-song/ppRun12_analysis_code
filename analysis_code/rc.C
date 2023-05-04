void rc()
{
    // set up log binning
    const int nbin = 10;
    const double t_min = 0.;
    const double t_max = 4.;
    double t_interval = (t_max - t_min)/nbin;
    double t_bin [nbin + 1] = {};
    for (int ibin = 0; ibin < nbin + 1; ibin++)
    {
        t_bin[ibin] = pow(10, t_min + ibin * t_interval);
    }    

    // set up histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(11);

    TH1F *h_t = new TH1F("h_t", ";t [GeV^{-1}]", nbin, t_bin);

    TFile *f = new TFile("Results/pythia/pythia_pid.root");
    TTree *t = (TTree *)f->Get("ResultTree");

    Double_t weight;
    Int_t njets;
    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("njets", &njets);

    Double_t pt[100], pt2[100], epair[100], z[100], dr[100];
    Int_t b[100];
    t->SetBranchAddress("pt", pt);
    t->SetBranchAddress("pt2", pt2);
    t->SetBranchAddress("b", b);
    t->SetBranchAddress("epair", epair);
    t->SetBranchAddress("z", z);
    t->SetBranchAddress("dr", dr);

    int nev = t->GetEntries();
    cout << "number of events: " << nev << endl;

    vector<Double_t> tau_vector;
    vector<Int_t> b_vector;
    vector<Double_t> w_vector;

    for (int iev = 0; iev < nev; iev++)
    {
        t->GetEntry(iev);
        for (int ijet = 0; ijet < njets; ijet++)
        {
            if (pt[ijet] < 15 || pt[ijet] > 20) // jet pT cut
                continue;
            if (b[ijet] == 0 || b[ijet] == -9)
                continue;
            if (pt2[ijet] < 0) // track pT cut
                continue;
            Double_t tau = 1 / ( 2 * epair[ijet] * z[ijet] * (1-z[ijet]) * (1-cos(dr[ijet]))); 
            tau_vector.push_back(tau);
            b_vector.push_back(b[ijet]);
            w_vector.push_back(weight);

            h_t->Fill(tau, weight);
        }
    }

    double tau [nbin] = {};
    double tau_err [nbin] = {};
    double rc [nbin] = {};
    double rc_err [nbin] = {0,0,0,0,0,0,0,0,0,0};
    double nss [nbin] = {};
    double nos [nbin] = {};

    for (int ijet = 0; ijet < tau_vector.size(); ijet++)
    {
        int ibin = h_t->FindFixBin(tau_vector.at(ijet));
        if (ibin == 0 || ibin == nbin + 1)
            continue;
        if (b_vector.at(ijet) == 1)
            nss[ibin-1] += w_vector.at(ijet);
        else
            nos[ibin-1] += w_vector.at(ijet);
    }

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        tau[ibin] = h_t->GetBinCenter(ibin+1);
        tau_err[ibin] = h_t->GetBinWidth(ibin+1)/2;
        rc[ibin] = (nss[ibin] - nos[ibin])/(nss[ibin] + nos[ibin]);
        cout << rc[ibin] << endl;
    }    

    auto g = new TGraphErrors(nbin, tau, rc, tau_err, rc_err);

    TCanvas *c = new TCanvas("c", "c");
    c->SetLogx(1);
    /*c->SetLogy(1);
    
    Double_t integral = h_t->Integral();
    h_t->Scale(1. / integral);
    h_t->Draw();*/

    g->Draw();
    g->GetXaxis()->SetTitle("t[GeV^{-1}]");
    g->GetYaxis()->SetTitle("r_{c}");

    c->SaveAs("plots/rc.png", "png");
}
