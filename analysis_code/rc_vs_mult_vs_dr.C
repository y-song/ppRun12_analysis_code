void rc_vs_mult_vs_dr()
{
    const double pt_min = 10.;
    const double pt_max = 100.;
    
    // set up log binning
    const int nbin = 10;
    const double dr_min = -2.;
    const double dr_max = 0.;
    double dr_interval = (dr_max - dr_min) / nbin;
    double dr_bin[nbin + 1] = {};
    for (int ibin = 0; ibin < nbin + 1; ibin++)
    {
        dr_bin[ibin] = pow(10, dr_min + ibin * dr_interval);
    }

    const int n_mult_bin = 3;
    const double mult_min = 1.5;
    const double mult_max = 13.5;
    double mult_bin[n_mult_bin + 1] = {1.5, 4.5, 9.5, 13.5}; // 234, 5678, 9/10/11/12/13

    // set up histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(11);

    TH1F *h_dr = new TH1F("h_dr", ";#DeltaR", nbin, dr_bin);
    TH1F *h_mult = new TH1F("h_mult", ";jet multiplicity", n_mult_bin, mult_bin);
    TFile *f = new TFile("Results/pythia/pythia_pid.root");
    TTree *t = (TTree *)f->Get("ResultTree");

    Double_t weight;
    Int_t njets;
    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("njets", &njets);

    Double_t pt[100], pt2[100], epair[100], z[100], dr[100];
    Int_t b[100], n[100], pid1[100], pid2[100];
    t->SetBranchAddress("pt", pt);
    t->SetBranchAddress("pt2", pt2);
    t->SetBranchAddress("pid1", pid1);
    t->SetBranchAddress("pid2", pid2);
    t->SetBranchAddress("b", b);
    t->SetBranchAddress("epair", epair);
    t->SetBranchAddress("z", z);
    t->SetBranchAddress("dr", dr);
    t->SetBranchAddress("n", n);

    int nev = t->GetEntries();
    cout << "number of events: " << nev << endl;

    double x[nbin] = {};
    double x_err[nbin] = {};
    double rc[n_mult_bin][nbin] = {};
    double rc_err[n_mult_bin][nbin];
    double nss[n_mult_bin][nbin] = {};
    double nos[n_mult_bin][nbin] = {};
    int n_unweighted[n_mult_bin][nbin] = {};
    
    for (int iev = 0; iev < nev; iev++)
    {
        t->GetEntry(iev);
        for (int ijet = 0; ijet < njets; ijet++)
        {
            if (pt[ijet] < pt_min || pt[ijet] > pt_max) // jet pT cut
                continue;
            if (b[ijet] == 0 || b[ijet] == -9)
                continue;
            if (pt2[ijet] < 0) // track pT cut
                continue;
            //if (abs(pid1[ijet]*pid2[ijet]) != 4892944) // pair PID check: 44521, 103041, 4892944
            //    continue;
            //Double_t tau = 1 / (2 * epair[ijet] * z[ijet] * (1 - z[ijet]) * (1 - cos(dr[ijet])));

            int i_mult_bin = h_mult->FindFixBin(n[ijet]*1.0);
            int i_dr_bin = h_dr->FindFixBin(dr[ijet]);

            if (i_dr_bin == 0 || i_dr_bin == nbin + 1 || i_mult_bin == 0 || i_mult_bin == n_mult_bin + 1)
                continue;
            if (b[ijet] == 1)
                nss[i_mult_bin - 1][i_dr_bin - 1] += weight;
            else
                nos[i_mult_bin - 1][i_dr_bin - 1] += weight;
            n_unweighted[i_mult_bin - 1][i_dr_bin - 1] += 1;
        }
    }

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        x[ibin] = h_dr->GetBinCenter(ibin + 1);
        x_err[ibin] = h_dr->GetBinWidth(ibin + 1) / 2;
        for (int i_mult_bin = 0; i_mult_bin < n_mult_bin; i_mult_bin++)
        {
            rc[i_mult_bin][ibin] = (nss[i_mult_bin][ibin] - nos[i_mult_bin][ibin]) / (nss[i_mult_bin][ibin] + nos[i_mult_bin][ibin]);
            rc_err[i_mult_bin][ibin] = 2 * sqrt(nss[i_mult_bin][ibin] * nos[i_mult_bin][ibin] / (pow((nss[i_mult_bin][ibin] + nos[i_mult_bin][ibin]), 2) * n_unweighted[i_mult_bin][ibin]));
            cout << rc[i_mult_bin][ibin] << endl;
        }
        //cout << "\n" << endl;
    }

    auto g0 = new TGraphErrors(nbin, x, rc[0], x_err, rc_err[0]);
    auto g1 = new TGraphErrors(nbin, x, rc[1], x_err, rc_err[1]);
    auto g2 = new TGraphErrors(nbin, x, rc[2], x_err, rc_err[2]);

    TCanvas *c = new TCanvas("c", "c");
    c->SetLogx(1);
    /*c->SetLogy(1);

    Double_t integral = h_t->Integral();
    h_t->Scale(1. / integral);
    h_t->Draw();*/

    g0->Draw();
    g1->Draw("same");
    g2->Draw("same");

    g0->GetXaxis()->SetTitle("#DeltaR");
    g0->GetYaxis()->SetTitle("r_{c}");
    g0->GetYaxis()->SetRangeUser(-1., 0.3);
    g0->SetTitle("r_{c} vs #DeltaR for different multiplicity ranges");
    g1->SetLineColor(4);
    g2->SetLineColor(2);

    c->SaveAs("plots/pythia_rc_mult_dr.png", "png");
}