void pid_pair_check()
{

    // set up histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(11);

    TH1F *h_pid_index_pair = new TH1F("h_pid_index_pair", ";pair PID product index", 16, -9, 7);
    TH1F *h_pid_absindex_pair = new TH1F("h_pid_absindex_pair", ";pair PID product index absolute value", 9, 0, 10);

    TFile *f = new TFile("Results/pythia/pythia_pid_new.root");
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
            if (pt[ijet] < 10) // jet pT cut
                continue;
            if (b[ijet] == 0 || b[ijet] == -9)
                continue;
            if (pt2[ijet] < 0.2) // track pT cut
                continue;
            if (pt2[ijet] > pt1[ijet])
            {
                cout << "subleading track has higher pT than leading!" << endl;
                break;
            }
            int pid_product = pid1[ijet] * pid2[ijet];
            switch (pid_product)
            {
            case 44521: // 211^2
                h_pid_index_pair->Fill(1., weight);
                h_pid_absindex_pair->Fill(1., weight);
                break;
            case -44521:
                h_pid_index_pair->Fill(-1., weight);
                h_pid_absindex_pair->Fill(1., weight);
                break;
            case 103041: // 321^2
                h_pid_index_pair->Fill(2., weight);
                h_pid_absindex_pair->Fill(2., weight);
                break;
            case -103041:
                h_pid_index_pair->Fill(-2., weight);
                h_pid_absindex_pair->Fill(2., weight);
                break;
            case 4892944: // 2212^2
                h_pid_index_pair->Fill(3., weight);
                h_pid_absindex_pair->Fill(3., weight);
                break;
            case -4892944:
                h_pid_index_pair->Fill(-3., weight);
                h_pid_absindex_pair->Fill(3., weight);
                break;
            case 67731: // 211*321
                h_pid_index_pair->Fill(4., weight);
                h_pid_absindex_pair->Fill(4., weight);
                break;
            case -67731:
                h_pid_index_pair->Fill(-4., weight);
                h_pid_absindex_pair->Fill(4., weight);
                break;
            case 466732: // 211*2212
                h_pid_index_pair->Fill(5., weight);
                h_pid_absindex_pair->Fill(5., weight);
                break;
            case -466732:
                h_pid_index_pair->Fill(-5., weight);
                h_pid_absindex_pair->Fill(5., weight);
                break;
            case 710052: // 321*2212
                h_pid_index_pair->Fill(6., weight);
                h_pid_absindex_pair->Fill(6., weight);
                break;
            case -710052:
                h_pid_index_pair->Fill(-6., weight);
                h_pid_absindex_pair->Fill(6., weight);
                break;
            default:
                // cout << pid_product << endl;
                h_pid_index_pair->Fill(-9., weight);
                h_pid_absindex_pair->Fill(9., weight);
                break;
            }
        }
    }

    TCanvas *c = new TCanvas("c", "c");
    c->SetGrid(0);

    /*Double_t integral = h_pid_index_pair->Integral();
    h_pid_index_pair->Scale(1. / integral);
    h_pid_index_pair->Draw();
    h_pid_index_pair->GetYaxis()->SetRangeUser(0., 1.);*/

    Double_t integral = h_pid_absindex_pair->Integral();
    h_pid_absindex_pair->Scale(1. / integral);
    h_pid_absindex_pair->Draw("hist");

    h_pid_absindex_pair->GetYaxis()->SetRangeUser(0., 1.);
    h_pid_absindex_pair->GetYaxis()->SetTitle("fraction");

    h_pid_absindex_pair->GetXaxis()->SetTitle("Leading & subleading track species");
    h_pid_absindex_pair->GetXaxis()->SetBinLabel(1, "#pi#pi");
    h_pid_absindex_pair->GetXaxis()->SetBinLabel(2, "KK");
    h_pid_absindex_pair->GetXaxis()->SetBinLabel(3, "pp");
    h_pid_absindex_pair->GetXaxis()->SetBinLabel(4, "#piK");
    h_pid_absindex_pair->GetXaxis()->SetBinLabel(5, "#pip");
    h_pid_absindex_pair->GetXaxis()->SetBinLabel(6, "Kp");
    h_pid_absindex_pair->GetXaxis()->SetBinLabel(9, "other");

    c->SaveAs("plots/leading_subleading_pid.png", "png");

}