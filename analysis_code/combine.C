void combine(){

    TFile *fin = new TFile("Results/MatchedHadded328.root");
    TFile *fout = new TFile("Results/MatchedHadded328_pt.root", "RECREATE");
    fout->cd();

    TH1D * mcpt = (TH1D*) fin->Get("mcpt");
    TH1D * rcpt = (TH1D*) fin->Get("rcpt");
    TH1D * rcpt_matched = (TH1D*) fin->Get("rcpt_matched");

    TCanvas * c1 = new TCanvas();
    mcpt->SetTitle("pythia");
    mcpt->SetLineColor(1);
    mcpt->Draw();
    rcpt->SetTitle("embedding");
    rcpt->SetLineColor(2);
    rcpt->Draw("Same");
    rcpt_matched->SetTitle("embedding matched to pythia");
    rcpt_matched->SetLineColor(4);
    rcpt_matched->Draw("Same");
    TLegend * l1 = new TLegend();
    l1->AddEntry(mcpt);
    l1->AddEntry(rcpt);
    l1->AddEntry(rcpt_matched);
    c1->SetLogy();
    c1->Print();
   
    mcpt->Write("mcpt");
    rcpt->Write("rcpt");
    rcpt_matched->Write("rcpt_matched");

    fout->Write();
    fout->Close();

}
