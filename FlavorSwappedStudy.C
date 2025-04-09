
void FlavorSwappedStudy() {

    TFile *fm = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg14-minervame1D1M1NWeightedAve.root", "READ");
    TFile *fe = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg12-minervame1D1M1NWeightedAve.root", "READ");

    TFile *fmbar = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg-14-minervame1D1M1NWeightedAve.root", "READ");
    TFile *febar = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg-12-minervame1D1M1NWeightedAve.root", "READ");


    PlotUtils::MnvH1D* hm = (PlotUtils::MnvH1D*)fm->Get("flux_E_unweighted");
    PlotUtils::MnvH1D* he = (PlotUtils::MnvH1D*)fe->Get("flux_E_unweighted");

    PlotUtils::MnvH1D* hmbar = (PlotUtils::MnvH1D*)fmbar->Get("flux_E_unweighted");
    PlotUtils::MnvH1D* hebar = (PlotUtils::MnvH1D*)febar->Get("flux_E_unweighted");

    for(int i=0; i<hm->GetNbinsX()+1; i++){
        std::cout << hm->GetXaxis()->GetBinLowEdge(i+1) << "\n";
    }
    //std::cout << hm->GetNbinsX();
/*
    he->Divide(he, hm, 1., 1.);
    hebar->Divide(hebar, hmbar, 1., 1.);

    TCanvas *c = new TCanvas("c","",800,600);
    c->SetLogx();
    //he->SetAxisRange()
    he->SetTitle("FHC nue/numu");
    he->Draw();
    c->SaveAs("fluxRatio_FHC.png");
    hebar->SetTitle("FHC anti nue/anti numu");
    hebar->Draw();
    c->SaveAs("fluxRatio_FHCbar.png");
*/
}