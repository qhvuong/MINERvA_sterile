
void GetFileNameToText() {
/*
    std::ofstream outFile("file_option/playlist_datale1Ap6.txt");

    const char folder[] = "/pnfs/minerva/persistent/DataPreservation/Low_Energy_Era/p6Prime/data";
    
    for(int j=20;j<24;j++){
    for(int k=0;k<100;k++){
    for (int i = 1; i <= 10; i++) {
        TString filename = TString::Format("%s/minerva1/grid/minerva/ana/numibeam/v10r8p9/00/00/%d/%02d/MV_0000%d%02d_Subruns_%04d_MasterAnaDev_AnaData_Tuple_v10r8p9.root", folder, j,k,j,k,i);
        //std::cout << filename << "\n";
        if (gSystem->AccessPathName(filename)) {  // Returns true if file DOES NOT exist
            std::cout << k << "\t" << i << "\n";
        }
        else {
            outFile << Form("root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr/minerva/persistent/DataPreservation/Low_Energy_Era/p6Prime/data/minerva1/grid/minerva/ana/numibeam/v10r8p9/00/00/%d/%02d/MV_0000%d%02d_Subruns_%04d_MasterAnaDev_AnaData_Tuple_v10r8p9.root",j,k,j,k,i) << "\n";
        }
    }}}
*/

    std::ofstream outFile("file_option/playlist_mcle5Ap6.txt");

    const char xrootd[] = "root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr";
    const char dataPath[] = "minerva/persistent/DataPreservation/Low_Energy_Era/p6Prime/mc/minerva5/grid/central_value/minerva/ana/v10r8p9/00/05/02";

    for(int j=0;j<50;j++){
    for (int i = 1; i <= 200; i++) {
        TString filename = TString::Format("/pnfs/%s/%02d/SIM_minerva_000502%02d_Subruns_%04d_MasterAnaDev_Ana_Tuple_v10r8p9.root", dataPath, j,j,i);
        //std::cout << filename << "\n";
        if (gSystem->AccessPathName(filename)) {  // Returns true if file DOES NOT exist
            std::cout << j << "\t" << i << "\n";
        }
        else {
            outFile << Form("%s/%s/%02d/SIM_minerva_000502%02d_Subruns_%04d_MasterAnaDev_Ana_Tuple_v10r8p9.root", xrootd, dataPath, j,j,i) << "\n";
        }
    }}

}

/*
    TFile *fm = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg14-minervame1D1M1NWeightedAve.root", "READ");
    TFile *fe = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg12-minervame1D1M1NWeightedAve.root", "READ");

    TFile *fmbar = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg-14-minervame1D1M1NWeightedAve.root", "READ");
    TFile *febar = new TFile("$PLOTUTILSROOT/data/flux/flux-g4numiv6-pdg-12-minervame1D1M1NWeightedAve.root", "READ");


    PlotUtils::MnvH1D* hm = (PlotUtils::MnvH1D*)fm->Get("flux_E_unweighted");
    PlotUtils::MnvH1D* he = (PlotUtils::MnvH1D*)fe->Get("flux_E_unweighted");

    PlotUtils::MnvH1D* hmbar = (PlotUtils::MnvH1D*)fmbar->Get("flux_E_unweighted");
    PlotUtils::MnvH1D* hebar = (PlotUtils::MnvH1D*)febar->Get("flux_E_unweighted");

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
