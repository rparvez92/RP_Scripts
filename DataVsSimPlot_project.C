#include <iostream>
#include <TStopwatch.h>

void DataVsSimPlot_project() {

    // Declare the stopwatch variable.
    TStopwatch timer;

    // ==== FILES ====
    TFile *fSim  = new TFile("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    //TFile *fData = new TFile("./Rsidis_ROOTfiles/hms_coin_replay_production_24118_-1.root");
    TFile *fData = new TFile("./Rsidis_ROOTfiles/coin_replay_production_23909_-1.root");
    TFile *fDummy = new TFile("./Rsidis_ROOTfiles/hms_coin_replay_production_23925_-1.root");
    std::cout << "Found the files" << std::endl;


    // Grab the trees
    TTree *simTree  = (TTree*) fSim->Get("h10");
    TTree *dataTree = (TTree*) fData->Get("T");
    TTree *dummyTree = (TTree*) fDummy->Get("T");

    // ==== HISTOGRAMS ====
    int nbins = 100;
    double xmin = -35.0, xmax = 35.0;
    TH1D *hSimDelta  = new TH1D("hSimDelta",  "HMS delta;Delta;Counts", nbins, xmin, xmax);
    TH1D *hDataDelta = new TH1D("hDataDelta", "HMS delta;Delta;Counts", nbins, xmin, xmax);
    TH1D *hDummyDelta = new TH1D("hDummyDelta", "HMS delta;Delta;Counts", nbins, xmin, xmax);

    TH1D *hDataSubDummyDelta = new TH1D("hDataSubDummyDelta", "HMS DataSubDummy delta;DataSubDummy Delta;Counts", nbins, xmin, xmax);



    // ==== NORMALIZATION FACTORS ====
    double normfac = 10.96;      // sim normalization
    double charge_data  = 15.210096; // mC
    double hms_eff_data = 0.9986;
    double charge_dummy  = 457.099476; // mC
    double hms_eff_dummy = 0.9988;

    //Cuts
    TCut mc_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8) && (stop_id==0))";
    TCut mc_norm_cuts = Form("weight * %f", normfac);
    TCut data_delta_cuts = "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && H.cal.etottracknorm>0.7 && H.cer.npeSum>2.0)";

    //Project()
    simTree->Project("hSimDelta", "hsdelta", mc_delta_cuts * mc_norm_cuts);
    dataTree->Project("hDataDelta", "H.gtr.dp", data_delta_cuts * "1.0/(18.585*0.9995*0.9985)");
    dummyTree->Project("hDummyDelta", "H.gtr.dp", data_delta_cuts * "2.0/(28.597)");

    //Data - Dummy
    hDataSubDummyDelta -> Add(hDataDelta, hDummyDelta, 1.0, -1.0/3.550);


    // ==== STYLE ====
    hSimDelta->SetLineColor(kBlue+2);
    hSimDelta->SetLineWidth(2);

    //hDataDelta->SetLineColor(kBlack);
    //hDataDelta->SetLineWidth(2);

    hDataSubDummyDelta->SetLineColor(kRed);
    hDataSubDummyDelta->SetLineWidth(2);

    //Stop the timer
    timer.Stop();
    std::cout << "Real Time: " << timer.RealTime() << " seconds" << std::endl;
    std::cout << "CPU Time: " << timer.CpuTime() << " seconds" << std::endl;


    // ==== DRAW ====
    TCanvas *c1 = new TCanvas("c1","Data-Dummy vs Simulation: Delta",900,700);
    hSimDelta->Draw("hist");
    //hDataDelta->Draw("hist same");
    hDataSubDummyDelta->Draw("hist same");

    auto legend = new TLegend(0.6,0.7,0.85,0.85);
    legend->AddEntry(hSimDelta,"Sim","l");
    //legend->AddEntry(hDataDelta,"Data","l");
    legend->AddEntry(hDataSubDummyDelta,"Data-Dummy","l");
    legend->Draw();

    c1->Update(); // Show the plot interactively
}
