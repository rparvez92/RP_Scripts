#include <iostream>
#include <TStopwatch.h>

void DataVsSimPlot_ratio() {

    // Declare the stopwatch variable.
    TStopwatch timer;

    // ==== FILES ====
    TFile *fSim  = new TFile("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    //TFile *fData = new TFile("./Rsidis_ROOTfiles/hms_coin_replay_production_24118_-1.root");
    TFile *fData = new TFile("./Rsidis_ROOTfiles/coin_replay_production_23909_-1.root");
    TFile *fDummy = new TFile("./Rsidis_ROOTfiles/hms_coin_replay_production_23925_-1.root");
    //TFile *fData = new TFile("./Rsidis_ROOTfiles/hms_coin_replay_production_24544_-1.root");
    //TFile *fDummy = new TFile("./Rsidis_ROOTfiles/hms_coin_replay_production_24546_-1.root");
    std::cout << "Found the files" << std::endl;


    // Grab the trees
    TTree *simTree  = (TTree*) fSim->Get("h10");
    TTree *dataTree = (TTree*) fData->Get("T");
    TTree *dummyTree = (TTree*) fDummy->Get("T");

    // ==== HISTOGRAMS ====
    int nbins = 100;
    double xmin = -9.0, xmax = 9.0;
    TH1D *hSimDelta  = new TH1D("hSimDelta",  "HMS delta;Delta;Counts", nbins, xmin, xmax);
    TH1D *hDataDelta = new TH1D("hDataDelta", "HMS delta;Delta;Counts", nbins, xmin, xmax);
    TH1D *hDummyDelta = new TH1D("hDummyDelta", "HMS delta;Delta;Counts", nbins, xmin, xmax);

    TH1D *hDataSubDummyDelta = new TH1D("hDataSubDummyDelta", "HMS DataSubDummy delta;DataSubDummy Delta;Counts", nbins, xmin, xmax);



    // ==== NORMALIZATION FACTORS ====
    double normfac = 10.96;      // sim normalization
    double charge_data  = 18.585; //447.088; // mC
    //double charge_data  = 49.262; // for HMSDIS run 24544 //mC
    double hms_eff_data = 0.9957;
    double charge_dummy  = 28.597; //457.099; // mC
    //double charge_dummy  = 28.948; // for HMSDIS run 24546 // mC
    double hms_eff_dummy = 0.9964;

    //Cuts
    TCut mc_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8) && (stop_id==0))";
    TCut mc_norm_cuts = Form("weight * %f", normfac);
    TCut data_delta_cuts = "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && H.cal.etottracknorm>0.7 && H.cer.npeSum>2.0)";
    TCut data_scale = Form("1.0 / (%f * %f)", charge_data, hms_eff_data);
    TCut dummy_scale = Form("2.0 / (%f * %f)", charge_dummy, hms_eff_dummy);


    //Project()
    simTree->Project("hSimDelta", "hsdelta", mc_delta_cuts * mc_norm_cuts);
    //simTree->Project("hSimDelta", "hsytar", mc_delta_cuts * mc_norm_cuts);
    dataTree->Project("hDataDelta", "H.gtr.dp", data_delta_cuts * data_scale);
    //dataTree->Project("hDataDelta", "H.gtr.y", data_delta_cuts * data_scale);
    //dataTree->Project("hDataDelta", "H.gtr.dp", data_delta_cuts * "1.0/(18.585*0.9995*0.9985)");
    dummyTree->Project("hDummyDelta", "H.gtr.dp", data_delta_cuts * dummy_scale);
    //dummyTree->Project("hDummyDelta", "H.gtr.y", data_delta_cuts * dummy_scale);
    //dummyTree->Project("hDummyDelta", "H.gtr.dp", data_delta_cuts * "2.0/(28.597)");


    //Data - Dummy
    hDataSubDummyDelta -> Add(hDataDelta, hDummyDelta, 1.0, -1.0/3.550);
    //hDataSubDummyDelta -> Add(hDataDelta, hDummyDelta, 1.0, -1.0/3.550);


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


    // Get the max values from both histograms
    double maxSim   = hSimDelta->GetMaximum();
    double maxDataSubDummy  = hDataSubDummyDelta->GetMaximum();

    // Find the higher of the two, and scale up slightly for margin
    double ymax = std::max(maxSim, maxDataSubDummy) * 1.1;  // add 10% margin


    // ==== DRAW ====
    TCanvas *c1 = new TCanvas("c1","Simulation Vs Data-Dummy: Delta",900,700);
    c1->Divide(1,2); // Divide into 2 vertical pads

    TPad *pad_u = (TPad*)c1->cd(1);
    pad_u->SetPad(0.0, 0.3, 1.0, 1.0);  // top 70% for histograms
    pad_u->SetTopMargin(0.1);
    pad_u->SetBottomMargin(0.0);
    pad_u->SetRightMargin(0.05);
    pad_u->SetLeftMargin(0.1);

    TPad *pad_d = (TPad*)c1->cd(2);
    pad_d->SetPad(0.0, 0.0, 1.0, 0.3); // bottom 30% for ratio
    pad_d->SetTopMargin(0.0);
    pad_d->SetBottomMargin(0.35); // leave space for x-axis labels
    pad_d->SetRightMargin(0.05);
    pad_d->SetLeftMargin(0.1);


    pad_u->cd();
    hSimDelta->SetMaximum(ymax); // Set the y-axis range of the first-drawn histogram
    //hSimDelta->GetXaxis()->SetRangeUser(xmin, xmax);
    hSimDelta->SetLineColor(kBlue);
    hSimDelta->GetXaxis()->SetLabelSize(0);
    hSimDelta->GetXaxis()->SetTitleSize(0);
    hSimDelta->Draw("HIST");
    //hDataSubDummyDelta->GetXaxis()->SetRangeUser(xmin, xmax);
    hDataSubDummyDelta->SetLineColor(kRed);
    hDataSubDummyDelta->GetXaxis()->SetLabelSize(0);
    hDataSubDummyDelta->GetXaxis()->SetTitleSize(0);
    hDataSubDummyDelta->Draw("HIST SAME");

    auto legend = new TLegend(0.6,0.7,0.8,0.8);
    legend->AddEntry(hSimDelta,"Sim","l");
    legend->AddEntry(hDataSubDummyDelta,"Data-Dummy","l");
    legend->Draw();




    pad_d->cd();
    TH1D *hRatio = (TH1D*)hSimDelta->Clone("hRatio"); // Clone to avoid modifying original
    hRatio->Divide(hDataSubDummyDelta); // Sim / (Data - Dummy)

    hRatio->SetLineColor(kBlack);
    hRatio->SetStats(0);
    hRatio->SetTitle("");

    hRatio->GetYaxis()->SetTitle("Ratio");
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetNdivisions(500);
    //hRatio->GetXaxis()->SetRangeUser(-12, 12);
    hRatio->GetXaxis()->SetTitle("Delta");
    hRatio->GetXaxis()->SetTitleSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetNdivisions(505);

    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);
    hRatio->Draw("E1");  // Error bars

    // Fix for invisible y-axis tick labels
    gPad->Update(); // Forces ROOT to layout the frame
    hRatio->GetYaxis()->SetLabelOffset(0.01); // Small offset
    hRatio->GetYaxis()->SetNdivisions(505); // 5 major ticks, 5 minor
    gPad->RedrawAxis(); // Redraws ticks and labels properly


    c1->cd();
    c1->Update();
    //c1->SaveAs("SimVsData_Ratio.pdf");

}
