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
    double charge_data  = 447.088; // mC
    double hms_eff_data = 0.9985;
    double charge_dummy  = 457.099; // mC
    double hms_eff_dummy = 0.9988;

    //Cuts
    TCut mc_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8) && (stop_id==0))";
    TCut mc_norm_cuts = Form("weight * %f", normfac);
    TCut data_delta_cuts = "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && H.cal.etottracknorm>0.7 && H.cer.npeSum>2.0)";

    //Project()
    simTree->Project("hSimDelta", "hsdelta", mc_delta_cuts * mc_norm_cuts);
    dataTree->Project("hDataDelta", "H.gtr.dp", data_delta_cuts * "1.0/(18.585*0.9995*0.9985)");
    //dataTree->Project("hDataDelta", "H.gtr.dp", data_delta_cuts * "1.0/(447.088 * 0.9985)");
    //dataTree->Project("hDataDelta", "H.gtr.dp", data_delta_cuts);
    dummyTree->Project("hDummyDelta", "H.gtr.dp", data_delta_cuts * "2.0/(28.597)");
    //dummyTree->Project("hDummyDelta", "H.gtr.dp", data_delta_cuts * "1.0/(457.099*0.9988)");
    //dummyTree->Project("hDummyDelta", "H.gtr.dp", data_delta_cuts);

    // Scale Data
    //hDataDelta->Scale(1.0/(charge_data * hms_eff_data));
    //hDummyDelta->Scale(1.0/(charge_dummy * hms_eff_dummy));

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
    pad_u->SetBottomMargin(0.02);       // no x-axis labels

    TPad *pad_d = (TPad*)c1->cd(2);
    pad_d->SetPad(0.0, 0.0, 1.0, 0.3);  // bottom 30% for ratio
    pad_d->SetTopMargin(0.05);
    pad_d->SetBottomMargin(0.3);        // leave space for x-axis labels


    pad_u->cd();
    hSimDelta->SetMaximum(ymax);    // Set the y-axis range of the first-drawn histogram
    hSimDelta->GetXaxis()->SetRangeUser(-35, 35);
    hSimDelta->SetLineColor(kBlue);
    hSimDelta->Draw("HIST");
    hSimDelta->GetXaxis()->SetLabelSize(0);
    hSimDelta->GetXaxis()->SetTitleSize(0);
    hDataSubDummyDelta->GetXaxis()->SetRangeUser(-35, 35);
    hDataSubDummyDelta->GetXaxis()->SetLabelSize(0);
    hDataSubDummyDelta->GetXaxis()->SetTitleSize(0);
    hDataSubDummyDelta->SetLineColor(kRed);
    hDataSubDummyDelta->Draw("HIST SAME");

    auto legend = new TLegend(0.6,0.7,0.8,0.8);
    legend->AddEntry(hSimDelta,"Sim","l");
    //legend->AddEntry(hDataDelta,"Data","l");
    legend->AddEntry(hDataSubDummyDelta,"Data-Dummy","l");
    legend->Draw();




    pad_d->cd();

    // Clone to avoid modifying original
    TH1D *hRatio = (TH1D*)hSimDelta->Clone("hRatio");
    hRatio->Divide(hDataSubDummyDelta);  // Sim / (Data - Dummy)
    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(0.5);

    hRatio->GetXaxis()->SetRangeUser(-12, 12);
    hRatio->SetLineColor(kBlack);
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Ratio");
    hRatio->GetYaxis()->SetTitleSize(0.08);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetNdivisions(500);
    hRatio->GetXaxis()->SetLabelSize(0.08);
    hRatio->GetXaxis()->SetTitleSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetTitle("Delta");
    hRatio->SetStats(0);

    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);
    hRatio->Draw("E1");  // Error bars


    c1->cd();
    c1->Update();
    //c1->SaveAs("SimVsData_Ratio.pdf");

}
