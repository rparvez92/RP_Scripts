#include <iostream>
#include <TStopwatch.h>

void DataVsSimPlot() {

    // Declare the stopwatch variable.
    TStopwatch timer;

    // ==== FILES ====
    //TFile *fSim  = new TFile("hms_29p05deg_1p531_lh2_rsidis.root");
    //TFile *fSim  = new TFile("/w/hallc-scshelf2102/c-rsidis/rparvez/single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    //TFile *fData = new TFile("../hallc_replay_rsidis/Rsidis_ROOTfiles/hms_coin_replay_production_24118_-1.root");

    TFile *fSim  = new TFile("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    TFile *fData = new TFile("./Rsidis_ROOTfiles/hms_coin_replay_production_24118_-1.root");


    if (!fSim || !fSim->IsOpen()) { 
        Error("DataVsSimPlot_general","Simulation file not found!"); return; 
    }
    if (!fData || !fData->IsOpen()) { 
        Error("DataVsSimPlot_general","Data file not found!"); return; 
    }

    std::cout << "Found the files" << std::endl;

    TTree *simTree  = (TTree*) fSim->Get("h10");
    TTree *dataTree = (TTree*) fData->Get("T");

    // ==== SIMULATION VARIABLES ====
    double hsdelta, hsytar, hsxptar, hsyptar, hsxpfp, hsypfp, sweight, stop_id;

    // Turn OFF all the branches. Otherwise running the code will take too long time
    simTree->SetBranchStatus("*", 0);
    // Turn ON only those branches and then set branch address.
    simTree->SetBranchStatus("hsdelta", 1);
    simTree->SetBranchAddress("hsdelta", &hsdelta);
    simTree->SetBranchStatus("stop_id", 1);
    simTree->SetBranchAddress("stop_id", &stop_id);
    simTree->SetBranchStatus("hsytar", 1);
    simTree->SetBranchAddress("hsytar", &hsytar);
    simTree->SetBranchStatus("hsxptar", 1);
    simTree->SetBranchAddress("hsxptar", &hsxptar);
    simTree->SetBranchStatus("hsyptar", 1);
    simTree->SetBranchAddress("hsyptar", &hsyptar);
    simTree->SetBranchStatus("hsxpfp", 1);
    simTree->SetBranchAddress("hsxpfp", &hsxpfp);
    simTree->SetBranchStatus("hsypfp", 1);
    simTree->SetBranchAddress("hsypfp", &hsypfp);
    simTree->SetBranchStatus("weight", 1);
    simTree->SetBranchAddress("weight", &sweight);

    // ==== DATA VARIABLES ====
    double H_gtr_dp, H_gtr_x, H_gtr_y, H_gtr_th, H_gtr_ph;
    double H_cal_etottracknorm, H_dc_x_fp, H_dc_y_fp, H_dc_xp_fp, H_dc_yp_fp;
    double H_kin_Q2, H_cer_npeSum;


    // Turn OFF all the branches. Otherwise running the code will take too long time
    dataTree->SetBranchStatus("*", 0);
    // Turn ON only those branches and then set branch address.
    dataTree->SetBranchStatus("H.gtr.dp", 1);
    dataTree->SetBranchAddress("H.gtr.dp", &H_gtr_dp);
    dataTree->SetBranchStatus("H.gtr.x", 1);
    dataTree->SetBranchAddress("H.gtr.x", &H_gtr_x);
    dataTree->SetBranchStatus("H.gtr.y", 1);
    dataTree->SetBranchAddress("H.gtr.y", &H_gtr_y);
    dataTree->SetBranchStatus("H.gtr.th", 1);
    dataTree->SetBranchAddress("H.gtr.th", &H_gtr_th);
    dataTree->SetBranchStatus("H.gtr.ph", 1);
    dataTree->SetBranchAddress("H.gtr.ph", &H_gtr_ph);
    dataTree->SetBranchStatus("H.cal.etottracknorm", 1);
    dataTree->SetBranchAddress("H.cal.etottracknorm", &H_cal_etottracknorm);
    dataTree->SetBranchStatus("H.dc.x_fp", 1);
    dataTree->SetBranchAddress("H.dc.x_fp", &H_dc_x_fp);
    dataTree->SetBranchStatus("H.dc.y_fp", 1);
    dataTree->SetBranchAddress("H.dc.y_fp", &H_dc_y_fp);
    dataTree->SetBranchStatus("H.dc.xp_fp", 1);
    dataTree->SetBranchAddress("H.dc.xp_fp", &H_dc_xp_fp);
    dataTree->SetBranchStatus("H.dc.yp_fp", 1);
    dataTree->SetBranchAddress("H.dc.yp_fp", &H_dc_yp_fp);
    dataTree->SetBranchStatus("H.kin.Q2", 1);
    dataTree->SetBranchAddress("H.kin.Q2", &H_kin_Q2);
    dataTree->SetBranchStatus("H.cer.npeSum", 1);
    dataTree->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);

    // ==== NORMALIZATION FACTORS ====
    double charge  = 15.210096; // mC
    double normfac = 10.96;      // sim normalization
    double hms_eff = 0.9986;

    // ==== HISTOGRAMS ====
    int nbins = 100;
    double xmin = -35.0, xmax = 35.0;
    TH1D *hSimDelta  = new TH1D("hSimDelta",  "HMS delta;Delta;Counts", nbins, xmin, xmax);
    TH1D *hDataDelta = new TH1D("hDataDelta", "HMS delta;Delta;Counts", nbins, xmin, xmax);

    TH1D *hSimXfp  = new TH1D("hSimXfp",  "HMS Xfp;Xfp;Counts", nbins, xmin, xmax);
    TH1D *hDataXfp  = new TH1D("hDataXfp",  "HMS Xfp;Xfp;Counts", nbins, xmin, xmax);

    // If any value falls outside current bin range, extend the axes.
    hSimDelta->SetCanExtend(TH1::kAllAxes);
    hDataDelta->SetCanExtend(TH1::kAllAxes);

    hSimXfp->SetCanExtend(TH1::kAllAxes);
    hDataXfp->SetCanExtend(TH1::kAllAxes);

    // ==== LOOP OVER SIMULATION ====
    Long64_t nSim = simTree->GetEntries();
    std::cout << "Total number of Sim entries: " << nSim << std::endl;

    for (Long64_t i=0; i<nSim; i++) {
        simTree->GetEntry(i);

        // CUTS
        if (hsdelta > -8.0 && hsdelta < 8.0 && stop_id==0) {
            hSimDelta->Fill(hsdelta, normfac * sweight);  // weighted filling
        }

        // CUTS
        if (hsxpfp > -35.0 && hsxpfp < 35.0 && stop_id==0) {
            hSimXfp->Fill(hsxpfp, normfac * sweight);  // weighted filling
        }

    }

    // Number of filled Sim entries
    std::cout << "Number of filled Sim delta entries: " << hSimDelta -> GetEntries() << std::endl; 
    std::cout << "Number of filled Sim xfp entries: " << hSimXfp -> GetEntries() << std::endl; 


    // ==== LOOP OVER DATA ====
    Long64_t nData = dataTree->GetEntries();
    std::cout << "Total number of Data entries: " << nData << std::endl;

    for (Long64_t i=0; i<nData; i++) {
        dataTree->GetEntry(i);

        // CUTS
        if (H_gtr_dp > -8.0 && H_gtr_dp < 8.0 &&
            H_cal_etottracknorm > 0.7 &&
            H_cer_npeSum > 2) {
            hDataDelta->Fill(H_gtr_dp);
        }

        // CUTS
        if (H_dc_xp_fp > -35.0 && H_dc_xp_fp < 35.0 &&
            H_cal_etottracknorm > 0.7 &&
            H_cer_npeSum > 2) {
            hDataXfp->Fill(H_dc_xp_fp);
        }

    }

    // Number of filled Data entries 
    std::cout << "Number of filled Data delta entries: " << hDataDelta -> GetEntries() << std::endl; 
    std::cout << "Number of filled Data xfp entries: " << hDataXfp -> GetEntries() << std::endl;

    // Normalize data by charge*hms_eff
    hDataDelta->Scale(1.0/(charge * hms_eff));
    hDataXfp->Scale(1.0/(charge * hms_eff));

    // ==== STYLE ====
    hSimDelta->SetLineColor(kBlue+2);
    hSimDelta->SetLineWidth(2);

    hDataDelta->SetLineColor(kBlack);
    hDataDelta->SetLineWidth(2);

    hSimXfp->SetLineColor(kBlue+2);
    hSimXfp->SetLineWidth(2);

    hDataXfp->SetLineColor(kBlack);
    hDataXfp->SetLineWidth(2);




    //Stop the timer
    timer.Stop();
    std::cout << "Real Time: " << timer.RealTime() << " seconds" << std::endl;
    std::cout << "CPU Time: " << timer.CpuTime() << " seconds" << std::endl;


    // ==== DRAW ====
    TCanvas *c1 = new TCanvas("c1","Data vs Simulation: Delta",900,700);
    hSimDelta->Draw("hist");
    hDataDelta->Draw("hist same");

    auto legend = new TLegend(0.6,0.7,0.85,0.85);
    legend->AddEntry(hSimDelta,"Sim","l");
    legend->AddEntry(hDataDelta,"Data","l");
    legend->Draw();
    // For log scale
    // c1->SetLogy();

    // For saving to PNG
    // c1->SaveAs("Data_vs_Sim.png");



    TCanvas *c2 = new TCanvas("c2","Data vs Simulation: Xfp",900,700);
    hSimXfp->Draw("hist");
    hDataXfp->Draw("hist same");

    auto legend1 = new TLegend(0.6,0.7,0.85,0.85);
    legend1->AddEntry(hSimXfp,"Sim","l");
    legend1->AddEntry(hDataXfp,"Data","l");
    legend1->Draw();

    // For log scale
    // c2->SetLogy();

    // For saving to PNG
    // c2->SaveAs("Data_vs_Sim.png");

    c2->Update(); // Show the plot interactively




    auto legend2 = new TLegend(0.6,0.7,0.85,0.85);
    legend2->AddEntry(hSimDelta,"Sim","l");
    legend2->AddEntry(hDataDelta,"Data","l");
    legend2->Draw();

    // For log scale
    // c1->SetLogy();

    // For saving to PNG
    // c1->SaveAs("Data_vs_Sim.png");

    c1->Update(); // Show the plot interactively

}
