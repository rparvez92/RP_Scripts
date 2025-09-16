// ROOT macro to automate plotting of data vs simulation comparison

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo> //For typeid function
#include "SimToDataMap.h"
#include "ReportParser.h"
#include "PlotComparisonAndRatio.h"

// Creating an anonymous namespace to store unique_ptrs in a global vector, so that the objects
// pointed by these pointers do not get destroyed at the end of function call.
namespace {
  std::vector<std::unique_ptr<TH1>> g_keep_hists;
}

// Function for returning file path
static std::string DnDRootPath(int run) {
  return Form("./Rsidis_ROOTfiles/hms_coin_replay_production_%d_-1.root", run);
}
static std::string DnDReportPath(int run) {
  return Form("./REPORT_OUTPUT/HMS/PRODUCTION/replay_hms_coin_production_%d_-1.report", run);
}


//============START BUILDING HISTOGRAMS============\\


// Create and project a normalized histogram for a SINGLE data run
static std::unique_ptr<TH1D> ProjectOneDnDRun(int run,
						const std::string& dndVar,
						int nbins,
						double xmin,
						double xmax,
						const TCut& dnd_delta_cuts,
						double& Qsum_mC) {

    // Get the data or dummy file
    std::string fpath = DnDRootPath(run);
    std::unique_ptr<TFile> fDnD(TFile::Open(fpath.c_str(), "READ"));

    // Get the data or dummy tree
    TTree* tDnD = (TTree*)fDnD->Get("T");

    // Get the values from report file
    ReportValues V = ParseReportFile(DnDReportPath(run));
    // Add the charge values from all runs
    Qsum_mC += V.charge_mC;
    // cout some run constants for debug
    cout << "dndRun " << run << ": charge = " << V.charge_mC << ", hms_eff = " << V.hms_eff << ", ps_factor = " << V.ps_factor << endl;

    // Warn if BAD values are found
    if (V.charge_mC <= 0 || V.hms_eff <= 0 || V.ps_factor <= 0) {
      std::cerr << "[WARN] Bad/zero values in report for run " << run << ". Check report file.\n";
    }

    // Create the scale for dnd
    TCut dnd_scale = Form("%d / (%f)", V.ps_factor, V.hms_eff);

    // Create a smart pointer histogram for dnd
    auto h = std::make_unique<TH1D>(Form("hDnD_run_%d_%s", run, dndVar.c_str()), "", nbins, xmin, xmax);
    // Keep the error info and project respective variable
    h->Sumw2(true);
    tDnD->Project(h->GetName(), dndVar.c_str(), dnd_delta_cuts * dnd_scale);

    // Because ROOT attaches any newly created histogram to the current directory or file,
    // when that file gets closed, ROOT will delete everything that file owned. Therefore,
    // to detach histograms from opened files we detach them from current file by using
    // SetDirectory(nullptr) and smart pointer has now its exclusive ownership.
    h->SetDirectory(nullptr);

    return h;
}


// Build a SIM histogram
static std::unique_ptr<TH1D> BuildSim(const std::string& simVar,
					TTree* tSim,
					int nbins,
					double xmin,
					double xmax,
					const TCut& sim_delta_cuts,
					const TCut& sim_norm_cuts) {

    // Create an empty histogram and project the correct branch with cuts
    auto h = std::make_unique<TH1D>(Form("hSim_%s", simVar.c_str()), "", nbins, xmin, xmax);
    h->Sumw2(true);
    tSim->Project(h->GetName(), simVar.c_str(), sim_delta_cuts * sim_norm_cuts);

    //Detach ownership from current directory
    h->SetDirectory(nullptr);

    return h;
}


// Build an averaged Data histogram for many runs
static std::unique_ptr<TH1D> BuildDataAvg(const std::vector<int>& dataRuns,
						const std::string& dndVar,
						int nbins,
						double xmin,
						double xmax,
						const TCut& dnd_delta_cuts) {

    // Create a smart pointer for averaged histogram
    std::unique_ptr<TH1D> hAvgData;
    // Variable for total charge
    double QtotData = 0.0;

    // Loop over the Data Runs
    for (int run : dataRuns) {
      // Create a histogram to store single run histogram
      auto h = ProjectOneDnDRun(run, dndVar, nbins, xmin, xmax, dnd_delta_cuts, QtotData);

      // If single run histogram can't be made, skip this run
      if (!h) {cout << "skipped this run = " << run  << endl; continue;}

      // hAvgData is empty initially, so we clone the single run histogram
      if (!hAvgData) {
        hAvgData.reset((TH1D*)h->Clone(Form("hDataAvg_%s", dndVar.c_str())));
        // Detach ownership from curent directory
        hAvgData->SetDirectory(nullptr);
      }
      // else we add the single run histograms repetatively
      else {
        hAvgData->Add(h.get(), 1.0);
      }
    }

    // Average the histogram
    if (hAvgData && QtotData > 0){
	cout << "Total Data Charge : " << QtotData << endl;
	hAvgData->Scale(1.0 / QtotData);
    }
    return hAvgData;
}


// Build an averaged Dummy histogram for many runs
static std::unique_ptr<TH1D> BuildDummyAvg(const std::vector<int>& dummyRuns,
						const std::string& dndVar,
						int nbins,
						double xmin,
						double xmax,
						const TCut& dnd_delta_cuts) {

    // Create a smart pointer for averaged histogram
    std::unique_ptr<TH1D> hAvgDummy;
    // Variable for total charge
    double QtotDummy = 0.0;

    // Loop over the Dummy Runs
    for (int run : dummyRuns) {
      // Create a histogram to store single run histogram
      auto h = ProjectOneDnDRun(run, dndVar, nbins, xmin, xmax, dnd_delta_cuts, QtotDummy);

      // If single run histogram can't be made, skip this run
      if (!h) {cout << "skipped this run = " << run  << endl; continue;}

      // hAvgDummy is empty initially, so we clone the single run histogram
      if (!hAvgDummy) {
        hAvgDummy.reset((TH1D*)h->Clone(Form("hDummyAvg_%s", dndVar.c_str())));
        // Detach ownership from curent directory
        hAvgDummy->SetDirectory(nullptr);
      }
      // else we add the single run histograms repetatively
      else {
        hAvgDummy->Add(h.get(), 1.0);
      }
    }

    // Average the histogram
    if (hAvgDummy && QtotDummy > 0){
	cout << "Total Dummy Charge : " << QtotDummy << endl;
	hAvgDummy->Scale(1.0 / QtotDummy);
    }
    return hAvgDummy;
}


//============END BUILDING HISTOGRAMS============\\


// The multi-run plotting function
void PlotVariablesMultiRuns(const std::vector<int>& dataRuns,
				const std::vector<int>& dummyRuns,
				const std::string& simVar,
				TTree* tSim,
				int nbins,
				double xmin,
				double xmax,
				double wall_thickness_ratio,
				TCut sim_delta_cuts,
				TCut sim_norm_cuts,
				TCut dnd_delta_cuts) {

    // Get the data variable name
    std::string dndVar = SimToDataMap(simVar);

    // Build histograms
    auto hSim       = BuildSim(simVar, tSim, nbins, xmin, xmax, sim_delta_cuts, sim_norm_cuts);
    auto hDataAvg   = BuildDataAvg(dataRuns, dndVar, nbins, xmin, xmax, dnd_delta_cuts);
    auto hDummyAvg  = BuildDummyAvg(dummyRuns, dndVar, nbins, xmin, xmax, dnd_delta_cuts);

    // Error if building histograms are incorrect
    if (!hSim || !hDataAvg || !hDummyAvg) {
      std::cerr << "[ERROR] Missing hist (sim/data/dummy). Aborting " << simVar << std::endl;
      return;
    }

    // DataAvg - DummyAvg / wall_thickness_ratio
    std::unique_ptr<TH1D> hDataSubDummy(new TH1D(Form("hDataSubDummy_%s", dndVar.c_str()), "", nbins, xmin, xmax)); //creates an empty histogram and wraps it in a smart pointer

    // Detach the ownership from current directory
    hDataSubDummy->SetDirectory(nullptr);

    hDataSubDummy->Sumw2(true); //keep the error information
    hDataSubDummy->Add(hDataAvg.get(), 1.0); //add the summed data histogram first
    hDataSubDummy->Add(hDummyAvg.get(), -1.0 / wall_thickness_ratio); //subtract the dummy histogram

    // Plot comparison and ratio
    PlotComparisonAndRatio(hSim.get(), hDataSubDummy.get(), simVar);


    // Push back the histos to namespace smart pointer vector
    // so that these histos can stay alive when current fuction call ends.
    g_keep_hists.push_back(std::move(hSim));
    g_keep_hists.push_back(std::move(hDataAvg));
    g_keep_hists.push_back(std::move(hDummyAvg));
    g_keep_hists.push_back(std::move(hDataSubDummy));

}


// MAIN FUNCTION
void DataVsSimPlot_MultiDataMultiDummy() {

    // Files that don't depend on run numbers, i.e. sim files
    TFile* fSim = TFile::Open("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    TTree* tSim = (TTree*) fSim->Get("h10");

    // Data runs and Dummy runs
    std::vector<int> dataRuns = {24017, 24083, 24136, 24183};
    std::vector<int> dummyRuns = {24118, 24235, 24340, 24465};
    //int dummyRun = 24118; // single dummy for now

    // Constants
    int nbins = 100;
    double xmin = -8.0, xmax = 8.0;
    double wall_thickness_ratio = 3.55;

    // Cuts and normalizations
    double normfac = 10.96;
    TCut sim_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8) && (stop_id==0))";
    TCut sim_norm_cuts  = Form("weight * %f", normfac);
    TCut dnd_delta_cuts = "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && H.cal.etottracknorm>0.7 && H.cer.npeSum>2.0)";

    // Plot each variable
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsdelta", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsytar", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsxptar", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsyptar", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "xb", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "q2", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "w", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);

}

