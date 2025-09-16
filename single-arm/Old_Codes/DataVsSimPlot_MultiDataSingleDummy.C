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
static std::unique_ptr<TH1D> ProjectOneDataRun(int run,
						const std::string& dndVar,
						int nbins,
						double xmin,
						double xmax,
						const TCut& dnd_delta_cuts,
						double& Qsum_mC) {

    // Get the data file
    std::string fpath = DnDRootPath(run);
    std::unique_ptr<TFile> f(TFile::Open(fpath.c_str(), "READ"));

    // Get the data tree
    TTree* tData = (TTree*)f->Get("T");

    // Get the values from report file
    ReportValues V = ParseReportFile(DnDReportPath(run));
    // Add the charge values from all runs
    Qsum_mC += V.charge_mC;
    // cout some run constants for debug
    cout << "dataRun " << run << ": charge = " << V.charge_mC << ", hms_eff = " << V.hms_eff << ", ps_factor = " << V.ps_factor << endl;

    // Warn if BAD values are found
    if (V.charge_mC <= 0 || V.hms_eff <= 0 || V.ps_factor <= 0) {
      std::cerr << "[WARN] Bad/zero values in report for run " << run << ". Check report file.\n";
    }

    // Create the scale for data
    //TCut data_scale = Form("%d / (%f * %f)", V.ps_factor, V.charge_mC, V.hms_eff);
    TCut data_scale = Form("%d / (%f)", V.ps_factor, V.hms_eff);

    // Create a smart pointer histogram for data
    auto h = std::make_unique<TH1D>(Form("hData_run%d_%s", run, dndVar.c_str()), "", nbins, xmin, xmax);
    // Keep the error info and project respective variable
    h->Sumw2(true);
    tData->Project(h->GetName(), dndVar.c_str(), dnd_delta_cuts * data_scale);

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

/*
// Build a summed, already-normalized Data histogram over many runs
static std::unique_ptr<TH1D> BuildDataSum(const std::vector<int>& dataRuns,
						const std::string& dndVar,
						int nbins,
						double xmin,
						double xmax,
						const TCut& dnd_delta_cuts) {

    // Create a smart pointer for summed histogram
    std::unique_ptr<TH1D> hSum;
    // Loop over the Data Runs
    for (int run : dataRuns) {
      // Create a histogram to store single run histogram
      auto h = ProjectOneDataRun(run, dndVar, nbins, xmin, xmax, dnd_delta_cuts);

      // If single run histogram can't be made, skip this run
      if (!h) {cout << "skipped this run = " << run  << endl; continue;}

      // hSum is empty initially, so we clone the single run histogram
      if (!hSum) {
        hSum.reset((TH1D*)h->Clone(Form("hDataSum_%s", dndVar.c_str())));
        // Detach ownership from curent directory
        hSum->SetDirectory(nullptr);
      }
      // else we add the single run histograms repetatively
      else {
        hSum->Add(h.get(), 1.0);
      }
    }

    return hSum;
}
*/

// Build an averaged Data histogram for many runs
static std::unique_ptr<TH1D> BuildDataAvg(const std::vector<int>& dataRuns,
						const std::string& dndVar,
						int nbins,
						double xmin,
						double xmax,
						const TCut& dnd_delta_cuts) {

    // Create a smart pointer for averaged histogram
    std::unique_ptr<TH1D> hAvg;
    // Variable for total charge
    double Qtot = 0.0;

    // Loop over the Data Runs
    for (int run : dataRuns) {
      // Create a histogram to store single run histogram
      auto h = ProjectOneDataRun(run, dndVar, nbins, xmin, xmax, dnd_delta_cuts, Qtot);

      // If single run histogram can't be made, skip this run
      if (!h) {cout << "skipped this run = " << run  << endl; continue;}

      // hSum is empty initially, so we clone the single run histogram
      if (!hAvg) {
        hAvg.reset((TH1D*)h->Clone(Form("hDataSum_%s", dndVar.c_str())));
        // Detach ownership from curent directory
        hAvg->SetDirectory(nullptr);
      }
      // else we add the single run histograms repetatively
      else {
        hAvg->Add(h.get(), 1.0);
      }
    }

    // Average the histogram
    if (hAvg && Qtot > 0){
	cout << "Total Data Charge : " << Qtot << endl;
	hAvg->Scale(1.0 / Qtot);
    }
    return hAvg;
}




// Build a normalized Dummy histogram
static std::unique_ptr<TH1D> BuildDummy(int dummyRun,
					const std::string& dndVar,
					int nbins,
					double xmin,
					double xmax,
					const TCut& dnd_delta_cuts) {

    // Get the dummy file
    std::string fpath = DnDRootPath(dummyRun);
    std::unique_ptr<TFile> f(TFile::Open(fpath.c_str(), "READ"));

    // Get the dummy tree
    TTree* tDummy = (TTree*)f->Get("T");

    // Get the values from report file
    ReportValues V = ParseReportFile(DnDReportPath(dummyRun));
    cout << "dummyRun " << dummyRun << ": charge = " << V.charge_mC << ", hms_eff = " << V.hms_eff << ", ps_factor = " << V.ps_factor << endl;

    // Create the scale for dummy
    TCut dummy_scale = Form("%d / (%f * %f)", V.ps_factor, V.charge_mC, V.hms_eff);

    // Create a smart pointer histogram for dummy
    auto h = std::make_unique<TH1D>(Form("hDummy_run%d_%s", dummyRun, dndVar.c_str()), "", nbins, xmin, xmax);

    // Keep error info and project respective variable
    h->Sumw2(true);
    tDummy->Project(h->GetName(), dndVar.c_str(), dnd_delta_cuts * dummy_scale);

    // Detach ownership from current directory
    h->SetDirectory(nullptr);

    return h;
}


//============END BUILDING HISTOGRAMS============\\


// The multi-run plotting function
void PlotVariablesMultiRuns(const std::vector<int>& dataRuns,
				int dummyRun,
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
    auto hDataSum   = BuildDataAvg(dataRuns, dndVar, nbins, xmin, xmax, dnd_delta_cuts);
    auto hDummy     = BuildDummy(dummyRun, dndVar, nbins, xmin, xmax, dnd_delta_cuts);

    // Error if building histograms are incorrect
    if (!hSim || !hDataSum || !hDummy) {
      std::cerr << "[ERROR] Missing hist (sim/data/dummy). Aborting " << simVar << std::endl;
      return;
    }

    // DataSum - Dummy / wall_thickness_ratio
    std::unique_ptr<TH1D> hDataSubDummy(new TH1D(Form("hDataSubDummy_%s", dndVar.c_str()), "", nbins, xmin, xmax)); //creates an empty histogram and wraps it in a smart pointer

    // Detach the ownership from current directory
    hDataSubDummy->SetDirectory(nullptr);

    hDataSubDummy->Sumw2(true); //keep the error information
    hDataSubDummy->Add(hDataSum.get(), 1.0); //add the summed data histogram first
    hDataSubDummy->Add(hDummy.get(), -1.0 / wall_thickness_ratio); //then subtract the dummy histogram

    // Plot comparison and ratio
    PlotComparisonAndRatio(hSim.get(), hDataSubDummy.get(), simVar);


    // Push back the histos to namespace smart pointer vector
    // so that these histos can stay alive when current fuction call ends.
    g_keep_hists.push_back(std::move(hSim));
    g_keep_hists.push_back(std::move(hDataSum));
    g_keep_hists.push_back(std::move(hDummy));
    g_keep_hists.push_back(std::move(hDataSubDummy));

}


// MAIN FUNCTION
void DataVsSimPlot_MultiDataSingleDummy() {

    // Files that don't depend on run numbers, i.e. sim files
    TFile* fSim = TFile::Open("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    TTree* tSim = (TTree*) fSim->Get("h10");

    // Data runs and Dummy run
    std::vector<int> dataRuns = {24017, 24083, 24136, 24183};
    int dummyRun = 24118; // single dummy for now

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
    PlotVariablesMultiRuns(dataRuns, dummyRun, "hsdelta", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRun, "hsytar", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRun, "hsxptar", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRun, "hsyptar", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRun, "xb", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRun, "q2", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRun, "w", tSim, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);

}

