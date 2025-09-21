// ROOT macro to automate plotting of data vs simulation comparison

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo> //For typeid function
#include "Mapping.h"
#include "ReportParser.h"
#include "PlotComparisonAndRatio.h"
#include "CoincidenceRandomSubtraction.h" // For coincidence time and random subtraction


// Creating an anonymous namespace to store unique_ptrs in a global vector, so that the objects
// pointed by these pointers do not get destroyed at the end of function call.
namespace {
  std::vector<std::unique_ptr<TH1>> g_keep_hists;
}

// Function for returning file path
static std::string DnDRootPath(int run) {
  return Form("./Rsidis_ROOTfiles/coin_replay_production_%d_-1.root", run);
}
static std::string DnDReportPath(int run) {
  return Form("./REPORT_OUTPUT/COIN/PRODUCTION/replay_coin_production_%d_-1.report", run);
}


//============START BUILDING HISTOGRAMS============\\


// Create and project a normalized histogram for a SINGLE data or dummy run
static std::unique_ptr<TH1D> ProjectOneDnDRun(int run,
						const std::string& dndVar,
						int nbins,
						double xmin,
						double xmax,
						TCut& dnd_delta_cuts,
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
    // Keep the error info
    h->Sumw2(true);

    // When z variable is passed, dndVar is P.gtr.p/H.kin.primary.nu
    // We need to make sure that denominator is not 0
    if (dndVar.find("P.gtr.p/H.kin.primary.nu") != std::string::npos) {
        dnd_delta_cuts = dnd_delta_cuts && "(H.kin.primary.nu>0)";
    }

    // Apply Coincidence Time Configuration: (defaults: [20,80] ns, RF=4 ns, ±1 ns coin window)
    CoincidenceConfig ctCfg;
    // Fill random-subtracted histogram for this run
    FillRandomSubtractedHistogram(tDnD, TString(dnd_delta_cuts.GetTitle()), dndVar.c_str(), h.get(), ctCfg);

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

    // Scale by total generated events
    const Long64_t nGenSim = tSim->GetEntries();
    cout << "Total generated events for Simulation: " << nGenSim << endl;
    h->Scale(1.0 / double(nGenSim));

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
						TCut& dnd_delta_cuts) {

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
						TCut& dnd_delta_cuts) {

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
    TFile* fSim = TFile::Open("./simc_worksim/coin_7p87deg_3p632gev_hyd_rsidis.root");
    TTree* tSim = (TTree*) fSim->Get("h10");

    // Data runs and Dummy runs
    std::vector<int> dataRuns = {24329, 24330, 24331, 24332};
    std::vector<int> dummyRuns = {24335, 24336, 24337, 24338};

    // Cuts and normalizations
    double normfac = 0.842205E+11; //Found from Miscellaneous section in respective .hist file
    TCut sim_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8))";
    TCut sim_norm_cuts  = Form("Weight * %f", normfac);
    TCut dnd_delta_cuts = "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && H.cal.etottracknorm>0.7 && H.cer.npeSum>2.0)";

    // Constants
    //int nbins = 100;
    //double xmin = -8.0, xmax = 8.0;
    //double xmin = 0.0, xmax = 1.0;
    double wall_thickness_ratio = 3.82; //Dummy_thicknes / Data_thickness

    // Declare variable specific bin parameters
    struct VarBinning { int nbins; double xmin; double xmax; }; //Declaring a Structure
    // Creating a map of variable to bin paramtervalues
    std::unordered_map<std::string, VarBinning> binsFor = {
      {"hsdelta", {100, -12.0, 12.0}},	{"hsytar", {100, -12.0, 12.0}},
      {"hsxptar", {100, -1.0, 1.0}},	{"hsyptar", {100, -1.0, 1.0}},
      {"ssdelta", {100, -20.0, 20.0}},	{"ssytar", {100, -20.0, 20.0}},
      {"ssxptar", {100, -1.0, 1.0}},	{"ssyptar", {100, -1.0, 1.0}},
      {"z", {100, 0.0, 1.0}},		{"xbj", {100, 0.0, 1.0}},
      {"Q2", {100, 0.0, 8.0}},		{"W", {100, 0.0, 8.0}},
      {"nu", {100, 0.0, 8.0}},		{"epsilon", {100, 0.0, 1.0}},
      {"thetapq", {100, -8.0, 8.0}},	{"phipq", {100, -8.0, 8.0}},
    };

    // Plot each variable
    // HMS Variables
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsdelta", tSim, binsFor["hsdelta"].nbins, binsFor["hsdelta"].xmin, binsFor["hsdelta"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsytar", tSim, binsFor["hsytar"].nbins, binsFor["hsytar"].xmin, binsFor["hsytar"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsxptar", tSim, binsFor["hsxptar"].nbins, binsFor["hsxptar"].xmin, binsFor["hsxptar"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "hsyptar", tSim, binsFor["hsyptar"].nbins, binsFor["hsyptar"].xmin, binsFor["hsyptar"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    // SHMS Variables
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "ssdelta", tSim, binsFor["ssdelta"].nbins, binsFor["ssdelta"].xmin, binsFor["ssdelta"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "ssytar", tSim, binsFor["ssytar"].nbins, binsFor["ssytar"].xmin, binsFor["ssytar"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "ssxptar", tSim, binsFor["ssxptar"].nbins, binsFor["ssxptar"].xmin, binsFor["ssxptar"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "ssyptar", tSim, binsFor["ssyptar"].nbins, binsFor["ssyptar"].xmin, binsFor["ssyptar"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    // Kinematic Variables
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "z", tSim, binsFor["z"].nbins, binsFor["z"].xmin, binsFor["z"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "xbj", tSim, binsFor["xbj"].nbins, binsFor["xbj"].xmin, binsFor["xbj"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "Q2", tSim, binsFor["Q2"].nbins, binsFor["Q2"].xmin, binsFor["Q2"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "W", tSim, binsFor["W"].nbins, binsFor["W"].xmin, binsFor["W"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    PlotVariablesMultiRuns(dataRuns, dummyRuns, "nu", tSim, binsFor["nu"].nbins, binsFor["nu"].xmin, binsFor["nu"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "epsilon", tSim, binsFor["epsilon"].nbins, binsFor["epsilon"].xmin, binsFor["epsilon"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "thetapq", tSim, binsFor["thetapq"].nbins, binsFor["thetapq"].xmin, binsFor["thetapq"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);
    //PlotVariablesMultiRuns(dataRuns, dummyRuns, "phipq", tSim, binsFor["phipq"].nbins, binsFor["phipq"].xmin, binsFor["phipq"].xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts);

}

