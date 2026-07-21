// DataVsSimPlot_MultiDataMultiDummy.C
//
// Compare RSIDIS data vs SIMC simulation for multiple data+dummmy runs
// with positron subtraction and dummy subtraction, then plot
// data/sim ratios for physics variables.
//
// This version uses rsidis_bigtable_pass0.csv instead of REPORT_OUTPUT
// and ReportParser.h for run-by-run normalization information.

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>     // std::snprintf
#include <map>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCut.h"
#include "TString.h"
#include "TROOT.h"

// Your helper headers
#include "Mapping.h"                    // SimToDataMap(simVar)
#include "PlotComparisonAndRatio.h"     // PlotComparisonAndRatio(dataHist, simHist, varName)
#include "CoincidenceRandomSubtraction.h"  // FillRandomSubtractedHistogram(...)
// #include "ReportParser.h"           // <-- NO LONGER NEEDED

//----------------------------------------------------------------------------
// Anonymous namespace for internal helpers and globals
//----------------------------------------------------------------------------
namespace {

  // Keep histograms alive after files close
  std::vector<std::unique_ptr<TH1>> g_keep_hists;

  // ---------------------- Bigtable CSV reader ------------------------------

  // One entry from rsidis_bigtable_pass0.csv with only the fields we need.
  struct BigTableEntry {
    double charge_mC_from_BCM4A;  // from column BCM4A_Q
    double hms_efficiency;        // from column h_esing_Eff
    double prescale_factor_ps5;   // from column ps5
  };

  // Map: run number → BigTableEntry
  std::unordered_map<int, BigTableEntry> g_bigtable_map;
  bool g_bigtable_loaded = false;

  // Load the CSV into memory (no-op if already loaded)
  void LoadBigTableCSV(const std::string &csv_path = "rsidis_bigtable_pass0.csv")
  {
    if (g_bigtable_loaded) return;

    std::ifstream infile(csv_path);
    if (!infile.is_open()) {
      std::cerr << "[ERROR] Could not open bigtable CSV: " << csv_path << std::endl;
      return;
    }

    std::string header_line;
    if (!std::getline(infile, header_line)) {
      std::cerr << "[ERROR] Bigtable CSV is empty: " << csv_path << std::endl;
      return;
    }

    // Parse header
    std::vector<std::string> column_names;
    {
      std::stringstream header_stream(header_line);
      std::string col;
      while (std::getline(header_stream, col, ',')) {
        column_names.push_back(col);
      }
    }

    auto find_col_index = [&](const std::string &name) -> int {
      for (size_t i = 0; i < column_names.size(); ++i) {
        if (column_names[i] == name) return static_cast<int>(i);
      }
      return -1;
    };

    // Column indices we care about
    const int idx_run            = find_col_index("run");
    const int idx_BCM4A_Q        = find_col_index("BCM4A_Q");
    const int idx_h_esing_Eff    = find_col_index("h_esing_Eff");
    const int idx_ps5            = find_col_index("ps5");

    if (idx_run < 0 || idx_BCM4A_Q < 0 || idx_h_esing_Eff < 0 || idx_ps5 < 0) {
      std::cerr << "[ERROR] One or more required columns are missing in bigtable CSV.\n"
                   "        Needed: run, BCM4A_Q, h_esing_Eff, ps5\n";
      return;
    }

    // Read each data row
    std::string line;
    while (std::getline(infile, line)) {
      if (line.empty()) continue;

      std::stringstream line_stream(line);
      std::vector<std::string> fields;
      std::string field;

      while (std::getline(line_stream, field, ',')) {
        fields.push_back(field);
      }

      int max_idx = std::max({idx_run, idx_BCM4A_Q, idx_h_esing_Eff, idx_ps5});
      if (static_cast<int>(fields.size()) <= max_idx) continue;

      int run_number = std::stoi(fields[idx_run]);

      BigTableEntry entry;
      entry.charge_mC_from_BCM4A = std::stod(fields[idx_BCM4A_Q]);
      entry.hms_efficiency       = std::stod(fields[idx_h_esing_Eff]);
      entry.prescale_factor_ps5  = std::stod(fields[idx_ps5]);

      g_bigtable_map[run_number] = entry;
    }

    g_bigtable_loaded = true;

    std::cout << "[INFO] Loaded " << g_bigtable_map.size()
              << " rows from bigtable CSV: " << csv_path << std::endl;
  }

  // Retrieve pointer to entry for given run (or nullptr if missing)
  const BigTableEntry* GetBigTableEntry(int run_number)
  {
    if (!g_bigtable_loaded) LoadBigTableCSV();

    auto it = g_bigtable_map.find(run_number);
    if (it == g_bigtable_map.end()) {
      std::cerr << "[ERROR] Run " << run_number
                << " not found in rsidis_bigtable_pass0.csv" << std::endl;
      return nullptr;
    }
    return &it->second;
  }

  // ---------------------- PATH HELPERS -------------------------------------

  // Pass0 coin replay ROOT file for a given run
  std::string DnDRootPath(int run_number)
  {
    char buffer[512];
    std::snprintf(buffer, sizeof(buffer),
                  "./Pass0_ROOTfiles/coin_replay_production_%d_-1.root",
                  run_number);
    return std::string(buffer);
  }

} // end anonymous namespace

//----------------------------------------------------------------------------
// Project a single data or dummy run into a histogram, with:
//   - reconstructed variable expression (dataVariableExpression)
//   - delta and PID cuts (deltaAndPidCuts)
//   - scale factor from bigtable: ps5 / hms_eff
//   - random-coincidence subtraction
//
// The function also updates accumulatedCharge_mC by adding charge_mC_from_BCM4A
// for this run.
//----------------------------------------------------------------------------
std::unique_ptr<TH1D> ProjectOneDnDRun(
    int run_number,
    const std::string &dataVariableExpression,
    const std::string &histogramName,
    int numberOfBins,
    double axisMin,
    double axisMax,
    const TCut &deltaAndPidCuts,
    double &accumulatedCharge_mC)
{
  // Open the Pass0 ROOT file
  std::string root_file_path = DnDRootPath(run_number);
  TFile *dataFile = TFile::Open(root_file_path.c_str(), "READ");
  if (!dataFile || dataFile->IsZombie()) {
    std::cerr << "[ERROR] Could not open data file: " << root_file_path << std::endl;
    if (dataFile) dataFile->Close();
    return nullptr;
  }

  TTree *dataTree = dynamic_cast<TTree *>(dataFile->Get("T"));
  if (!dataTree) {
    std::cerr << "[ERROR] Could not find tree 'T' in file: "
              << root_file_path << std::endl;
    dataFile->Close();
    return nullptr;
  }

  // Get normalization info for this run from bigtable
  const BigTableEntry *entry = GetBigTableEntry(run_number);
  if (!entry) {
    std::cerr << "[ERROR] Missing bigtable entry for run " << run_number << std::endl;
    dataFile->Close();
    return nullptr;
  }

  accumulatedCharge_mC += entry->charge_mC_from_BCM4A;

  std::cout << "[RUN " << run_number << "] "
            << "charge_mC (BCM4A_Q) = " << entry->charge_mC_from_BCM4A
            << "  hms_eff (h_esing_Eff) = " << entry->hms_efficiency
            << "  ps_factor (ps5) = " << entry->prescale_factor_ps5
            << std::endl;

  if (entry->charge_mC_from_BCM4A <= 0.0 ||
      entry->hms_efficiency       <= 0.0 ||
      entry->prescale_factor_ps5  <= 0.0) {
    std::cerr << "[WARNING] Non-positive normalization parameter for run "
              << run_number << std::endl;
  }

  // Scale factor for this run: prescale / efficiency
  // This multiplies event weights in the cut expression.
  TCut scaleFactorCut = Form("%f / (%f)",
                             entry->prescale_factor_ps5,
                             entry->hms_efficiency);

  // Combined cuts: physics cuts * scale factor
  TCut totalCutsWithScale = deltaAndPidCuts * scaleFactorCut;

  // Create the histogram for this run
  auto histForRun = std::make_unique<TH1D>(histogramName.c_str(),
                                           histogramName.c_str(),
                                           numberOfBins, axisMin, axisMax);
  histForRun->Sumw2(true);

  // Configure and perform random-coincidence subtraction
  CoincidenceConfig coincidenceConfig; // default settings from your header

  // FillRandomSubtractedHistogram(TTree*         coinTree,
  //                               TString        weightAndCutExpression,
  //                               const char*    variableExpression,
  //                               TH1*           outputHist,
  //                               CoincidenceConfig cfg);
  FillRandomSubtractedHistogram(dataTree,
                                TString(totalCutsWithScale.GetTitle()),
                                dataVariableExpression.c_str(),
                                histForRun.get(),
                                coincidenceConfig);

  histForRun->SetDirectory(nullptr);
  dataFile->Close();

  return histForRun;
}

//----------------------------------------------------------------------------
// Build simulation histogram (normalized per generated event)
//----------------------------------------------------------------------------
std::unique_ptr<TH1D> BuildSimHistogram(
    const std::string &simVariableExpression,
    TTree *simTree,
    const std::string &histogramName,
    int numberOfBins,
    double axisMin,
    double axisMax,
    const TCut &simDeltaCuts,
    const TCut &simNormCuts)
{
  if (!simTree) {
    std::cerr << "[ERROR] simTree is null in BuildSimHistogram" << std::endl;
    return nullptr;
  }

  auto simHist = std::make_unique<TH1D>(histogramName.c_str(),
                                        histogramName.c_str(),
                                        numberOfBins, axisMin, axisMax);
  simHist->Sumw2(true);

  TCut combinedSimCuts = simDeltaCuts * simNormCuts;

  TString drawCommand = TString::Format("%s>>%s(%d,%f,%f)",
                                        simVariableExpression.c_str(),
                                        histogramName.c_str(),
                                        numberOfBins, axisMin, axisMax);

  Long64_t numberOfGeneratedEvents = simTree->GetEntries();

  simTree->Draw(drawCommand, combinedSimCuts, "goff");

  if (numberOfGeneratedEvents > 0) {
    simHist->Scale(1.0 / static_cast<double>(numberOfGeneratedEvents));
  }

  simHist->SetDirectory(nullptr);
  return simHist;
}

//----------------------------------------------------------------------------
// Build charge-normalized average histogram over a list of runs.
//
// Each run is projected with ProjectOneDnDRun (which accumulates charge).
// After summing all runs, the histogram is scaled by 1 / totalCharge_mC.
//----------------------------------------------------------------------------
std::unique_ptr<TH1D> BuildChargeNormalizedAverageHistogram(
    const std::vector<int> &runList,
    const std::string &dataVariableExpression,
    const std::string &baseHistogramName,
    int numberOfBins,
    double axisMin,
    double axisMax,
    const TCut &deltaAndPidCuts,
    double &totalCharge_mC)
{
  totalCharge_mC = 0.0;

  std::unique_ptr<TH1D> averageHist;

  for (int run_number : runList) {
    std::string histNameForRun =
        baseHistogramName + Form("_run%d", run_number);

    auto histForRun = ProjectOneDnDRun(run_number,
                                       dataVariableExpression,
                                       histNameForRun,
                                       numberOfBins,
                                       axisMin,
                                       axisMax,
                                       deltaAndPidCuts,
                                       totalCharge_mC);
    if (!histForRun) continue;

    if (!averageHist) {
      averageHist.reset(
          dynamic_cast<TH1D *>(histForRun->Clone(
              (baseHistogramName + "_avg").c_str())));
      averageHist->SetDirectory(nullptr);
    } else {
      averageHist->Add(histForRun.get());
    }

    g_keep_hists.push_back(std::move(histForRun));
  }

  if (averageHist && totalCharge_mC > 0.0) {
    averageHist->Scale(1.0 / totalCharge_mC);
  }

  return averageHist;
}

//----------------------------------------------------------------------------
// Do everything for a single variable:
//  - map sim variable name → data expression via SimToDataMap
//  - build sim histogram
//  - build data/dummy for electrons and positrons
//  - perform positron subtraction and dummy subtraction
//  - call PlotComparisonAndRatio
//----------------------------------------------------------------------------
void PlotVariablesMultiRuns(
    const std::vector<int> &electronDataRuns,
    const std::vector<int> &electronDummyRuns,
    const std::vector<int> &positronDataRuns,
    const std::vector<int> &positronDummyRuns,
    const std::string &simVariableName,
    TTree *simTree,
    int numberOfBins,
    double axisMin,
    double axisMax,
    double wallThicknessRatio,
    const TCut &simDeltaCuts,
    const TCut &simNormCuts,
    const TCut &dataDeltaAndPidCuts)
{
  // Map SIMC branch name to corresponding data expression
  std::string dataVariableExpression = SimToDataMap(simVariableName);

  // --- Build SIM histogram ---
  auto simHist = BuildSimHistogram(simVariableName,
                                   simTree,
                                   "hSim",
                                   numberOfBins,
                                   axisMin,
                                   axisMax,
                                   simDeltaCuts,
                                   simNormCuts);

  // --- Build charge-normalized average histograms for electrons ---
  double totalCharge_electronData_mC  = 0.0;
  double totalCharge_electronDummy_mC = 0.0;

  auto hDataElectronAvg = BuildChargeNormalizedAverageHistogram(
      electronDataRuns, dataVariableExpression,
      "hDataElectron", numberOfBins, axisMin, axisMax,
      dataDeltaAndPidCuts, totalCharge_electronData_mC);

  auto hDummyElectronAvg = BuildChargeNormalizedAverageHistogram(
      electronDummyRuns, dataVariableExpression,
      "hDummyElectron", numberOfBins, axisMin, axisMax,
      dataDeltaAndPidCuts, totalCharge_electronDummy_mC);

  // --- Build charge-normalized average histograms for positrons ---
  double totalCharge_positronData_mC  = 0.0;
  double totalCharge_positronDummy_mC = 0.0;

  auto hDataPositronAvg = BuildChargeNormalizedAverageHistogram(
      positronDataRuns, dataVariableExpression,
      "hDataPositron", numberOfBins, axisMin, axisMax,
      dataDeltaAndPidCuts, totalCharge_positronData_mC);

  auto hDummyPositronAvg = BuildChargeNormalizedAverageHistogram(
      positronDummyRuns, dataVariableExpression,
      "hDummyPositron", numberOfBins, axisMin, axisMax,
      dataDeltaAndPidCuts, totalCharge_positronDummy_mC);

  if (!simHist || !hDataElectronAvg || !hDummyElectronAvg ||
      !hDataPositronAvg || !hDummyPositronAvg) {
    std::cerr << "[ERROR] Missing one or more histograms in PlotVariablesMultiRuns for "
              << simVariableName << std::endl;
    return;
  }

  // --- Data - Positron ---
  auto hDataMinusPositron = std::make_unique<TH1D>(
      *hDataElectronAvg);
  hDataMinusPositron->SetName(
      ("hDataMinusPositron_" + simVariableName).c_str());
  hDataMinusPositron->Add(hDataPositronAvg.get(), -1.0);

  // --- Dummy - Positron ---
  auto hDummyMinusPositron = std::make_unique<TH1D>(
      *hDummyElectronAvg);
  hDummyMinusPositron->SetName(
      ("hDummyMinusPositron_" + simVariableName).c_str());
  hDummyMinusPositron->Add(hDummyPositronAvg.get(), -1.0);

  // --- Final: (Data - Positron) - (1/wall)* (Dummy - Positron) ---
  auto hDataMinusDummy = std::make_unique<TH1D>(
      *hDataMinusPositron);
  hDataMinusDummy->SetName(
      ("hDataMinusDummy_" + simVariableName).c_str());
  double dummyScaleFactor = 1.0 / wallThicknessRatio;
  hDataMinusDummy->Add(hDummyMinusPositron.get(), -dummyScaleFactor);

  // Compare final data histogram to simulation
  PlotComparisonAndRatio(hDataMinusDummy.get(),
                         simHist.get(),
                         simVariableName.c_str());

  // Keep histograms alive
  g_keep_hists.push_back(std::move(simHist));
  g_keep_hists.push_back(std::move(hDataElectronAvg));
  g_keep_hists.push_back(std::move(hDummyElectronAvg));
  g_keep_hists.push_back(std::move(hDataPositronAvg));
  g_keep_hists.push_back(std::move(hDummyPositronAvg));
  g_keep_hists.push_back(std::move(hDataMinusPositron));
  g_keep_hists.push_back(std::move(hDummyMinusPositron));
  g_keep_hists.push_back(std::move(hDataMinusDummy));
}

//----------------------------------------------------------------------------
// Main macro
//----------------------------------------------------------------------------
void DataVsSimPlot_GPT()
{
  gROOT->SetBatch(kTRUE);

  // ---------------------- Open SIM file ------------------------
  TFile *fSim = TFile::Open(
      "./simc_worksim/coin_7p87deg_3p632gev_hyd_rsidis.root",
      "READ");
  if (!fSim || fSim->IsZombie()) {
    std::cerr << "[ERROR] Could not open SIMC file." << std::endl;
    if (fSim) fSim->Close();
    return;
  }

  TTree *tSim = dynamic_cast<TTree *>(fSim->Get("h10"));
  if (!tSim) {
    std::cerr << "[ERROR] Could not find tree 'h10' in SIMC file." << std::endl;
    fSim->Close();
    return;
  }

  // ---------------------- Run lists ----------------------------
  // (These are the sets you were using; adjust here as needed.)
  std::vector<int> electronDataRuns   = {24329, 24330};
  std::vector<int> electronDummyRuns  = {23839, 23840};  // <-- adjust to your actual list
  std::vector<int> positronDataRuns   = {24603, 24603};
  std::vector<int> positronDummyRuns  = {24601, 24601};

  // ---------------------- Normalization & cuts -----------------
  // SIM normalization factor (from SIMC .hist "Miscellaneous" section)
  double simNormalizationFactor = 0.842205e11;

  TCut sim_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8.0))";
  TCut sim_norm_cuts  = Form("Weight * %e", simNormalizationFactor);

  // Data delta and PID cuts (HMS)
  TCut data_delta_pid_cuts =
      "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && "
      " (H.cal.etottracknorm>0.7) && (H.cer.npeSum>2.0))";

  // Dummy wall thickness ratio (dummy thickness / LH2 thickness)
  double wall_thickness_ratio = 3.82;

  // ---------------------- Binning ------------------------------
  struct VarBinning {
    int    nbins;
    double xmin;
    double xmax;
  };

  std::map<std::string, VarBinning> binsFor;

  // Fill in the binning (you can extend/add variables as in your old macro)
  binsFor["Q2"] = {20, 1.0, 4.0};  // <-- set to your preferred Q2 binning

  // ---------------------- Call for one variable ----------------
  VarBinning q2bin = binsFor["Q2"];

  PlotVariablesMultiRuns(
      electronDataRuns,
      electronDummyRuns,
      positronDataRuns,
      positronDummyRuns,
      "Q2",                // sim variable name
      tSim,
      q2bin.nbins,
      q2bin.xmin,
      q2bin.xmax,
      wall_thickness_ratio,
      sim_delta_cuts,
      sim_norm_cuts,
      data_delta_pid_cuts);

  // You can re-enable additional variables by adding more calls here,
  // once you populate binsFor["xbj"], binsFor["z"], etc.
}
