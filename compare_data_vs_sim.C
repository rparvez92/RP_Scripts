void compare_data_vs_sim() {
  // ----- User Configuration -----
  const double normfac = 10.96; //24.8755; // MC normalization (e.g., LD2)
  const double data_charge = 15.210096; //3.0;
  const double data_scaler = 30.5; //18.103;
  const double dummy_charge = 457.099476; //1.0;
  const double dummy_scaler = 28.597; //18.975;
  const double lumi_ratio = data_scaler / dummy_scaler;

  // Variable setup: {varName, MC branch, Data branch, title, nBins, xMin, xMax}
  std::vector<std::tuple<TString, TString, TString, TString, int, double, double>> varList = {
    {"xbj", "xb", "H.kin.primary.x_bj", "x-Bjorken", 50, 0.1, 0.4},
    {"delta", "hsdelta", "H.gtr.dp", "HMS Delta", 50, -8.0, 8.0},
    {"ytar", "hsytar", "H.gtr.y", "HMS ytar", 50, -4.0, 4.0}
  };

  // ----- Open files -----
  //TFile* fmc = TFile::Open("/group/c-rsidis/gaskelld/mc-single-arm/worksim/hms_29p05deg_1p531gev_deut_rsidis.root");
  TFile* fmc = TFile::Open("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
  TTree* tmc = (TTree*)fmc->Get("h10");

  //TFile* fdata = TFile::Open("/home/cdaq/rsidis-2025/hallc_replay_rsidis/ROOTfiles/coin_replay_production_23919_-1.root");
  TFile* fdata = TFile::Open("./Rsidis_ROOTfiles/hms_coin_replay_production_24118_-1.root");
  TTree* tdata = (TTree*)fdata->Get("T");

  //TFile* fdummy = TFile::Open("/home/cdaq/rsidis-2025/hallc_replay_rsidis/ROOTfiles/coin_replay_production_23968_-1.root");
  TFile* fdummy = TFile::Open("./Rsidis_ROOTfiles/hms_coin_replay_production_23925_-1.root");
  TTree* tdummy = (TTree*)fdummy->Get("T");

  // ----- Define cuts -----
  TCut mc_cut = "(hsdelta > -8.0 && hsdelta < 8.0 && stop_id == 0)";
  TCut mc_weight = Form("weight * %f", normfac);
  TCut data_cut = "(H.gtr.dp > -8.0 && H.gtr.dp < 8.0 && H.cal.etottracknorm > 0.6 && H.cer.npeSum > 1.0)";
  TCut data_weight = Form("%f/(%f)", data_charge, data_scaler);
  TCut dummy_weight = Form("%f/(%f)", dummy_charge, dummy_scaler);

  // ----- Create canvas -----
  TCanvas* c = new TCanvas("c", "Data - Dummy vs Sim", 1400, 500);
  c->Divide(3, 1);

  // ----- Loop over variables -----
  for (size_t i = 0; i < varList.size(); ++i) {
    const auto& [varName, mcBranch, dataBranch, title, nBins, xMin, xMax] = varList[i];

    // Histograms
    TString hname_mc = "h_mc_" + varName;
    TString hname_data = "h_data_" + varName;
    TString hname_dummy = "h_dummy_" + varName;
    TString hname_sub = "h_sub_" + varName;

    TH1F* h_mc = new TH1F(hname_mc, title, nBins, xMin, xMax);
    TH1F* h_data = new TH1F(hname_data, title, nBins, xMin, xMax);
    TH1F* h_dummy = new TH1F(hname_dummy, title, nBins, xMin, xMax);
    TH1F* h_sub = new TH1F(hname_sub, title + " Data - Dummy", nBins, xMin, xMax);

    h_mc->Sumw2(); h_data->Sumw2(); h_dummy->Sumw2(); h_sub->Sumw2();

    // Fill MC
    tmc->Project(h_mc->GetName(), mcBranch, mc_cut * mc_weight);

    // Fill data and dummy
    tdata->Project(h_data->GetName(), dataBranch, data_cut * data_weight);
    tdummy->Project(h_dummy->GetName(), dataBranch, data_cut * dummy_weight);

    // Subtract dummy from data
    h_sub->Add(h_data, h_dummy, 1.0, -1.0 / lumi_ratio);

    // Draw
    c->cd(i + 1);
    h_sub->SetLineColor(kRed);
    h_sub->SetLineWidth(2);
    h_mc->SetLineColor(kBlue);
    h_mc->SetLineStyle(2);
    h_mc->SetLineWidth(2);

    h_sub->SetMinimum(0);
    h_sub->Draw("hist");
    h_mc->Draw("hist same");

    // Legend
    auto leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(h_sub, "Data - Dummy", "l");
    leg->AddEntry(h_mc, "MC Simulation", "l");
    leg->Draw();
  }
}
