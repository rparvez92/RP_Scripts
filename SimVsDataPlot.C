#include <iostream>

void SimVsDataPlot() {

	// Grab and Open the root files
	TFile *fSim = TFile::Open("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
	TFile *fData = TFile::Open("./Rsidis_ROOTfiles/hms_coin_replay_production_24118_-1.root");

	// Show error if not found
	if (!fSim || !fSim->IsOpen()) { 
		Error("SimVsDataPlot","Simulation file not found!"); return; 
	}
	if (!fData || !fData->IsOpen()) { 
		Error("SimVsDataPlot","Data file not found!"); return; 
	}
	cout << "Found the root files." << endl;

	// Grab the Trees
	TTree *tree_s = (TTree*) fSim -> Get("h10");
	TTree *tree_d = (TTree*) fData -> Get("T");

	// Print nEntries
	cout << "Sim nEntries: " << tree_s -> GetEntries() << endl;
	cout << "Data nEntries: " << tree_d -> GetEntries() << endl;

	// Declare Histograms
	TH1F *hSimXbj = new TH1F("hSimXbj", "SimXbj", 100, 0, 1);
	TH1F *hSimDelta = new TH1F("hSimDelta", "SimDelta", 100, -35.0, 35.0);
	TH1F *hSimYtar = new TH1F("hSimYtar", "SimYtar", 100, -4.0, 4.0);

	TH1F *hDataDelta = new TH1F("hDataDelta", "DataDelta", 100, -35.0, 35.0);


	// Keep the sum of squered weights for error calculation
	hSimXbj -> Sumw2();
	hSimDelta -> Sumw2();
	hSimYtar -> Sumw2();

	hDataDelta -> Sumw2();

	// Declare constants
	double normfac = 10.96; // Sim normalization
	double hms_eff = 0.9986; // Data normalization
	double charge = 15.210096; // mC; Data normalization

	// Declare Cuts
	TCut sim_delta_cut = "(hsdelta > -8.0) && (hsdelta < 8.0) && (stop_id == 0)";
	TCut sim_norm_cut = Form("weight * %f", normfac);

	TCut data_delta_cut = "(H.gtr.dp > -8.0) && (H.gtr.dp < 8.0) && (H.cal.etottracknorm > 0.7) && (H.cer.npesum > 2.0)";

	// Using Project().
	// Note: Project() is not good for nested conditional filling histogram.
	tree_s -> Project("hSimDelta", "hsdelta", sim_delta_cut * sim_norm_cut);
	tree_d -> Project("hDataDelta", "H.gtr.dp", data_delta_cut);


	// Normalize Data
	hDataDelta -> Scale(1.0/(charge * hms_eff));

	// ==== STYLE ====
	hSimDelta->SetLineColor(kBlue+2);
	hSimDelta->SetLineWidth(2);

	hDataDelta->SetLineColor(kBlack);
	hDataDelta->SetLineWidth(2);

	//hSimXfp->SetLineColor(kBlue+2);
	//hSimXfp->SetLineWidth(2);

	//hDataXfp->SetLineColor(kBlack);
	//hDataXfp->SetLineWidth(2);




	// Draw
	TCanvas *c1 = new TCanvas("c1", "Sim Vs Data: Delta", 900, 700);
	hSimDelta -> Draw("hist");
	hDataDelta -> Draw("hist same");


	auto legend1 = new TLegend(0.6,0.7,0.85,0.85);
	legend1 -> AddEntry(hSimDelta,"Sim","l");
	legend1 -> AddEntry(hDataDelta,"Data","l");
	legend1 -> Draw();
	// For log scale
	// c1->SetLogy();

	// For saving to PNG
	// c1->SaveAs("Data_vs_Sim.png");




	// Update canvas
	c1 -> Update();
}
