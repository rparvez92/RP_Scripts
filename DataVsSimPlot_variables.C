// ROOT macro to automate plotting of simulation vs data comparison

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo> //For typeid function
#include "SimToDataMap.h"
#include "PlotComparisonAndRatio.h"



// MAIN FUNCTION
void DataVsSimPlot_variables() {
    //int DataRunNmbr, DummyRunNmbr;
    //DataRunNmbr = 24544; DummyRunNmbr = 24546;

    //std::cout << "Enter simulation run number: "; std::cin >> run_sim;
    //std::cout << "Enter data run number: "; std::cin >> DataRunNmbr;
    //std::cout << "Enter dummy run number: "; std::cin >> DummyRunNmbr;

    // Get the root files
    //TFile* fSim = TFile::Open(Form("./single_arm_worksim/hms_29p05deg_1p531_hyd_rsidis.root", SimRunNmbr));
    TFile* fSim = TFile::Open("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    TFile* fData = TFile::Open("./Rsidis_ROOTfiles/coin_replay_production_23909_-1.root");
    //TFile* fData = TFile::Open(Form("./Rsidis_ROOTfiles/hms_coin_replay_production_%d_-1.root", DataRunNmbr));
    TFile* fDummy = TFile::Open("./Rsidis_ROOTfiles/hms_coin_replay_production_23925_-1.root");
    //TFile* fDummy = TFile::Open(Form("./Rsidis_ROOTfiles/hms_coin_replay_production_%d_-1.root", DummyRunNmbr));
    cout <<"Found the files. " << endl;


    // Get the trees
    TTree* tSim = (TTree*) fSim -> Get("h10");
    TTree* tData = (TTree*) fData -> Get("T");
    TTree* tDummy = (TTree*) fDummy -> Get("T");


    // Define Constants
    int nbins = 100;
    double xmin = -9.0, xmax = 9.0;

    double normfac = 10.96;
    double charge_data = 18.585;
    double hms_eff_data = 0.9957;
    double charge_dummy = 28.597;
    double hms_eff_dummy = 0.9964;

    //Cuts
    TCut sim_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8) && (stop_id==0))";
    TCut sim_norm_cuts = Form("weight * %f", normfac);
    TCut dnd_delta_cuts = "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && H.cal.etottracknorm>0.7 && H.cer.npeSum>2.0)"; //dnd = Data And Dummy
    TCut data_scale = Form("1.0 / (%f * %f)", charge_data, hms_eff_data);
    TCut dummy_scale = Form("2.0 / (%f * %f)", charge_dummy, hms_eff_dummy);


    // Create variable lists
    std::vector<std::string> SimVarList = {"hsdelta", "hsytar"};
    //std::vector<std::string> DnDVarList = {"H.gtr.dp", "H.gtr.y"}; //DnD = Data And Dummy


    // Iterate over variable list
    for (auto var : SimVarList) {
	//Debugging
	cout << "Current sim var = " << var.c_str() << " of type " << typeid(var).name() << endl;
	cout << "Current dnd var = " << SimToDataMap(var).c_str() << " of type " << typeid(SimToDataMap(var)).name() << endl;

	// Declare the histograms
        TH1D* hSim  = new TH1D(Form("hSim_%s", var.c_str()),  "", nbins, xmin, xmax);
        TH1D* hData  = new TH1D(Form("hData_%s", SimToDataMap(var).c_str()),  "", nbins, xmin, xmax);
        TH1D* hDummy = new TH1D(Form("hDummy_%s", SimToDataMap(var).c_str()), "", nbins, xmin, xmax);
        TH1D* hDataSubDummy  = new TH1D(Form("hDataSubDummy_%s", SimToDataMap(var).c_str()),  "", nbins, xmin, xmax);

	hSim -> SetDirectory(0);
	hData -> SetDirectory(0);
	hDummy -> SetDirectory(0);
	hDataSubDummy -> SetDirectory(0);

	// Project the histograms
	tSim->Project(hSim->GetName(),  var.c_str(), sim_delta_cuts * sim_norm_cuts);
	tData->Project(hData->GetName(), SimToDataMap(var).c_str(), dnd_delta_cuts * data_scale);
	tDummy->Project(hDummy->GetName(), SimToDataMap(var).c_str(), dnd_delta_cuts * dummy_scale);

	// Data - Dummy
	hDataSubDummy -> Add(hData, hDummy, 1.0, -1.0/3.550);

	// Create comparison and ratio plot
        PlotComparisonAndRatio(hSim, hDataSubDummy, var);

	// Delete the histos so that next iteration can use them
        delete hSim;
        delete hData;
        delete hDummy;
        delete hDataSubDummy;

    }

}
