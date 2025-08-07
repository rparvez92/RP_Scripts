// ROOT macro to automate plotting of simulation vs data comparison

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo> //For typeid function
#include "SimToDataMap.h"
#include "PlotComparisonAndRatio.h"

// This function plots the variables
void PlotVariables(std::string simVar, TTree* tSim, TTree* tData, TTree* tDummy, int nbins, double xmin, double xmax, double wall_thickness_ratio, TCut sim_delta_cuts, TCut sim_norm_cuts, TCut dnd_delta_cuts, TCut data_scale, TCut dummy_scale) {

    // Get Data Variables
    std::string dndVar = SimToDataMap(simVar);

    // Declare histograms
    TH1D* hSim  = new TH1D(Form("hSim_%s", simVar.c_str()), "", nbins, xmin, xmax);
    TH1D* hData = new TH1D(Form("hData_%s", dndVar.c_str()), "", nbins, xmin, xmax);
    TH1D* hDummy = new TH1D(Form("hDummy_%s", dndVar.c_str()), "", nbins, xmin, xmax);
    TH1D* hDataSubDummy = new TH1D(Form("hDataSubDummy_%s", dndVar.c_str()), "", nbins, xmin, xmax);

    // Project the branches
    tSim->Project(hSim->GetName(), simVar.c_str(), sim_delta_cuts * sim_norm_cuts);
    tData->Project(hData->GetName(), dndVar.c_str(), dnd_delta_cuts * data_scale);
    tDummy->Project(hDummy->GetName(), dndVar.c_str(), dnd_delta_cuts * dummy_scale);

    // Data - Dummy
    hDataSubDummy->Add(hData, hDummy, 1.0, -1.0 / wall_thickness_ratio);

    // Plot the comparison and ratio
    PlotComparisonAndRatio(hSim, hDataSubDummy, simVar);

}


// MAIN FUNCTION
void DataVsSimPlot_variables() {

    // Get the root files
    TFile* fSim = TFile::Open("./single_arm_worksim/hms_29p05deg_1p531gev_hyd_rsidis.root");
    //TFile* fData = TFile::Open("./Rsidis_ROOTfiles/hms_coin_replay_production_24575_-1.root");
    TFile* fData = TFile::Open("./Rsidis_ROOTfiles/coin_replay_production_23909_-1.root");
    //TFile* fData = TFile::Open(Form("./Rsidis_ROOTfiles/hms_coin_replay_production_%d_-1.root", DataRunNmbr));
    //TFile* fDummy = TFile::Open("./Rsidis_ROOTfiles/hms_coin_replay_production_24577_-1.root");
    TFile* fDummy = TFile::Open("./Rsidis_ROOTfiles/hms_coin_replay_production_23925_-1.root");
    //TFile* fDummy = TFile::Open(Form("./Rsidis_ROOTfiles/hms_coin_replay_production_%d_-1.root", DummyRunNmbr));
    cout <<"Found the files. " << endl;


    // Get the trees
    TTree* tSim = (TTree*) fSim -> Get("h10");
    TTree* tData = (TTree*) fData -> Get("T");
    TTree* tDummy = (TTree*) fDummy -> Get("T");


    // Define Constants
    int nbins = 100;
    double xmin = -8.0, xmax = 8.0;

    double normfac = 10.96;
    double charge_data = 18.585; //26.646; //18.585;
    double hms_eff_data = 0.9957; //0.9955; //0.9957;
    double ps_factor_data = 1.0; //2.0; //1.0;
    double charge_dummy = 28.597; //30.000; //28.597;
    double hms_eff_dummy = 0.9964; //0.9964; //0.9964;
    double ps_factor_dummy = 2.0; //2.0; //2.0;
    double wall_thickness_ratio = 3.55;

    //Cuts
    TCut sim_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8) && (stop_id==0))";
    TCut sim_norm_cuts = Form("weight * %f", normfac);
    TCut dnd_delta_cuts = "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && H.cal.etottracknorm>0.7 && H.cer.npeSum>2.0)"; //dnd = Data And Dummy
    TCut data_scale = Form("%f / (%f * %f)", ps_factor_data, charge_data, hms_eff_data);
    TCut dummy_scale = Form("%f / (%f * %f)", ps_factor_dummy, charge_dummy, hms_eff_dummy);


    // Plot Variables
    PlotVariables("hsdelta", tSim, tData, tDummy, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts, data_scale, dummy_scale);
    PlotVariables("hsytar", tSim, tData, tDummy, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts, data_scale, dummy_scale);
    PlotVariables("hsxptar", tSim, tData, tDummy, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts, data_scale, dummy_scale);
    PlotVariables("hsyptar", tSim, tData, tDummy, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts, data_scale, dummy_scale);
    PlotVariables("xb", tSim, tData, tDummy, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts, data_scale, dummy_scale);
    PlotVariables("q2", tSim, tData, tDummy, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts, data_scale, dummy_scale);
    PlotVariables("w", tSim, tData, tDummy, nbins, xmin, xmax, wall_thickness_ratio, sim_delta_cuts, sim_norm_cuts, dnd_delta_cuts, data_scale, dummy_scale);


}

