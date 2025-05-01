/**
 *@file FlatVertexProcessor.cxx
 *@brief Generating flat tuples for preselection and tight selection analysis
 *@author Sarah, SLAC
 */

#include "FlatVertexProcessor.h"
#include <iostream>
#include <fstream>
#include <map>

FlatVertexProcessor::FlatVertexProcessor(const std::string& name, Process& process) : Processor(name, process) {

}

//TODO Check this destructor

FlatVertexProcessor::~FlatVertexProcessor(){}

void FlatVertexProcessor::configure(const ParameterSet& parameters) {
    std::cout << "Configuring FlatVertexProcessor" << std::endl;
    try
    {
        debug_   = parameters.getInteger("debug", debug_);
        anaName_ = parameters.getString("anaName", anaName_);
        tsColl_  = parameters.getString("tsColl", tsColl_);
        vtxColl_ = parameters.getString("vtxColl", vtxColl_);
        hitColl_ = parameters.getString("hitColl", hitColl_);
        mcColl_  = parameters.getString("mcColl", mcColl_);
        ecalColl_ = parameters.getString("ecalColl", ecalColl_);
        isRadPDG_ = parameters.getInteger("isRadPDG", isRadPDG_);

        selectionCfg_ = parameters.getString("vtxSelectionjson", selectionCfg_);
        // histoCfg_ = parameters.getString("histoCfg", histoCfg_);
        // mcHistoCfg_ = parameters.getString("mcHistoCfg", mcHistoCfg_);
        timeOffset_ = parameters.getDouble("CalTimeOffset", timeOffset_);
        beamE_  = parameters.getDouble("beamE", beamE_);
        isData_  = parameters.getInteger("isData", isData_);
        analysis_ = parameters.getString("analysis");

        pSmearingFile_ = parameters.getString("pSmearingFile", pSmearingFile_);
        pBiasingFile_  = parameters.getString("pBiasingFile", pBiasingFile_);

        //region definitions
        regionSelections_ = parameters.getVString("regionDefinitions", regionSelections_);

        //v0 projection fits
        v0ProjectionFitsCfg_ = parameters.getString("v0ProjectionFitsCfg", v0ProjectionFitsCfg_);

        //beamspot positions
        beamPosCfg_ = parameters.getString("beamPosCfg", beamPosCfg_);
        //track time bias corrections
        eleTrackTimeBias_ = parameters.getDouble("eleTrackTimeBias", eleTrackTimeBias_);
        posTrackTimeBias_ = parameters.getDouble("posTrackTimeBias", posTrackTimeBias_);

        //bField scale factor (ONLY to correct for mistakes in ntuplizing). If < 0  is not used
        bFieldScaleFactor_ = parameters.getDouble("bFieldScaleFactor", -1);
    }
    catch (std::runtime_error& error)
    {
        std::cout << error.what() << std::endl;
    }
}

void FlatVertexProcessor::initialize(TTree* tree) {
    tree_ = tree;
    _ah = std::make_shared<AnaHelpers>();

    vtxSelector = std::make_shared<BaseSelector>(anaName_ + "_" + "vtxSelection", selectionCfg_);
    vtxSelector->setDebug(debug_);
    vtxSelector->LoadSelection();

    // Load Run Dependent V0 target projection fits from json
    if (!v0ProjectionFitsCfg_.empty()){
        std::ifstream v0proj_file(v0ProjectionFitsCfg_);
        v0proj_file >> v0proj_fits_;
        v0proj_file.close();
    }

    // Run Dependent Corrections
    // Beam Position
    if (!beamPosCfg_.empty()){
        std::ifstream bpc_file(beamPosCfg_);
        bpc_file >> bpc_configs_;
        bpc_file.close();
    }

    //For each region initialize plots

    for (unsigned int i_reg = 0; i_reg < regionSelections_.size(); i_reg++) {
        std::string regname = AnaHelpers::getFileName(regionSelections_[i_reg], false);
        std::cout << "Setting up region:: " << regname << std::endl;
        _reg_vtx_selectors[regname] = std::make_shared<BaseSelector>(anaName_ + "_" + regname, regionSelections_[i_reg]);
        _reg_vtx_selectors[regname]->setDebug(debug_);
        _reg_vtx_selectors[regname]->LoadSelection();

        // Build a flat tuple for vertex and track params
        _reg_tuples[regname] = std::make_shared<FlatTupleMaker>(anaName_ + "_" + regname + "_tree");

        // vtx vars
        _reg_tuples[regname]->addVariable("unc_vtx_mass");
        _reg_tuples[regname]->addVariable("unc_vtx_z");
        _reg_tuples[regname]->addVariable("unc_vtx_chi2");
        _reg_tuples[regname]->addVariable("unc_vtx_psum");
        _reg_tuples[regname]->addVariable("unc_vtx_px");
        _reg_tuples[regname]->addVariable("unc_vtx_py");
        _reg_tuples[regname]->addVariable("unc_vtx_pz");
        _reg_tuples[regname]->addVariable("unc_vtx_x");
        _reg_tuples[regname]->addVariable("unc_vtx_y");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_pos_clus_dt");
        _reg_tuples[regname]->addVariable("run_number");
	_reg_tuples[regname]->addVariable("event_number");
	_reg_tuples[regname]->addVariable("singles0trigger");
	_reg_tuples[regname]->addVariable("singles1trigger");
	_reg_tuples[regname]->addVariable("singles2trigger");
	_reg_tuples[regname]->addVariable("singles3trigger");
        _reg_tuples[regname]->addVariable("unc_vtx_cxx");
        _reg_tuples[regname]->addVariable("unc_vtx_cyy");
        _reg_tuples[regname]->addVariable("unc_vtx_czz");
        _reg_tuples[regname]->addVariable("unc_vtx_cyx");
        _reg_tuples[regname]->addVariable("unc_vtx_czy");
        _reg_tuples[regname]->addVariable("unc_vtx_czx");
        _reg_tuples[regname]->addVariable("unc_vtx_proj_x");
        _reg_tuples[regname]->addVariable("unc_vtx_proj_y");
        _reg_tuples[regname]->addVariable("unc_vtx_proj_x_sig");
        _reg_tuples[regname]->addVariable("unc_vtx_proj_y_sig");
        _reg_tuples[regname]->addVariable("unc_vtx_proj_sig");
        _reg_tuples[regname]->addVariable("unc_vtx_deltaZ");
        _reg_tuples[regname]->addVariable("nGoodVtx");


        // track vars
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_p");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_t");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_d0");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_phi0");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_omega");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_tanLambda");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_z0");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_chi2ndf");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_clust_dt");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_z0Err");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_d0Err");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_tanLambdaErr");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_PhiErr");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_OmegaErr");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_L1_isolation");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_nhits");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_lastlayer");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_si0");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_si1");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_ecal_x");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_ecal_y");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_z");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_px");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_py");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_track_pz");

        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_clust_dt");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_p");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_t");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_d0");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_phi0");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_omega");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_tanLambda");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_z0");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_chi2ndf");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_z0Err");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_d0Err");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_tanLambdaErr");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_PhiErr");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_OmegaErr");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_L1_isolation");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_nhits");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_lastlayer");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_si0");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_si1");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_ecal_x");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_ecal_y");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_z");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_px");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_py");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_track_pz");

        // clust vars
        _reg_tuples[regname]->addVariable("unc_vtx_ele_clust_E");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_clust_x");
        _reg_tuples[regname]->addVariable("unc_vtx_ele_clust_corr_t");

        _reg_tuples[regname]->addVariable("unc_vtx_pos_clust_E");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_clust_x");
        _reg_tuples[regname]->addVariable("unc_vtx_pos_clust_corr_t");

        if(!isData_)
        {
            // _reg_tuples[regname]->addVariable("true_vtx_z");
            // _reg_tuples[regname]->addVariable("true_vtx_mass");
            _reg_tuples[regname]->addVariable("ap_true_vtx_z");
            _reg_tuples[regname]->addVariable("ap_true_vtx_mass");
            _reg_tuples[regname]->addVariable("ap_true_vtx_energy");
            // _reg_tuples[regname]->addVariable("vd_true_vtx_z");
            // _reg_tuples[regname]->addVariable("vd_true_vtx_mass");
            // _reg_tuples[regname]->addVariable("vd_true_vtx_energy");
            _reg_tuples[regname]->addVariable("hitCode");
            _reg_tuples[regname]->addVariable("L1hitCode");
            _reg_tuples[regname]->addVariable("L2hitCode");
            _reg_tuples[regname]->addVariable("momPDG_ele");
            _reg_tuples[regname]->addVariable("momPDG_pos");
            _reg_tuples[regname]->addVariable("originPDG_ele");
            _reg_tuples[regname]->addVariable("originPDG_pos");
        }

        _regions.push_back(regname);
    }

    // Get list of branches in tree to help protect accessing them
    int nBr = tree_->GetListOfBranches()->GetEntries();
    if(debug_) {
        std::cout << "Tree has " << nBr << " branches" << std::endl;
    }
    for (int iBr = 0; iBr < nBr; iBr++) {
        TBranch *br = dynamic_cast<TBranch*>(tree_->GetListOfBranches()->At(iBr));
        brMap_.insert(std::map<const char *, int, char_cmp>::value_type(br->GetName(), 1));
        if (debug_) {
            std::cout << br->GetName() << ": " << brMap_[br->GetName()] << std::endl;
        }
    }

    // init Reading Tree
    tree_->SetBranchAddress("EventHeader", &evth_ , &bevth_);
    if (brMap_.find(tsColl_.c_str()) != brMap_.end()) {
        tree_->SetBranchAddress(tsColl_.c_str(), &ts_, &bts_);
    }
    tree_->SetBranchAddress(vtxColl_.c_str(), &vtxs_, &bvtxs_);

    if (brMap_.find(hitColl_.c_str()) != brMap_.end()) { 
        tree_->SetBranchAddress(hitColl_.c_str(), &hits_, &bhits_);
    }
    if (!isData_ && !mcColl_.empty()) {
        tree_->SetBranchAddress(mcColl_.c_str(), &mcParts_, &bmcParts_);
    }

    if (!trkColl_.empty()) {
        tree_->SetBranchAddress(trkColl_.c_str(),&trks_, &btrks_);
    }
    if (not pSmearingFile_.empty()) {
      // just using the same seed=42 for now
      smearingTool_ = std::make_shared<TrackSmearingTool>(pSmearingFile_, true);
    }

    if (not pBiasingFile_.empty()) {
      biasingTool_ = std::make_shared<TrackBiasingTool>(pBiasingFile_);
    }
}

bool FlatVertexProcessor::process(IEvent* ievent) {
    if (debug_) {
        std:: cout << "----------------- Event " << evth_->getEventNumber() << " -----------------" << std::endl;
    }
    HpsEvent* hps_evt = (HpsEvent*) ievent;
    double weight = 1.;
    int run_number = evth_->getRunNumber();
    int closest_run;
    if (debug_) std::cout << "Check pbc_configs" << std::endl;
    if (!bpc_configs_.empty()) {
        for (auto run : bpc_configs_.items()){
            int check_run = std::stoi(run.key());
            if (check_run > run_number) {
                break;
            } else {
                closest_run = check_run;
            }
        }
        beamPosCorrections_ = {bpc_configs_[std::to_string(closest_run)]["beamspot_x"], 
                               bpc_configs_[std::to_string(closest_run)]["beamspot_y"],
                               bpc_configs_[std::to_string(closest_run)]["beamspot_z"]};
    }


    //Get "true" values
    //AP
    double apMass = -0.9;
    double apZ = -0.9;
    double apEnergy = -0.9;
    //Simp
    double vdMass = -0.9;
    double vdZ = -0.9;
    double vdEnergy = -0.9;
    
    if (mcParts_) {
        for (int i = 0; i < mcParts_->size(); i++)
        {
            if (mcParts_->at(i)->getPDG() == 622)
            {
                apMass = mcParts_->at(i)->getMass();
                apZ = mcParts_->at(i)->getVertexPosition().at(2);
                apEnergy = mcParts_->at(i)->getEnergy();
            }
            if (mcParts_->at(i)->getPDG() == 625)
            {
                vdMass = mcParts_->at(i)->getMass();
                vdZ = mcParts_->at(i)->getVertexPosition().at(2);
                vdEnergy = mcParts_->at(i)->getEnergy();
            }
        }
    }

    //Store processed number of events
    std::vector<Vertex*> selected_vtxs;
    bool passVtxPresel = false;

    if(debug_){
        std::cout << "Number of vertices found in event: " << vtxs_->size() << std::endl;
    }

    // Loop over vertices in event and make selections
    for (int i_vtx = 0; i_vtx < vtxs_->size(); i_vtx++) {
        vtxSelector->getCutFlowHisto()->Fill(0., weight);
	
        Vertex* vtx = vtxs_->at(i_vtx);
        Particle* ele = nullptr;
        Particle* pos = nullptr;

        if (isData_) {
            if (!vtxSelector->passCutEq("Pair1_eq", (int)evth_->isPair1Trigger(),weight))
                break;
	    if (!vtxSelector->passCutEq("Singles0_eq", (int)ts_->isSingles0Trigger(), weight))
		break;
	    if (!vtxSelector->passCutEq("Singles1_eq", (int)ts_->isSingles1Trigger(), weight))
                break;
	    if (!vtxSelector->passCutEq("Singles2_eq", (int)ts_->isSingles2Trigger(), weight))
                break;
	    if (!vtxSelector->passCutEq("Singles3_eq", (int)ts_->isSingles3Trigger(), weight))
                break;
	    if (!vtxSelector->passCutEq("Singles2or3_eq", (int)(ts_->isSingles2Trigger() || ts_->isSingles3Trigger()), weight))
                break;
        }

        bool foundParts = _ah->GetParticlesFromVtx(vtx, ele, pos);
        if (!foundParts) {
            if (debug_) {
                std::cout << "FlatVertexProcessor::WARNING::Found vtx without ele/pos. Skip." << std::endl;
            }
            continue;
        }

        if (debug_) std::cout << "got parts" << std::endl;
        // Track ele_trk = ele->getTrack();
        // Track pos_trk = pos->getTrack();
        Track ele_trk;
        Track pos_trk;
        Track* ele_trk_ptr;
        Track* pos_trk_ptr;
	
        if (!trkColl_.empty()) {
            bool foundTracks = _ah->MatchToGBLTracks((ele->getTrack()).getID(),(pos->getTrack()).getID(),
                    ele_trk_ptr, pos_trk_ptr, *trks_);

            if (!foundTracks) {
                if(debug_) std::cout<<"VertexAnaProcessor::ERROR couldn't find ele/pos in the GBLTracks collection"<<std::endl;
                continue;
            }
            ele_trk = *ele_trk_ptr;
            pos_trk = *pos_trk_ptr;
        }
        else {
            ele_trk = ele->getTrack();
            pos_trk = pos->getTrack();
        }

        if (debug_) {
            std::cout << "Check Ele/Pos Track momenta" << std::endl;
            std::cout << ele_trk.getP() << " " << pos_trk.getP() << std::endl;
            std::cout << ele_trk.getOmega() << " " << pos_trk.getOmega() << std::endl;
        }

	    // Beam Position Corrections
        ele_trk.applyCorrection("z0", beamPosCorrections_.at(1));
        pos_trk.applyCorrection("z0", beamPosCorrections_.at(1));
        // Track Time Corrections
        double corr_eleTrackTime = ele_trk.getTrackTime() + eleTrackTimeBias_;
	double corr_posTrackTime = pos_trk.getTrackTime() + posTrackTimeBias_;

	    // Correct for the momentum bias
        if (biasingTool_) {
            // Correct for wrong track momentum - Bug Fix
            // In case there was mis-configuration during reco/hpstr-ntuple step, correct
            // the momentum magnitude here using the right bField for the data taking year
            if (bFieldScaleFactor_ > 0) {
                biasingTool_->updateWithBiasP(ele_trk,bFieldScaleFactor_);
                biasingTool_->updateWithBiasP(pos_trk,bFieldScaleFactor_);
            }

            biasingTool_->updateWithBiasP(ele_trk);
            biasingTool_->updateWithBiasP(pos_trk);
        }

        if (debug_) {
            std::cout << "Corrected Ele/Pos Track momenta" << std::endl;
            std::cout << ele_trk.getP() << " " << pos_trk.getP() << std::endl;
        }

        double invm_smear = 1.;
        if (smearingTool_) {
            double unsmeared_prod = ele_trk.getP()*pos_trk.getP();
            smearingTool_->updateWithSmearP(ele_trk);
            smearingTool_->updateWithSmearP(pos_trk);
            double smeared_prod = ele_trk.getP()*pos_trk.getP();
            invm_smear = sqrt(smeared_prod/unsmeared_prod);
        }

        double ele_E = ele->getEnergy();
        double pos_E = pos->getEnergy();
        if (debug_) std::cout << "got tracks" << std::endl;

        CalCluster eleClus = ele->getCluster();
        CalCluster posClus = pos->getCluster();

        // Compute analysis variables here.
        TLorentzVector p_ele;
        p_ele.SetPxPyPzE(ele_trk.getMomentum()[0], ele_trk.getMomentum()[1],
                         ele_trk.getMomentum()[2], ele->getEnergy());
        TLorentzVector p_pos;
        p_pos.SetPxPyPzE(pos_trk.getMomentum()[0], pos_trk.getMomentum()[1],
                         pos_trk.getMomentum()[2], ele->getEnergy());


        if (debug_) std::cout << "start selection" << std::endl;
        // Ele Track Time
        if (!vtxSelector->passCutLt("eleTrkTime_lt", fabs(corr_eleTrackTime), weight))
            continue;

        // Pos Track Time
        if (!vtxSelector->passCutLt("posTrkTime_lt", fabs(corr_posTrackTime), weight))
            continue;

        // Ele Track-cluster match
        if (!vtxSelector->passCutLt("eleTrkCluMatch_lt", ele->getGoodnessOfPID(), weight))
            continue;

        // Pos Track-cluster match
        if (!vtxSelector->passCutLt("posTrkCluMatch_lt", pos->getGoodnessOfPID(),weight))
            continue;

        //Require Positron Cluster exists
        if (!vtxSelector->passCutGt("posClusE_gt", posClus.getEnergy(), weight))
            continue;

        // Require Positron Cluster does NOT exists
        if (!vtxSelector->passCutLt("posClusE_lt", posClus.getEnergy(), weight))
            continue;


        double corr_eleClusterTime = ele->getCluster().getTime() - timeOffset_;
        double corr_posClusterTime = pos->getCluster().getTime() - timeOffset_;

        double botClusTime = 0.0;
        if (ele->getCluster().getPosition().at(1) < 0.0) botClusTime = ele->getCluster().getTime();
        else botClusTime = pos->getCluster().getTime();

        // Bottom Cluster Time
        if (!vtxSelector->passCutLt("botCluTime_lt", botClusTime, weight))
            continue;

        if (!vtxSelector->passCutGt("botCluTime_gt", botClusTime, weight))
            continue;

        // Ele Pos Cluster Time Difference
        if (!vtxSelector->passCutLt("eleposCluTimeDiff_lt", fabs(corr_eleClusterTime - corr_posClusterTime), weight))
            continue;

        //Ele Track-Cluster Time Difference
        if (!vtxSelector->passCutLt("eleTrkCluTimeDiff_lt", fabs(corr_eleClusterTime - corr_posClusterTime), weight))
            continue;

        //Pos Track-Cluster Time Difference
        if (!vtxSelector->passCutLt("posTrkCluTimeDiff_lt", fabs(corr_posClusterTime - corr_posClusterTime), weight))
            continue;

	// Ele Pos Track Time Difference
        if (!vtxSelector->passCutLt("eleposTrkTimeDiff_lt", fabs(corr_eleTrackTime - corr_posTrackTime), weight))
            continue;

        TVector3 ele_mom;
        ele_mom.SetX(ele_trk.getMomentum()[0]);
        ele_mom.SetY(ele_trk.getMomentum()[1]);
        ele_mom.SetZ(ele_trk.getMomentum()[2]);

        TVector3 pos_mom;
        pos_mom.SetX(pos_trk.getMomentum()[0]);
        pos_mom.SetY(pos_trk.getMomentum()[1]);
        pos_mom.SetZ(pos_trk.getMomentum()[2]);

        // Ele Track Quality - Chi2
        if (!vtxSelector->passCutLt("eleTrkChi2_lt", ele_trk.getChi2(), weight))
            continue;

        //Pos Track Quality - Chi2
        if (!vtxSelector->passCutLt("posTrkChi2_lt", pos_trk.getChi2(), weight))
            continue;

        //Ele Track Quality - Chi2Ndf
        if (!vtxSelector->passCutLt("eleTrkChi2Ndf_lt", ele_trk.getChi2Ndf(), weight))
            continue;

        //Pos Track Quality - Chi2Ndf
        if (!vtxSelector->passCutLt("posTrkChi2Ndf_lt", pos_trk.getChi2Ndf(), weight))
            continue;

        // Beam Electron cut
        if (!vtxSelector->passCutLt("eleMom_lt", ele_mom.Mag(), weight))
            continue;

        // Ele min momentum cut
        if (!vtxSelector->passCutGt("eleMom_gt", ele_mom.Mag(), weight))
            continue;

        // Pos min momentum cut
        if (!vtxSelector->passCutGt("posMom_gt", pos_mom.Mag(), weight))
            continue;

        // Ele nHits
        int ele2dHits = ele_trk.getTrackerHitCount();
        if (!ele_trk.isKalmanTrack())
            ele2dHits*=2;

        if (!vtxSelector->passCutGt("eleN2Dhits_gt", ele2dHits, weight))  {
            continue;
        }

        // Pos nHits
        int pos2dHits = pos_trk.getTrackerHitCount();
        if (!pos_trk.isKalmanTrack())
            pos2dHits*=2;

        if (!vtxSelector->passCutGt("posN2Dhits_gt", pos2dHits, weight))  {
            continue;
        }

        // Less than N shared hits for ele/pos track
        if (!vtxSelector->passCutLt("eleNshared_lt", ele_trk.getNShared(), weight)) {
            continue;
        }

        if (!vtxSelector->passCutLt("posNshared_lt", pos_trk.getNShared(), weight)) {
            continue;
        }

        // Vertex Quality
        if (!vtxSelector->passCutLt("chi2unc_lt", vtx->getChi2(), weight))
            continue;

        // Max vtx momentum
        if (!vtxSelector->passCutLt("maxVtxMom_lt", (ele_mom+pos_mom).Mag(), weight))
            continue;

        // Min vtx momentum
        if (!vtxSelector->passCutGt("minVtxMom_gt", (ele_mom+pos_mom).Mag(), weight))
            continue;

        passVtxPresel = true;

        selected_vtxs.push_back(vtx);
        vtxSelector->clearSelector();
    }

    if (debug_) std::cout << "start regions" << std::endl;
    for (auto region : _regions) {
        int nGoodVtx = 0;
        Vertex* goodVtx = nullptr;
        std::vector<Vertex*> goodVtxs;

        float truePsum = -1;
        float trueEsum = -1;

        for (auto vtx : selected_vtxs) {

            // No cuts.
            _reg_vtx_selectors[region]->getCutFlowHisto()->Fill(0., weight);

            Particle* ele = nullptr;
            Particle* pos = nullptr;

            _ah->GetParticlesFromVtx(vtx, ele, pos);

            CalCluster eleClus = ele->getCluster();
            CalCluster posClus = pos->getCluster();

            // vtx X position
            if (!_reg_vtx_selectors[region]->passCutLt("uncVtxX_lt", fabs(vtx->getX()), weight))
                continue;

            // vtx Y position
            if (!_reg_vtx_selectors[region]->passCutLt("uncVtxY_lt", fabs(vtx->getY()), weight))
                continue;

            // vtx Z position
            if (!_reg_vtx_selectors[region]->passCutGt("uncVtxZ_gt", vtx->getZ(), weight))
                continue;

            double ele_E = ele->getEnergy();
            double pos_E = pos->getEnergy();

            //Compute analysis variables here.

            // Track ele_trk = ele->getTrack();
            // Track pos_trk = pos->getTrack();
            Track* ele_trk_ptr;
            Track* pos_trk_ptr;
            Track ele_trk;
            Track pos_trk;
      
            if (!trkColl_.empty()) {
                bool foundTracks = _ah->MatchToGBLTracks((ele->getTrack()).getID(),(pos->getTrack()).getID(),
                        ele_trk_ptr, pos_trk_ptr, *trks_);

                if (!foundTracks) {
                    if(debug_) std::cout<<"VertexAnaProcessor::ERROR couldn't find ele/pos in the GBLTracks collection"<<std::endl;
                    continue;
                }
                ele_trk = *ele_trk_ptr;
                pos_trk = *pos_trk_ptr;
            }
            else {
                ele_trk = ele->getTrack();
                pos_trk = pos->getTrack();
            }

            //Beam Position Corrections
            ele_trk.applyCorrection("z0", beamPosCorrections_.at(1));
            pos_trk.applyCorrection("z0", beamPosCorrections_.at(1));
            //Track Time Corrections
            double corr_eleTrackTime = ele_trk.getTrackTime() + eleTrackTimeBias_;
	        double corr_posTrackTime = pos_trk.getTrackTime() + posTrackTimeBias_;
	    
	        if (biasingTool_) {
	            //Correct the wrong Bfield first
	            if (bFieldScaleFactor_ > 0) {
                    biasingTool_->updateWithBiasP(ele_trk, bFieldScaleFactor_);
                    biasingTool_->updateWithBiasP(pos_trk, bFieldScaleFactor_);
	            }
	      
                biasingTool_->updateWithBiasP(ele_trk);
                biasingTool_->updateWithBiasP(pos_trk);
            }
	    	    
            double invm_smear = 1.;
            if (smearingTool_) {
                double unsmeared_prod = ele_trk.getP()*pos_trk.getP();
                smearingTool_->updateWithSmearP(ele_trk);
                smearingTool_->updateWithSmearP(pos_trk);
                double smeared_prod = ele_trk.getP()*pos_trk.getP();
                invm_smear = sqrt(smeared_prod/unsmeared_prod);
            }

            TVector3 recEleP(ele->getMomentum()[0], ele->getMomentum()[1], ele->getMomentum()[2]);
            TLorentzVector p_ele;
            p_ele.SetPxPyPzE(ele_trk.getMomentum()[0], ele_trk.getMomentum()[1], ele_trk.getMomentum()[2], ele_E);
            TLorentzVector p_pos;
            p_pos.SetPxPyPzE(pos_trk.getMomentum()[0], pos_trk.getMomentum()[1], pos_trk.getMomentum()[2], pos_E);

            // Get the layers hit on each track
            std::vector<int> ele_hit_layers = ele_trk.getHitLayers();
            int ele_Si0 = 0;
            int ele_Si1 = 0;
            int ele_lastlayer = 0;
            for(int i=0; i < ele_hit_layers.size(); i++)
            {
                int layer = ele_hit_layers.at(i);
                ele_lastlayer = layer;
                if (layer == 0) ele_Si0++;
                if (layer == 1) ele_Si1++;
            }

            std::vector<int> pos_hit_layers = pos_trk.getHitLayers();
            int pos_Si0 = 0;
            int pos_Si1 = 0;
            int pos_lastlayer = 0;
            for(int i=0; i < pos_hit_layers.size(); i++)
            {
                int layer = pos_hit_layers.at(i);
                pos_lastlayer = layer;
                if (layer == 0) pos_Si0++;
                if (layer == 1) pos_Si1++;
            }

            // Defining these here so they are in scope elsewhere
            TVector3 trueEleP;
            TVector3 truePosP;

            if (debug_) {
                std::cout << "Check on ele_Track" << std::endl;
                std::cout << "Number of hits:" << ele_trk.getTrackerHitCount() << std::endl;
            }

            bool foundL1ele = false;
            bool foundL2ele = false;
            _ah->InnermostLayerCheck(&ele_trk, foundL1ele, foundL2ele);
            if (debug_) {
                std::cout << "Check on ele_Track" << std::endl;
                std::cout << "Number of hits:" << ele_trk.getTrackerHitCount() << std::endl;
            }

            bool foundL1pos = false;
            bool foundL2pos = false;
            _ah->InnermostLayerCheck(&pos_trk, foundL1pos, foundL2pos);
            if (debug_) {
                std::cout << "Check on pos_Track" << std::endl;
                std::cout << "Innermost:" << foundL1pos << " Second Innermost:" << foundL2pos << std::endl;
            }

            // PRESELECTION CUTS
            if (isData_) {
                if (!_reg_vtx_selectors[region]->passCutEq("Pair1_eq", (int)evth_->isPair1Trigger(), weight))
                    break;
		if (!_reg_vtx_selectors[region]->passCutEq("Singles0_eq", (int)ts_->isSingles0Trigger(), weight))
                    break;
            	if (!_reg_vtx_selectors[region]->passCutEq("Singles1_eq", (int)ts_->isSingles1Trigger(), weight))
                    break;
            	if (!_reg_vtx_selectors[region]->passCutEq("Singles2_eq", (int)ts_->isSingles2Trigger(), weight))
                    break;
            	if (!_reg_vtx_selectors[region]->passCutEq("Singles3_eq", (int)ts_->isSingles3Trigger(), weight))
		    break;
		if (!_reg_vtx_selectors[region]->passCutEq("Singles2or3_eq", (int)(ts_->isSingles2Trigger() || ts_->isSingles3Trigger()), weight))
                    break;
            }

            // Ele Track Time
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkTime_lt", fabs(corr_eleTrackTime), weight))
                continue;

            // Pos Track Time
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkTime_lt",fabs(corr_posTrackTime), weight))
                continue;

            // Ele Track-cluster match
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkCluMatch_lt", ele->getGoodnessOfPID(), weight))
                continue;

            // Pos Track-cluster match
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkCluMatch_lt", pos->getGoodnessOfPID(), weight))
                continue;

            // Require Positron Cluster exists
            if (!_reg_vtx_selectors[region]->passCutGt("posClusE_gt", posClus.getEnergy(), weight))
                continue;

            // Require Positron Cluster does NOT exists
            if (!_reg_vtx_selectors[region]->passCutLt("posClusE_lt", posClus.getEnergy(), weight))
                continue;

            double corr_eleClusterTime = ele->getCluster().getTime() - timeOffset_;
            double corr_posClusterTime = pos->getCluster().getTime() - timeOffset_;

            double botClusTime = 0.0;
            if (ele->getCluster().getPosition().at(1) < 0.0) {
                botClusTime = ele->getCluster().getTime();
            }
            else botClusTime = pos->getCluster().getTime();

            // Bottom Cluster Time
            if (!_reg_vtx_selectors[region]->passCutLt("botCluTime_lt", botClusTime, weight))
                continue;

            if (!_reg_vtx_selectors[region]->passCutGt("botCluTime_gt", botClusTime, weight))
                continue;

            // Ele Pos Cluster Time Difference
            if (!_reg_vtx_selectors[region]->passCutLt("eleposCluTimeDiff_lt", fabs(corr_eleClusterTime - corr_posClusterTime), weight))
                continue;

            // Ele Track-Cluster Time Difference
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkCluTimeDiff_lt", fabs(corr_eleTrackTime - corr_posClusterTime), weight))
                continue;

            // Pos Track-Cluster Time Difference
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkCluTimeDiff_lt", fabs(corr_posTrackTime - corr_posClusterTime), weight))
                continue;

	    // Ele Pos Track Time Difference
            if (!_reg_vtx_selectors[region]->passCutLt("eleposTrkTimeDiff_lt", fabs(corr_eleTrackTime - corr_posTrackTime), weight))
                continue;

            TVector3 ele_mom;
            ele_mom.SetX(ele_trk.getMomentum()[0]);
            ele_mom.SetY(ele_trk.getMomentum()[1]);
            ele_mom.SetZ(ele_trk.getMomentum()[2]);

            TVector3 pos_mom;
            pos_mom.SetX(pos_trk.getMomentum()[0]);
            pos_mom.SetY(pos_trk.getMomentum()[1]);
            pos_mom.SetZ(pos_trk.getMomentum()[2]);

            // Ele Track Quality - Chi2
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkChi2_lt", ele_trk.getChi2(), weight))
                continue;

            // Pos Track Quality - Chi2
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkChi2_lt", pos_trk.getChi2(), weight))
                continue;

            // Ele Track Quality - Chi2Ndf
            if (!_reg_vtx_selectors[region]->passCutLt("eleTrkChi2Ndf_lt", ele_trk.getChi2Ndf(), weight))
                continue;

            // Pos Track Quality - Chi2Ndf
            if (!_reg_vtx_selectors[region]->passCutLt("posTrkChi2Ndf_lt", pos_trk.getChi2Ndf(), weight))
                continue;

            // Beam Electron cut
            if (!_reg_vtx_selectors[region]->passCutLt("eleMom_lt", ele_mom.Mag(), weight))
                continue;

            // Ele min momentum cut
            if (!_reg_vtx_selectors[region]->passCutGt("eleMom_gt", ele_mom.Mag(), weight))
                continue;

            // Pos min momentum cut
            if (!_reg_vtx_selectors[region]->passCutGt("posMom_gt", pos_mom.Mag(), weight))
                continue;

            // Ele nHits
            int ele2dHits = ele_trk.getTrackerHitCount();
            if (!ele_trk.isKalmanTrack())
                ele2dHits*=2;

            if (!_reg_vtx_selectors[region]->passCutGt("eleN2Dhits_gt", ele2dHits, weight))  {
                continue;
            }

            // Pos nHits
            int pos2dHits = pos_trk.getTrackerHitCount();
            if (!pos_trk.isKalmanTrack())
                pos2dHits*=2;

            if (!_reg_vtx_selectors[region]->passCutGt("posN2Dhits_gt", pos2dHits, weight))  {
                continue;
            }

            //Less than N shared hits for ele/pos track
            if (!_reg_vtx_selectors[region]->passCutLt("eleNshared_lt", ele_trk.getNShared(), weight)) {
                continue;
            }

            if (!_reg_vtx_selectors[region]->passCutLt("posNshared_lt", pos_trk.getNShared(), weight)) {
                continue;
            }

            // Vertex Quality
            if (!_reg_vtx_selectors[region]->passCutLt("chi2unc_lt", vtx->getChi2(), weight))
                continue;

            // Max vtx momentum
            if (!_reg_vtx_selectors[region]->passCutLt("maxVtxMom_lt", (ele_mom + pos_mom).Mag(), weight))
                continue;

            // Min vtx momentum
            if (!_reg_vtx_selectors[region]->passCutGt("minVtxMom_gt", (ele_mom + pos_mom).Mag(), weight))
                continue;

            // END PRESELECTION CUTS

            // L1 requirement
            if (!_reg_vtx_selectors[region]->passCutEq("L1Requirement_eq", (int)(foundL1ele && foundL1pos), weight))
                continue;

            // L2 requirement
            if (!_reg_vtx_selectors[region]->passCutEq("L2Requirement_eq", (int)(foundL2ele && foundL2pos), weight))
                continue;

            // L1 requirement for positron
            if (!_reg_vtx_selectors[region]->passCutEq("L1PosReq_eq", (int)(foundL1pos), weight))
                continue;

            // ESum low cut
            if (!_reg_vtx_selectors[region]->passCutLt("eSum_lt", (ele_E + pos_E), weight))
                continue;

            // ESum high cut
            if (!_reg_vtx_selectors[region]->passCutGt("eSum_gt", (ele_E + pos_E), weight))
                continue;

            // PSum low cut
            if (!_reg_vtx_selectors[region]->passCutLt("pSum_lt", (p_ele.P() + p_pos.P()), weight))
                continue;

            // PSum high cut
            if (!_reg_vtx_selectors[region]->passCutGt("pSum_gt", (p_ele.P() + p_pos.P()), weight))
                continue;

            // Require Electron Cluster exists
            if (!_reg_vtx_selectors[region]->passCutGt("eleClusE_gt", eleClus.getEnergy(), weight))
                continue;


            // Require Electron Cluster does NOT exists
            if (!_reg_vtx_selectors[region]->passCutLt("eleClusE_lt", eleClus.getEnergy(), weight))
                continue;

            // No shared hits requirement
            if (!_reg_vtx_selectors[region]->passCutEq("ele_sharedL0_eq", (int)ele_trk.getSharedLy0(), weight))
                continue;
            if (!_reg_vtx_selectors[region]->passCutEq("pos_sharedL0_eq", (int)pos_trk.getSharedLy0(), weight))
                continue;
            if (!_reg_vtx_selectors[region]->passCutEq("ele_sharedL1_eq", (int)ele_trk.getSharedLy1(), weight))
                continue;
            if (!_reg_vtx_selectors[region]->passCutEq("pos_sharedL1_eq", (int)pos_trk.getSharedLy1(), weight))
                continue;

            // Min vtx Y pos
            if (!_reg_vtx_selectors[region]->passCutGt("VtxYPos_gt", vtx->getY(), weight))
                continue;

            // Max vtx Y pos
            if (!_reg_vtx_selectors[region]->passCutLt("VtxYPos_lt", vtx->getY(), weight))
                continue;

            // Tracking Volume for positron
            if (!_reg_vtx_selectors[region]->passCutGt("volPos_top", p_pos.Py(), weight))
                continue;

            if (!_reg_vtx_selectors[region]->passCutLt("volPos_bot", p_pos.Py(), weight))
                continue;

            if (!_reg_vtx_selectors[region]->passCutLt("deltaZ_lt", std::abs((ele_trk.getZ0()/ele_trk.getTanLambda()) - (pos_trk.getZ0()/pos_trk.getTanLambda())), weight))
                continue;

            // If this is MC check if MCParticle matched to the electron track is from rad or recoil
            if (!isData_)
            {
                // Count the number of hits per part on the ele track
                std::map<int, int> nHits4part;
                for(int i =0; i < ele_trk.getMcpHits().size(); i++)
                {
                    int partID = ele_trk.getMcpHits().at(i).second;
                    if ( nHits4part.find(partID) == nHits4part.end() )
                    {
                        // not found
                        nHits4part[partID] = 1;
                    }
                    else
                    {
                        // found
                        nHits4part[partID]++;
                    }
                }

                // Determine the MC part with the most hits on the track
                int maxNHits = 0;
                int maxID = 0;
                for (std::map<int,int>::iterator it = nHits4part.begin(); it != nHits4part.end(); ++it)
                {
                    if (it->second > maxNHits)
                    {
                        maxNHits = it->second;
                        maxID = it->first;
                    }
                }

                // Find the correct mc part and grab mother id
                int isRadEle = -999;
                int isRecEle = -999;


                trueEleP.SetXYZ(-999,-999,-999);
                truePosP.SetXYZ(-999,-999,-999);

                if (mcParts_) {
                    float trueEleE = -1;
                    float truePosE = -1;
                    for (int i = 0; i < mcParts_->size(); i++)
                    {
                        int momPDG = mcParts_->at(i)->getOriginPDG();
                        if (mcParts_->at(i)->getPDG() == 11 && momPDG == isRadPDG_)
                        {
                            std::vector<double> lP = mcParts_->at(i)->getMomentum();
                            trueEleP.SetXYZ(lP[0],lP[1],lP[2]);
                            trueEleE = mcParts_->at(i)->getEnergy();
                        }
                        if (mcParts_->at(i)->getPDG() == -11 && momPDG == isRadPDG_)
                        {
                            std::vector<double> lP = mcParts_->at(i)->getMomentum();
                            truePosP.SetXYZ(lP[0],lP[1],lP[2]);
                            truePosE = mcParts_->at(i)->getEnergy();
                        }
                        if (trueEleP.X() != -999 && truePosP.X() != -999){
                            truePsum =  trueEleP.Mag() + trueEleP.Mag();
                            trueEsum = trueEleE + truePosE;
                        }

                        if (mcParts_->at(i)->getID() != maxID) continue;
                        // Default isRadPDG = 622
                        if (momPDG == isRadPDG_) isRadEle = 1;
                        if (momPDG == 623) isRecEle = 1;
                    }
                }
                double momRatio = recEleP.Mag() / trueEleP.Mag();
                double momAngle = trueEleP.Angle(recEleP) * TMath::RadToDeg();
                if (!_reg_vtx_selectors[region]->passCutLt("momRatio_lt", momRatio, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutGt("momRatio_gt", momRatio, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutLt("momAngle_lt", momAngle, weight)) continue;

                if (!_reg_vtx_selectors[region]->passCutEq("isRadEle_eq", isRadEle, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutEq("isNotRadEle_eq", isRadEle, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutEq("isRecEle_eq", isRecEle, weight)) continue;
            }

            goodVtx = vtx;
            nGoodVtx++;
            goodVtxs.push_back(vtx);
        } // selected vertices

        // N selected vertices - this is quite a silly cut to make at the end. But okay. that's how we decided atm.
        if (!_reg_vtx_selectors[region]->passCutEq("nVtxs_eq", nGoodVtx, weight))
            continue;

        _reg_tuples[region]->setVariableValue("nGoodVtx", nGoodVtx);

        // Loop over all selected vertices in the region
        for (std::vector<Vertex*>::iterator it = goodVtxs.begin(); it != goodVtxs.end(); it++) {

            Vertex* vtx = *it;

            Particle* ele = nullptr;
            Particle* pos = nullptr;

            if (!vtx || !_ah->GetParticlesFromVtx(vtx, ele, pos))
                continue;

            CalCluster eleClus = ele->getCluster();
            CalCluster posClus = pos->getCluster();

            double corr_eleClusterTime = ele->getCluster().getTime() - timeOffset_;
            double corr_posClusterTime = pos->getCluster().getTime() - timeOffset_;

            double ele_E = ele->getEnergy();
            double pos_E = pos->getEnergy();

            //Compute analysis variables here.
            // Track ele_trk = ele->getTrack();
            // Track pos_trk = pos->getTrack();
            //Get the shared info - TODO change and improve

            Track ele_trk;
            Track pos_trk;
            Track* ele_trk_ptr;
            Track* pos_trk_ptr;
        
            if (!trkColl_.empty()) {
                bool foundTracks = _ah->MatchToGBLTracks((ele->getTrack()).getID(),(pos->getTrack()).getID(),
                        ele_trk_ptr, pos_trk_ptr, *trks_);

                if (!foundTracks) {
                    if(debug_) std::cout<<"VertexAnaProcessor::ERROR couldn't find ele/pos in the GBLTracks collection"<<std::endl;
                    continue;
                }
                ele_trk = *ele_trk_ptr;
                pos_trk = *pos_trk_ptr;
            }
            else {
                ele_trk = ele->getTrack();
                pos_trk = pos->getTrack();
            }
            
            //Track Time Corrections
            double corr_eleTrackTime = ele_trk.getTrackTime() + eleTrackTimeBias_;
            double corr_posTrackTime = pos_trk.getTrackTime() + posTrackTimeBias_;

            // Track Momentum bias
            if (biasingTool_) {
                // Correct for wrong track momentum - Bug Fix
                // In case there was mis-configuration during reco/hpstr-ntuple step, correct
                // the momentum magnitude here using the right bField for the data taking year
                if (bFieldScaleFactor_ > 0) {
                    biasingTool_->updateWithBiasP(ele_trk, bFieldScaleFactor_);
                    biasingTool_->updateWithBiasP(pos_trk, bFieldScaleFactor_);
                }
                biasingTool_->updateWithBiasP(ele_trk);
                biasingTool_->updateWithBiasP(pos_trk);
            }

            double invm_smear = 1.;
            if (smearingTool_) {
                double unsmeared_prod = ele_trk.getP()*pos_trk.getP();
                smearingTool_->updateWithSmearP(ele_trk);
                smearingTool_->updateWithSmearP(pos_trk);
                double smeared_prod = ele_trk.getP()*pos_trk.getP();
                invm_smear = sqrt(smeared_prod/unsmeared_prod);
            }

            // Get the layers hit on each track
            std::vector<int> ele_hit_layers = ele_trk.getHitLayers();
            int ele_Si0 = 0;
            int ele_Si1 = 0;
            int ele_lastlayer = 0;
            for (int i=0; i<ele_hit_layers.size(); i++)
            {
                int layer = ele_hit_layers.at(i);
                ele_lastlayer = layer;
                if (layer == 0) ele_Si0++;
                if (layer == 1) ele_Si1++;
            }

            std::vector<int> pos_hit_layers = pos_trk.getHitLayers();
            int pos_Si0 = 0;
            int pos_Si1 = 0;
            int pos_lastlayer = 0;
            for (int i=0; i<pos_hit_layers.size(); i++)
            {
                int layer = pos_hit_layers.at(i);
                pos_lastlayer = layer;
                if (layer == 0) pos_Si0++;
                if (layer == 1) pos_Si1++;
            }

            // Vertex Covariance
            std::vector<float> vtx_cov = vtx->getCovariance();
            float cxx = vtx_cov.at(0);
            float cyx = vtx_cov.at(1);
            float cyy = vtx_cov.at(2);
            float czx = vtx_cov.at(3);
            float czy = vtx_cov.at(4);
            float czz = vtx_cov.at(5);

            // MC Truth hits in first 4 sensors
            int L1L2hitCode = 0; //hit code '1111' means truth ax+ster hits in L1_ele, L1_pos, L2_ele, L2_pos
            int L1hitCode = 0; //hit code '1111' means truth in L1_ele_ax, L1_ele_ster, L1_pos_ax, L1_pos_ster
            int L2hitCode = 0; // hit code '1111' means truth in L2_ele_ax, L2_ele_ster, L2_pos_ax, L2_pos_ster
            if (!isData_){
                // Get hit codes. Only sure this works for 2016 KF as is.
                utils::get2016KFMCTruthHitCodes(&ele_trk, &pos_trk, L1L2hitCode, L1hitCode, L2hitCode);
                // L1L2 truth hit selection
                if (!_reg_vtx_selectors[region]->passCutLt("hitCode_lt", ((double)L1L2hitCode)-0.5, weight)) continue;
                if (!_reg_vtx_selectors[region]->passCutGt("hitCode_gt", ((double)L1L2hitCode)+0.5, weight)) continue;
            }

            // track isolations
            // Only calculate isolations if both track L1 and L2 hits exist
            bool hasL1ele = false;
            bool hasL2ele = false;
            _ah->InnermostLayerCheck(&ele_trk, hasL1ele, hasL2ele);

            bool hasL1pos = false;
            bool hasL2pos = false;
            _ah->InnermostLayerCheck(&pos_trk, hasL1pos, hasL2pos);

            double ele_trk_iso_L1 = 99999.9;
            double pos_trk_iso_L1 = 99999.9;
            if(hasL1ele && hasL2ele && hasL1pos && hasL2pos){
                if (ele_trk.isKalmanTrack()){
                    ele_trk_iso_L1 = utils::getKalmanTrackL1Isolations(&ele_trk, hits_);
                    pos_trk_iso_L1 = utils::getKalmanTrackL1Isolations(&pos_trk, hits_);
                }
            }

            TVector3 ele_mom;
            ele_mom.SetX(ele_trk.getMomentum()[0]);
            ele_mom.SetY(ele_trk.getMomentum()[1]);
            ele_mom.SetZ(ele_trk.getMomentum()[2]);

            TVector3 pos_mom;
            pos_mom.SetX(pos_trk.getMomentum()[0]);
            pos_mom.SetY(pos_trk.getMomentum()[1]);
            pos_mom.SetZ(pos_trk.getMomentum()[2]);

            double psum = ele_mom.Mag()+pos_mom.Mag();

            // Ele nHits
            int ele2dHits = ele_trk.getTrackerHitCount();
            if (!ele_trk.isKalmanTrack())
                ele2dHits*=2;

            // pos nHits
            int pos2dHits = pos_trk.getTrackerHitCount();

            TVector3 recEleP(ele->getMomentum()[0], ele->getMomentum()[1], ele->getMomentum()[2]);
            TLorentzVector p_ele;
            p_ele.SetPxPyPzE(ele_trk.getMomentum()[0], ele_trk.getMomentum()[1], ele_trk.getMomentum()[2], ele_E);
            TLorentzVector p_pos;
            p_pos.SetPxPyPzE(pos_trk.getMomentum()[0], pos_trk.getMomentum()[1], pos_trk.getMomentum()[2], pos_E);

            double reconz = vtx->getZ(); 
            double ele_trk_z0 = ele_trk.getZ0();
            double ele_trk_z0err = ele_trk.getZ0Err();
            double pos_trk_z0 = pos_trk.getZ0();
            double pos_trk_z0err = pos_trk.getZ0Err();

            //DeltaZ
            double deltaZ = std::abs( (ele_trk_z0/ele_trk.getTanLambda()) - (pos_trk_z0/pos_trk.getTanLambda()) );
            
            //Project vertex to target
            double vtx_proj_x = -999.9;
            double vtx_proj_y = -999.9;
            double vtx_proj_x_sig = -999.9;
            double vtx_proj_y_sig = -999.9;
            double vtx_proj_sig = -999.9;
            if (!v0ProjectionFitsCfg_.empty()) {
                vtx_proj_sig = utils::v0_projection_to_target_significance(v0proj_fits_, evth_->getRunNumber(),
                        vtx_proj_x, vtx_proj_y, vtx_proj_x_sig, vtx_proj_y_sig, vtx->getX(), vtx->getY(),
                        reconz, vtx->getP().X(), vtx->getP().Y(), vtx->getP().Z());
            }

            // TODO put this in the Vertex!
            TVector3 vtxPosSvt;
            vtxPosSvt.SetX(vtx->getX());
            vtxPosSvt.SetY(vtx->getY());
            vtxPosSvt.SetZ(vtx->getZ());
            vtxPosSvt.RotateY(-0.0305);

            int momPDG_ele = -999;
            int momPDG_pos = -999;
            int originPDG_ele = -999;
            int originPDG_pos = -999;

            if (mcParts_) {    
                for (int i = 0; i < mcParts_->size(); i++) {    
                    int momPDG = mcParts_->at(i)->getMomPDG();
                    int originPDG = mcParts_->at(i)->getOriginPDG();
                    if (mcParts_->at(i)->getPDG() == 11) {
                        momPDG_ele = momPDG;
                        originPDG_ele = originPDG;
                    }
                    if (mcParts_->at(i)->getPDG() == -11) {
                        momPDG_pos = momPDG;
                        originPDG_pos = originPDG;
                    }
                }
            }
            
            // Just for the selected vertex
            if (!isData_){
                _reg_tuples[region]->setVariableValue("ap_true_vtx_z", apZ);
                _reg_tuples[region]->setVariableValue("ap_true_vtx_mass", apMass);
                _reg_tuples[region]->setVariableValue("ap_true_vtx_energy", apEnergy);
                // _reg_tuples[region]->setVariableValue("vd_true_vtx_z", vdZ);
                // _reg_tuples[region]->setVariableValue("vd_true_vtx_mass", vdMass);
                // _reg_tuples[region]->setVariableValue("vd_true_vtx_energy", vdEnergy);
                _reg_tuples[region]->setVariableValue("hitCode", float(L1L2hitCode));
                _reg_tuples[region]->setVariableValue("L1hitCode", float(L1hitCode));
                _reg_tuples[region]->setVariableValue("L2hitCode", float(L2hitCode));
                _reg_tuples[region]->setVariableValue("momPDG_ele", int(momPDG_ele));
                _reg_tuples[region]->setVariableValue("momPDG_pos", int(momPDG_pos));
                _reg_tuples[region]->setVariableValue("originPDG_ele", int(originPDG_ele));
                _reg_tuples[region]->setVariableValue("originPDG_pos", int(originPDG_pos));
            }

            _reg_tuples[region]->setVariableValue("unc_vtx_mass", vtx->getInvMass());
            _reg_tuples[region]->setVariableValue("unc_vtx_z"   , vtxPosSvt.Z());
            _reg_tuples[region]->setVariableValue("unc_vtx_chi2", vtx->getChi2());
            _reg_tuples[region]->setVariableValue("unc_vtx_psum", p_ele.P() + p_pos.P());
	    _reg_tuples[region]->setVariableValue("unc_vtx_px", vtx->getP().X());
            _reg_tuples[region]->setVariableValue("unc_vtx_py", vtx->getP().Y());
            _reg_tuples[region]->setVariableValue("unc_vtx_pz", vtx->getP().Z());
            _reg_tuples[region]->setVariableValue("unc_vtx_x", vtx->getX());
            _reg_tuples[region]->setVariableValue("unc_vtx_y", vtx->getY());
            _reg_tuples[region]->setVariableValue("unc_vtx_proj_x", vtx_proj_x);
            _reg_tuples[region]->setVariableValue("unc_vtx_proj_y", vtx_proj_y);
            _reg_tuples[region]->setVariableValue("unc_vtx_proj_x_sig", vtx_proj_x_sig);
            _reg_tuples[region]->setVariableValue("unc_vtx_proj_y_sig", vtx_proj_y_sig);
            _reg_tuples[region]->setVariableValue("unc_vtx_proj_sig", vtx_proj_sig);
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_pos_clust_dt", corr_eleClusterTime - corr_posClusterTime);

            _reg_tuples[region]->setVariableValue("unc_vtx_cxx", cxx);
            _reg_tuples[region]->setVariableValue("unc_vtx_cyy", cyy);
            _reg_tuples[region]->setVariableValue("unc_vtx_czz", czz);
            _reg_tuples[region]->setVariableValue("unc_vtx_cyx", cyx);
            _reg_tuples[region]->setVariableValue("unc_vtx_czy", czy);
            _reg_tuples[region]->setVariableValue("unc_vtx_czx", czx);
            _reg_tuples[region]->setVariableValue("unc_vtx_deltaZ", deltaZ);

            //track vars
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_p", ele_trk.getP());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_t", corr_eleTrackTime);
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_d0", ele_trk.getD0());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_phi0", ele_trk.getPhi());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_omega", ele_trk.getOmega());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_tanLambda", ele_trk.getTanLambda());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_z0", ele_trk.getZ0());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_chi2ndf", ele_trk.getChi2Ndf());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_clust_dt", corr_eleTrackTime - corr_posClusterTime);
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_z0Err",ele_trk.getZ0Err());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_d0Err", ele_trk.getD0Err());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_tanLambdaErr", ele_trk.getTanLambdaErr());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_PhiErr", ele_trk.getPhiErr());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_OmegaErr", ele_trk.getOmegaErr());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_nhits",ele2dHits);
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_lastlayer",ele_lastlayer);
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_si0",ele_Si0);
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_si1",ele_Si1);
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_L1_isolation", ele_trk_iso_L1);

            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_p", pos_trk.getP());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_t", corr_posTrackTime);
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_d0", pos_trk.getD0());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_phi0", pos_trk.getPhi());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_omega", pos_trk.getOmega());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_tanLambda", pos_trk.getTanLambda());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_z0", pos_trk.getZ0());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_chi2ndf", pos_trk.getChi2Ndf());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_clust_dt", corr_posTrackTime - corr_posClusterTime);
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_z0Err",pos_trk.getZ0Err());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_d0Err", pos_trk.getD0Err());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_tanLambdaErr", pos_trk.getTanLambdaErr());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_PhiErr", pos_trk.getPhiErr());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_OmegaErr", pos_trk.getOmegaErr());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_nhits",pos2dHits);
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_lastlayer",pos_lastlayer);
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_si0",pos_Si0);
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_si1",pos_Si1);
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_L1_isolation", pos_trk_iso_L1);

            //clust vars
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_clust_E", eleClus.getEnergy());
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_clust_x", eleClus.getPosition().at(0));
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_clust_corr_t",corr_eleClusterTime);

            _reg_tuples[region]->setVariableValue("unc_vtx_pos_clust_E", posClus.getEnergy());
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_clust_x", posClus.getPosition().at(0));
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_clust_corr_t",corr_posClusterTime);
            _reg_tuples[region]->setVariableValue("run_number", evth_->getRunNumber());
	    _reg_tuples[region]->setVariableValue("event_number", evth_->getEventNumber());
	    _reg_tuples[region]->setVariableValue("singles0trigger", (int)ts_->isSingles0Trigger());
	    _reg_tuples[region]->setVariableValue("singles1trigger", (int)ts_->isSingles1Trigger());
	    _reg_tuples[region]->setVariableValue("singles2trigger", (int)ts_->isSingles2Trigger());
	    _reg_tuples[region]->setVariableValue("singles3trigger", (int)ts_->isSingles3Trigger());

            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_ecal_x", ele_trk.getPositionAtEcal().at(0));
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_ecal_y", ele_trk.getPositionAtEcal().at(1));
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_z", ele_trk.getPosition().at(2));
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_ecal_x", pos_trk.getPositionAtEcal().at(0));
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_ecal_y", pos_trk.getPositionAtEcal().at(1));
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_z", pos_trk.getPosition().at(2));
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_px", ele_trk.getMomentum().at(0));
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_py", ele_trk.getMomentum().at(1));
            _reg_tuples[region]->setVariableValue("unc_vtx_ele_track_pz", ele_trk.getMomentum().at(2));
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_px", pos_trk.getMomentum().at(0));
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_py", pos_trk.getMomentum().at(1));
            _reg_tuples[region]->setVariableValue("unc_vtx_pos_track_pz", pos_trk.getMomentum().at(2));

            _reg_tuples[region]->fill();   
        }

    }// regions

    return true;
}

void FlatVertexProcessor::finalize() {
    outF_->cd();
    vtxSelector->getCutFlowHisto()->Write();

    for (auto region : _regions) {
        // std::string dirName = anaName_ + "_" + region;
        // outF_->cd(dirName.c_str());
        _reg_vtx_selectors[region]->getCutFlowHisto()->Write();
        _reg_tuples[region]->writeTree();
    }

    outF_->Close();
}

DECLARE_PROCESSOR(FlatVertexProcessor);
