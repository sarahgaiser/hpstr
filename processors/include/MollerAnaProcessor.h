#ifndef __MOLLER_ANAPROCESSOR_H__
#define __MOLLER_ANAPROCESSOR_H__


//HPSTR
#include "HpsEvent.h"
#include "TSData.h"
#include "Collections.h"
#include "Track.h"
#include "CalCluster.h"
#include "CalHit.h"
#include "Particle.h"
#include "Vertex.h"
#include "Processor.h"
#include "HistoManager.h"
#include "FlatTupleMaker.h"
#include "AnaHelpers.h"
#include "BaseSelector.h"

#include "MollerAnaHistos.h"

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TRefArray.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TF1.h"

//C++
#include <memory>
#include "math.h"


struct char_cmp {
    bool operator () (const char *a,const char *b) const
    {
        return strcmp(a,b)<0;
    }
};


class MollerAnaProcessor : public Processor {

    public:
		MollerAnaProcessor(const std::string& name, Process& process);
        ~MollerAnaProcessor();
        virtual bool process(IEvent* ievent);

        virtual void initialize(TTree* tree);

        virtual void finalize();

        virtual void configure(const ParameterSet& parameters);

    private:
        std::shared_ptr<BaseSelector> trackSelector;
        std::string trackSelectionCfg_;

        std::shared_ptr<BaseSelector> vtxSelector;
        std::string vertexSelectionCfg_;

        std::vector<std::string> regionSelections_; //!< description

        //Containers to hold histogrammer info
        MollerAnaHistos* trackHistos{nullptr};
        std::string  histTrackCfgFilename_;

        MollerAnaHistos* vertexHistos{nullptr};
        std::string  histVertexCfgFilename_;

        TTree* tree_{nullptr};
        TBranch* btsData_{nullptr}; //!< TS Bank branch
        TBranch* btrks_{nullptr}; //!< tracks branch 
        TBranch* bvtxs_{nullptr}; //!< vertices branch

        TSData* tsData_{};
        std::vector<Track*>* trks_{};
        std::vector<Vertex*>* vtxs_{};

        std::string anaName_{"mollerAna"};
        std::string tsColl_{"TSBank"};
        std::string trkColl_{"KalmanFullTracks"};
        std::string vtxColl_{"UnconstrainedMollerVertices"};

        std::map<std::string, std::shared_ptr<BaseSelector>> _reg_selectors; //!< description
        std::map<std::string, std::shared_ptr<TrackHistos>> _reg_histos; //!< description
        std::map<std::string, std::shared_ptr<FlatTupleMaker>> _reg_tuples; //!< description

        std::vector<std::string> _regions; //!< description

        typedef std::map<std::string, std::shared_ptr<TrackHistos>>::iterator reg_it; //!< description

        double beamE_{3.7};
        int isData_{1};
        int makeFlatTuple_{0}; //!< make true in config to save flat tuple
        std::string analysis_{"MollerAnalysis"};

        //Debug level
        int debug_{0};

        std::shared_ptr<AnaHelpers> _ah;

        // Kinematic equations
        TF1* func_E_vs_theta_after_roation;
        TF1* func_theta1_vs_theta2_after_roation;

        // save a tree for information of tracks from vertices
        // std::shared_ptr<FlatTupleMaker> _reg_tracks_from_vertices;

};

#endif
