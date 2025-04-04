#include "MCAnaHistos.h"
#include <math.h>

void MCAnaHistos::Define2DHistos() {
    std::string h_name = "";
    for (auto hist : _h_configs.items()) {
        if (hist.key() == "pos_pxpy_hh")
        {
            for (int pxz = hist.value().at("lowPxz");
                 pxz < hist.value().at("highPxz");
                 pxz += (int) hist.value().at("stepPxz"))
            {
                h_name = m_name + "_pos_pxpy_" + std::to_string(pxz) + "_hh";
                histos2d[h_name] = plot2D(h_name, hist.value().at("xtitle"),
                                          hist.value().at("binsX"), hist.value().at("minX"),
                                          hist.value().at("maxX"),  hist.value().at("ytitle"),
                                          hist.value().at("binsY"), hist.value().at("minY"),
                                          hist.value().at("maxY"));
            }
        }
        if (hist.key() == "ele_pxpy_hh")
        {
            for (int pxz = hist.value().at("lowPxz");
                 pxz < hist.value().at("highPxz");
                 pxz += (int) hist.value().at("stepPxz"))
            {
                h_name = m_name + "_ele_pxpy_" + std::to_string(pxz) + "_hh";
                histos2d[h_name] = plot2D(h_name, hist.value().at("xtitle"),
                                          hist.value().at("binsX"), hist.value().at("minX"),
                                          hist.value().at("maxX"),  hist.value().at("ytitle"),
                                          hist.value().at("binsY"), hist.value().at("minY"),
                                          hist.value().at("maxY"));
            }
        }
    }
}

void MCAnaHistos::FillMCParticles(std::vector<MCParticle*> *mcParts, std::string analysis, float weight) {
    if(mcParts == nullptr)
        std::cout << "MCPARTS IS NULL" << std::endl;
    int nParts = mcParts->size();
    Fill1DHisto("numMCparts_h", (float)nParts, weight);
    int nMuons = 0;
    int nElec = 0;
    int nPos = 0;
    int nGamma = 0;
    double minMuonE = -99.9;

    TLorentzVector ele;
    TLorentzVector pos;

    for (int i = 0; i < nParts; i++)
    {
        MCParticle *part = mcParts->at(i);
        int pdg = part->getPDG();
        int momPdg = part->getMomPDG();
	//	if ( pdg > 600)
	//  std::cout<<"Found particle with momPDG = "<<momPdg<<" part = " << pdg << " mass " << part->getMass() << std::endl;

        double energy = part->getEnergy();
        double massMeV = 1000.0*part->getMass();
        double zPos = part->getVertexPosition().at(2);
        std::vector<double> partP = part->getMomentum();
        TLorentzVector part4P(partP.at(0), partP.at(1), partP.at(2), energy);
        part4P.RotateY(-0.0305);
        double momentum = part4P.P();

        if (momPdg == 623)
        {
            Fill1DHisto("recoil_ele_p_h",momentum, weight);
            Fill1DHisto("recoil_ele_px_h",part4P.X(), weight);
            Fill1DHisto("recoil_ele_py_h",part4P.Y(), weight);
            Fill1DHisto("recoil_ele_pz_h",part4P.Z(), weight);
        }

        if (pdg == 622)
        {
            Fill1DHisto("mc622Mass_h", massMeV, weight);
            Fill1DHisto("mc622Z_h", zPos, weight);
            Fill1DHisto("mc622Energy_h", energy, weight);
            Fill1DHisto("mc622P_h", momentum, weight);
        }

        if (pdg == 625)
        {
            Fill1DHisto("mc625Mass_h", massMeV, weight);
            Fill1DHisto("mc625Z_h", zPos, weight);
            Fill1DHisto("mc625Energy_h", energy, weight);
            Fill1DHisto("mc625P_h", momentum, weight);
        }

        if (pdg == 624)
        {
            Fill1DHisto("mc624Mass_h", massMeV, weight);
            Fill1DHisto("mc624Z_h", zPos, weight);
            Fill1DHisto("mc625Energy_h", energy, weight);
            Fill1DHisto("mc625P_h", momentum, weight);
        }


        if (fabs(pdg) == 13)
        {
            nMuons++;
            if (energy < minMuonE || minMuonE < 0.0)
            {
                minMuonE = energy;
            }
        }

        bool partOfInt = false;
        //	std::cout<<analysis<<std::endl;
        if (analysis == "simps"){
            if (fabs(pdg) == 11 && momPdg == 625)
                partOfInt = true;
        } else {
            if ((momPdg == 623 || momPdg == 622) && (fabs(pdg) == 11))
                partOfInt = true;
        }

        if (partOfInt == true)
        {
            double PperpB = 1000.0*sqrt( (part4P.Px()*part4P.Px()) + (part4P.Pz()*part4P.Pz()) );
            int Pxz = int(floor(PperpB));
            int round = Pxz%100;
            Pxz = Pxz - round;
            if (pdg == 11)
            {
                Fill1DHisto("ele_pxz_h", PperpB, weight);
                Fill2DHisto("ele_pxpy_" + std::to_string(Pxz) + "_hh", part4P.Px(), part4P.Py(), weight);
                
                Fill1DHisto("truthRadElecE_h", energy, weight);
                Fill1DHisto("truthRadEleczPos_h", zPos, weight);
                Fill1DHisto("truthRadElecPt_h", part4P.Pt(), weight);
                Fill1DHisto("truthRadElecPz_h", part4P.Pz(), weight);

                ele = part4P;
            }
            if (pdg == -11)
            {
                Fill1DHisto("pos_pxz_h", PperpB, weight);
                Fill2DHisto("pos_pxpy_" + std::to_string(Pxz) + "_hh", part4P.Px(), part4P.Py(), weight);

                Fill1DHisto("truthRadPosE_h", energy, weight);
                Fill1DHisto("truthRadPoszPos_h", zPos, weight);
                Fill1DHisto("truthRadPosPt_h", part4P.Pt(), weight);
                Fill1DHisto("truthRadPosPz_h", part4P.Pz(), weight);

                pos = part4P;
            }
        }

        if (analysis == "beam") {
            if (pdg == 11) {
                nElec++;
                Fill1DHisto("truthElecE_h", energy, weight);
                Fill1DHisto("truthElecPt_h", part4P.Pt(), weight);
                Fill1DHisto("truthElecPz_h", part4P.Pz(), weight);
            }

            if (pdg == -11) {
                nPos++;
                Fill1DHisto("truthPosE_h", energy, weight);
            
            }

            if (pdg == 22) {
                nGamma++;
                Fill1DHisto("truthGammaE_h", energy, weight);
                Fill1DHisto("truthGammaELow_h", energy*1000.0, weight);// Scaled to MeV
            }
        }

        Fill1DHisto("MCpartsEnergy_h", energy, weight);
        Fill1DHisto("MCpartsEnergyLow_h", energy*1000.0, weight);// Scaled to MeV
    }

    //TLorentzVector res = ele + pos;
    //std::cout<<" My resonance mass is "<< res.M()<< std::endl;

    Fill1DHisto("numMuons_h", nMuons, weight);
    Fill1DHisto("minMuonE_h", minMuonE, weight);
    Fill1DHisto("minMuonEhigh_h", minMuonE, weight);

    Fill1DHisto("numElectrons_h", nElec, weight);
    Fill1DHisto("numPositrons_h", nPos, weight);
    Fill1DHisto("numGammas_h", nGamma, weight);
}

void MCAnaHistos::FillMCTrackerHits(std::vector<MCTrackerHit*> *mcTrkrHits, float weight ) {
    int nHits = mcTrkrHits->size();
    Fill1DHisto("numMCTrkrHit_h", nHits, weight);
    for (int i=0; i < nHits; i++)
    {
        MCTrackerHit *hit = mcTrkrHits->at(i);
        int pdg = hit->getPDG();
        Fill1DHisto("mcTrkrHitEdep_h", hit->getEdep()*1000.0, weight); // Scaled to MeV
        Fill1DHisto("mcTrkrHitPdgId_h", (float)hit->getPDG(), weight);
    }
}

void MCAnaHistos::FillMCEcalHits(std::vector<MCEcalHit*> *mcEcalHits, float weight ) {
    int nHits = mcEcalHits->size();
    Fill1DHisto("numMCEcalHit_h", nHits, weight);
    for (int i=0; i < nHits; i++)
    {
        MCEcalHit *hit = mcEcalHits->at(i);
        Fill1DHisto("mcEcalHitEnergy_h", hit->getEnergy()*1000.0, weight); // Scaled to MeV
    }
}

void MCAnaHistos::FillMCParticleHistos(MCParticle* mcpart, std::string label, double weight){
  double px=mcpart->getMomentum().at(0);
  double py=mcpart->getMomentum().at(1);
  double pz=mcpart->getMomentum().at(2);

  double p=sqrt(px*px+py*py+pz*pz);

  Fill1DHisto(label+"_px_h",px,weight);
  Fill1DHisto(label+"_py_h",py,weight);
  Fill1DHisto(label+"_pz_h",pz,weight);
  Fill1DHisto(label+"_p_h",p,weight);

  double thetaX=px/p;  
  double thetaY=py/p;
  Fill1DHisto(label+"_thetax_h",thetaX,weight);
  Fill1DHisto(label+"_thetay_h",thetaY,weight);

  return;
}

void MCAnaHistos::FillMCPairHistos(MCParticle* ele,MCParticle* pos,std::string label, double weight){
  double px1=ele->getMomentum().at(0);
  double py1=ele->getMomentum().at(1);
  double pz1=ele->getMomentum().at(2);
  double px2=pos->getMomentum().at(0);
  double py2=pos->getMomentum().at(1);
  double pz2=pos->getMomentum().at(2);

  double E1=sqrt(px1*px1+py1*py1+pz1*pz1);
  double E2=sqrt(px2*px2+py2*py2+pz2*pz2);

  double pxTot=px1+px2;
  double pyTot=py1+py2;
  double pzTot=pz1+pz2;
  double pTotSq=pxTot*pxTot+pyTot*pyTot+pzTot*pzTot;
  double Etot=E1+E2;

  double mass=sqrt(Etot*Etot-pTotSq);

  Fill1DHisto(label+"_pxV0_h",pxTot,weight);
  Fill1DHisto(label+"_pyV0_h",pyTot,weight);
  Fill1DHisto(label+"_pzV0_h",pzTot,weight);
  Fill1DHisto(label+"_pV0_h",sqrt(pTotSq),weight);
  Fill1DHisto(label+"_mass_h",mass,weight);
  
  double thetaX=pxTot/sqrt(pTotSq);  
  double thetaY=pyTot/sqrt(pTotSq);
  Fill1DHisto(label+"_thetaxV0_h",thetaX,weight);
  Fill1DHisto(label+"_thetayV0_h",thetaY,weight);

}
