#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <stdlib.h>

void CountEvents_Data(){
    
    Int_t totalSize = 0;
    
    TFile *f = TFile::Open("../datasets/ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016G-23Sep2016-v1_1of5.root");
    if (f == 0) {
        // if we cannot open the file, print an error message and return immediatly
        std::cout << "error: ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016G-23Sep2016-v1_1of5.root cannot be extracted." << std::endl;
        return;
    }
    
    gStyle->SetPalette(kOcean);
    
    const static int maxNumJets  = 130;
    
    const static int maxNumElecs = 100;
    
    int   m_nLeps;
    float m_lep_pt    [maxNumElecs];
    float m_lep_eta   [maxNumElecs];
    float m_lep_phi   [maxNumElecs];
    int   m_lep_charge[maxNumElecs];
    float m_lep_iso   [maxNumElecs];
    
    int m_nJets;
    float m_Vtype ;
    float m_Z_mass, m_Z_pt, m_Z_eta, m_Z_phi;
    float m_jet_pt[maxNumJets], m_jet_eta[maxNumJets], m_jet_phi[maxNumJets], m_jet_msv[maxNumJets];;
    
    float m_MET_et;
    float m_MET_phi;
    float m_MET_sumet;
    
    int m_event;
    float m_nPVs;
    unsigned int m_run;
    unsigned int m_lumi;
    
    TTreeReader myReader("EventTree", f);
    TTree *tree = (TTree*)f->Get("tree");
    
    Int_t nEntry = tree->GetEntries();
    
    tree->SetBranchAddress( "nvLeptons"           , &m_nLeps      );
    tree->SetBranchAddress( "vLeptons_pt"         ,  m_lep_pt     );
    tree->SetBranchAddress( "vLeptons_eta"        ,  m_lep_eta    );
    tree->SetBranchAddress( "vLeptons_phi"        ,  m_lep_phi    );
    tree->SetBranchAddress( "vLeptons_charge"     ,  m_lep_charge );
    tree->SetBranchAddress( "vLeptons_pfRelIso04" ,  m_lep_iso    );
    
    tree->SetBranchAddress( "V_mass", &m_Z_mass      );
    tree->SetBranchAddress( "V_pt"  , &m_Z_pt        );
    tree->SetBranchAddress( "V_eta" , &m_Z_eta       );
    tree->SetBranchAddress( "V_phi" , &m_Z_phi       );
    
    tree->SetBranchAddress( "nJet"           , &m_nJets       );
    tree->SetBranchAddress( "Jet_pt"         ,  m_jet_pt      );
    tree->SetBranchAddress( "Jet_eta"        ,  m_jet_eta     );
    tree->SetBranchAddress( "Jet_phi"        ,  m_jet_phi     );
    tree->SetBranchAddress( "Vtype"          , &m_Vtype       );
    tree->SetBranchAddress( "Jet_vtxMass"    ,  m_jet_msv     );

    tree->SetBranchAddress( "met_pt"     , &m_MET_et    );
    tree->SetBranchAddress( "met_phi"    , &m_MET_phi   );
    tree->SetBranchAddress( "met_sumEt"  , &m_MET_sumet );
    
    tree->SetBranchAddress( "evt" , &m_event);
    tree->SetBranchAddress( "nPVs", &m_nPVs );
    tree->SetBranchAddress( "run" , &m_run  );
    tree->SetBranchAddress( "lumi", &m_lumi );
    
    
    
    TCanvas *c1 = new TCanvas("c1","multipad1",1000,600);
    TCanvas *c2 = new TCanvas("c2","multipad2",1000,600);
    TCanvas *c3 = new TCanvas("c3","multipad3",1000,600);
    gStyle->SetOptStat(1);
    c1->Divide(2,2,0.02,0.02);
    c2->Divide(2,1,0.02,0.02);
    c3->Divide(3,2,0.02,0.02);
    
    //c1->SetLogy();
    //c2->SetLogy();
    
    //create histogram
    /*
     THStack *ptBf           = new THStack("ptbf","pt distribution of muons (before cut)");
     THStack *etaBf          = new THStack("etabf","#eta distribution of muons (before cut)");
     THStack *pt             = new THStack("pt","pt distribution of muons");
     THStack *eta            = new THStack("eta","#eta distribution of muons");
     */
    
    TH1F *zbosonBf          = new TH1F("zmbf","Dilepton Mass (w/o Lepton cut)",100,0,300.);
    TH1F *zboson            = new TH1F("zm","Dilepton Mass (w/ Lepton cut)",100,0,300.);
    TH1F *nPV               = new TH1F("nPVs","Primary Vertices (Z+jet w/MET cut)",100,0,100.);
    TH1F *MuonPtBf1         = new TH1F("muonPtbf1","Muon p_T (w/o cut)",100,0,300.);
    TH1F *MuonEtaBf1        = new TH1F("muonEtabf1","Muon #eta (w/o cut)",100,-4,4);
    TH1F *MuonPt1           = new TH1F("muonPt1","Muon p_T (Z+jet w/MET cut)",100,0,300.);
    TH1F *MuonEta1          = new TH1F("muonEta1","Muon #eta (Z+jet w/MET cut)",100,-4,4);
    TH1F *MuonPtBf2         = new TH1F("muonPtbf2","Muon p_T (w/o cut)",100,0,300.);
    TH1F *MuonEtaBf2        = new TH1F("muonEtabf2","Muon #eta (w/o cut)",100,-4,4);
    TH1F *MuonPt2           = new TH1F("muonPt2","Muon p_T (Z+jet w/MET cut)",100,0,300.);
    TH1F *MuonEta2          = new TH1F("muonEta2","Muon #eta (Z+jet w/MET cut)",100,-4,4);
	TH1F *JetMult           = new TH1F("jetMult","Jet Multiplicity (Valid Z-boson)",10,0,10.);
    TH1F *JetPt             = new TH1F("jetPt","Jet p_T (Z+jet w/MET cut)",100,0,200.);
    TH1F *JetEta            = new TH1F("jetEta","Jet #eta (Z+jet w/MET cut)",100,-4,4);
    TH1F *nMET				= new TH1F("nmet","MET (Valid Z-boson)",100,0,200.);
	TH1F *MSV				= new TH1F("msv","M_{SV} (Z+jet w/MET cut)",100,0,20.);
    
	TH1F *DY_Pt             = new TH1F("DYPt","pt distribution of DY process",100,0,300.);
    TH1F *DY_Eta            = new TH1F("DYEta","#eta distribution of DY process",100,-4,4);
    
    //setting color scheme
    MuonPtBf1->SetFillColorAlpha(kRed,0.35);
    MuonPtBf1->SetFillStyle(3002);
    MuonPtBf2->SetFillColorAlpha(kBlue,0.35);
    MuonPtBf2->SetFillStyle(3001);
    
	MuonEtaBf1->SetFillColorAlpha(kRed,0.35);
    MuonEtaBf1->SetFillStyle(3002);
    MuonEtaBf2->SetFillColorAlpha(kBlue,0.35);
    MuonEtaBf2->SetFillStyle(3001);
    
	MuonPt1->SetFillColorAlpha(kRed,0.35);
    MuonPt1->SetFillStyle(3002);
    MuonPt2->SetFillColorAlpha(kBlue,0.35);
    MuonPt2->SetFillStyle(3001);
    
	MuonEta1->SetFillColorAlpha(kRed,0.35);
    MuonEta1->SetFillStyle(3002);
    MuonEta2->SetFillColorAlpha(kBlue,0.35);
    MuonEta2->SetFillStyle(3001);
    
	zbosonBf->SetFillColorAlpha(kRed,0.35);
    zbosonBf->SetFillStyle(3002);
    zboson->SetFillColorAlpha(kBlue,0.35);
    zboson->SetFillStyle(3001);

	JetMult->SetFillColorAlpha(kRed,0.35);
	JetMult->SetFillStyle(3001);
	JetPt->SetFillColorAlpha(kRed,0.35);
	JetPt->SetFillStyle(3001);
	JetEta->SetFillColorAlpha(kRed,0.35);
	JetEta->SetFillStyle(3001);
	nPV->SetFillColorAlpha(kRed,0.35);
	nPV->SetFillStyle(3001);
	nMET->SetFillColorAlpha(kRed,0.35);
	nMET->SetFillStyle(3001);
	MSV->SetFillColorAlpha(kRed,0.35);
	MSV->SetFillStyle(3001);

    
    std::cout << "Total Events: " << nEntry << std::endl;
    
    //loop over events
    for(int iEntry=0; iEntry<nEntry ;++iEntry)
    {
        tree->GetEntry(iEntry);
        
        //if(iEntry%1000 == 0) std::cout << "Event: " << iEntry << std::endl;
        
        
        if(m_Vtype == 0 || m_nLeps == 2)
        {
            bool Zmass_Check = false;
            bool MET_Check = false;
            bool Lep0_Check = false;
            bool Lep1_Check = false;
            double dif = 0;
            int max = 0;
            
            dif = m_lep_pt[0] - m_lep_pt[1];
            
            
            zbosonBf->Fill(m_Z_mass);
            if(dif > 0){
                MuonPtBf1->Fill(m_lep_pt[0]);
                MuonPtBf2->Fill(m_lep_pt[1]);
                MuonEtaBf1->Fill(m_lep_eta[0]);
                MuonEtaBf2->Fill(m_lep_eta[1]);
            }else if(dif < 0){
                MuonPtBf1->Fill(m_lep_pt[1]);
                MuonPtBf2->Fill(m_lep_pt[0]);
                MuonEtaBf1->Fill(m_lep_eta[1]);
                MuonEtaBf2->Fill(m_lep_eta[0]);
            }
            
            if( m_Z_mass>70. && m_Z_mass<110. ) Zmass_Check = true;
            if( m_MET_et<=40. ) MET_Check = true;
            if( m_lep_pt[0] > 20 && abs(m_lep_eta[0]) < 2.4 && m_lep_iso[0] < 0.25) Lep0_Check = true;
            if( m_lep_pt[1] > 20 && abs(m_lep_eta[1]) < 2.4 && m_lep_iso[1] < 0.25) Lep1_Check = true;
            
            if( Zmass_Check==true ){

				zboson->Fill(m_Z_mass);
                nMET->Fill(m_MET_et);
                JetMult->Fill(m_nJets);

                if(Lep0_Check==true && Lep1_Check==true && MET_Check==true){
                    
                    nPV->Fill(m_nPVs);
                    
                    if(dif > 0){
                        MuonPt1->Fill(m_lep_pt[0]);
                        MuonPt2->Fill(m_lep_pt[1]);
                        MuonEta1->Fill(m_lep_eta[0]);
                        MuonEta2->Fill(m_lep_eta[1]);
                    }else if(dif < 0){
                        MuonPt1->Fill(m_lep_pt[1]);
                        MuonPt2->Fill(m_lep_pt[0]);
                        MuonEta1->Fill(m_lep_eta[1]);
                        MuonEta2->Fill(m_lep_eta[0]);
                    }
                    
                    for(int j=0; j<m_nJets; j++){
                        
                        if( m_jet_pt[j] < 30. || abs(m_jet_eta[j]) > 2.4) continue;
                        
                        JetPt->Fill(m_jet_pt[j]);
                        JetEta->Fill(m_jet_eta[j]);
                        MSV->Fill(m_jet_msv[j]);

                    }
                    
                }
            }
        }
        
    }
    
    c1->cd(1);
    gPad->SetLogy();
    MuonPtBf1->Draw();
    MuonPtBf2->Draw("same");
    c1->cd(2);
    gPad->SetLogy();
    MuonPt1->Draw();
    MuonPt2->Draw("same");
    c1->cd(3);
    gPad->SetLogy();
    MuonEtaBf1->Draw();
    MuonEtaBf2->Draw("same");
    c1->cd(4);
    gPad->SetLogy();
    MuonEta1->Draw();
    MuonEta2->Draw("same");
    
    c2->cd(1);
    gPad->SetLogy();
    zbosonBf->Draw("");
    c2->cd(2);
    gPad->SetLogy();
    zboson->Draw("");
    
	c3->cd(1);
	gPad->SetLogy();
	nMET->Draw("");
	c3->cd(2);
	gPad->SetLogy();
	JetMult->Draw("");
	c3->cd(3);
	gPad->SetLogy();
	MSV->Draw("");
	c3->cd(4);
	gPad->SetLogy();
	nPV->Draw("");
	c3->cd(5);
	gPad->SetLogy();
	JetPt->Draw("");
	c3->cd(6);
	gPad->SetLogy();
	JetEta->Draw("");

	c1->Print("Lepton.pdf");
	c2->Print("Z-boson.pdf");
	c3->Print("Jet.pdf");

	gSystem->Exit(0);


    
}
