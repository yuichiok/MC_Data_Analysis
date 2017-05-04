#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <stdlib.h>

void CountEvents_MC(){
    
    Int_t totalSize = 0;
    
    TFile *f = TFile::Open("root://cmseos.fnal.gov//store/group/leptonjets/noreplica/godshalk/ZC2017_03_ZJNtuples2016/MC_SUMMER2016_PRIMARY_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v2_01of18.root");
    if (f == 0) {
        // if we cannot open the file, print an error message and return immediatly
        std::cout << "error: MC_SUMMER2016_PRIMARY_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v2_01of18.root cannot be extracted." << std::endl;
        return;
    }
    
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
    float m_jet_pt[maxNumJets], m_jet_eta[maxNumJets], m_jet_phi[maxNumJets];
    
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
    
    tree->SetBranchAddress( "met_pt"     , &m_MET_et    );
    tree->SetBranchAddress( "met_phi"    , &m_MET_phi   );
    tree->SetBranchAddress( "met_sumEt"  , &m_MET_sumet );
    
    tree->SetBranchAddress( "evt" , &m_event);
    tree->SetBranchAddress( "nPVs", &m_nPVs );
    tree->SetBranchAddress( "run" , &m_run  );
    tree->SetBranchAddress( "lumi", &m_lumi );
    
    
    
    TCanvas *c1 = new TCanvas("c1","multipad1",1000,500);
    TCanvas *c2 = new TCanvas("c2","multipad2",1000,500);
    gStyle->SetOptStat(1);
    c1->Divide(3,2,0.02,0.02);
    c2->Divide(2,1,0.02,0.02);
    
    
    
    //create histogram
    TH1F *zboson            = new TH1F("zm","Z0 invariant mass distribution",100,0,300.);
    TH1F *nPV               = new TH1F("nPVs","Number of Primary Vertices",50,0,50.);
    TH1F *MuonPt            = new TH1F("muonPt","pt distribution of muons",100,0,300.);
    TH1F *MuonEta           = new TH1F("muonEta","#eta distribution of muons",100,-2.5,2.5);
    TH1F *JetPt             = new TH1F("jetPt","pt distribution of jets",100,0,300.);
    TH1F *JetEta            = new TH1F("jetEta","#eta distribution of jets",100,-2.5,2.5);
    
    TH1F *DY_Pt             = new TH1F("DYPt","pt distribution of DY process",100,0,300.);
    TH1F *DY_Eta            = new TH1F("DYEta","#eta distribution of DY process",100,-2.5,2.5);
    
    
    
    std::cout << "Total Events: " << nEntry << std::endl;
    
    //loop over events
    for(int iEntry=0; iEntry<100 ;++iEntry)
    {
        tree->GetEntry(iEntry);
        
        //if(iEntry%1000 == 0) std::cout << "Event: " << iEntry << std::endl;
        
        
        if(m_Vtype == 0)
        {
            for(int i=0; i<m_nLeps; i++)
            {
                if( m_lep_pt[i] < 20 || abs(m_lep_eta[i]) > 2.4 || m_lep_iso[i] > 0.25) continue;
                
                MuonPt->Fill(m_lep_pt[i]);
                MuonEta->Fill(m_lep_eta[i]);
                
            }
            
        }
        if(m_Vtype == 0 || m_Vtype == 1)
        {
            for(int i=0; i<m_nLeps; i++)
            {
                if( m_lep_pt[i] < 20 || abs(m_lep_eta[i]) > 2.4 || m_lep_iso[i] > 0.25) continue;
                
                DY_Pt->Fill(m_lep_pt[i]);
                DY_Eta->Fill(m_lep_eta[i]);
                
            }
            
        }
        
        for(int j=0; j<m_nJets; j++){
            
            if( m_jet_pt[j] < 30. || abs(m_jet_eta[j]) > 2.4) continue;
            
            JetPt->Fill(m_jet_pt[j]);
            JetEta->Fill(m_jet_eta[j]);
            
        }
        
        
        nPV->Fill(m_nPVs);
        
        if( m_Z_mass>70. && m_Z_mass<110. ) zboson->Fill(m_Z_mass);
        
    }
    
    
    c1->cd(1);
    zboson->Draw("");
    c1->cd(2);
    nPV->Draw("");
    c1->cd(3);
    MuonPt->Draw("");
    c1->cd(4);
    MuonEta->Draw("");
    c1->cd(5);
    JetPt->Draw("");
    c1->cd(6);
    JetEta->Draw("");
    
    c2->cd(1);
    DY_Pt->Draw("");
    c2->cd(2);
    DY_Eta->Draw("");
    
    c1->Print("c1.pdf");
    
    gSystem->Exit(0);
    
}
