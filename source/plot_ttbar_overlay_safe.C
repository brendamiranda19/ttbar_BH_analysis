#include "TFile.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFrame.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

// =====================================================================
//  Normaliza histograma para uma seção de choque total (pb)
// =====================================================================
void NormalizeToSigma(TH1* h, double sigma_pb)
{
    if (!h) return;
    double I = h->Integral("width");
    if (I <= 0.0) return;
    h->Scale(sigma_pb / I);
}
// =====================================================================
//  Verificação
// =====================================================================
void CheckNormalization(TH1* h, double sigma_pb)
{
    if (!h) return;
    double I = h->Integral("width");
    cout << "Integral (pb): " << I << " | Seção de choque total (pb): " << sigma_pb << endl;
}
// =====================================================================
// MACRO PRINCIPAL
// =====================================================================
void plot_ttbar_overlay()
{
    // ================================================================
    // ATLAS-LIKE CLEAN STYLE
    // ================================================================
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetFrameFillColor(0);

    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");

    gStyle->SetTitleSize(0.048, "XYZ");
    gStyle->SetLabelSize(0.038, "XYZ");

    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.7, "Y");

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetEndErrorSize(0);
    gStyle->SetLineWidth(2);
    
    // ================================================================
    // >>> AJUSTE AQUI: coloque as seções de choque (pb) que você quer usar <<<
    const double sigmaLO  = 18847; // pb
    const double sigmaNLO = 28112; // pb
    const double sigmaBH4 = 6.8e2; // pb
    const double sigmaBH6 = 4.56e2; // pb
    const double sigmaBHX = 3.16e2; // pb 

    // arquivos
    TFile* fLO  = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/pqcd_LO_100tev_5M_output.root");
    TFile* fNLO = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/pqcd_NLO_100tev_5M_output.root");
    TFile* fBH4 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh8.root");
    TFile* fBH6 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh9.root");
    TFile* fBHX = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh10.root");

    if (!fLO || !fNLO || !fBH4 || !fBH6 || !fBHX) {
        cout << "ERRO abrindo arquivos .root" << endl;
        return;
    }

    // histos
    TH1F* hEta_LO  = (TH1F*)fLO ->Get("h_top_eta");
    TH1F* hPt_LO  = (TH1F*)fLO ->Get("h_top_pt");
    TH1F* hEta_NLO = (TH1F*)fNLO->Get("h_top_eta");
    TH1F* hPt_NLO  = (TH1F*)fNLO->Get("h_top_pt");

    TH1F* hEta_BH4 = (TH1F*)fBH4->Get("h_top_eta");
    TH1F* hPt_BH4  = (TH1F*)fBH4->Get("h_top_pt");
    TH1F* hEta_BH6 = (TH1F*)fBH6->Get("h_top_eta");
    TH1F* hPt_BH6  = (TH1F*)fBH6->Get("h_top_pt");

    TH1F* hEta_BHX = (TH1F*)fBHX->Get("h_top_eta");
    TH1F* hPt_BHX  = (TH1F*)fBHX->Get("h_top_pt");

    if (!hEta_LO || !hPt_LO || !hEta_NLO || !hPt_NLO ||
        !hEta_BH4 || !hPt_BH4 || !hEta_BH6 || !hPt_BH6 ||
        !hEta_BHX || !hPt_BHX) {
        cout << "ERRO: algum histograma 'eta' ou 'pt' nao foi encontrado." << endl;
        return;
    }

    // ================================================================
    // CLONES:
    //  - ETA: vamos produzir dσ/dη (pb) -> normaliza
    //  - PT : vamos produzir dσ/dpT (pb/GeV) -> normaliza
    // ================================================================
    TH1F* hEta_LO_p  = (TH1F*)hEta_LO ->Clone("hEta_LO_p");
    TH1F* hEta_NLO_p = (TH1F*)hEta_NLO->Clone("hEta_NLO_p");
    TH1F* hEta_BH4_p = (TH1F*)hEta_BH4->Clone("hEta_BH4_p");
    TH1F* hEta_BH6_p = (TH1F*)hEta_BH6->Clone("hEta_BH6_p");
    TH1F* hEta_BHX_p = (TH1F*)hEta_BHX->Clone("hEta_BHX_p");

    TH1F* hPt_LO_p  = (TH1F*)hPt_LO ->Clone("hPt_LO_p");
    TH1F* hPt_NLO_p = (TH1F*)hPt_NLO->Clone("hPt_NLO_p");
    TH1F* hPt_BH4_p = (TH1F*)hPt_BH4->Clone("hPt_BH4_p");
    TH1F* hPt_BH6_p = (TH1F*)hPt_BH6->Clone("hPt_BH6_p");
    TH1F* hPt_BHX_p = (TH1F*)hPt_BHX->Clone("hPt_BHX_p");

    // ================================================================
    // NOVO: NORMALIZA PARA σ_total (pb) 
    // 1) normaliza usando *Integral(width) -> soma a largura dos bins X os conteúdos, dando a integral total em pb para eta e pt 
    // ================================================================
    if (sigmaLO  > 0) NormalizeToSigma(hEta_LO_p,  sigmaLO);
    if (sigmaNLO > 0) NormalizeToSigma(hEta_NLO_p, sigmaNLO);
    if (sigmaBH4 > 0) NormalizeToSigma(hEta_BH4_p, sigmaBH4);
    if (sigmaBH6 > 0) NormalizeToSigma(hEta_BH6_p, sigmaBH6);
    if (sigmaBHX > 0) NormalizeToSigma(hEta_BHX_p, sigmaBHX);

    if (sigmaLO  > 0) NormalizeToSigma(hPt_LO_p,  sigmaLO);
    if (sigmaNLO > 0) NormalizeToSigma(hPt_NLO_p, sigmaNLO);
    if (sigmaBH4 > 0) NormalizeToSigma(hPt_BH4_p, sigmaBH4);
    if (sigmaBH6 > 0) NormalizeToSigma(hPt_BH6_p, sigmaBH6);
    if (sigmaBHX > 0) NormalizeToSigma(hPt_BHX_p, sigmaBHX);              

    //==============================================================
    // Checar a normalização
    //==============================================================
     cout << "\n===== CHECK NORMALIZATION: ETA =====" << endl;
     cout << "LO   ";  CheckNormalization(hEta_LO_p,  sigmaLO);
     cout << "NLO  ";  CheckNormalization(hEta_NLO_p, sigmaNLO);
     cout << "BH4  ";  CheckNormalization(hEta_BH4_p, sigmaBH4);
     cout << "BH6  ";  CheckNormalization(hEta_BH6_p, sigmaBH6);
     cout << "BHX  ";  CheckNormalization(hEta_BHX_p, sigmaBHX);
     cout << "----------------------------------" << endl;

      cout << "\n===== CHECK NORMALIZATION: PT =====" << endl;
     cout << "LO   ";  CheckNormalization(hPt_LO_p,  sigmaLO);
     cout << "NLO  ";  CheckNormalization(hPt_NLO_p, sigmaNLO);
     cout << "BH4  ";  CheckNormalization(hPt_BH4_p, sigmaBH4);
     cout << "BH6  ";  CheckNormalization(hPt_BH6_p, sigmaBH6);
     cout << "BHX  ";  CheckNormalization(hPt_BHX_p, sigmaBHX);
     cout << "----------------------------------" << endl;  

    // canvas
    TCanvas* c = new TCanvas("c", "ttbar pQCD + BH", 1200, 500);
    // Canvas principal
    TPad* pad1 = new TPad("pad1","",0.00,0.25,0.50,1.00);
    TPad* pad2 = new TPad("pad2","",0.50,0.25,1.00,1.00);
    TPad* padLeg = new TPad("padLeg","",0.00,0.00,1.00,0.25);

    pad1->SetLeftMargin(0.18);
    pad1->SetRightMargin(0.05);
    pad1->SetBottomMargin(0.15);
    pad1->SetTopMargin(0.08);

    pad2->SetLeftMargin(0.18);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.15);
    pad2->SetTopMargin(0.08);

    pad1->Draw();
    pad2->Draw();
    padLeg->Draw();
  
    // =================================================================
    // PAINEL 1 — η  (dσ/dη)
    // =================================================================
    pad1->cd();
    gPad->SetLogy();

    hEta_LO_p->SetLineColor(kBlack);
    hEta_LO_p->SetLineWidth(4);

    hEta_NLO_p->SetLineColor(kGray+1);
    hEta_NLO_p->SetLineStyle(2);
    hEta_NLO_p->SetLineWidth(4);

    hEta_BH4_p->SetLineColor(kRed);
    hEta_BH4_p->SetLineStyle(3);
    hEta_BH4_p->SetLineWidth(4);

    hEta_BH6_p->SetLineColor(kGreen+2);
    hEta_BH6_p->SetLineStyle(4);
    hEta_BH6_p->SetLineWidth(4);

    hEta_BHX_p->SetLineColor(kBlue);
    hEta_BHX_p->SetLineStyle(5);
    hEta_BHX_p->SetLineWidth(4);

    hEta_LO_p->SetTitle(";#eta;d#sigma/d#eta [pb]");
    hEta_LO_p->GetXaxis()->CenterTitle(true);
    hEta_LO_p->GetYaxis()->CenterTitle(true);

    hEta_LO_p->Draw("HIST");
    hEta_NLO_p->Draw("HIST SAME");
    hEta_BH4_p->Draw("HIST SAME");
    hEta_BH6_p->Draw("HIST SAME");
    hEta_BHX_p->Draw("HIST SAME");

  
   // =================================================================
   // PAINEL 2 — pT  (dσ/dpT)
   // =================================================================
    pad2->cd();
    gPad->SetLogy();

    if (hPt_LO_p) {
        hPt_LO_p->SetLineColor(kBlack);
        hPt_LO_p->SetLineWidth(4);

        hPt_LO_p->GetXaxis()->SetRangeUser(0, 1200);
        hPt_LO_p->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");
        hPt_LO_p->GetXaxis()->CenterTitle(true);
        hPt_LO_p->GetYaxis()->CenterTitle(true);
        hPt_LO_p->GetYaxis()->SetTitleOffset(1.7);

        hPt_LO_p->Draw("HIST");
    }

    if (hPt_NLO_p) {
        hPt_NLO_p->SetLineColor(kGray+1);
        hPt_NLO_p->SetLineWidth(4);
        hPt_NLO_p->SetLineStyle(2);
        hPt_NLO_p->Draw("HIST SAME");
    }

    if (hPt_BH4_p) { hPt_BH4_p->SetLineColor(kRed);     hPt_BH4_p->SetLineStyle(3); hPt_BH4_p->SetLineWidth(4); hPt_BH4_p->Draw("HIST SAME"); }
    if (hPt_BH6_p) { hPt_BH6_p->SetLineColor(kGreen+2); hPt_BH6_p->SetLineStyle(4); hPt_BH6_p->SetLineWidth(4); hPt_BH6_p->Draw("HIST SAME"); }
    if (hPt_BHX_p) { hPt_BHX_p->SetLineColor(kBlue);    hPt_BHX_p->SetLineStyle(5); hPt_BHX_p->SetLineWidth(4); hPt_BHX_p->Draw("HIST SAME"); }

    
   pad2->Modified();
   pad2->Update();
  
   // Legenda centralizada no rodapé
    padLeg->cd();
    TLegend* leg = new TLegend(0.10,0.15,0.90,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetNColumns(2);
    leg->SetTextSize(0.15);

    leg->AddEntry(hPt_LO_p,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l");
    leg->AddEntry(hPt_NLO_p, "pp #rightarrow t#bar{t} NLO (pQCD)", "l");
    leg->AddEntry(hPt_BH4_p, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",    "l");
    leg->AddEntry(hPt_BH6_p, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",    "l");
    leg->AddEntry(hPt_BHX_p, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",   "l");

    leg->Draw();

    c->cd();   // volta para o canvas principal
    TLatex header;
    header.SetNDC(true);
    header.SetTextFont(42);
    header.SetTextSize(0.035);
    header.SetTextAlign(13);  // topo-esquerda

    header.DrawLatex(0.02, 0.97, "FCC  #sqrt{s} = 100 TeV");

    // salvar
     c->Update();
     c->SaveAs("ttbar_overlay_2pads_legenda.png");
     c->SaveAs("ttbar_overlay_2pads_legenda.pdf");
 
    // ===============================
    // PAINEL 3 - pT (dσ/dpT) - NOVO: gráfico separado só para pT, usando os histos "hPt_X_p" normalizados
    // ===============================
    TCanvas* c2 = new TCanvas("c2", "Figura 2", 1200, 800);  // tamanho livre

    c2->cd();   // entra no canvas
    gPad->SetLeftMargin(0.16);           
    gPad->SetLogy();       // escala log

    if (hPt_LO_p) {
        hPt_LO_p->SetLineColor(kBlack);
        hPt_LO_p->SetLineWidth(4);

        hPt_LO_p->GetXaxis()->SetRangeUser(0, 3500);
        hPt_LO_p->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");
        hPt_LO_p->GetXaxis()->CenterTitle(true);
        hPt_LO_p->GetYaxis()->CenterTitle(true);
        hPt_LO_p->GetYaxis()->SetTitleOffset(1.7);

        hPt_LO_p->Draw("HIST");
    }

    if (hPt_NLO_p) {
        hPt_NLO_p->SetLineColor(kGray+1);
        hPt_NLO_p->SetLineWidth(4);
        hPt_NLO_p->SetLineStyle(2);
        hPt_NLO_p->Draw("HIST SAME");
    }

    if (hPt_BH4_p) { hPt_BH4_p->SetLineColor(kRed);     hPt_BH4_p->SetLineStyle(3); hPt_BH4_p->SetLineWidth(4); hPt_BH4_p->Draw("HIST SAME"); }
    if (hPt_BH6_p) { hPt_BH6_p->SetLineColor(kGreen+2); hPt_BH6_p->SetLineStyle(4); hPt_BH6_p->SetLineWidth(4); hPt_BH6_p->Draw("HIST SAME"); }
    if (hPt_BHX_p) { hPt_BHX_p->SetLineColor(kBlue);    hPt_BHX_p->SetLineStyle(5); hPt_BHX_p->SetLineWidth(4); hPt_BHX_p->Draw("HIST SAME"); }


    c2->Modified();
    c2->Update();

    TLegend* legPt2 = new TLegend(0.55, 0.65, 0.88, 0.88);

    legPt2->SetBorderSize(0);
    legPt2->SetFillStyle(0);
    legPt2->SetTextFont(42);
    legPt2->SetTextSize(0.03);

    legPt2->AddEntry(hPt_LO_p,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l");
    legPt2->AddEntry(hPt_NLO_p, "pp #rightarrow t#bar{t} NLO (pQCD)", "l");
    legPt2->AddEntry(hPt_BH4_p, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l");
    legPt2->AddEntry(hPt_BH6_p, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l");
    legPt2->AddEntry(hPt_BHX_p, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l");

    legPt2->Draw();

    // salvar
    c2->Update();
    c2->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh.pdf");
    c2->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh.eps");
    c2->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh.png");
   

    cout << "Figuras salvas em ttbar_overlay_100TeV_5M_n6_md4_mbh.[pdf,eps,png]" << endl;
}