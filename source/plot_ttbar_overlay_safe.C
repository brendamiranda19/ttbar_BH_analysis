
#include "TFile.h"                    //Classe TFile da framework do ROOT, serve para manipular arquivos .root (abrir, criar, ler, escrever)
#include "TLatex.h"                   //Classe TLatex da framework do ROOT, serve para escrever texto no canvas usando sintaxe LaTeX (fórmulas matemáticas, símbolos, etc.)
#include "TH1F.h"                     //Classe TH1F da framework do ROOT, serve para criar e manipular histogramas unidimensionais de tipo float
#include "TCanvas.h"                  //Classe TCanvas da framework do ROOT, serve para criar e manipular o canvas onde os histogramas serão desenhados
#include "TLegend.h"                  //Classe TLegend da framework do ROOT, serve para criar e manipular legendas para os histogramas
#include "TStyle.h"                   //Classe TStyle da framework do ROOT, serve para definir o estilo visual dos gráficos (cores, fontes, margens, etc.)
#include "TROOT.h"                    //Classe TROOT da framework do ROOT, serve para inicializar e configurar o ambiente do ROOT (definir estilos, carregar bibliotecas, etc.)
#include "TMath.h"                    //Classe TMath da framework do ROOT, serve para funções matemáticas comuns (funções trigonométricas, exponenciais, etc.)
#include "TFrame.h"                   //Classe TFrame da framework do ROOT, manipula os frames dos gráficos (eixos, títulos, etc.)
#include <vector>                     //Classe vector da biblioteca standard do C++, serve para criar e manipular vetores
#include <fstream>                    //Classe fstream da biblioteca standard do C++, serve para manipular arquivos de texto (abrir, criar, ler, escrever)
#include <iostream>                   //Classe iostream da biblioteca standard do C++, serve para manipulação de entrada e saída de dados (cout, cin, etc.) 
#include <algorithm>                  //Classe algorithm da biblioteca standard do C++, serve para funções de algoritmos de estrutura de dados (sort, find, copy, etc.)

using namespace std;                  //Evita ter que escrever "std::" antes de cada classe ou função da biblioteca standard do C++

// =====================================================================
// MACRO PRINCIPAL
// =====================================================================
void plot_ttbar_overlay_safe()                                            
{
    // ================================================================
    // ATLAS-LIKE CLEAN STYLE
    // ================================================================
    gROOT->SetStyle("Plain");                                           
    gStyle->SetOptStat(0);                                              
    gStyle->SetOptTitle(0);                                             //Função do ROOT para desabilitar o título dos histogramas

    gStyle->SetCanvasColor(0);                                          //Define a cor do canvas (fundo do gráfico) como 0, que é branco        
    gStyle->SetPadColor(0);                                             //Define a cor da área de plotagem (pad) como 0, que é branco
    gStyle->SetFrameFillColor(0);                                       //Define a cor de preenchimento do frame (a borda do gráfico) como 0, que é branco

    gStyle->SetPadBorderMode(0);                                        //Define como a borda do pad (área de plotagem) deve ser desenhada, 0 significa sem borda
    gStyle->SetFrameBorderMode(0);                                      //Define como a borda do frame (borda do gráfico) deve ser desenhada, 0 significa sem borda

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
    // >>> Coloque as seções de choque (pb) que você quer usar <<<
    //const double sigmattbar_tW = 3.0466e4; // pb tW 2354+ ttbar NLO 28112
    //const double sigmaLO  = 18847; // pb                                // Seção de choque total para o processo pp -> ttbar no nível de Leading Order (LO) em pb         
    //const double sigmaNLO = 28112; // pb                                // Seção de choque total para o processo pp -> ttbar no nível de Next-to-Leading Order (NLO) em pb  
    //const double sigmaBH4 = 6.8e2; // pb                                // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV em pb
    //const double sigmaBH6 = 4.56e2; // pb                               // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV em pb
    //const double sigmaBHX = 3.16e2; // pb                               // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV em pb

    // arquivos
    TFile* fTotal = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_ttbar_tw_inclusive.root");
    //TFile* ftW_antitop = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_single_top_tWchannel_antitop.root");
    //TFile* ftW_top = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_single_top_tWchannel_top.root");
    //TFile* fLO  = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_pqcd_LO_100tev_5M.root");     //entrada de dados para o processo pp -> ttbar no nível de Leading Order (LO)
    //TFile* fNLO = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_pqcd_NLO_100tev_5M.root");      //entrada de dados para o processo pp -> ttbar no nível de Next-to-Leading Order (NLO)
    TFile* fBH4 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh8.root");    //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV
    TFile* fBH6 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh9.root");    //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV
    TFile* fBHX = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh10.root");   //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV

    if (!fTotal || !fBH4 || !fBH6 || !fBHX) {                    // Add fLO se usar o arquivo de LO, o mesmo p !ftW_antitop || !ftW_top || !fNLO 
        cout << "ERRO abrindo arquivos .root" << endl;                 //Verifica se os arquivos .root foram abertos corretamente, se algum deles for um ponteiro nulo, imprime uma mensagem de erro e retorna sem fazer nada
        return;
    }

    // histos
    TH1F* hPt_total = (TH1F*)fTotal->Get("h_top_pt");
    //TH1F* hPt_tW_antitop  = (TH1F*)ftW_antitop->Get("h_top_pt"); 
    //TH1F* hPt_tW_top  = (TH1F*)ftW_top->Get("h_top_pt"); 
    //TH1F* hPt_LO  = (TH1F*)fLO ->Get("h_top_pt");                     
    //TH1F* hPt_NLO  = (TH1F*)fNLO->Get("h_top_pt");                    
    TH1F* hPt_BH4  = (TH1F*)fBH4->Get("h_top_pt");                   
    TH1F* hPt_BH6  = (TH1F*)fBH6->Get("h_top_pt");                  
    TH1F* hPt_BHX  = (TH1F*)fBHX->Get("h_top_pt");                 

    TH1F* hEta_total = (TH1F*)fTotal->Get("h_top_eta");
    //TH1F* hEta_tW_antitop  = (TH1F*)ftW_antitop->Get("h_top_eta");   
    //TH1F* hEta_tW_top  = (TH1F*)ftW_top->Get("h_top_eta");
    //TH1F* hEta_LO  = (TH1F*)fLO ->Get("h_top_eta");             
    //TH1F* hEta_NLO = (TH1F*)fNLO->Get("h_top_eta");                  
    TH1F* hEta_BH4 = (TH1F*)fBH4->Get("h_top_eta");                  
    TH1F* hEta_BH6 = (TH1F*)fBH6->Get("h_top_eta");                 
    TH1F* hEta_BHX = (TH1F*)fBHX->Get("h_top_eta");      

    if ( !hEta_total || !hPt_total ||                                                                  //!hEta_tW_antitop || !hPt_tW_antitop || !hEta_tW_top || !hPt_tW_top || !hEta_NLO || !hPt_NLO ||
        !hEta_BH4 || !hPt_BH4 || !hEta_BH6 || !hPt_BH6 ||             //Verifica se os histogramas de eta e pt foram obtidos corretamente dos arquivos .root, se algum deles for um ponteiro nulo, imprime uma mensagem de erro e retorna sem fazer nada
        !hEta_BHX || !hPt_BHX) {
        cout << "ERRO: algum histograma 'eta' ou 'pt' nao foi encontrado." << endl;
        return;
    }

    // aumenta a largura dos bins de Pt
    hPt_total->Rebin(2);
    //hPt_tW_antitop->Rebin(2);
    //hPt_tW_top->Rebin(2);
    //hPt_LO->Rebin(2);
    //hPt_NLO->Rebin(2);
    hPt_BH4->Rebin(2);
    hPt_BH6->Rebin(2);
    hPt_BHX->Rebin(2);

//Somar as curvas de tbarW e tW
//TH1F* hPt_tW = (TH1F*)hPt_tW_top->Clone("hPt_tW");
//hPt_tW->Add(hPt_tW_antitop);

//TH1F* hEta_tW = (TH1F*)hEta_tW_top->Clone("hEta_tW");
//hEta_tW->Add(hEta_tW_antitop);

    // =================================================================
    // PAINEL 1 — η  (dσ/dη)
    // =================================================================
    // --- CANVAS 1: Pseudorrapidez (Eta) ---
    TCanvas *c1 = new TCanvas("c1", "Distribuicao de Eta", 800, 600);
    c1->cd();                                                        //Desenha o canvas c1, onde o gráfico de pt será desenhado                                         
    gPad->SetLeftMargin(0.16);                                      //Margem esquerda do pad para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y 
    gPad->SetLogy();  

    if (hEta_total) {
        hEta_total->SetLineColor(kBlack);    
        hEta_total->SetTitle("Pseudorrapidez do Top; #eta; d#sigma/d#eta [pb]");                             
        hEta_total->SetLineStyle(2);                                   
        hEta_total->SetLineWidth(4);  
        hEta_total->Draw("HIST");
    }

    //if (hEta_LO) {
    //hEta_LO->Draw("HIST");
    //hEta_LO->SetLineColor(kBlack);                                   
    //hEta_LO->SetLineWidth(4);     
    //}

    //if (hEta_NLO) {
        //hEta_NLO->SetLineColor(kBlack);    
        //hEta_NLO->SetTitle("Pseudorrapidez do Top; #eta; d#sigma/d#eta [pb]");                             
        //hEta_NLO->SetLineStyle(2);                                   
        //hEta_NLO->SetLineWidth(4);  
        //hEta_NLO->Draw("HIST");
    //}

    //if (hEta_tW_antitop) {
        //hEta_tW_antitop->SetLineColor(kGray+1);                                   
        //hEta_tW_antitop->SetLineWidth(4);    
        //hEta_tW_antitop->Draw("HIST SAME"); 
    //}

    //if (hEta_tW_top) {
        //hEta_tW_top->SetLineColor(kOrange);                                   
        //hEta_tW_top->SetLineWidth(4);
        //hEta_tW_top->SetLineStyle(4);
        //hEta_tW_top->Draw("HIST SAME");     
    //}

    //if (hEta_tW) {
        //hEta_tW->SetLineColor(kGray+1);                                   
        //hEta_tW->SetLineWidth(4);
        //hEta_tW->SetLineStyle(4);
        //hEta_tW->Draw("HIST SAME");     
    //}

    if (hEta_BH4) { 
        hEta_BH4->SetLineColor(kRed);                                    
        hEta_BH4->SetLineStyle(3);                                      
        hEta_BH4->SetLineWidth(4);
        hEta_BH4->Draw("HIST SAME");  
    }

    if (hEta_BH6) {
        hEta_BH6->SetLineColor(kGreen+2);                               
        hEta_BH6->SetLineStyle(4);                                       
        hEta_BH6->SetLineWidth(4); 
        hEta_BH6->Draw("HIST SAME");
    }

    if (hEta_BHX) {
        hEta_BHX->SetLineColor(kBlue);                                   
        hEta_BHX->SetLineStyle(5);                                      
        hEta_BHX->SetLineWidth(4);
        hEta_BHX->Draw("HIST SAME");
    }

    hEta_total->SetMaximum(hEta_total->GetMaximum() * 500);
    hEta_total->SetTitle(";#eta;d#sigma/d#eta [pb]");                   
    hEta_total->GetXaxis()->CenterTitle(true);                         
    hEta_total->GetYaxis()->CenterTitle(true);   
    
    c1->Modified();
    c1->Update();

    TLegend* legPt1 = new TLegend(0.53, 0.68, 0.88, 0.88); 

    legPt1->SetBorderSize(0);                                       
    legPt1->SetFillStyle(0);                                        
    legPt1->SetTextFont(42);                                        
    legPt1->SetTextSize(0.03);            

    legPt1->AddEntry(hEta_total, "pp #rightarrow t#bar{t}+tW", "l");
    //legPt1->AddEntry(hEta_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l"); 
    //legPt1->AddEntry(hEta_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l");
    //legPt1->AddEntry(hEta_tW, "pp #rightarrow tW  (pQCD)", "l"); 
    //legPt1->AddEntry(hEta_tW_antitop, "pp #rightarrow #bar{t}W^{+} (pQCD)", "l"); 
    //legPt1->AddEntry(hEta_tW_top, "pp #rightarrow tW  (pQCD)", "l"); 
    legPt1->AddEntry(hEta_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l"); 
    legPt1->AddEntry(hEta_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l"); 
    legPt1->AddEntry(hEta_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l");
    legPt1->Draw();


    c1->Update();
    c1->SaveAs("ttbar_tW_inclusive_100TeV_5M_n6_md4_mbh_fig1.eps");  
    c1->SaveAs("ttbar_tW_inclusive_100TeV_5M_n6_md4_mbh_fig1.png");  

   // =================================================================
   // PAINEL 2 — pT  (dσ/dpT)
   // ================================================================= 
    TCanvas* c2 = new TCanvas("c2", "Figura 2", 1200, 800);         //Tamanho do canvas para o gráfico de pt

    c2->cd();                                                       //Desenha o canvas c2, onde o gráfico de pt será desenhado 
    gPad->SetLeftMargin(0.16);                                      //Margem esquerda do pad para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y 
    gPad->SetLogy();                                                //Define o eixo Y do canvas c2 como logarítmico

    //if (hPt_LO) {                                                  //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe
        //hPt_LO->SetLineColor(kBlack);                              //Define a cor da linha do histograma 
        //hPt_LO->SetLineWidth(4);                                   //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível
        //hPt_LO->GetXaxis()->SetRangeUser(0, 1200);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 1200 GeV
        //hPt_LO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        //hPt_LO->GetXaxis()->CenterTitle(true);                     
        //hPt_LO->GetYaxis()->CenterTitle(true);                  
        //hPt_LO->GetYaxis()->SetTitleOffset(1.7);                   //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO)
        //hPt_LO->DrawCopy("HIST");                                      //Desenha o histograma de pt para o processo pp -> ttbar (LO)
    //}

    if (hPt_total) {                                                 //Verifica se o histograma de pt para o processo pp -> ttbar (NLO) existe
        hPt_total->SetLineColor(kBlack);   
        hPt_total->GetXaxis()->SetRangeUser(0, 1200);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 1200 GeV
        hPt_total->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        hPt_total->GetXaxis()->CenterTitle(true);                     
        hPt_total->GetYaxis()->CenterTitle(true);                  
        hPt_total->GetYaxis()->SetTitleOffset(1.7);    
        hPt_total->SetMaximum(hPt_total->GetMaximum() * 300);                         
        hPt_total->SetLineWidth(4);                                                                    
        hPt_total->DrawCopy("HIST");                           
    }

    //if (hPt_NLO) {                                                 //Verifica se o histograma de pt para o processo pp -> ttbar (NLO) existe
        //hPt_NLO->SetLineColor(kBlack);   
        //hPt_NLO->GetXaxis()->SetRangeUser(0, 1200);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 1200 GeV
        //hPt_NLO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        //hPt_NLO->GetXaxis()->CenterTitle(true);                     
        //hPt_NLO->GetYaxis()->CenterTitle(true);                  
        //hPt_NLO->GetYaxis()->SetTitleOffset(1.7);    
        //hPt_NLO->SetMaximum(hPt_NLO->GetMaximum() * 300);                         
        //hPt_NLO->SetLineWidth(4);                                                                    
        //hPt_NLO->DrawCopy("HIST");                           
    //}

    //if (hPt_tW) {                                                  //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe
        //hPt_tW->SetLineColor(kGray+1);    
        //hPt_tW->SetLineStyle(2);                          //Define a cor da linha do histograma 
        //hPt_tW->SetLineWidth(4);                                   //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível
        //hPt_tW->DrawCopy("HIST SAME");                                      //Desenha o histograma de pt para o processo pp -> ttbar (LO)
    //}

    //if (hPt_tW_antitop) {                                                  //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe
        //hPt_tW_antitop->SetLineColor(kGray+1);    
        //hPt_tW_antitop->SetLineStyle(2);                          //Define a cor da linha do histograma 
        //hPt_tW_antitop->SetLineWidth(4);                                   //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível
        //hPt_tW_antitop->DrawCopy("HIST SAME");                                      //Desenha o histograma de pt para o processo pp -> ttbar (LO)
    //}

    //if (hPt_tW_top) {
        //hPt_tW_top->SetLineColor(kOrange);  
        //hPt_tW_top->SetLineStyle(4);                                 
        //hPt_tW_top->SetLineWidth(4); 
        //hPt_tW_top->Draw("HIST SAME");    
    //}

    if (hPt_BH4) { 
        hPt_BH4->SetLineColor(kRed);                                    
        hPt_BH4->SetLineStyle(4);                                      
        hPt_BH4->SetLineWidth(4); 
        hPt_BH4->Draw("HIST SAME"); 
    }

    if (hPt_BH6) {
        hPt_BH6->SetLineColor(kGreen+2);                               
        hPt_BH6->SetLineStyle(5);                                       
        hPt_BH6->SetLineWidth(4);
        hPt_BH6->Draw("HIST SAME"); 
    }

    if (hPt_BHX) {
        hPt_BHX->SetLineColor(kBlue);                                   
        hPt_BHX->SetLineStyle(6);                                      
        hPt_BHX->SetLineWidth(4);
        hPt_BHX->Draw("HIST SAME");
    }

    c2->Modified();                                                 
    c2->Update();                                                       //Atualiza o canvas c2 para garantir que todas as modificações feitas, como o desenho dos histogramas, a configuração dos pads e a adição da legenda, sejam refletidas corretamente no canvas c2 antes de salvar a figura

    TLegend* legPt2 = new TLegend(0.53, 0.68, 0.88,  0.88); 

    legPt2->SetBorderSize(0);                                       
    legPt2->SetFillStyle(0);                                        
    legPt2->SetTextFont(42);                                        
    legPt2->SetTextSize(0.03);            

    legPt2->AddEntry(hPt_total,  "pp #rightarrow t#bar{t}+tW",  "l"); 
    //legPt2->AddEntry(hPt_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l"); 
    //legPt2->AddEntry(hPt_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l");
    //legPt2->AddEntry(hPt_tW, "pp #rightarrow tW (pQCD)", "l"); 
    //legPt2->AddEntry(hPt_tW_antitop, "pp #rightarrow #bar{t}W^{+} (pQCD)", "l"); 
    //legPt2->AddEntry(hPt_tW_top, "pp #rightarrow tW^{-} (pQCD)", "l"); 
    legPt2->AddEntry(hPt_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l"); 
    legPt2->AddEntry(hPt_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l"); 
    legPt2->AddEntry(hPt_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l");

    legPt2->Draw();     
    c2->SaveAs("ttbar_tW_inclusive_100TeV_5M_n6_md4_mbh_fig2.eps");            //Salva o canvas "c2" como um arquivo de imagem no formato EPS com o nome "ttbar_overlay_100TeV_5M_n6_md4_mbh.eps", para que a figura possa ser visualizada e compartilhada em formato de documento EPS, mantendo a qualidade vetorial da figura
    c2->SaveAs("ttbar_tW_inclusive_100TeV_5M_n6_md4_mbh_fig2.png");            //Salva o canvas "c2" como um arquivo de imagem no formato PNG com o nome "ttbar_overlay_100TeV_5M_n6_md4_mbh.png", para que a figura possa ser visualizada e compartilhada em formato de imagem

    // ===============================
    // PAINEL 3 - pT (dσ/dpT) - Gráfico separado só para pT até 3500 GeV
    // ===============================
    TCanvas* c3 = new TCanvas("c3", "Figura 3", 1200, 800);         //Tamanho do canvas para o gráfico de pt

    c3->cd();                                                       //Desenha o canvas c3, onde o gráfico de pt será desenhado 
    gPad->SetLeftMargin(0.16);                                      //Margem esquerda do pad para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y 
    gPad->SetLogy();                                                //Define o eixo Y do canvas c3 como logarítmico

    //if (hPt_LO) {                                                   //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe
        //hPt_LO->SetLineColor(kBlack);                               //Define a cor da linha do histograma 
        //hPt_LO->SetLineWidth(4);                                    //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível

        //hPt_LO->GetXaxis()->SetRangeUser(800, 3500);    
        //hPt_LO->SetMaximum(hEta_NLO->GetMaximum() * 200);             //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 3500 GeV
        //hPt_LO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        //hPt_LO->GetXaxis()->CenterTitle(true);                     
        //hPt_LO->GetYaxis()->CenterTitle(true);                  
        //hPt_LO->GetYaxis()->SetTitleOffset(1.7);                   //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO)

        //hPt_LO->DrawCopy("HIST");                                  //Desenha o histograma de pt para o processo pp -> ttbar (LO) usando DrawCopy para garantir que as configurações de estilo sejam aplicadas corretamente
    //}

    if (hPt_total) {                                                 //Verifica se o histograma de pt para o processo pp -> ttbar (NLO) existe
        hPt_total->SetLineColor(kBlack);
        hPt_total->GetXaxis()->SetRangeUser(0, 3500);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 1200 GeV
        hPt_total->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        hPt_total->GetXaxis()->CenterTitle(true);                     
        hPt_total->GetYaxis()->CenterTitle(true);                  
        hPt_total->GetYaxis()->SetTitleOffset(1.7);                   //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO)
        hPt_total->SetMaximum(hPt_total->GetMaximum() * 300);                            
        hPt_total->SetLineWidth(4);                                                                   
        hPt_total->DrawCopy("HIST");                           
    }

    //if (hPt_NLO) {                                                 //Verifica se o histograma de pt para o processo pp -> ttbar (NLO) existe
        //hPt_NLO->SetLineColor(kBlack);
        //hPt_NLO->GetXaxis()->SetRangeUser(0, 3500);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 1200 GeV
        //hPt_NLO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        //hPt_NLO->GetXaxis()->CenterTitle(true);                     
        //hPt_NLO->GetYaxis()->CenterTitle(true);                  
        //hPt_NLO->GetYaxis()->SetTitleOffset(1.7);                   //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO)
        //hPt_NLO->SetMaximum(hPt_NLO->GetMaximum() * 300);                            
        //hPt_NLO->SetLineWidth(4);                                                                   
        //hPt_NLO->DrawCopy("HIST");                           
    //}

    //if (hPt_tW) {
        //hPt_tW->SetLineColor(kGray+1);                                   
        //hPt_tW->SetLineWidth(4);  
        //hPt_tW->SetLineStyle(4);   
        //hPt_tW->Draw("HIST SAME");   
    //}

    //if (hPt_tW_antitop) {                                                  //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe
        //hPt_tW_antitop->SetLineColor(kGray+1);                              //Define a cor da linha do histograma 
        //hPt_tW_antitop->SetLineWidth(4);
        //hPt_tW_antitop->SetLineStyle(2);                                  
        //hPt_tW_antitop->DrawCopy("HIST SAME");                                      //Desenha o histograma de pt para o processo pp -> ttbar (LO)
    //}

    //if (hPt_tW_top) {
        //hPt_tW_top->SetLineColor(kOrange);                                   
        //hPt_tW_top->SetLineWidth(4);  
        //hPt_tW_top->SetLineStyle(4);   
        //hPt_tW_top->Draw("HIST SAME");   
    //}

    if (hPt_BH4) { 
        hPt_BH4->SetLineColor(kRed);                                    
        hPt_BH4->SetLineStyle(4);                                      
        hPt_BH4->SetLineWidth(4);
        hPt_BH4->Draw("HIST SAME");  
    }

    if (hPt_BH6) {
        hPt_BH6->SetLineColor(kGreen+2);                               
        hPt_BH6->SetLineStyle(5);                                       
        hPt_BH6->SetLineWidth(4);
        hPt_BH6->Draw("HIST SAME"); 
    }

    if (hPt_BHX) {
        hPt_BHX->SetLineColor(kBlue);                                   
        hPt_BHX->SetLineStyle(6);                                      
        hPt_BHX->SetLineWidth(4);
        hPt_BHX->Draw("HIST SAME");
    }

    c3->Modified();                                                 
    c3->Update();                                                 

    TLegend* legPt3 = new TLegend(0.53, 0.68, 0.88,  0.88); 

    legPt3->SetBorderSize(0);                                       
    legPt3->SetFillStyle(0);                                        
    legPt3->SetTextFont(42);                                        
    legPt3->SetTextSize(0.03);            

    legPt3->AddEntry(hPt_total,  "pp #rightarrow t#bar{t}+tW",  "l"); 
    //legPt3->AddEntry(hPt_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l"); 
    //legPt3->AddEntry(hPt_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l");
    //legPt3->AddEntry(hPt_tW, "pp #rightarrow tW (pQCD)", "l"); 
    //legPt3->AddEntry(hPt_tW_antitop, "pp #rightarrow #bar{t}W (pQCD)", "l"); 
    //legPt3->AddEntry(hPt_tW_top, "pp #rightarrow tW (pQCD)", "l");
    legPt3->AddEntry(hPt_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l"); 
    legPt3->AddEntry(hPt_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l"); 
    legPt3->AddEntry(hPt_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l");

    legPt3->Draw();
    c3->SaveAs("ttbar_tW_inclusive_100TeV_5M_n6_md4_mbh_fig3.eps");            
    c3->SaveAs("ttbar_tW_inclusive_100TeV_5M_n6_md4_mbh_fig3.png");  

    cout << string(70, '=') << "\n" << endl;
    cout << "Figuras salvas em ttbar_tW_inclusive_100TeV_5M_n6_md4_mbh.[eps,png]" << endl;   //Imprime no console a mensagem indicando que as figuras foram salvas com os nomes "ttbar_tW_100TeV_5M_n6_md4_mbh.eps" e "ttbar_tW_100TeV_5M_n6_md4_mbh.png", para informar ao usuário que as figuras foram geradas e estão disponíveis para visualização e compartilhamento
}