
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
void plot_ttbar_overlay_safe()                                            //função void é a função principal que será executada quando rodar o código,os arquivos
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
    const double sigmaLO  = 18847; // pb                                // Seção de choque total para o processo pp -> ttbar no nível de Leading Order (LO) em pb         
    const double sigmaNLO = 28112; // pb                                // Seção de choque total para o processo pp -> ttbar no nível de Next-to-Leading Order (NLO) em pb  
    const double sigmaBH4 = 6.8e2; // pb                                // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV em pb
    const double sigmaBH6 = 4.56e2; // pb                               // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV em pb
    const double sigmaBHX = 3.16e2; // pb                               // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV em pb

    // arquivos
    TFile* fLO  = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/cut_pt1000/analysis_pqcd_LO_100tev_5M_cut_pt1000.root");     //entrada de dados para o processo pp -> ttbar no nível de Leading Order (LO)
    TFile* fNLO = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/cut_pt1000/analysis_pqcd_NLO_100tev_5M_cut_pt1000.root");      //entrada de dados para o processo pp -> ttbar no nível de Next-to-Leading Order (NLO)
    TFile* fBH4 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/cut_pt1000/analysis_100tev_n6_md4_mbh8_cut_pt1000.root");    //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV
    TFile* fBH6 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/cut_pt1000/analysis_100tev_n6_md4_mbh9_cut_pt1000.root");    //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV
    TFile* fBHX = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/cut_pt1000/analysis_100tev_n6_md4_mbh10_cut_pt1000.root");   //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV

    if (!fLO || !fNLO || !fBH4 || !fBH6 || !fBHX) {
        cout << "ERRO abrindo arquivos .root" << endl;                 //Verifica se os arquivos .root foram abertos corretamente, se algum deles for um ponteiro nulo, imprime uma mensagem de erro e retorna sem fazer nada
        return;
    }

    // histos
    TH1F* hPt_LO  = (TH1F*)fLO ->Get("h_top_pt");                     
    TH1F* hPt_NLO  = (TH1F*)fNLO->Get("h_top_pt");                    
    TH1F* hPt_BH4  = (TH1F*)fBH4->Get("h_top_pt");                   
    TH1F* hPt_BH6  = (TH1F*)fBH6->Get("h_top_pt");                  
    TH1F* hPt_BHX  = (TH1F*)fBHX->Get("h_top_pt");                 


    TH1F* hEta_LO  = (TH1F*)fLO ->Get("h_top_eta");                
    TH1F* hEta_NLO = (TH1F*)fNLO->Get("h_top_eta");                  
    TH1F* hEta_BH4 = (TH1F*)fBH4->Get("h_top_eta");                  
    TH1F* hEta_BH6 = (TH1F*)fBH6->Get("h_top_eta");                 
    TH1F* hEta_BHX = (TH1F*)fBHX->Get("h_top_eta");                   


    

    if (!hEta_LO || !hPt_LO || !hEta_NLO || !hPt_NLO ||
        !hEta_BH4 || !hPt_BH4 || !hEta_BH6 || !hPt_BH6 ||             //Verifica se os histogramas de eta e pt foram obtidos corretamente dos arquivos .root, se algum deles for um ponteiro nulo, imprime uma mensagem de erro e retorna sem fazer nada
        !hEta_BHX || !hPt_BHX) {
        cout << "ERRO: algum histograma 'eta' ou 'pt' nao foi encontrado." << endl;
        return;
    }
    // =================================================================
    // PAINEL 1 — η  (dσ/dη)
    // =================================================================
    // --- CANVAS 1: Pseudorrapidez (Eta) ---
    TCanvas *c1 = new TCanvas("c1", "Distribuicao de Eta", 800, 600);
    c1->cd();                                                        //Desenha o canvas c1, onde o gráfico de pt será desenhado                                         
    gPad->SetLeftMargin(0.16);                                      //Margem esquerda do pad para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y 
    gPad->SetLogy();  

    if (hEta_LO) {
    hEta_LO->SetTitle("Pseudorrapidez do Top; #eta; d#sigma/d#eta [pb]");
    hEta_LO->Draw("HIST");
    hEta_LO->SetLineColor(kBlack);                                   
    hEta_LO->SetLineWidth(4);     
    }
    if (hEta_NLO) {
    hEta_NLO->Draw("HIST SAME");
    hEta_NLO->SetLineColor(kGray+1);                                 
    hEta_NLO->SetLineStyle(2);                                   
    hEta_NLO->SetLineWidth(4);  
    }

    if (hEta_BH4) { 
    hEta_BH4->Draw("HIST SAME");
    hEta_BH4->SetLineColor(kRed);                                    
    hEta_BH4->SetLineStyle(3);                                      
    hEta_BH4->SetLineWidth(4);  
    }

    if (hEta_BH6) {
    hEta_BH6->Draw("HIST SAME");
    hEta_BH6->SetLineColor(kGreen+2);                               
    hEta_BH6->SetLineStyle(4);                                       
    hEta_BH6->SetLineWidth(4); 
    }

    if (hEta_BHX) {
    hEta_BHX->Draw("HIST SAME");
    hEta_BHX->SetLineColor(kBlue);                                   
    hEta_BHX->SetLineStyle(5);                                      
    hEta_BHX->SetLineWidth(4);
    }

    hEta_LO->SetMaximum(hEta_NLO->GetMaximum() * 300);
    hEta_LO->SetTitle(";#eta;d#sigma/d#eta [pb]");                   
    hEta_LO->GetXaxis()->CenterTitle(true);                         
    hEta_LO->GetYaxis()->CenterTitle(true);   
    
    c1->Modified();
    c1->Update();

    TLegend* legPt1 = new TLegend(0.53, 0.68, 0.88, 0.88); 

    legPt1->SetBorderSize(0);                                       
    legPt1->SetFillStyle(0);                                        
    legPt1->SetTextFont(42);                                        
    legPt1->SetTextSize(0.03);            

    legPt1->AddEntry(hEta_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l"); 
    legPt1->AddEntry(hEta_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l"); 
    legPt1->AddEntry(hEta_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l"); 
    legPt1->AddEntry(hEta_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l"); 
    legPt1->AddEntry(hEta_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l");
    legPt1->Draw();


    c1->Update();
    c1->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh_fig1.eps");  
    c1->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh_fig1.png");  

   // =================================================================
   // PAINEL 2 — pT  (dσ/dpT)
   // ================================================================= 
    TCanvas* c2 = new TCanvas("c2", "Figura 2", 1200, 800);         //Tamanho do canvas para o gráfico de pt

    c2->cd();                                                       //Desenha o canvas c2, onde o gráfico de pt será desenhado 
    gPad->SetLeftMargin(0.16);                                      //Margem esquerda do pad para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y 
    gPad->SetLogy();                                                //Define o eixo Y do canvas c2 como logarítmico

    if (hPt_LO) {                                                  //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe
        hPt_LO->SetLineColor(kBlack);                              //Define a cor da linha do histograma 
        hPt_LO->SetLineWidth(4);                                   //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível

        hPt_LO->GetXaxis()->SetRangeUser(0, 1200);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 1200 GeV
        hPt_LO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        hPt_LO->GetXaxis()->CenterTitle(true);                     
        hPt_LO->GetYaxis()->CenterTitle(true);                  
        hPt_LO->GetYaxis()->SetTitleOffset(1.7);                   //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO)

        hPt_LO->DrawCopy("HIST");                                      //Desenha o histograma de pt para o processo pp -> ttbar (LO)
    }

    if (hPt_NLO) {                                                 //Verifica se o histograma de pt para o processo pp -> ttbar (NLO) existe
        hPt_NLO->SetLineColor(kGray+1);                            
        hPt_NLO->SetLineWidth(4);                                  
        hPt_NLO->SetLineStyle(2);                                  
        hPt_NLO->Draw("HIST SAME");                           
    }

    if (hPt_BH4) { hPt_BH4->SetLineColor(kRed);     hPt_BH4->SetLineStyle(3); hPt_BH4->SetLineWidth(4); hPt_BH4->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV 
    if (hPt_BH6) { hPt_BH6->SetLineColor(kGreen+2); hPt_BH6->SetLineStyle(4); hPt_BH6->SetLineWidth(4); hPt_BH6->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp...
    if (hPt_BHX) { hPt_BHX->SetLineColor(kBlue);    hPt_BHX->SetLineStyle(5); hPt_BHX->SetLineWidth(4); hPt_BHX->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp...


    c2->Modified();                                                 
    c2->Update();                                                       //Atualiza o canvas c2 para garantir que todas as modificações feitas, como o desenho dos histogramas, a configuração dos pads e a adição da legenda, sejam refletidas corretamente no canvas c2 antes de salvar a figura

    TLegend* legPt2 = new TLegend(0.53, 0.68, 0.88,  0.88); 

    legPt2->SetBorderSize(0);                                       
    legPt2->SetFillStyle(0);                                        
    legPt2->SetTextFont(42);                                        
    legPt2->SetTextSize(0.03);            

    legPt2->AddEntry(hPt_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l"); 
    legPt2->AddEntry(hPt_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l"); 
    legPt2->AddEntry(hPt_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l"); 
    legPt2->AddEntry(hPt_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l"); 
    legPt2->AddEntry(hPt_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l");

    legPt2->Draw();     
    c2->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh_fig2.eps");            //Salva o canvas "c2" como um arquivo de imagem no formato EPS com o nome "ttbar_overlay_100TeV_5M_n6_md4_mbh.eps", para que a figura possa ser visualizada e compartilhada em formato de documento EPS, mantendo a qualidade vetorial da figura
    c2->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh_fig2.png");            //Salva o canvas "c2" como um arquivo de imagem no formato PNG com o nome "ttbar_overlay_100TeV_5M_n6_md4_mbh.png", para que a figura possa ser visualizada e compartilhada em formato de imagem

    // ===============================
    // PAINEL 3 - pT (dσ/dpT) - Gráfico separado só para pT até 3500 GeV
    // ===============================
    TCanvas* c3 = new TCanvas("c3", "Figura 3", 1200, 800);         //Tamanho do canvas para o gráfico de pt

    c3->cd();                                                       //Desenha o canvas c3, onde o gráfico de pt será desenhado 
    gPad->SetLeftMargin(0.16);                                      //Margem esquerda do pad para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y 
    gPad->SetLogy();                                                //Define o eixo Y do canvas c3 como logarítmico

    if (hPt_LO) {                                                   //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe
        hPt_LO->SetLineColor(kBlack);                               //Define a cor da linha do histograma 
        hPt_LO->SetLineWidth(4);                                    //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível

        hPt_LO->GetXaxis()->SetRangeUser(0, 3500);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 3500 GeV
        hPt_LO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO)
        hPt_LO->GetXaxis()->CenterTitle(true);                     
        hPt_LO->GetYaxis()->CenterTitle(true);                  
        hPt_LO->GetYaxis()->SetTitleOffset(1.7);                   //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO)

        hPt_LO->DrawCopy("HIST");                                  //Desenha o histograma de pt para o processo pp -> ttbar (LO) usando DrawCopy para garantir que as configurações de estilo sejam aplicadas corretamente
    }

    if (hPt_NLO) {                                                 //Verifica se o histograma de pt para o processo pp -> ttbar (NLO) existe
        hPt_NLO->SetLineColor(kGray+1);                            
        hPt_NLO->SetLineWidth(4);                                  
        hPt_NLO->SetLineStyle(2);                                  
        hPt_NLO->Draw("HIST SAME");                           
    }

    if (hPt_BH4) { hPt_BH4->SetLineColor(kRed);     hPt_BH4->SetLineStyle(3); hPt_BH4->SetLineWidth(4); hPt_BH4->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV 
    if (hPt_BH6) { hPt_BH6->SetLineColor(kGreen+2); hPt_BH6->SetLineStyle(4); hPt_BH6->SetLineWidth(4); hPt_BH6->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp...
    if (hPt_BHX) { hPt_BHX->SetLineColor(kBlue);    hPt_BHX->SetLineStyle(5); hPt_BHX->SetLineWidth(4); hPt_BHX->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp...


    c3->Modified();                                                 
    c3->Update();                                                 

    TLegend* legPt3 = new TLegend(0.53, 0.68, 0.88,  0.88); 

    legPt3->SetBorderSize(0);                                       
    legPt3->SetFillStyle(0);                                        
    legPt3->SetTextFont(42);                                        
    legPt3->SetTextSize(0.03);            

    legPt3->AddEntry(hPt_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l"); 
    legPt3->AddEntry(hPt_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l"); 
    legPt3->AddEntry(hPt_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l"); 
    legPt3->AddEntry(hPt_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l"); 
    legPt3->AddEntry(hPt_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l");

    legPt3->Draw();
    c3->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh_fig3.eps");            
    c3->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh_fig3.png");  

    cout << string(70, '=') << "\n" << endl;
    cout << "Figuras salvas em ttbar_overlay_100TeV_5M_n6_md4_mbh.[eps,png]" << endl;   //Imprime no console a mensagem indicando que as figuras foram salvas com os nomes "ttbar_overlay_100TeV_5M_n6_md4_mbh.eps" e "ttbar_overlay_100TeV_5M_n6_md4_mbh.png", para informar ao usuário que as figuras foram geradas e estão disponíveis para visualização e compartilhamento
}