//COMENTAR TUDO DO CODIGO PARA QUE FIQUE CLARO O QUE CADA PARTE FAZ - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB)
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
//  Normaliza histograma para uma seção de choque total (pb)
// =====================================================================
void NormalizeToSigma(TH1* h, double sigma_pb)                            //Normaliza o histograma "h" para a seção de choque total "sigma_pb" em pb    
{
    if (!h) return;                                                      //Verifica se h é um ponteiro nulo, ou seja, se o histograma existe. Se não existir, retorna sem fazer nada
    double I = h->Integral("width");                                    //Calcula a integral do histograma usando a opção "width", que faz a multiplicação do conteúdo de cada bin (cada bin guarda um número) X largura do bin, dando a integral total em pb para eta e pt
    if (I <= 0.0) return;                                               //Se a integral for menor ou igual a zero, retorna sem fazer nada (para evitar divisão por zero ou normalização errada)
    h->Scale(sigma_pb / I);                                             //Escala do histograma multiplicando cada bin por sigma_pb / I, ou seja, normaliza o histograma para que a integral total seja igual a sigma_pb (a seção de choque total em pb)
}
// =====================================================================
// MACRO PRINCIPAL
// =====================================================================
void plot_ttbar_overlay()                                                //função void é a função principal que será executada quando rodar o código, e é onde todo o código de leitura dos arquivos, manipulação dos histogramas e plotagem dos gráficos será escrito
{
    // ================================================================
    // ATLAS-LIKE CLEAN STYLE
    // ================================================================
    gROOT->SetStyle("Plain");                                           //Definição do estilo "Plain" do ROOT, que é um estilo básico e limpo para os gráficos (sem fundo colorido, sem grades, etc.)
    gStyle->SetOptStat(0);                                              //Função do ROOT para desabilitar a caixa de estatísticas que aparece no canto dos histogramas (que mostra média, desvio padrão, etc.)
    gStyle->SetOptTitle(0);                                             //Função do ROOT para desabilitar o título dos histogramas

    gStyle->SetCanvasColor(0);                                          //Define a cor do canvas (fundo do gráfico) como 0, que é branco        
    gStyle->SetPadColor(0);                                             //Define a cor da área de plotagem (pad) como 0, que é branco
    gStyle->SetFrameFillColor(0);                                       //Define a cor de preenchimento do frame (a borda do gráfico) como 0, que é branco

    gStyle->SetPadBorderMode(0);                                        //Define como a borda do pad (área de plotagem) deve ser desenhada, 0 significa sem borda
    gStyle->SetFrameBorderMode(0);                                      //Define como a borda do frame (borda do gráfico) deve ser desenhada, 0 significa sem borda

    gStyle->SetTextFont(42);                                            //Define a fonte do texto em 42, que é uma fonte sem serifa (Arial) para os títulos e rótulos dos eixos
    gStyle->SetLabelFont(42, "XYZ");                                    //Define a fonte dos eixos X, Y e Z em 42, que é uma fonte sem serifa (Arial) para os rótulos dos eixos
    gStyle->SetTitleFont(42, "XYZ");                                    //Define a fonte dos títulos em 42 para os eixos xyz, que é uma fonte sem serifa (Arial) para os títulos dos eixos

    gStyle->SetTitleSize(0.048, "XYZ");                                 //Define o tamanho dos títulos em 0.048 para os eixos X, Y e Z, que é um tamanho grande para os títulos dos eixos
    gStyle->SetLabelSize(0.038, "XYZ");                                 //Define o tamanho dos rótulos dos eixos em 0.038 para os eixos X, Y e Z, que é um tamanho médio para os rótulos dos eixos    

    gStyle->SetTitleOffset(1.2, "X");                                   //Define o espaço entre o título (referente ao eixo X) e o eixo X em 1.2, que é um valor padrão para deixar o título do eixo X bem posicionado
    gStyle->SetTitleOffset(1.7, "Y");                                   //Define o espaço entre o título (referente ao eixo Y) e o eixo Y em 1.7, que é um valor maior para deixar o título do eixo Y mais afastado do eixo, já que o título do eixo Y é rotacionado e fica mais próximo do eixo

    gStyle->SetPadTickX(1);                                             //Define o estilo dos marcadores dos ticks (as linhas pequenas que indicam os valores nos eixos) no eixo X como 1, que é um estilo de marcador simples                      
    gStyle->SetPadTickY(1);                                             //Define o estilo dos marcadores dos ticks (as linhas pequenas que indicam os valores nos eixos) no eixo Y como 1, que é um estilo de marcador simples

    gStyle->SetEndErrorSize(0);                                         //Define o tamanho das linhas horizontais no final dos erros (as barras de erro) como 0, ou seja, sem linhas horizontais no final das barras de erro
    gStyle->SetLineWidth(2);                                            //Define a largura das curvas, histogramas e linhas em 2, que é uma largura média para deixar as curvas e histogramas bem visíveis
    
    // ================================================================
    // >>> Coloque as seções de choque (pb) que você quer usar <<<
    const double sigmaLO  = 18847; // pb                                // Seção de choque total para o processo pp -> ttbar no nível de Leading Order (LO) em pb         
    const double sigmaNLO = 28112; // pb                                // Seção de choque total para o processo pp -> ttbar no nível de Next-to-Leading Order (NLO) em pb  
    const double sigmaBH4 = 6.8e2; // pb                                // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV em pb
    const double sigmaBH6 = 4.56e2; // pb                               // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV em pb
    const double sigmaBHX = 3.16e2; // pb                               // Seção de choque total para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV em pb

    // arquivos
    TFile* fLO  = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/pqcd_LO_100tev_5M_output.root");       //entrada de dados para o processo pp -> ttbar no nível de Leading Order (LO)
    TFile* fNLO = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/pqcd_NLO_100tev_5M_output.root");      //entrada de dados para o processo pp -> ttbar no nível de Next-to-Leading Order (NLO)
    TFile* fBH4 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh8.root");    //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV
    TFile* fBH6 = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh9.root");    //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV
    TFile* fBHX = TFile::Open("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/analysis_outputs/analysis_100tev_n6_md4_mbh10.root");   //entrada de dados para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV

    if (!fLO || !fNLO || !fBH4 || !fBH6 || !fBHX) {
        cout << "ERRO abrindo arquivos .root" << endl;                 //Verifica se os arquivos .root foram abertos corretamente, se algum deles for um ponteiro nulo, imprime uma mensagem de erro e retorna sem fazer nada
        return;
    }

    // histos
    TH1F* hEta_LO  = (TH1F*)fLO ->Get("h_top_eta");                   //Obtém o histograma "h_top_eta" do arquivo fLO, que é o histograma de distribuição de pseudorrapidez (eta) dos tops para o processo pp -> ttbar (LO)
    TH1F* hEta_NLO = (TH1F*)fNLO->Get("h_top_eta");                   //Obtém o histograma "h_top_eta" do arquivo fNLO, que é o histograma de distribuição de pseudorrapidez (eta) dos tops para o processo pp -> ttbar (NLO)
    
    TH1F* hPt_LO  = (TH1F*)fLO ->Get("h_top_pt");                     //Obtém o histograma "h_top_pt" do arquivo fLO, que é o histograma de distribuição de momento transverso (pT) dos tops para o processo pp -> ttbar (LO)
    TH1F* hPt_NLO  = (TH1F*)fNLO->Get("h_top_pt");                    //Obtém o histograma "h_top_pt" do arquivo fNLO, que é o histograma de distribuição de momento transverso (pT) dos tops para o processo pp -> ttbar (NLO)

    TH1F* hEta_BH4 = (TH1F*)fBH4->Get("h_top_eta");                   //Obtém o histograma "h_top_eta" do arquivo fBH4, que é o histograma de distribuição de pseudorrapidez (eta) dos tops para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV
    TH1F* hEta_BH6 = (TH1F*)fBH6->Get("h_top_eta");                   //Obtém o histograma "h_top_eta" do arquivo fBH6, que é o histograma de distribuição de pseudorrapidez (eta) dos tops para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV
    TH1F* hEta_BHX = (TH1F*)fBHX->Get("h_top_eta");                   //Obtém o histograma "h_top_eta" do arquivo fBHX, que é o histograma de distribuição de pseudorrapidez (eta) dos tops para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV

    TH1F* hPt_BH4  = (TH1F*)fBH4->Get("h_top_pt");                    //Obtém o histograma "h_top_pt" do arquivo fBH4, que é o histograma de distribuição de momento transverso (pT) dos tops para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV
    TH1F* hPt_BH6  = (TH1F*)fBH6->Get("h_top_pt");                    //Obtém o histograma "h_top_pt" do arquivo fBH6, que é o histograma de distribuição de momento transverso (pT) dos tops para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV
    TH1F* hPt_BHX  = (TH1F*)fBHX->Get("h_top_pt");                    //Obtém o histograma "h_top_pt" do arquivo fBHX, que é o histograma de distribuição de momento transverso (pT) dos tops para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV

    if (!hEta_LO || !hPt_LO || !hEta_NLO || !hPt_NLO ||
        !hEta_BH4 || !hPt_BH4 || !hEta_BH6 || !hPt_BH6 ||             //Verifica se os histogramas de eta e pt foram obtidos corretamente dos arquivos .root, se algum deles for um ponteiro nulo, imprime uma mensagem de erro e retorna sem fazer nada
        !hEta_BHX || !hPt_BHX) {
        cout << "ERRO: algum histograma 'eta' ou 'pt' nao foi encontrado." << endl;
        return;
    }

    // ================================================================
    // NORMALIZA PARA σ_total (pb) 
    // 1) normaliza usando *Integral(width) -> soma a largura dos bins X os conteúdos, dando a integral total em pb para eta e pt 
    // ================================================================
    if (sigmaLO  > 0) NormalizeToSigma(hEta_LO,  sigmaLO);             //Normaliza o histograma de eta para o processo pp -> ttbar (LO) usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaLO em pb
    if (sigmaNLO > 0) NormalizeToSigma(hEta_NLO, sigmaNLO);            //Normaliza o histograma de eta para o processo pp -> ttbar (NLO) usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaNLO em pb
    if (sigmaBH4 > 0) NormalizeToSigma(hEta_BH4, sigmaBH4);            //Normaliza o histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaBH4 em pb
    if (sigmaBH6 > 0) NormalizeToSigma(hEta_BH6, sigmaBH6);            //Normaliza o histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaBH6 em pb
    if (sigmaBHX > 0) NormalizeToSigma(hEta_BHX, sigmaBHX);            //Normaliza o histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaBHX em pb

    if (sigmaLO  > 0) NormalizeToSigma(hPt_LO,  sigmaLO);              //Normaliza o histograma de pt para o processo pp -> ttbar (LO) usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaLO em pb
    if (sigmaNLO > 0) NormalizeToSigma(hPt_NLO, sigmaNLO);             //Normaliza o histograma de pt para o processo pp -> ttbar (NLO) usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaNLO em pb
    if (sigmaBH4 > 0) NormalizeToSigma(hPt_BH4, sigmaBH4);             //Normaliza o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaBH4 em pb
    if (sigmaBH6 > 0) NormalizeToSigma(hPt_BH6, sigmaBH6);             //Normaliza o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaBH6 em pb
    if (sigmaBHX > 0) NormalizeToSigma(hPt_BHX, sigmaBHX);             //Normaliza o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV usando a função NormalizeToSigma, retornando o histograma normalizado para a seção de choque total sigmaBHX em pb

    // canvas
    TCanvas* c = new TCanvas("c", "ttbar pQCD + BH", 1200, 500);       //Cria o canvas "c" com o título "ttbar pQCD + BH" e dimensões 1200x500 pixels, onde os histogramas serão desenhados
    // Canvas principal
    TPad* pad1 = new TPad("pad1","",0.00,0.25,0.50,1.00);              //Cria o pad "pad1" dentro do canvas "c" com as coordenadas (0.00, 0.25) para o canto inferior esquerdo e (0.50, 1.00) para o canto superior direito, ou seja, ocupando a metade esquerda do canvas para desenhar o histograma de eta
    TPad* pad2 = new TPad("pad2","",0.50,0.25,1.00,1.00);              //Cria o pad "pad2" dentro do canvas "c" com as coordenadas (0.50, 0.25) para o canto inferior esquerdo e (1.00, 1.00) para o canto superior direito, ou seja, ocupando a metade direita do canvas para desenhar o histograma de pt
    TPad* padLeg = new TPad("padLeg","",0.00,0.00,1.00,0.25);          //Cria o pad "padLeg" dentro do canvas "c" com as coordenadas (0.00, 0.00) para o canto inferior esquerdo e (1.00, 0.25) para o canto superior direito, ou seja, ocupando a parte inferior do canvas para desenhar a legenda

    pad1->SetLeftMargin(0.18);                                         //Definição da margem esquerda do pad1 como 0.18, para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y
    pad1->SetRightMargin(0.05);                                        //Definição da margem direita do pad1 como 0.05, para deixar um pouco de espaço entre o final do gráfico e a borda do pad
    pad1->SetBottomMargin(0.15);                                       //Definição da margem inferior do pad1 como 0.15, para deixar espaço para o título do eixo X e os rótulos dos ticks do eixo X
    pad1->SetTopMargin(0.08);                                          //Definição da margem superior do pad1 como 0.08, para deixar um pouco de espaço entre o início do gráfico e a borda do pad

    pad2->SetLeftMargin(0.18);                                         //Definição da margem esquerda do pad2 como 0.18, para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y
    pad2->SetRightMargin(0.05);                                        //Definição da margem direita do pad2 como 0.05, para deixar um pouco de espaço entre o final do gráfico e a borda do pad
    pad2->SetBottomMargin(0.15);                                       //Definição da margem inferior do pad2 como 0.15, para deixar espaço para o título do eixo X e os rótulos dos ticks do eixo X
    pad2->SetTopMargin(0.08);                                          //Definição da margem superior do pad2 como 0.08, para deixar um pouco de espaço entre o início do gráfico e a borda do pad

    pad1->Draw();                                                      //Desenha o pad1 no canvas, sendo o pad onde o histograma de eta será desenhado
    pad2->Draw();                                                      //Desenha o pad2 no canvas, sendo o pad onde o histograma de pt será desenhado
    padLeg->Draw();                                                    //Desenha o padLeg no canvas, sendo o pad onde a legenda será desenhada
  
    // =================================================================
    // PAINEL 1 — η  (dσ/dη)
    // =================================================================
    pad1->cd();                                                        //pad1 é o pad onde o histograma de eta será desenhado, então chamamos pad1->cd() para definir que as próximas operações de desenho serão feitas nesse pad
    gPad->SetLogy();                                                   //Define o eixo Y do pad1 como logarítmico, para melhor visualização das distribuições de eta, que podem variar em várias ordens de magnitude

    hEta_LO->SetLineColor(kBlack);                                   //Define a cor da linha do histograma de eta para o processo pp -> ttbar (LO) como preto (kBlack)
    hEta_LO->SetLineWidth(4);                                        //Define a largura da linha do histograma de eta para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível

    hEta_NLO->SetLineColor(kGray+1);                                 //Define a cor da linha do histograma de eta para o processo pp -> ttbar (NLO) como cinza claro (kGray+1)
    hEta_NLO->SetLineStyle(2);                                       //Define o estilo da linha do histograma de eta para o processo pp -> ttbar (NLO) como 2, que é uma linha tracejada
    hEta_NLO->SetLineWidth(4);                                       //Define a largura da linha do histograma de eta para o processo pp -> ttbar (NLO) como 4, para deixar a linha bem visível

    hEta_BH4->SetLineColor(kRed);                                    //Define a cor da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV como vermelho (kRed)
    hEta_BH4->SetLineStyle(3);                                       //Define o estilo da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV como 3, que é uma linha pontilhada
    hEta_BH4->SetLineWidth(4);                                       //Define a largura da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV como 4, para deixar a linha bem visível

    hEta_BH6->SetLineColor(kGreen+2);                                //Define a cor da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV como verde (kGreen+2)
    hEta_BH6->SetLineStyle(4);                                       //Define o estilo da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV como 4, que é uma linha com traços e pontos
    hEta_BH6->SetLineWidth(4);                                       //Define a largura da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV como 4, para deixar a linha bem visível

    hEta_BHX->SetLineColor(kBlue);                                   //Define a cor da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV como azul (kBlue)
    hEta_BHX->SetLineStyle(5);                                       //Define o estilo da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV como 5, que é uma linha com traços longos e pontos
    hEta_BHX->SetLineWidth(4);                                       //Define a largura da linha do histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV como 4, para deixar a linha bem visível

    hEta_LO->SetTitle(";#eta;d#sigma/d#eta [pb]");                   //Define o título do histograma de eta para o processo pp -> ttbar (LO) como ";#eta;d#sigma/d#eta [pb]", onde o primeiro ";" indica que não há título para o eixo Y, "#eta" é o título do eixo X, e "d#sigma/d#eta [pb]" é o título do eixo Y, indicando que a distribuição é a seção de choque diferencial em função de eta em pb
    hEta_LO->GetXaxis()->CenterTitle(true);                          //Centraliza o título do eixo X do histograma de eta para o processo pp -> ttbar (LO), para deixar o título do eixo X bem posicionado
    hEta_LO->GetYaxis()->CenterTitle(true);                          //Centraliza o título do eixo Y do histograma de eta para o processo pp -> ttbar (LO), para deixar o título do eixo Y bem posicionado

    hEta_LO->Draw("HIST");                                           //Desenha o histograma de eta para o processo pp -> ttbar (LO) usando a opção "HIST", que desenha apenas a linha do histograma sem os erros
    hEta_NLO->Draw("HIST SAME");                                     //Desenha o histograma de eta para o processo pp -> ttbar (NLO) usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe ao histograma anterior (LO)
    hEta_BH4->Draw("HIST SAME");                                     //Desenha o histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO e NLO)
    hEta_BH6->Draw("HIST SAME");                                     //Desenha o histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO, NLO e BH4)
    hEta_BHX->Draw("HIST SAME");                                     //Desenha o histograma de eta para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO, NLO, BH4 e BH6)

  
   // =================================================================
   // PAINEL 2 — pT  (dσ/dpT)
   // =================================================================
    pad2->cd();                                                        //pad2 é o pad onde o histograma de pt será desenhado, então chamamos pad2->cd() para definir que as próximas operações de desenho serão feitas nesse pad
    gPad->SetLogy();                                                   //Define o eixo Y do pad2 como logarítmico, para melhor visualização das distribuições de pt, que podem variar em várias ordens de magnitude 

    if (hPt_LO) {
        hPt_LO->SetLineColor(kBlack);                                //Define a cor da linha do histograma de pt para o processo pp -> ttbar (LO) como preto (kBlack)
        hPt_LO->SetLineWidth(4);                                     //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível

        hPt_LO->GetXaxis()->SetRangeUser(0, 1200);                   //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 1200 GeV, para focar na região de interesse do gráfico de pt
        hPt_LO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");    //Define o título do histograma de pt para o processo pp -> ttbar (LO) como ";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]", onde o primeiro ";" indica que não há título para o eixo Y, "p_{T} [GeV]" é o título do eixo X, e "d#sigma/dp_{T} [pb/GeV]" é o título do eixo Y, indicando que a distribuição é a seção de choque diferencial em função de pt em pb/GeV
        hPt_LO->GetXaxis()->CenterTitle(true);                       //Centraliza o título do eixo X do histograma de pt para o processo pp -> ttbar (LO), para deixar o título do eixo X bem posicionado
        hPt_LO->GetYaxis()->CenterTitle(true);                       //Centraliza o título do eixo Y do histograma de pt para o processo pp -> ttbar (LO), para deixar o título do eixo Y bem posicionado
        hPt_LO->GetYaxis()->SetTitleOffset(1.7);                     //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO) como 1.7, para deixar o título do eixo Y mais afastado do eixo, já que o título do eixo Y é rotacionado e fica mais próximo do eixo

        hPt_LO->Draw("HIST");                                        //Desenha o histograma de pt para o processo pp -> ttbar (LO) usando a opção "HIST", que desenha apenas a linha do histograma sem os erros
    }

    if (hPt_NLO) {
        hPt_NLO->SetLineColor(kGray+1);                              //Define a cor da linha do histograma de pt para o processo pp -> ttbar (NLO) como cinza claro (kGray+1)
        hPt_NLO->SetLineWidth(4);                                    //Define a largura da linha do histograma de pt para o processo pp -> ttbar (NLO) como 4, para deixar a linha bem visível
        hPt_NLO->SetLineStyle(2);                                    //Define o estilo da linha do histograma de pt para o processo pp -> ttbar (NLO) como 2, que é uma linha tracejada
        hPt_NLO->Draw("HIST SAME");                                  //Desenha o histograma de pt para o processo pp -> ttbar (NLO) usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe ao histograma anterior (LO)
    }

    if (hPt_BH4) { hPt_BH4->SetLineColor(kRed);     hPt_BH4->SetLineStyle(3); hPt_BH4->SetLineWidth(4); hPt_BH4->Draw("HIST SAME"); }       //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO e NLO)
    if (hPt_BH6) { hPt_BH6->SetLineColor(kGreen+2); hPt_BH6->SetLineStyle(4); hPt_BH6->SetLineWidth(4); hPt_BH6->Draw("HIST SAME"); }       //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO, NLO e BH4)
    if (hPt_BHX) { hPt_BHX->SetLineColor(kBlue);    hPt_BHX->SetLineStyle(5); hPt_BHX->SetLineWidth(4); hPt_BHX->Draw("HIST SAME"); }       //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO, NLO, BH4 e BH6)

    
   pad2->Modified();                                                  //Atualiza o pad2 com as modificações feitas, para garantir que as mudanças nas margens, escala logarítmica e desenho dos histogramas sejam refletidas corretamente no pad2
   pad2->Update();                                                    //Atualiza o pad2 para garantir que as mudanças feitas sejam feitas corretamente no pad2, como o desenho dos histogramas e a configuração do eixo Y em escala logarítmica
  
   // =================================================================
   // Legenda centralizada no rodapé
   // =================================================================
    padLeg->cd();                                                     //padLeg é o pad onde a legenda será desenhada, então chamamos padLeg->cd() para definir que as próximas operações de desenho serão feitas nesse pad
    TLegend* leg = new TLegend(0.10,0.15,0.90,0.85);                  //Cria a legenda com as coordenadas (0.10, 0.15) para o canto inferior esquerdo e (0.90, 0.85) para o canto superior direito, ou seja, ocupando a maior parte do padLeg para acomodar a legenda dos processos desenhados nos pads anteriores
    leg->SetBorderSize(0);                                            //Define o tamanho da borda da legenda como 0, ou seja, sem borda, para deixar a legenda mais limpa e integrada ao gráfico
    leg->SetFillStyle(0);                                             //Define o estilo de preenchimento da legenda como 0, ou seja, sem preenchimento, para deixar a legenda transparente e não cobrir os histogramas desenhados nos pads anteriores
    leg->SetTextFont(42);                                             //Define a fonte do texto da legenda como 42, que é uma fonte sem serifa e legível, para deixar a legenda clara e fácil de ler
    leg->SetNColumns(2);                                              //Define o número de colunas da legenda como 2, para organizar os itens da legenda em duas colunas e aproveitar melhor o espaço disponível no padLeg
    leg->SetTextSize(0.15);

    leg->AddEntry(hPt_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l");   //Legenda para o histograma de pt do processo pp -> ttbar (LO), com a descrição "pp #rightarrow t#bar{t} LO (pQCD)" e o estilo de linha "l" (linha)
    leg->AddEntry(hPt_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l");   //Legenda para o histograma de pt do processo pp -> ttbar (NLO), com a descrição "pp #rightarrow t#bar{t} NLO (pQCD)" e o estilo de linha "l" (linha)
    leg->AddEntry(hPt_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",    "l");    //Legenda para o histograma de pt do processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV, com a descrição "n=6, M_{D}=4 TeV, M_{BH}=8 TeV" e o estilo de linha "l" (linha)
    leg->AddEntry(hPt_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",    "l");    //Legenda para o histograma de pt do processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV, com a descrição "n=6, M_{D}=4 TeV, M_{BH}=9 TeV" e o estilo de linha "l" (linha)
    leg->AddEntry(hPt_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",   "l");    //Legenda para o histograma de pt do processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV, com a descrição "n=6, M_{D}=4 TeV, M_{BH}=10 TeV" e o estilo de linha "l" (linha)

    leg->Draw();                                                      //Desenha a legenda no padLeg, para que a legenda seja exibida na parte inferior do canvas, centralizada e organizada em duas colunas, associando cada item da legenda aos histogramas desenhados nos pads anteriores

    c->cd();                                                          //Volta para o canvas principal, para que as próximas operações de desenho sejam feitas no canvas principal, como o título geral do gráfico
    TLatex header;                                                    //Cria um objeto TLatex chamado header, que será usado para desenhar o título geral do gráfico no canvas principal 
    header.SetNDC(true);                                              //Define o sistema de coordenadas do header como NDC (Normalized Device Coordinates), para que as posições do texto sejam definidas em relação ao tamanho do canvas, facilitando a centralização e posicionamento do título
    header.SetTextFont(42);                                           //Define a fonte do texto do header como 42, que é uma fonte sem serifa e legível, para deixar o título claro e fácil de ler
    header.SetTextSize(0.035);                                        //Define o tamanho do texto do header como 0.035, para deixar o título em um tamanho adequado e proporcional ao canvas
    header.SetTextAlign(13);                                          //Define o alinhamento do texto do header como 13, que é alinhamento topo-esquerda, para posicionar o título no canto superior esquerdo do canvas

    header.DrawLatex(0.02, 0.97, "FCC  #sqrt{s} = 100 TeV");          //Desenha o texto "FCC  #sqrt{s} = 100 TeV" no canvas principal usando o objeto header, com as coordenadas (0.02, 0.97) em NDC, ou seja, próximo ao canto superior esquerdo do canvas, indicando que a figura é para o FCC com energia de 100 TeV
    
    // ===============================
    // Salvar
    // ===============================
     c->Update();                                                    //Atualiza o canvas para garantir que todas as modificações feitas, como o desenho dos histogramas, a configuração dos pads e a adição do título, sejam refletidas corretamente no canvas antes de salvar a figura
     c->SaveAs("ttbar_overlay_2pads_legenda.png");                   //Salva o canvas "c" como um arquivo de imagem no formato PNG com o nome "ttbar_overlay_2pads_legenda.png", para que a figura possa ser visualizada e compartilhada em formato de imagem
     c->SaveAs("ttbar_overlay_2pads_legenda.pdf");                   //Salva o canvas "c" como um arquivo de imagem no formato PDF com o nome "ttbar_overlay_2pads_legenda.pdf", para que a figura possa ser visualizada e compartilhada em formato de documento PDF, mantendo a qualidade vetorial da figura
 
    // ===============================
    // PAINEL 3 - pT (dσ/dpT) - Gráfico separado só para pT até 3500 GeV, para mostrar melhor a região de alta pT onde os efeitos de BH podem ser mais evidentes, já que no gráfico anterior o eixo X foi limitado a 1200 GeV para focar na região de interesse, mas isso pode esconder as diferenças entre os modelos em pT mais altos, então criamos um gráfico separado para mostrar essa região
    // ===============================
    TCanvas* c2 = new TCanvas("c2", "Figura 2", 1200, 800);         //Tamanho do canvas para o gráfico de pt, para ter mais espaço para mostrar a região de alta pT e os detalhes das distribuições nessa região

    c2->cd();                                                       //Desenha o canvas c2, onde o gráfico de pt será desenhado 
    gPad->SetLeftMargin(0.16);                                      //Margem esquerda do pad para deixar espaço para o título do eixo Y e os rótulos dos ticks do eixo Y 
    gPad->SetLogy();                                                //Define o eixo Y do canvas c2 como logarítmico, para melhor visualização das distribuições de pt, que podem variar em várias ordens de magnitude, especialmente na região de alta pT onde os efeitos de BH podem ser mais evidentes

    if (hPt_LO) {                                                  //Verifica se o histograma de pt para o processo pp -> ttbar (LO) existe, para evitar erros caso o histograma não tenha sido obtido corretamente dos arquivos .root
        hPt_LO->SetLineColor(kBlack);                              //Define a cor da linha do histograma de pt para o processo pp -> ttbar (LO) como preto (kBlack)
        hPt_LO->SetLineWidth(4);                                   //Define a largura da linha do histograma de pt para o processo pp -> ttbar (LO) como 4, para deixar a linha bem visível

        hPt_LO->GetXaxis()->SetRangeUser(0, 3500);                 //Define o intervalo do eixo X do histograma de pt para o processo pp -> ttbar (LO) de 0 a 3500 GeV, para mostrar melhor a região de alta pT onde os efeitos de BH podem ser mais evidentes, já que no gráfico anterior o eixo X foi limitado a 1200 GeV para focar na região de interesse, mas isso pode esconder as diferenças entre os modelos em pT mais altos, então criamos um gráfico separado para mostrar essa região
        hPt_LO->SetTitle(";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]");  //Titulo do histograma de pt para o processo pp -> ttbar (LO), com a descrição ";p_{T} [GeV];d#sigma/dp_{T} [pb/GeV]", onde o primeiro ";" indica que não há título para o eixo Y, "p_{T} [GeV]" é o título do eixo X, e "d#sigma/dp_{T} [pb/GeV]" é o título do eixo Y, indicando que a distribuição é a seção de choque diferencial em função de pt em pb/GeV
        hPt_LO->GetXaxis()->CenterTitle(true);                     //Centraliza o título do eixo X do histograma de pt para o processo pp -> ttbar (LO), para deixar o título do eixo X bem posicionado
        hPt_LO->GetYaxis()->CenterTitle(true);                     //Centraliza o título do eixo Y do histograma de pt para o processo pp -> ttbar (LO), para deixar o título do eixo Y bem posicionado
        hPt_LO->GetYaxis()->SetTitleOffset(1.7);                   //Define o deslocamento do título do eixo Y do histograma de pt para o processo pp -> ttbar (LO) como 1.7, para deixar o título do eixo Y mais afastado do eixo, já que o título do eixo Y é rotacionado e fica mais próximo do eixo

        hPt_LO->Draw("HIST");                                      //Desenha o histograma de pt para o processo pp -> ttbar (LO) usando a opção "HIST", que desenha apenas a linha do histograma sem os erros, para mostrar a distribuição de pt para o processo pp -> ttbar (LO) no gráfico separado para a região de alta pT
    }

    if (hPt_NLO) {                                                 //Verifica se o histograma de pt para o processo pp -> ttbar (NLO) existe, para evitar erros caso o histograma não tenha sido obtido corretamente dos arquivos .root
        hPt_NLO->SetLineColor(kGray+1);                            //Define a cor da linha do histograma de pt para o processo pp -> ttbar (NLO) como cinza claro (kGray+1)
        hPt_NLO->SetLineWidth(4);                                  //Define a largura da linha do histograma de pt para o processo pp -> ttbar (NLO) como 4, para deixar a linha bem visível
        hPt_NLO->SetLineStyle(2);                                  //Define o estilo da linha do histograma de pt para o processo pp -> ttbar (NLO) como 2, que é uma linha tracejada, para diferenciar visualmente do histograma de pt para o processo pp -> ttbar (LO) e destacar a contribuição das correções de NLO
        hPt_NLO->Draw("HIST SAME");                                //Desenha o histograma de pt para o processo pp -> ttbar (NLO) usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe ao histograma anterior (LO), para mostrar a distribuição de pt para o processo pp -> ttbar (NLO) no gráfico separado para a região de alta pT, permitindo comparar visualmente com o histograma de pt para o processo pp -> ttbar (LO) e destacar as diferenças entre os dois processos
    }

    if (hPt_BH4) { hPt_BH4->SetLineColor(kRed);     hPt_BH4->SetLineStyle(3); hPt_BH4->SetLineWidth(4); hPt_BH4->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO e NLO), para mostrar a distribuição de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV no gráfico separado para a região de alta pT, permitindo comparar visualmente com os histogramas de pt para os processos pp -> ttbar (LO) e pp -> ttbar (NLO) e destacar as diferenças entre os modelos
    if (hPt_BH6) { hPt_BH6->SetLineColor(kGreen+2); hPt_BH6->SetLineStyle(4); hPt_BH6->SetLineWidth(4); hPt_BH6->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO, NLO e BH4), para mostrar a distribuição de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV no gráfico separado para a região de alta pT, permitindo comparar visualmente com os histogramas de pt para os processos pp -> ttbar (LO), pp -> ttbar (NLO) e pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV e destacar as diferenças entre os modelos
    if (hPt_BHX) { hPt_BHX->SetLineColor(kBlue);    hPt_BHX->SetLineStyle(5); hPt_BHX->SetLineWidth(4); hPt_BHX->Draw("HIST SAME"); }      //Desenha o histograma de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV usando a opção "HIST SAME", que desenha apenas a linha do histograma e sobrepõe aos histogramas anteriores (LO, NLO, BH4 e BH6), para mostrar a distribuição de pt para o processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV no gráfico separado para a região de alta pT, permitindo comparar visualmente com os histogramas de pt para os processos pp -> ttbar (LO), pp -> ttbar (NLO), pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV e pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV e destacar as diferenças entre os modelos


    c2->Modified();                                                 //Atualiza o canvas c2 com as modificações feitas, para garantir que as mudanças nas margens, escala logarítmica e desenho dos histogramas sejam refletidas corretamente no canvas c2
    c2->Update();                                                   //Atualiza o canvas c2 para garantir que as mudanças feitas sejam feitas corretamente no canvas c2, como o desenho dos histogramas e a configuração do eixo Y em escala logarítmica

    TLegend* legPt2 = new TLegend(0.55, 0.65, 0.88, 0.88);          //Cria a legenda para o gráfico de pt no canvas c2, com as coordenadas (0.55, 0.65) para o canto inferior esquerdo e (0.88, 0.88) para o canto superior direito, ou seja, posicionada no canto superior direito do canvas c2 para não sobrepor os histogramas desenhados e permitir uma visualização clara da legenda associada ao gráfico de pt

    legPt2->SetBorderSize(0);                                       //Define o tamanho da borda da legenda como 0, ou seja, sem borda, para deixar a legenda mais limpa e integrada ao gráfico
    legPt2->SetFillStyle(0);                                        //Define o estilo de preenchimento da legenda como 0, ou seja, sem preenchimento, para deixar a legenda transparente e não cobrir os histogramas desenhados no canvas c2
    legPt2->SetTextFont(42);                                        //Define a fonte do texto da legenda como 42, que é uma fonte sem serifa e legível, para deixar a legenda clara e fácil de ler
    legPt2->SetTextSize(0.03);                                      //Define o tamanho do texto da legenda como 0.03, para deixar a legenda em um tamanho adequado e proporcional ao canvas c2, garantindo que a legenda seja legível sem ocupar muito espaço no canvas

    legPt2->AddEntry(hPt_LO,  "pp #rightarrow t#bar{t} LO (pQCD)",  "l"); //Legenda para o histograma de pt do processo pp -> ttbar (LO), com a descrição "pp #rightarrow t#bar{t} LO (pQCD)" e o estilo de linha "l" (linha)
    legPt2->AddEntry(hPt_NLO, "pp #rightarrow t#bar{t} NLO (pQCD)", "l"); //Legenda para o histograma de pt do processo pp -> ttbar (NLO), com a descrição "pp #rightarrow t#bar{t} NLO (pQCD)" e o estilo de linha "l" (linha)
    legPt2->AddEntry(hPt_BH4, "n=6, M_{D}=4 TeV, M_{BH}=8 TeV",     "l"); //Legenda para o histograma de pt do processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=8 TeV, com a descrição "n=6, M_{D}=4 TeV, M_{BH}=8 TeV" e o estilo de linha "l" (linha)
    legPt2->AddEntry(hPt_BH6, "n=6, M_{D}=4 TeV, M_{BH}=9 TeV",     "l"); //Legenda para o histograma de pt do processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=9 TeV, com a descrição "n=6, M_{D}=4 TeV, M_{BH}=9 TeV" e o estilo de linha "l" (linha)
    legPt2->AddEntry(hPt_BHX, "n=6, M_{D}=4 TeV, M_{BH}=10 TeV",    "l"); //Legenda para o histograma de pt do processo pp -> ttbar via BH com n=6, MD=4 TeV, MBH=10 TeV, com a descrição "n=6, M_{D}=4 TeV, M_{BH}=10 TeV" e o estilo de linha "l" (linha)

    legPt2->Draw();                                                  //Desenha a legenda para o gráfico de pt no canvas c2, para que a legenda seja exibida no canto superior direito do canvas c2, associando cada item da legenda aos histogramas desenhados no gráfico de pt e permitindo uma visualização clara da legenda sem sobrepor os histogramas
    
    // =================================================================
    // Salvar
    // =================================================================
    c2->Update();                                                    //Atualiza o canvas c2 para garantir que todas as modificações feitas, como o desenho dos histogramas, a configuração dos pads e a adição da legenda, sejam refletidas corretamente no canvas c2 antes de salvar a figura
    c2->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh.eps");            //Salva o canvas "c2" como um arquivo de imagem no formato EPS com o nome "ttbar_overlay_100TeV_5M_n6_md4_mbh.eps", para que a figura possa ser visualizada e compartilhada em formato de documento EPS, mantendo a qualidade vetorial da figura
    c2->SaveAs("ttbar_overlay_100TeV_5M_n6_md4_mbh.png");            //Salva o canvas "c2" como um arquivo de imagem no formato PNG com o nome "ttbar_overlay_100TeV_5M_n6_md4_mbh.png", para que a figura possa ser visualizada e compartilhada em formato de imagem
   

    cout << "Figuras salvas em ttbar_overlay_100TeV_5M_n6_md4_mbh.[eps,png]" << endl;   //Imprime no console a mensagem indicando que as figuras foram salvas com os nomes "ttbar_overlay_100TeV_5M_n6_md4_mbh.eps" e "ttbar_overlay_100TeV_5M_n6_md4_mbh.png", para informar ao usuário que as figuras foram geradas e estão disponíveis para visualização e compartilhamento
}