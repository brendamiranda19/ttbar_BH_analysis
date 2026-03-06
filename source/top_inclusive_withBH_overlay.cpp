/*
 * Análise de produção de bb̄ em colisões O-O
 * Identifica quarks b ou hádrons B e analisa o sistema BB
 *
 * Compilar: g++ -o analyze_bb analyze_bb.cpp `root-config --cflags --libs`
 * Executar: ./analyze_bb input.root output.root [cross_section_nb]
 */

//PARA SEGUNDA-FEIRA 
//INCLUIR O ETA (HIST) -10 ATÉ 10
//ADICIONAR NO GITHUB ESSE CODIGO COM ETA E PT, PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB) - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB) - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB) - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB) - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB) - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB) - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB) - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB)
//COMENTAR CADA PARTE DO CÓDIGO PARA QUE FIQUE CLARO O QUE CADA PARTE FAZ - PARA QUE SEJA USADO COMO BASE PARA OS PROXIMOS PLOTS DE SOBREPOSIÇÃO (ETA E PT) COM OS HISTOGRAMAS NORMALIZADOS PARA SEÇÃO DE CHOQUE TOTAL (PB)
//COMENTAR **TUDO**

#include <TFile.h>                      //Classe TFile da framework do ROOT, serve para manipular arquivos .root (abrir, criar, ler, escrever)
#include <TH1F.h>                       //Classe TH1F da framework do ROOT, serve para criar e manipular histogramas unidimensionais de tipo float
#include <TH2F.h>                       //Classe TH2F da framework do ROOT, serve para criar e manipular histogramas bidimensionais de tipo float
#include <TLorentzVector.h>             //Classe TLorentzVector da framework do ROOT, serve para criar e manipular vetores de quatro componentes (px, py, pz, E) usados para representar partículas em física de altas energias
#include <TTree.h>                      //Classe TTree da framework do ROOT, serve para criar e manipular árvores de dados, que são estruturas hierárquicas usadas para armazenar grandes conjuntos de dados de forma eficiente
#include <algorithm>                    //Classe algorithm da biblioteca standard do C++, serve para funções de algoritmos de estrutura de dados (sort, find, copy, etc.)
#include <cmath>                        //Classe cmath da biblioteca standard do C++, serve para funções matemáticas comuns (funções trigonométricas, exponenciais, etc.)
#include <iostream>                     //Classe iostream da biblioteca standard do C++, serve para manipulação de entrada e saída de dados (cout, cin, etc.)
#include <string>                       //Classe string da biblioteca standard do C++, serve para criar e manipular strings (cadeias de caracteres)
#include <vector>                       //Classe vector da biblioteca standard do C++, serve para criar e manipular vetores

// Contadores 
double xs = 0.0;                        //nb isospin nPDFSet0  
double lumi = 1.0;                      //nb-1

int main(int argc, char **argv)         //Função principal do programa, onde a execução começa. Recebe os argumentos da linha de comando (argc é o número de argumentos, argv é um array de strings com os argumentos)
{

    if (argc < 3)                       //Verifica se o número de argumentos é menor que 3 (o nome do programa, o arquivo de entrada e o arquivo de saída), ou seja, se os argumentos necessários não foram fornecidos. Se for o caso, exibe uma mensagem de uso e retorna 1 para indicar um erro na execução
    {
        std::cout << "Uso: " << argv[0] << " <input.root> <output.root> [cross_section_nb]" << std::endl;  //Exibe a mensagem de uso do programa, indicando que o programa deve ser executado com o nome do arquivo de entrada (input.root), o nome do arquivo de saída (output.root) e opcionalmente a seção de choque em nb (cross_section_nb)
        std::cout << "  cross_section_nb: seção de choque em nb (opcional)" << std::endl;                  //Exibe uma mensagem adicional explicando que o argumento cross_section_nb é opcional e representa a seção de choque em nb, caso o usuário queira fornecer um valor diferente do padrão
        return 1;                                                                                          //Retorna 1 para indicar que o programa terminou com um erro devido à falta de argumentos necessários
    }

    std::string output_name = argv[2];

    // Seção de choque fornecida como argumento (opcional)
    if (argc >= 4)
    {
        xs = std::atof(argv[3]);
        std::cout << "Seção de choque fornecida: " << xs << " nb" << std::endl;
    }
    else
    {
        std::cout << "Seção de choque padrão: " << xs << " nb" << std::endl;
    }

    // Abrir arquivo de entrada
    TFile *input_file = TFile::Open(argv[1]);
    if (!input_file || input_file->IsZombie())
    {
        std::cerr << "Erro ao abrir arquivo: " << argv[1] << std::endl;
        return 1;
    }

    // Obter a árvore
    TTree *tree = (TTree *)input_file->Get("lheTree");
    if (!tree)
    {
        std::cerr << "Erro: árvore 'lheTree' não encontrada!" << std::endl;
        input_file->Close();
        return 1;
    }

    // Variáveis para leitura
    std::vector<int> *pdgID = nullptr;
    std::vector<int> *status = nullptr;
    std::vector<float> *px = nullptr;
    std::vector<float> *py = nullptr;
    std::vector<float> *pz = nullptr;
    std::vector<float> *e = nullptr;
    std::vector<float> *m = nullptr;

    tree->SetBranchAddress("pdgID", &pdgID);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("m", &m);

    // Histogramas
    TH1F *h_top_pt = new TH1F("h_top_pt", "pT dos tops;pT_{top} [GeV];Eventos", 100, 0, 4000);
    TH1F *h_ntop = new TH1F("h_ntop", "NTops;N_{top};Eventos", 10, 0, 10);

    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processando " << nEntries << " eventos..." << std::endl;
    TLorentzVector top;

    // Loop sobre eventos
    for (Long64_t iEvent = 0; iEvent < nEntries; iEvent++)
    {
        int n_top = 0;
        tree->GetEntry(iEvent);

        

        for (size_t iPart = 0; iPart < pdgID->size(); iPart++)
        {
            int pdg = fabs(pdgID->at(iPart));
            int stat = status->at(iPart);

            if (pdg == 6 && stat == 1)
            {
                top.SetPxPyPzE(px->at(iPart), py->at(iPart), pz->at(iPart), e->at(iPart));
                h_top_pt->Fill(top.Pt());
                n_top++;
            }
        }
        h_ntop->Fill(n_top);

        // Progresso
        if (iEvent % 1000 == 0)
        {
            std::cout << "Processado: " << iEvent << " / " << nEntries
                      << " (" << (100.0 * iEvent / nEntries) << "%)" << std::endl;
        }
    }

    // Salvar histogramas
    TFile *output_file = new TFile(output_name.c_str(), "RECREATE");

    h_top_pt->Write();
    h_ntop->Write();

    output_file->Close();

    std::cout << "\nHistogramas salvos em: " << output_name << std::endl;
    std::cout << "Arquivo de análise salvo em: out_bbbar.txt" << std::endl;

    // Fechar arquivo de entrada
    input_file->Close();
    delete input_file;
    delete output_file;

    return 0;
}