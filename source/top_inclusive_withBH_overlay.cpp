/*
 * Análise de produção de bb̄ em colisões O-O
 * Identifica quarks b ou hádrons B e analisa o sistema BB
 *
 * Compilar: g++ -o analyze_bb analyze_bb.cpp `root-config --cflags --libs`
 * Executar: ./analyze_bb input.root output.root [cross_section_nb]
 */
#include <TFile.h>                      
#include <TH1F.h>                       
#include <TH2F.h>                       
#include <TLorentzVector.h>  
#include <cstdlib>           
#include <TTree.h> 
#include <TChain.h>                     
#include <algorithm>                    
#include <cmath>                        
#include <iostream>                     
#include <string>                       
#include <vector>                       

// Contadores 
double xs = 0.0;                        //nb isospin nPDFSet0  
double lumi = 1.0;                      //nb-1

int main(int argc, char **argv)         //Função principal do programa
{

    if (argc < 3)                       
    {
        std::cout << "Uso: " << argv[0] << " <input.root> <output.root> [cross_section_nb]" << std::endl;  
        std::cout << "  cross_section_nb: seção de choque em nb (opcional)" << std::endl;                  
        return 1;                                                                                          
    }

    std::string output_name = argv[2];                                                                     
    
    // ==== Configuração da seção de choque ====
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
    //Para um unico arquivo de entrada, usaríamos:
    //TChain *tree = new TChain("lheTree");
    //tree->Add(argv[1]);

    //Para o caso específico de arquivos separados para top e antitop (caso tW), podemos adicionar ambos à TChain
    TChain *tree = new TChain("lheTree");
    tree->Add("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/single-top/results/root_files/single_top_tWchannel_top_converted.root");                      //sigmatW_top = 1177 pb
    tree->Add("/home/brenda_rolin/Documentos/programas/BlackMax-2.02.0/BlackMax/ttbar_BH_analysis/single-top/results/root_files/single_top_tWchannel_antitop_converted.root");                  //sigmatW_antitop = 1177 pb

    // Obter a árvore do arquivo de entrada
    if (!tree)                                                                                         //Verifica se a tree foi obtida com sucesso. Se o ponteiro tree for nulo
    {
        std::cerr << "Erro: árvore 'lheTree' não encontrada!" << std::endl;                            //Mensagem de erro indicando que a árvore "lheTree" não foi encontrada no arquivo de entrada, mostrando o nome da árvore que causou o erro
        return 1;                                                                                      //Retorna 1 para indicar que o programa terminou com um erro devido à falha ao encontrar a árvore "lheTree" no arquivo de entrada 
    }
    // Variáveis para leitura dos dados da árvore
    std::vector<int> *pdgID = nullptr;                                                                 //Ponteiro para um vetor de inteiros que armazenará os códigos PDG (Particle Data Group) das partículas presentes em cada evento
    std::vector<int> *status = nullptr;                                                                //Ponteiro para um vetor de inteiros que armazenará os status das partículas presentes em cada evento
    std::vector<float> *px = nullptr;                                                                  //Ponteiro para um vetor de floats que armazenará as componentes de momento px das partículas presentes em cada evento
    std::vector<float> *py = nullptr;                                                                  //Ponteiro para um vetor de floats que armazenará as componentes de momento py das partículas presentes em cada evento
    std::vector<float> *pz = nullptr;                                                                  //Ponteiro para um vetor de floats que armazenará as componentes de momento pz das partículas presentes em cada evento
    std::vector<float> *e = nullptr;                                                                   //Ponteiro para um vetor de floats que armazenará as energias das partículas presentes em cada evento
    std::vector<float> *m = nullptr;                                                                   //Ponteiro para um vetor de floats que armazenará as massas das partículas presentes em cada evento

    tree->SetBranchAddress("pdgID", &pdgID);                                                           //Endereço das branches da árvore para ler os dados. A função SetBranchAddress do ROOT é usada para associar as variáveis do programa (pdgID, status, px, py, pz, e, m)
    tree->SetBranchAddress("status", &status);                                                         //Associa a variável status ao endereço da branch "status" 
    tree->SetBranchAddress("px", &px);                                                                 //Associa a variável px ao endereço da branch "px" 
    tree->SetBranchAddress("py", &py);                                                                 //Associa a variável py ao endereço da branch "py" 
    tree->SetBranchAddress("pz", &pz);                                                                 //Associa a variável pz ao endereço da branch "pz" 
    tree->SetBranchAddress("e", &e);                                                                   //Associa a variável e ao endereço da branch "e" 
    tree->SetBranchAddress("m", &m);                                                                   //Associa a variável m ao endereço da branch "m"


    // Histogramas
    TH1F *h_top_pt = new TH1F("h_top_pt", "pT dos tops;pT_{top} [GeV];Eventos", 100, 0, 4000);         //Cria um histograma unidimensional do tipo float chamado h_top_pt
    TH1F *h_ntop = new TH1F("h_ntop", "NTops;N_{top};Eventos", 10, 0, 10);                             //Cria um histograma unidimensional do tipo float chamado h_ntop
    TH1F *h_top_eta = new TH1F("h_top_eta", "Eta dos tops;#eta_{top};Eventos", 100, -10, 10); 
    //TH1F *h_top_pt_selection = new TH1F("h_top_pt_selection", "pT dos tops (seleção);pT_{top} [GeV];Eventos", 100, 0, 4000); 
    //TH1F *h_top_eta_selection = new TH1F("h_top_eta_selection","eta dos tops (seleção);#eta_{top};Eventos", 100, -10, 10);
    Long64_t nEntries = tree->GetEntries();                                                           //Numero total de eventos na árvore "lheTree"
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

                n_top++;
                h_top_pt->Fill(top.Pt());                                                                                                
                h_top_eta->Fill(top.Eta());   
            
                //if (top.Pt() > 1000)
                //{
                //h_top_pt_selection->Fill(top.Pt());        //para pt
                //h_top_eta_selection->Fill(top.Eta());      //para eta                                  
                //}
            } 
        } //fim do loop sobre partículas em um evento

        h_ntop->Fill(n_top); //multiplicidade
        // Progresso 
        if (iEvent % 1000 == 0)                                                                    
        {
            std::cout << "Processado: " << iEvent << " / " << nEntries                              
                      << " (" << (100.0 * iEvent / nEntries) << "%)" << std::endl;                  
        }
    } //fim do loop sobre eventos

    h_top_pt->Scale(xs/nEntries); 
    h_top_eta->Scale(xs/nEntries);
    //h_top_pt_selection->Scale(xs/nEntries);
    //h_top_eta_selection->Scale(xs/nEntries);

    // Salvar histogramas
    TFile *output_file = new TFile(output_name.c_str(), "RECREATE");              

    h_top_pt->Write();                                                                             
    h_top_eta->Write();  
    //h_top_pt_selection->Write();
    h_ntop->Write();    
    //h_top_eta_selection->Write();                                                                        

    output_file->Close();                                                                          

    std::cout << "\nHistogramas salvos em: " << output_name << std::endl;                         
    std::cout << "Arquivo de análise salvo em: out_bbbar.txt" << std::endl;                       

    // Fechar arquivo de entrada                                                                                                                                            
    delete output_file;                                                                            

    return 0;                                                                       
}