# PROJETO DE PROGRAMAÇÃO: COMPARAÇÃO DE GENES
Este projeto tem como aplicação a comparação de 2 ou 3 variantes do mesmo gene, de forma a perceber a percentagem de semelhança (identity) entre as variantes. Adicionalmente são obtidos os alinhamentos entre as sequências.

### Features deste projeto:
- Leitura de sequências e ficheiros FASTA: públicos e providenciados pelo professor
- Calcular LCS para 2 e 3 Sequências (Algoritmo LCS)
- Alinhamento de 2 e 3 Sequências (Algoritmo Needleman-Wunsch)
- Estimar a identity
- Componente gráfica e visual:
    - Interface de controlo do programa
    - Mapa de cores sobre sequências alinhadas
    - Gráfico de dimensões de sequências

## Requirements

Projeto testado com versão `Python 3.13.3`

- Instalar seguinte package:
```powershell
pip install matplotlib
```

## Usage
Para correr este projeto, pode utilizar um dos seguintes comandos exemplo na pasta dos ficheiros python:
```powershell
python .\main.py                    # Comando default
python .\main.py .\data\filename    # Dar ficheiro de entrada para leitura durante inicialização
```

## Dados
### Ficheiros providenciados pelo docente:
- gene_example.fasta
- prot_example.fasta

### Ficheiros obtidos em repositorios [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
- butterfly_genes.fasta - Combinação das sequências:
    - https://www.ncbi.nlm.nih.gov/nuccore/XR_011173596.1 
    - https://www.ncbi.nlm.nih.gov/nuccore/XR_011173595.1
    - https://www.ncbi.nlm.nih.gov/nuccore/XR_005288263.1
- turtle_genes.fasta - Combinação das sequências: 
    - https://www.ncbi.nlm.nih.gov/nuccore/NM_001299614.1 
    - https://www.ncbi.nlm.nih.gov/nuccore/NM_001103726.3  
    - https://www.ncbi.nlm.nih.gov/nuccore/NM_001299384.1