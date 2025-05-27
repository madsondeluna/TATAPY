# TATA.py - Analisador de Sequências FASTA

Uma ferramenta de linha de comando, escrita em Python, para realizar diversas análises genéticas e de bioinformática em arquivos no formato FASTA. Este script é ideal para estudantes, pesquisadores e qualquer pessoa que precise de uma análise rápida e automatizada de sequências de DNA.

## Principais Funcionalidades

* **Análise de Composição**: Calcula o conteúdo GC de cada sequência e o conteúdo GC geral do arquivo.

* **Busca por Fatores de Transcrição**: Identifica a localização de sítios de ligação para diversos fatores de transcrição (SP1, TATA-box, CREB, etc.) usando expressões regulares (Regex).

* **Busca por ORFs (Open Reading Frames)**: Encontra potenciais regiões codificantes de proteínas (genes) acima de um tamanho mínimo especificado.

* **Análise de Proteínas**: Para cada ORF encontrado, calcula propriedades físico-químicas da proteína predita, como:

  * Peso Molecular (Da)
  * Ponto Isoelétrico (pI)
  * Índice de Hidrofobicidade (GRAVY)

* **Análise de Sítios de Restrição**: Simula a digesão das sequências com enzimas de restrição especificadas, mostrando onde os cortes ocorreriam.

## Requisitos

* Python 3.7 ou superior
* Biblioteca Biopython

## Instalação

Clone este repositório (ou baixe os arquivos):

```bash
git clone https://github.com/seu-usuario/seu-repositorio.git
cd seu-repositorio
```

Se não estiver usando Git, apenas baixe o arquivo `TATA.py` e coloque-o em uma pasta no seu computador.

Crie e ative um ambiente virtual (recomendado):

```bash
python -m venv venv
source venv/bin/activate  # No Windows: venv\Scripts\activate
```

Instale as dependências:

```bash
pip install biopython
```

## Uso

O script é executado a partir do terminal. A estrutura básica do comando é:

```bash
python TATA.py [ARQUIVO_FASTA] [OPÇÕES...]
```

### Argumentos e Opções

* `arquivo` (obrigatório): O caminho para o arquivo FASTA que você deseja analisar.
* `--composicao`: Ativa a análise de conteúdo GC.
* `--orfs <TAMANHO_MIN_AA>`: Ativa a busca por ORFs. Você deve fornecer o tamanho mínimo do ORF em aminoácidos.
* `--restricao <ENZ1,ENZ2,...>`: Ativa a análise de restrição. Forneça uma lista de enzimas separadas por vírgula, sem espaços (ex: EcoRI,BamHI).
* `--<fator>`: Ativa a busca por um sítio de fator de transcrição (ex: --tata, --sp1, --nfkb). Essas opções são mutuamente exclusivas.
* `-h` ou `--help`: Mostra a mensagem de ajuda com todas as opções disponíveis.

## Exemplos de Comandos

Imagine que você tem um arquivo chamado `meu_gene.fasta`.

1. Para obter uma análise de composição (GC content):

```bash
python TATA.py meu_gene.fasta --composicao
```

2. Para encontrar todos os ORFs com mais de 50 aminoácidos:

```bash
python TATA.py meu_gene.fasta --orfs 50
```

3. Para encontrar sítios da TATA-box e da SP1 (executar separadamente):

```bash
python TATA.py meu_gene.fasta --tata
python TATA.py meu_gene.fasta --sp1
```

4. Para simular uma digesão com EcoRI e HindIII:

```bash
python TATA.py meu_gene.fasta --restricao EcoRI,HindIII
```

5. Combinando Análises:

```bash
python TATA.py meu_gene.fasta --composicao --orfs 70 --restricao BamHI
```

## Exemplo de Saída

```text
================================================================================
RELATÓRIO DE ANÁLISE DO ARQUIVO: meu_gene.fasta
TOTAL DE SEQUÊNCIAS: 1
================================================================================

--- Análise de Composição (GC Content) ---
Conteúdo GC Geral do Arquivo: 55.45%
  -> Gene_X_humano: 55.45%

--- Análise de ORFs (Tamanho Mínimo: 70 aa) ---
  Nenhum ORF correspondente aos critérios foi encontrado.

--- Análise de Sítios de Restrição ---
  -> Sequência ID: Gene_X_humano
    - Enzima 'BamHI' corta 1 vez(es) na(s) posição(oes): [123]
================================================================================
```

## Como Contribuir

Contribuições são sempre bem-vindas! Se você tem uma ideia para uma nova análise ou encontrou um bug, sinta-se à vontade para:

1. Fazer um *fork* deste repositório.
2. Criar uma nova *branch* (`git checkout -b feature/NovaAnalise`).
3. Fazer o *commit* das suas mudanças (`git commit -am 'Adiciona nova análise X'`).
4. Fazer o *push* para a *branch* (`git push origin feature/NovaAnalise`).
5. Abrir um *pull request*.
