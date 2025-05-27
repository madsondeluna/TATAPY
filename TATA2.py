# =============================================================================
# ANALISADOR GENÉTICO DE SEQUÊNCIAS FASTA
# =============================================================================
# Autor: Seu Nome (com ajuda da IA Gemini)
# Versão: 2.0
# Descrição: Uma ferramenta de linha de comando para realizar análises
#            básicas e avançadas em arquivos FASTA.
# =============================================================================

import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Restriction import Analysis as RestrictionAnalysis, RestrictionBatch

# --- Dicionário de Fatores de Transcrição (expansível) ---
SITIOS_REGEX = {
    "SP1": r'GGGCGG', "TATA": r'TATA[AT]A[AT]', "CREB": r'TGACGTCA',
    "NFKB": r'GGGA[CT]TTCC', "AP1": r'TGA[CG]TCA', "CAAT": r'GGCCAATCT',
    "EBOX": r'CA[ACGT]{2}TG'
}

# =============================================================================
# FUNÇÕES DE ANÁLISE INDIVIDUAIS
# =============================================================================

def analisar_composicao(sequencias):
    """Calcula o conteúdo GC de cada sequência e o geral."""
    resultados = {'por_sequencia': {}}
    todas_as_bases = ""
    for seq_record in sequencias:
        gc_content = GC(seq_record.seq)
        resultados['por_sequencia'][seq_record.id] = f"{gc_content:.2f}%"
        todas_as_bases += str(seq_record.seq)
    
    resultados['geral'] = f"{GC(Seq(todas_as_bases)):.2f}%"
    return resultados

def buscar_sitios_transcricao(sequencias, fator_pesquisa):
    """Busca por sítios de fatores de transcrição usando Regex."""
    padrao_regex = SITIOS_REGEX.get(fator_pesquisa)
    if not padrao_regex:
        return {"erro": f"Fator '{fator_pesquisa}' desconhecido."}

    resultados = {"fator": fator_pesquisa, "encontrados": {}}
    for seq_record in sequencias:
        posicoes = []
        for match in re.finditer(padrao_regex, str(seq_record.seq), re.IGNORECASE):
            posicoes.append((match.start(), match.end()))
        if posicoes:
            resultados["encontrados"][seq_record.id] = posicoes
    return resultados

def analisar_orfs_e_proteinas(sequencias, tamanho_min_orf_aa):
    """Encontra ORFs e analisa a proteína mais longa de cada ORF encontrado."""
    resultados = {}
    for seq_record in sequencias:
        orfs_encontrados = []
        # Testa os 6 frames de leitura (3 foward, 3 reverse-complement)
        for strand, nuc_seq in [ (1, seq_record.seq), (-1, seq_record.seq.reverse_complement()) ]:
            for frame in range(3):
                seq_traduzida = nuc_seq[frame:].translate(to_stop=True) # Para no primeiro Stop Codon
                if len(seq_traduzida) >= tamanho_min_orf_aa:
                    try:
                        analise_proteina = ProteinAnalysis(str(seq_traduzida))
                        orf_info = {
                            'proteina': str(seq_traduzida),
                            'tamanho_aa': len(seq_traduzida),
                            'frame': frame + 1,
                            'direcao': strand,
                            'peso_molecular_Da': f"{analise_proteina.molecular_weight():.2f}",
                            'ponto_isoeletrico': f"{analise_proteina.isoelectric_point():.2f}",
                            'gravy': f"{analise_proteina.gravy():.2f}"
                        }
                        orfs_encontrados.append(orf_info)
                    except Exception:
                        continue # Ignora sequências que causam erro na análise de proteína
        if orfs_encontrados:
            # Ordena os orfs pelo tamanho para facilmente encontrar o maior
            orfs_encontrados.sort(key=lambda x: x['tamanho_aa'], reverse=True)
            resultados[seq_record.id] = orfs_encontrados
    return resultados

def analisar_restricao(sequencias, enzimas_str):
    """Analisa os sítios de corte para as enzimas de restrição fornecidas."""
    try:
        lista_enzimas = RestrictionBatch(enzimas_str.split(','))
    except ValueError as e:
        return {"erro": f"Nome de enzima inválido: {e}"}

    resultados = {}
    for seq_record in sequencias:
        analise = RestrictionAnalysis(lista_enzimas, seq_record.seq)
        resultados_enzimas = analise.run()
        # Filtra para mostrar apenas enzimas que realmente cortam a sequência
        cortes_reais = {str(k): v for k, v in resultados_enzimas.items() if len(v) > 0}
        if cortes_reais:
            resultados[seq_record.id] = cortes_reais
    return resultados

# =============================================================================
# FUNÇÃO PRINCIPAL DE ORQUESTRAÇÃO
# =============================================================================

def executar_analises(args):
    """Função central que lê o arquivo e chama as análises selecionadas."""
    try:
        sequencias = list(SeqIO.parse(args.arquivo, "fasta"))
    except FileNotFoundError:
        print(f"Erro: O arquivo '{args.arquivo}' não foi encontrado.")
        return
    if not sequencias:
        print("Nenhuma sequência encontrada no arquivo.")
        return

    # --- Coleta de Resultados ---
    resultados_finais = {'arquivo': args.arquivo, 'total_sequencias': len(sequencias)}

    if args.fator_escolhido:
        resultados_finais['busca_tf'] = buscar_sitios_transcricao(sequencias, args.fator_escolhido)
    
    if args.composicao:
        resultados_finais['composicao'] = analisar_composicao(sequencias)

    if args.orfs:
        resultados_finais['orfs'] = analisar_orfs_e_proteinas(sequencias, args.orfs)

    if args.restricao:
        resultados_finais['restricao'] = analisar_restricao(sequencias, args.restricao)
    
    imprimir_resultados(resultados_finais)

# =============================================================================
# FUNÇÃO DE IMPRESSÃO DOS RESULTADOS
# =============================================================================

def imprimir_resultados(resultados):
    """Imprime todos os resultados coletados de forma organizada."""
    print("\n" + "="*80)
    print(f"📄 RELATÓRIO DE ANÁLISE DO ARQUIVO: {resultados['arquivo']}")
    print(f"🧬 TOTAL DE SEQUÊNCIAS: {resultados['total_sequencias']}")
    print("="*80)

    if 'composicao' in resultados:
        print("\n--- 🔬 Análise de Composição (GC Content) ---")
        print(f"Conteúdo GC Geral do Arquivo: {resultados['composicao']['geral']}")
        for seq_id, gc in resultados['composicao']['por_sequencia'].items():
            print(f"  -> {seq_id}: {gc}")
    
    if 'busca_tf' in resultados:
        busca = resultados['busca_tf']
        print(f"\n--- 🔎 Resultados da Busca pelo Sítio '{busca['fator']}' ---")
        if 'erro' in busca:
            print(f"  Erro: {busca['erro']}")
        elif not busca['encontrados']:
            print("  Nenhuma ocorrência encontrada.")
        else:
            for seq_id, posicoes in busca['encontrados'].items():
                pos_str = ", ".join([f"[{p[0]}-{p[1]}]" for p in posicoes])
                print(f"  -> Em '{seq_id}': encontrado na(s) posição(ões) {pos_str}")

    if 'orfs' in resultados:
        print(f"\n--- 🧬 Análise de ORFs (Tamanho Mínimo: {args.orfs} aa) ---")
        if not resultados['orfs']:
            print("  Nenhum ORF correspondente aos critérios foi encontrado.")
        else:
            for seq_id, orfs in resultados['orfs'].items():
                print(f"\n  -> Sequência ID: {seq_id}")
                for i, orf in enumerate(orfs):
                    print(f"    - ORF #{i+1} (Maior Primeiro):")
                    print(f"      Tamanho: {orf['tamanho_aa']} aa | Direção: {'Forward' if orf['direcao']==1 else 'Reverse'} | Frame: {orf['frame']}")
                    print(f"      Peso Mol.: {orf['peso_molecular_Da']} Da | pI: {orf['ponto_isoeletrico']} | GRAVY: {orf['gravy']}")
    
    if 'restricao' in resultados:
        print(f"\n--- ✂️ Análise de Sítios de Restrição ---")
        if 'erro' in resultados['restricao']:
            print(f"  Erro: {resultados['restricao']['erro']}")
        elif not resultados['restricao']:
            print("  Nenhuma das enzimas especificadas corta as sequências.")
        else:
            for seq_id, enzimas in resultados['restricao'].items():
                print(f"\n  -> Sequência ID: {seq_id}")
                for enzima, posicoes in enzimas.items():
                    print(f"    - Enzima '{enzima}' corta {len(posicoes)} vez(es) na(s) posição(ões): {posicoes}")
    
    print("\n" + "="*80 + "\n")

# =============================================================================
# CONFIGURAÇÃO DA INTERFACE DE LINHA DE COMANDO E EXECUÇÃO
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Ferramenta avançada para análise de arquivos FASTA.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Argumento obrigatório
    parser.add_argument("arquivo", help="Caminho para o arquivo FASTA de entrada.")
    
    # Grupo de flags para busca de Fatores de Transcrição (mutuamente exclusivo)
    grupo_tf = parser.add_mutually_exclusive_group()
    for fator, regex in SITIOS_REGEX.items():
        grupo_tf.add_argument(f"--{fator.lower()}", dest="fator_escolhido", action="store_const", const=fator, help=f"Busca pelo sítio {fator} (Regex: {regex})")

    # Grupo de flags para análises adicionais (não são exclusivas)
    grupo_analises = parser.add_argument_group('Análises Adicionais')
    grupo_analises.add_argument("--composicao", action="store_true", help="Calcula o conteúdo GC de cada sequência e o geral.")
    grupo_analises.add_argument("--orfs", type=int, metavar="TAMANHO_MIN_AA", help="Encontra ORFs com um tamanho mínimo de aminoácidos.")
    grupo_analises.add_argument("--restricao", type=str, metavar="ENZ1,ENZ2", help="Analisa sítios de restrição. Forneça enzimas separadas por vírgula (ex: 'EcoRI,BamHI').")
    
    args = parser.parse_args()
    
    executar_analises(args)