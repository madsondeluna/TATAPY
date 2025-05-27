# Importa as bibliotecas necessárias
import argparse
import re # Biblioteca para Expressões Regulares (Regex)
from Bio import SeqIO
# A linha abaixo foi removida pois não estava sendo usada na sua versão final.
# from Bio.SeqUtils import GC 

# =============================================================
# SEÇÃO 1: ADICIONE OS NOVOS FATORES AQUI
# =============================================================
# Dicionário que mapeia o nome do fator ao seu padrão de Regex
SITIOS_REGEX = {
    "SP1": r'GGGCGG',
    "TATA": r'TATA[AT]A[AT]',
    "CREB": r'TGACGTCA',
    "NFKB": r'GGGA[CT]TTCC',      # Sítio para NF-κB (Nuclear Factor kappa B)
    "AP1": r'TGA[CG]TCA',         # Sítio para AP-1 (Activator Protein 1)
    "CAAT": r'GGCCAATCT',         # CAAT Box
    "EBOX": r'CA[ACGT]{2}TG'      # E-Box (CANNTG, onde N é qualquer base)
}

def buscar_sitios_regex(sequencia, padrao_regex):
    """
    Busca todas as ocorrências de um padrão Regex em uma sequência.

    Retorna:
        list: Uma lista de tuplas com a posição (início, fim) de cada match.
    """
    posicoes = []
    # Usamos re.finditer para encontrar todas as ocorrências não sobrepostas
    for match in re.finditer(padrao_regex, str(sequencia), re.IGNORECASE):
        # re.IGNORECASE faz a busca ser case-insensitive (ignora se é 'A' ou 'a')
        posicoes.append((match.start(), match.end()))
    return posicoes

def analisar_fasta(caminho_arquivo, fator_pesquisa=None):
    """
    Função principal para analisar um arquivo FASTA.

    Argumentos:
        caminho_arquivo (str): O caminho para o arquivo FASTA.
        fator_pesquisa (str, opcional): O nome do fator a ser pesquisado.
    """
    try:
        sequencias = list(SeqIO.parse(caminho_arquivo, "fasta"))
    except FileNotFoundError:
        print(f"Erro: O arquivo '{caminho_arquivo}' não foi encontrado.")
        return None
    except Exception as e:
        print(f"Ocorreu um erro ao ler o arquivo: {e}")
        return None

    if not sequencias:
        print("Nenhuma sequência encontrada no arquivo.")
        return None

    # --- Análise Básica (como antes) ---
    comprimentos = [len(seq.seq) for seq in sequencias]
    
    resultados = {
        "total_sequencias": len(sequencias),
        "comprimento_maximo": max(comprimentos),
        "comprimento_minimo": min(comprimentos),
        "comprimento_medio": f"{sum(comprimentos) / len(sequencias):.2f}",
        "resultados_busca": None # Inicializa os resultados da busca
    }

    # --- Análise de Sítios de Reconhecimento (nova parte) ---
    if fator_pesquisa:
        padrao_regex = SITIOS_REGEX.get(fator_pesquisa)
        if not padrao_regex:
            print(f"Erro: Fator '{fator_pesquisa}' desconhecido.")
            return resultados

        resultados_busca = {"fator": fator_pesquisa, "encontrados": {}}
        for seq_record in sequencias:
            posicoes = buscar_sitios_regex(seq_record.seq, padrao_regex)
            if posicoes: # Se encontrou alguma posição
                resultados_busca["encontrados"][seq_record.id] = posicoes
        
        resultados["resultados_busca"] = resultados_busca

    return resultados

def imprimir_resultados(resultados):
    """
    Imprime os resultados da análise de forma organizada no terminal.
    """
    if resultados is None:
        return

    print("\n--- Resultados da Análise do Arquivo FASTA ---")
    print(f"🧬 Total de Sequências: {resultados['total_sequencias']}")
    print(f"📏 Comprimento Máximo: {resultados['comprimento_maximo']} bp")
    print(f"📏 Comprimento Mínimo: {resultados['comprimento_minimo']} bp")
    print(f"📏 Comprimento Médio: {resultados['comprimento_medio']} bp")
    
    # Imprime os resultados da busca por sítios, se houver
    if resultados.get("resultados_busca"):
        busca = resultados["resultados_busca"]
        print(f"\n--- 🔎 Resultados da Busca pelo Sítio '{busca['fator']}' ---")
        if not busca['encontrados']:
            print("Nenhuma ocorrência encontrada nas sequências.")
        else:
            for seq_id, posicoes in busca['encontrados'].items():
                pos_str = ", ".join([f"[{p[0]}-{p[1]}]" for p in posicoes])
                print(f"-> Em '{seq_id}': encontrado na(s) posição(ões) {pos_str}")

    print("\n---------------------------------------------\n")

def main():
    """
    Função principal que configura o 'argparse' e executa o programa.
    """
    parser = argparse.ArgumentParser(
        description="Analisa um arquivo FASTA e opcionalmente busca por sítios de reconhecimento de fatores de transcrição.",
        formatter_class=argparse.RawTextHelpFormatter # Melhora a formatação da ajuda
    )
    
    parser.add_argument("arquivo", help="O caminho para o arquivo FASTA de entrada.")
    
    # Cria um grupo de argumentos mutuamente exclusivos.
    # O usuário só pode escolher UMA das flags abaixo.
    grupo_fatores = parser.add_mutually_exclusive_group()
    
    # Flags já existentes
    grupo_fatores.add_argument("--sp1", action="store_const", const="SP1", dest="fator_escolhido", help=f"Busca pelo sítio de reconhecimento SP1 (Regex: {SITIOS_REGEX['SP1']})")
    grupo_fatores.add_argument("--tata", action="store_const", const="TATA", dest="fator_escolhido", help=f"Busca pela TATA box (Regex: {SITIOS_REGEX['TATA']})")
    grupo_fatores.add_argument("--creb", action="store_const", const="CREB", dest="fator_escolhido", help=f"Busca pelo sítio de reconhecimento CREB (Regex: {SITIOS_REGEX['CREB']})")

    # =============================================================
    # SEÇÃO 2: ADICIONE AS NOVAS FLAGS AQUI
    # =============================================================
    grupo_fatores.add_argument("--nfkb", action="store_const", const="NFKB", dest="fator_escolhido", help=f"Busca pelo sítio de NFKB (Regex: {SITIOS_REGEX['NFKB']})")
    grupo_fatores.add_argument("--ap1", action="store_const", const="AP1", dest="fator_escolhido", help=f"Busca pelo sítio de AP-1 (Regex: {SITIOS_REGEX['AP1']})")
    grupo_fatores.add_argument("--caat", action="store_const", const="CAAT", dest="fator_escolhido", help=f"Busca pela CAAT Box (Regex: {SITIOS_REGEX['CAAT']})")
    grupo_fatores.add_argument("--ebox", action="store_const", const="EBOX", dest="fator_escolhido", help=f"Busca pela E-Box (Regex: {SITIOS_REGEX['EBOX']})")

    args = parser.parse_args()
    
    resultados_analise = analisar_fasta(args.arquivo, args.fator_escolhido)
    imprimir_resultados(resultados_analise)

if __name__ == "__main__":
    main()