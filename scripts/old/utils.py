def hamming_distance  (string_1, string_2):
  if len(string_1) > len(string_2):
    for i in range(len(string_1) - len(string_2)):
      string_2 += ' '
  if len(string_2) > len(string_1):
    for i in range(len(string_2) - len(string_1)):
      string_1 += ' '
  diferencas = sum(x != y for x, y in zip(string_1, string_2))
  return diferencas / len(string_1)
#   diferencas = 0
#   for i in range(len(string_1)):
#     if (string_1[i] != string_2[i]):
#       diferencas += 1
#   return diferencas

# Exemplo
# distancia = hamming_distance('CTAG', 'CTCA')
# print(f"Distancia de Hamming: {distancia:.4f}")


def calcular_tm(seq):
    # Contagem de bases
    A = seq.count("A")
    T = seq.count("T")
    C = seq.count("C")
    G = seq.count("G")

    # Fórmula de Wallace: Tm = 2*(A+T) + 4*(C+G)
    tm = 2 * (A + T) + 4 * (C + G)
    return tm

# sequencia = "ACGTACGTACGT"
# tm_calculado = calcular_tm(sequencia)
# print(f"A temperatura de melting da sequência {sequencia} é aproximadamente {tm_calculado}°C")


def dna_to_one_hot_flat(sequence):
    mapping = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}
    encoded = [mapping[nucleotide] for nucleotide in sequence]
    # Achatar a lista de listas para um único vetor
    flattened = [item for sublist in encoded for item in sublist]
    return flattened

def dna_to_one_hot_flat_df(sequence):
    import pandas as pd
    flattened = dna_to_one_hot_flat(sequence)
    # Criar as colunas com base nos nucleotídeos e posições
    nucleotides = ['A', 'T', 'C', 'G']
    columns = [f"{nucleotide}_{i//4 + 1}" for i, nucleotide in enumerate(nucleotides * (len(flattened) // 4))]

    # Criar DataFrame com a linha única representando todo o registro
    df_encoded = pd.DataFrame([flattened], columns=columns)
    return df_encoded

# def one_hot_flat_to_df2(sequence_flat):
#     # flattened = dna_to_one_hot_flat(sequence_flat)
#     # Criar as colunas com base nos nucleotídeos e posições
#     nucleotides = ['A', 'T', 'C', 'G']
#     columns = [f"{nucleotide}_{i//4 + 1}" for i, nucleotide in enumerate(nucleotides * (len(sequence_flat) // 4))]

#     # Criar DataFrame com a linha única representando todo o registro
#     df_encoded = pd.DataFrame([sequence_flat], columns=columns)
#     return df_encoded

# def encodeToArrayOfChar(sequence):
#     return  np.array(list(sequence))
#     # return  np.array(nucleotide for nucleotide in sequence)

# def encodToArrayOfBol(sequence):
#     mapping = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1],
#                'a': [1, 0, 0, 0], 't': [0, 1, 0, 0], 'c': [0, 0, 1, 0], 'g': [0, 0, 0, 1]}
#     return np.array([mapping[nucleotide] for nucleotide in sequence])

# def encodToSmallArrayOfBol(sequence):
#     mapping = {'A': [0, 0], 'T': [0, 1], 'C': [1, 0], 'G': [1, 1],
#                'a': [0, 0], 't': [0, 1], 'c': [1, 0], 'g': [1, 1]}
#     return np.array([mapping[nucleotide] for nucleotide in sequence])

# def encodeToStringOfBol(sequence):
#     mapping = {'A': '1000', 'T': '0100', 'C': '0010', 'G': '0001',
#                'a': '1000', 't': '0100', 'c': '0010', 'g': '0001'}
#     # return  ''.join(map(str, np.array([mapping[nucleotide] for nucleotide in sequence])))
#     return  np.array([mapping[nucleotide] for nucleotide in sequence])

# def encodeToSmallStringOfBol(sequence):
#     mapping = {'A': '00', 'T': '01', 'C': '10', 'G': '11',
#                'a': '00', 't': '01', 'c': '10', 'g': '11'}
#     return  np.array([mapping[nucleotide] for nucleotide in sequence])


def flank_str_amostra(dna_sequence, inicio, fim, tamanho):
    flank_length = min(tamanho, (fim - inicio))
    # Calcular as regiões flanqueadoras
    fragment = dna_sequence[inicio:fim]
    forward_primer = dna_sequence[max(0, inicio - flank_length):inicio]
    reverse_primer = dna_sequence[inicio + (fim - inicio):inicio + (fim - inicio) + flank_length]

    return {
        "fragment": fragment,
        "forward_primer": forward_primer,
        "reverse_primer": reverse_primer,
        "full_primer": forward_primer + fragment + reverse_primer,
    }

def create_primers(dna_sequence, fragment, flank_length=20):
    # Encontrar a posição do fragmento na sequência de DNA
    start_index = dna_sequence.find(fragment)
    if start_index == -1:
        raise ValueError("O fragmento não foi encontrado na sequência de DNA.")

    # Calcular as regiões flanqueadoras
    forward_primer = dna_sequence[max(0, start_index - flank_length):start_index]
    reverse_primer = dna_sequence[start_index + len(fragment):start_index + len(fragment) + flank_length]

    print('start_index: ', start_index)
    return {
        "forward_primer": forward_primer,
        "reverse_primer": reverse_primer
    }

# # Exemplo de uso
# dna_sequence = "ATCGTACCGTACGATCGTACCGGATCGTAGCGTACGATCG"
# fragment = "TACCGGATCG"
# primers = create_primers(dna_sequence, fragment)

# # Exibir os primers
# print("Forward Primer:", len(primers["forward_primer"]))
# print("Reverse Primer:", len(primers["reverse_primer"]))



def obtem_dna_complementar(dna_sequence):
  from Bio.Seq import Seq
  seq = Seq(str(dna_sequence))
  return str(seq.complement())



def calcular_percentual_bases(amostra):
  from Bio.Seq import Seq
  dna_sequence = Seq(amostra)
  total_bases = len(amostra)

  percentual = {
      "A": 100 * Seq.count(dna_sequence, 'A') / total_bases,
      "T": 100 * Seq.count(dna_sequence, 'T') / total_bases,
      "C": 100 * Seq.count(dna_sequence, 'C') / total_bases,
      "G": 100 * Seq.count(dna_sequence, 'G') / total_bases
  }
  return percentual



# cria_pasta_analise()
# Criar uma nova pasta: pasta_analise
def cria_pasta_analise(pasta_saida):
    import pandas as pd
    import os
    import datetime
    if len(pasta_analise) <= 0:
        data_hora = datetime.datetime.now()
        data_hora = data_hora.strftime('%Y-%m-%d_%I-%M-%S_%M%S')
        # data_hora = datetime.datetime.now().strftime('%Y-%m-%d_%I-%M-%S_%M%S')

        pasta_analise = pasta_saida + data_hora + '/'

        # Criar a pasta se não existir
        if not os.path.exists(pasta_analise):
            os.makedirs(pasta_analise)
            print(f"Pasta criada: {pasta_analise}")
        else:
            print(f"A pasta '{pasta_analise}' já existe!")




# busca_arquivos
def busca_arquivos(tipo, diretorio_raiz):
    from pathlib import Path
    origem = Path(diretorio_raiz).expanduser()
    arquivos = []
    for arquivo in origem.rglob(tipo):
        arquivos.append(arquivo)
    return arquivos 

def gravar_arquivo_saida(nome_amostra, tamanho_str, tipo, dados, pasta_analise):
  import os
  url_arquivo = pasta_analise + nome_amostra + '_' + tipo+ '_' + str(tamanho_str) + '.txt'

  if os.path.exists(url_arquivo):
    os.remove(url_arquivo)

  dados.to_csv(url_arquivo, sep=';', index=False)

  print(url_arquivo)
  return url_arquivo



# buscar aquivos na pasta de analises. tipos de arquivos: STR
def busca_arquivos(tipo, pasta_analise):
  import os
  arquivos_localizados = []
  # Percorrendo todos os diretórios e subdiretórios
  for raiz, _, arquivos in os.walk(pasta_analise):
      for nome_arquivo in arquivos:
          # Aplicar filtros
          if nome_arquivo.endswith(".txt") and tipo in nome_arquivo:
              caminho_arquivo = os.path.join(raiz, nome_arquivo)  # Caminho completo do arquivo
              # nome_pasta = os.path.basename(raiz)[0:5]  # Nome da pasta onde o arquivo está
              arquivos_localizados.append(caminho_arquivo)

  return sorted(arquivos_localizados)

# busca_arquivos('STR')



def quebra_nome_arquivo(arquivo):
    import os
    nome_arquivo = os.path.basename(arquivo) # Obter apenas o nome do arquivo sem o diretório
    parte_principal = os.path.splitext(nome_arquivo)[0] # Remover a extensão do arquivo
    return parte_principal.split('_')






# le arquivos de amostras 
def ler_amostras(arquivos_amostras):
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    amostras_originais = []
    for idx_amostras, arquivo in enumerate(arquivos_amostras):
        print(f'lendo {idx_amostras +1} de {len(arquivos_amostras)} amostra(s): {arquivo}')
        sequencias_arquivo_sfp = ''
        with open(arquivo) as handle:
            for seq in SimpleFastaParser(handle):
                sequencias_arquivo_sfp += str(seq).upper()
        amostras_originais.append(sequencias_arquivo_sfp)
    print('     ...fim')
    return amostras_originais








# testando STRs em tandem
def find_str(sequencia_dna, tamanho_str):
  lista = []
  fim_da_fila = 0
  str_localizado = False

  for i in range(len(sequencia_dna)):
    if str_localizado and fim_da_fila > 0 and i < fim_da_fila:
      continue

    amostra_str = ''
    fim_da_fila = 0
    if i < len(sequencia_dna)-tamanho_str+1:
      inicio_amostra_str = i
      fim_amostra_str = inicio_amostra_str+tamanho_str
      amostra_str = sequencia_dna[inicio_amostra_str:fim_amostra_str]

      quantidade_sequencias_repetidas = 0
      is_testar_proxima_Sequencia = True
      while is_testar_proxima_Sequencia:
        inicio_proxima_sequencia = i+(tamanho_str*(quantidade_sequencias_repetidas+1))
        fim_proxima_sequencia = inicio_proxima_sequencia+tamanho_str
        proxima_sequencia = ''
        if fim_proxima_sequencia <= len(sequencia_dna):
          proxima_sequencia = sequencia_dna[inicio_proxima_sequencia:fim_proxima_sequencia]

        if amostra_str == proxima_sequencia:
          quantidade_sequencias_repetidas = quantidade_sequencias_repetidas+1
          fim_da_fila = fim_proxima_sequencia
        else:
          is_testar_proxima_Sequencia = False
      # fim while

      if quantidade_sequencias_repetidas > 0:
        str_localizado = True
        lista.append({
              'posicao': i+1,
              'unidade': amostra_str,
              'inicio_loci': inicio_amostra_str,
              'fim_unidade': fim_amostra_str,
              'fim_loci': fim_da_fila,
              'copias': quantidade_sequencias_repetidas +1
        })
      else:
        str_localizado = False
  # fim for
  df = pd.DataFrame.from_records(lista,index=['posicao'])
  return df
