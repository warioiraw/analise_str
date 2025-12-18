import utils as ut

analises = [1, 2, 3, 4, 5, 6, 7, 8, 9]
pasta_saida = '/Users/wariol/Library/Mobile Documents/com~apple~CloudDocs/UFU/pesquisa/py/saida/'
pasta_amostras = '/Users/wariol/Library/Mobile Documents/com~apple~CloudDocs/UFU/pesquisa/py/data_raw/ncbi/'
pasta_analise = '/Users/wariol/Library/Mobile Documents/com~apple~CloudDocs/UFU/pesquisa/py/saida/2025-05-11_04-50-31_5031/'
labels_amostras = ['ASM1928827v1', 'ASM2978390v1', 'ASM2978392v1']
arquivos_amostras = [
    pasta_amostras + 'ASM1928827v1/ncbi_dataset/data/GCA_019288275.1/GCA_019288275.1_ASM1928827v1_genomic.fna',
    pasta_amostras + 'ASM2978390v1/ncbi_dataset/data/GCA_029783905.1/GCA_029783905.1_ASM2978390v1_genomic.fna',
    pasta_amostras + 'ASM2978392v1/ncbi_dataset/data/GCA_029783925.1/GCA_029783925.1_ASM2978392v1_genomic.fna'
]

amostras = ut.ler_amostras(arquivos_amostras)

distancia_ab = ut.hamming_distance(arquivos_amostras[0], arquivos_amostras[1])
print(f'distancia_ab: {distancia_ab}')

distancia_bc = ut.hamming_distance(arquivos_amostras[1], arquivos_amostras[2])
print(f'distancia_bc: {distancia_bc}')

distancia_ac = ut.hamming_distance(arquivos_amostras[0], arquivos_amostras[2])
print(f'distancia_ac: {distancia_ac}')


# testes = ut.hamming_distance('ACTG', 'ACTG')
# print(f'testes: {testes * 100}')