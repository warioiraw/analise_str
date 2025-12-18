import pandas as pd
import subprocess
import os
import tempfile


# dados = [
#     ["ASM285v2",     "ASM285v2 *",   "GCA_000002855.2", "PRJNA19275", "CBS 513.88",   "DSM, The Netherlands",    "Scaffold",   34, './Users/wariol/projetos/blast_local/ASM285v2/ncbi_dataset/data/GCA_000002855.2/GCA_000002855.2_ASM285v2_genomic.fna'],
#     ["ASM1928827v1", "ASM1928827v1", "GCA_019288275.1", "JAGRPH01",   "CBS 554.65",   "TU Wien",                 "Chromosome", 40, './Users/wariol/projetos/blast_local/ASM1928827v1/ncbi_dataset/data/GCA_019288275.1/GCA_019288275.1_ASM1928827v1_genomic.fna'],
#     ["ASM2978390v1", "ASM2978390v1", "GCA_029783905.1", "JAPVRD01",   "KJC3",         "Soongsil University",     "Chromosome", 40, './Users/wariol/projetos/blast_local/ASM2978390v1/ncbi_dataset/data/GCA_029783905.1/GCA_029783905.1_ASM2978390v1_genomic.fna'],
#     ["ASM2978392v1", "ASM2978392v1", "GCA_029783925.1", "JAPVRE01",   "KYF3",         "Soongsil University",     "Chromosome", 37, './Users/wariol/projetos/blast_local/ASM2978392v1/ncbi_dataset/data/GCA_029783925.1/GCA_029783925.1_ASM2978392v1_genomic.fna'],
#     ["ASM4765177v1", "ASM4765177v1", "GCA_047651775.1", "JBKZXA01",   "CCTCC 206047", "Zhejiang Uni.Technology", "Complete",   35, './Users/wariol/projetos/blast_local/ASM4765177v1/ncbi_dataset/data/GCA_047651775.1/GCA_047651775.1_ASM4765177v1_genomic.fna']
# ]
# colunas = [
#     "id",
#     "Nome da Montagem",
#     "Referência GenBank",
#     "WGS Accession",
#     "Cepa (strain)",
#     "Submissor",
#     "Nível de montagem",
#     "Tamanho (Mb)",
#     "Dados"
# ]
# df_amostras = pd.DataFrame(dados, columns=colunas)
# amostra = df_amostras.query("id == 'ASM285v2'")
# # print(amostra)
# arquivo = amostra.get('Dados').iloc[0]
# print(f"arquivo_amostra: '{arquivo}'")

# # 1 gerar banco de dados da amostra de referencia
# # arquivos de entrada e saída
# diretorio = os.path.dirname(arquivo)
# banco_dados = os.path.join(diretorio, "ASM285v2_db")
# print(f">>>>>> banco_dados: {banco_dados}")

# # >>> comando como lista (boa prática, evita problemas de espaços)
# comando = [
#     "makeblastdb",
#     "-in", arquivo,
#     "-dbtype", "nucl",
#     "-out", banco_dados
# ]

# # executa o comando
# subprocess.run(comando, check=True)

# print(f"Banco BLAST criado: {banco_dados}")

# 2 rodar o blastn local
# Arquivos
query_seq = "TCCCTCACAGACACAGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCCCTCTCCCTCTTGCTATCTTTTTTTACATATAT"
db_file = "/Users/wariol/projetos/blast_local/ASM285v2/ncbi_dataset/data/GCA_000002855.2/ASM285v2_db"
out_file = "/Users/wariol/projetos/blast_local/blast_results.tsv"
query_file = ""

# Criar um arquivo FASTA temporário para a query
with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_fasta:
    tmp_fasta.write(">query_sequence\n")
    tmp_fasta.write(query_seq + "\n")
    query_file = tmp_fasta.name  # caminho do arquivo temporário
# 
subprocess.run([
    "blastn",
    "-query", query_file,
    "-db", db_file,
    "-out", out_file,
    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    # "-word_size", "7",
    "-dust", "no"
])

# Carregar resultado no pandas
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_csv(out_file, sep="\t", names=cols)
print(df.head())