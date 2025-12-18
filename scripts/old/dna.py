import sys
# import sklearn
import numpy as np
import csv
import pandas as pd
# import numpy as npc

# import tensorflow as tf

def one_hot_encode(sequence):
    mapping = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1],
               'a': [1, 0, 0, 0], 't': [0, 1, 0, 0], 'c': [0, 0, 1, 0], 'g': [0, 0, 0, 1]}
    return np.array([mapping[nucleotide] for nucleotide in sequence])


# import seaborn as sbn
# import scipy as spc

# https://www.kaggle.com/code/bulentsiyah/classifying-dna-sequences-markov-models-knn-svm
# https://github.com/ajitsingh98/DNA-Classification-Machine-Learning-Projecthttps://github.com/ajitsingh98/DNA-Classification-Machine-Learning-Project
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8285202/
# https://archive.ics.uci.edu/ml/machine-learning-databases/molecular-biology/promoter-gene-sequences/promoters.data
# https://www.hashtagtreinamentos.com/criar-dataframes-com-pandas?gad_source=1&gclid=CjwKCAjw5qC2BhB8EiwAvqa41qXmklIgxL-pb_6Qe4R5HXqp8igMbwg8b9ZTwc3L6jihsy7F-B4WQRoCibwQAvD_BwE
# https://adityamahes.medium.com/ml-applications-for-data-mining-gene-sequences-76fa2f70c82d

# how can machine learning car be used to find patters in dna sequencies

dados = [
    {'class': '+', 'id': 'S10', 'sequencia': 'tactagcaatacgcttgcgttcggtggttaagtatgtataatgcgcgggcttgtcgt'},
    {'class': '+', 'id': 'AMPC', 'sequencia': 'tgctatcctgacagttgtcacgctgattggtgtcgttacaatctaacgcatcgccaa'},
    {'class': '+', 'id': 'AROH', 'sequencia': 'gtactagagaactagtgcattagcttatttttttgttatcatgctaaccacccggcg'},
    {'class': '+', 'id': 'DEOP2', 'sequencia': 'aattgtgatgtgtatcgaagtgtgttgcggagtagatgttagaatactaacaaactc'},
    {'class': '+', 'id': 'LEU1_TRNA', 'sequencia': 'tcgataattaactattgacgaaaagctgaaaaccactagaatgcgcctccgtggtag'},
    {'class': '+', 'id': 'MALEFG', 'sequencia': 'aggggcaaggaggatggaaagaggttgccgtataaagaaactagagtccgtttaggt'},
    {'class': '+', 'id': 'MALK', 'sequencia': 'cagggggtggaggatttaagccatctcctgatgacgcatagtcagcccatcatgaat'},
    {'class': '+', 'id': 'RECA', 'sequencia': 'tttctacaaaacacttgatactgtatgagcatacagtataattgcttcaacagaaca'},
    {'class': '+', 'id': 'RPOB', 'sequencia': 'cgacttaatatactgcgacaggacgtccgttctgtgtaaatcgcaatgaaatggttt'},
    {'class': '+', 'id': 'RRNAB_P1', 'sequencia': 'ttttaaatttcctcttgtcaggccggaataactccctataatgcgccaccactgaca'},
    {'class': '+', 'id': 'RRNAB_P2', 'sequencia': 'gcaaaaataaatgcttgactctgtagcgggaaggcgtattatgcacaccccgcgccg'},
    {'class': '+', 'id': 'RRNDEX_P2', 'sequencia': 'cctgaaattcagggttgactctgaaagaggaaagcgtaatatacgccacctcgcgac'},
    {'class': '+', 'id': 'RRND_P1', 'sequencia': 'gatcaaaaaaatacttgtgcaaaaaattgggatccctataatgcgcctccgttgaga'},
    {'class': '+', 'id': 'RRNE_P1', 'sequencia': 'ctgcaatttttctattgcggcctgcggagaactccctataatgcgcctccatcgaca'},
    {'class': '+', 'id': 'RRNG_P1', 'sequencia': 'tttatatttttcgcttgtcaggccggaataactccctataatgcgccaccactgaca'},
    {'class': '+', 'id': 'RRNG_P2', 'sequencia': 'aagcaaagaaatgcttgactctgtagcgggaaggcgtattatgcacaccgccgcgcc'},
    {'class': '+', 'id': 'RRNX_P1', 'sequencia': 'atgcatttttccgcttgtcttcctgagccgactccctataatgcgcctccatcgaca'},
    {'class': '+', 'id': 'TNAA', 'sequencia': 'aaacaatttcagaatagacaaaaactctgagtgtaataatgtagcctcgtgtcttgc'},
    {'class': '+', 'id': 'TYRT', 'sequencia': 'tctcaacgtaacactttacagcggcgcgtcatttgatatgatgcgccccgcttcccg'},
    {'class': '+', 'id': 'ARAC', 'sequencia': 'gcaaataatcaatgtggacttttctgccgtgattatagacacttttgttacgcgttt'},
    {'class': '+', 'id': 'LACI', 'sequencia': 'gacaccatcgaatggcgcaaaacctttcgcggtatggcatgatagcgcccggaagag'},
    {'class': '+', 'id': 'MALT', 'sequencia': 'aaaaacgtcatcgcttgcattagaaaggtttctggccgaccttataaccattaatta'},
    {'class': '+', 'id': 'TRP', 'sequencia': 'tctgaaatgagctgttgacaattaatcatcgaactagttaactagtacgcaagttca'},
    {'class': '+', 'id': 'TRPP2', 'sequencia': 'accggaagaaaaccgtgacattttaacacgtttgttacaaggtaaaggcgacgccgc'},
    {'class': '+', 'id': 'THR', 'sequencia': 'aaattaaaattttattgacttaggtcactaaatactttaaccaatataggcatagcg'},
    {'class': '+', 'id': 'BIOB', 'sequencia': 'ttgtcataatcgacttgtaaaccaaattgaaaagatttaggtttacaagtctacacc'},
    {'class': '+', 'id': 'FOL', 'sequencia': 'catcctcgcaccagtcgacgacggtttacgctttacgtatagtggcgacaatttttt'},
    {'class': '+', 'id': 'UVRBP1', 'sequencia': 'tccagtataatttgttggcataattaagtacgacgagtaaaattacatacctgcccg'},
    {'class': '+', 'id': 'UVRBP3', 'sequencia': 'acagttatccactattcctgtggataaccatgtgtattagagttagaaaacacgagg'},
    {'class': '+', 'id': 'LEXA', 'sequencia': 'tgtgcagtttatggttccaaaatcgccttttgctgtatatactcacagcataactgt'},
    {'class': '+', 'id': 'PORI-L', 'sequencia': 'ctgttgttcagtttttgagttgtgtataacccctcattctgatcccagcttatacgg'},
    {'class': '+', 'id': 'SPOT42', 'sequencia': 'attacaaaaagtgctttctgaactgaacaaaaaagagtaaagttagtcgcgtagggt'},
    {'class': '+', 'id': 'M1RNA', 'sequencia': 'atgcgcaacgcggggtgacaagggcgcgcaaaccctctatactgcgcgccgaagctg'},
    {'class': '+', 'id': 'GLNS', 'sequencia': 'taaaaaactaacagttgtcagcctgtcccgcttataagatcatacgccgttatacgt'},
    {'class': '+', 'id': 'TUFB', 'sequencia': 'atgcaattttttagttgcatgaactcgcatgtctccatagaatgcgcgctacttgat'},
    {'class': '+', 'id': 'SUBB-E', 'sequencia': 'ccttgaaaaagaggttgacgctgcaaggctctatacgcataatgcgccccgcaacgc'},
    {'class': '+', 'id': 'STR', 'sequencia': 'tcgttgtatatttcttgacaccttttcggcatcgccctaaaattcggcgtcctcata'},
    {'class': '+', 'id': 'SPC', 'sequencia': 'ccgtttattttttctacccatatccttgaagcggtgttataatgccgcgccctcgat'},
    {'class': '+', 'id': 'RPOA', 'sequencia': 'ttcgcatatttttcttgcaaagttgggttgagctggctagattagccagccaatctt'},
    {'class': '+', 'id': 'RPLJ', 'sequencia': 'tgtaaactaatgcctttacgtgggcggtgattttgtctacaatcttacccccacgta'},
    {'class': '+', 'id': 'PORI-R', 'sequencia': 'gatcgcacgatctgtatacttatttgagtaaattaacccacgatcccagccattctt'},
    {'class': '+', 'id': 'ALAS', 'sequencia': 'aacgcatacggtattttaccttcccagtcaagaaaacttatcttattcccacttttc'},
    {'class': '+', 'id': 'ARABAD', 'sequencia': 'ttagcggatcctacctgacgctttttatcgcaactctctactgtttctccatacccg'},
    {'class': '+', 'id': 'BIOA', 'sequencia': 'gccttctccaaaacgtgttttttgttgttaattcggtgtagacttgtaaacctaaat'},
    {'class': '+', 'id': 'DEOP1', 'sequencia': 'cagaaacgttttattcgaacatcgatctcgtcttgtgttagaattctaacatacggt'},
    {'class': '+', 'id': 'GALP2', 'sequencia': 'cactaatttattccatgtcacacttttcgcatctttgttatgctatggttatttcat'},
    {'class': '+', 'id': 'HIS', 'sequencia': 'atataaaaaagttcttgctttctaacgtgaaagtggtttaggttaaaagacatcagt'},
    {'class': '+', 'id': 'HISJ', 'sequencia': 'caaggtagaatgctttgccttgtcggcctgattaatggcacgatagtcgcatcggat'},
    {'class': '+', 'id': 'ILVGEDA', 'sequencia': 'ggccaaaaaatatcttgtactatttacaaaacctatggtaactctttaggcattcct'},
    {'class': '+', 'id': 'LACP1', 'sequencia': 'taggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtg'},
    {'class': '+', 'id': 'LPP', 'sequencia': 'ccatcaaaaaaatattctcaacataaaaaactttgtgtaatacttgtaacgctacat'},
    {'class': '+', 'id': 'TRPR', 'sequencia': 'tggggacgtcgttactgatccgcacgtttatgatatgctatcgtactctttagcgag'},
    {'class': '+', 'id': 'UVRB_P2', 'sequencia': 'tcagaaatattatggtgatgaactgtttttttatccagtataatttgttggcataat'},
    {'class': '-', 'id': ' 867', 'sequencia': 'atatgaacgttgagactgccgctgagttatcagctgtgaacgacattctggcgtcta'},
    {'class': '-', 'id': '1169', 'sequencia': 'cgaacgagtcaatcagaccgctttgactctggtattactgtgaacattattcgtctc'},
    {'class': '-', 'id': ' 802', 'sequencia': 'caatggcctctaaacgggtcttgaggggttttttgctgaaaggaggaactatatgcg'},
    {'class': '-', 'id': ' 521', 'sequencia': 'ttgacctactacgccagcattttggcggtgtaagctaaccattccggttgactcaat'},
    {'class': '-', 'id': ' 918', 'sequencia': 'cgtctatcggtgaacctccggtatcaacgctggaaggtgacgctaacgcagatgcag'},
    {'class': '-', 'id': '1481', 'sequencia': 'gccaatcaatcaagaacttgaagggtggtatcagccaacagcctgacatccttcgtt'},
    {'class': '-', 'id': '1024', 'sequencia': 'tggatggacgttcaacattgaggaaggcataacgctactacctgatgtttactccaa'},
    {'class': '-', 'id': '1149', 'sequencia': 'gaggtggctatgtgtatgaccgaacgagtcaatcagaccgctttgactctggtatta'},
    {'class': '-', 'id': ' 313', 'sequencia': 'cgtagcgcatcagtgctttcttactgtgagtacgcaccagcgccagaggacgacgac'},
    {'class': '-', 'id': ' 780', 'sequencia': 'cgaccgaagcgagcctcgtcctcaatggcctctaaacgggtcttgaggggttttttg'},
    {'class': '-', 'id': '1384', 'sequencia': 'ctacggtgggtacaatatgctggatggagatgcgttcacttctggtctactgactcg'},
    {'class': '-', 'id': ' 507', 'sequencia': 'atagtctcagagtcttgacctactacgccagcattttggcggtgtaagctaaccatt'},
    {'class': '-', 'id': '  39', 'sequencia': 'aactcaaggctgatacggcgagacttgcgagccttgtccttgcggtacacagcagcg'},
    {'class': '-', 'id': '1203', 'sequencia': 'ttactgtgaacattattcgtctccgcgactacgatgagatgcctgagtgcttccgtt'},
    {'class': '-', 'id': ' 988', 'sequencia': 'tattctcaacaagattaaccgacagattcaatctcgtggatggacgttcaacattga'},
    {'class': '-', 'id': '1171', 'sequencia': 'aacgagtcaatcagaccgctttgactctggtattactgtgaacattattcgtctccg'},
    {'class': '-', 'id': ' 753', 'sequencia': 'aagtgcttagcttcaaggtcacggatacgaccgaagcgagcctcgtcctcaatggcc'},
    {'class': '-', 'id': ' 630', 'sequencia': 'gaagaccacgcctcgccaccgagtagacccttagagagcatgtcagcctcgacaact'},
    {'class': '-', 'id': ' 660', 'sequencia': 'ttagagagcatgtcagcctcgacaacttgcataaatgctttcttgtagacgtgccct'},
    {'class': '-', 'id': '1216', 'sequencia': 'tattcgtctccgcgactacgatgagatgcctgagtgcttccgttactggattgtcac'},
    {'class': '-', 'id': ' 835', 'sequencia': 'tgctgaaaggaggaactatatgcgctcatacgatatgaacgttgagactgccgctga'},
    {'class': '-', 'id': '  35', 'sequencia': 'catgaactcaaggctgatacggcgagacttgcgagccttgtccttgcggtacacagc'},
    {'class': '-', 'id': '1218', 'sequencia': 'ttcgtctccgcgactacgatgagatgcctgagtgcttccgttactggattgtcacca'},
    {'class': '-', 'id': ' 668', 'sequencia': 'catgtcagcctcgacaacttgcataaatgctttcttgtagacgtgccctacgcgctt'},
    {'class': '-', 'id': ' 413', 'sequencia': 'aggaggaactacgcaaggttggaacatcggagagatgccagccagcgcacctgcacg'},
    {'class': '-', 'id': ' 991', 'sequencia': 'tctcaacaagattaaccgacagattcaatctcgtggatggacgttcaacattgagga'},
    {'class': '-', 'id': ' 751', 'sequencia': 'tgaagtgcttagcttcaaggtcacggatacgaccgaagcgagcctcgtcctcaatgg'},
    {'class': '-', 'id': ' 850', 'sequencia': 'ctatatgcgctcatacgatatgaacgttgagactgccgctgagttatcagctgtgaa'},
    {'class': '-', 'id': '  93', 'sequencia': 'gcggcagcacgtttccacgcggtgagagcctcaggattcatgtcgatgtcttccggt'},
    {'class': '-', 'id': '1108', 'sequencia': 'atccctaatgtctacttccggtcaatccatctacgttaaccgaggtggctatgtgta'},
    {'class': '-', 'id': ' 915', 'sequencia': 'tggcgtctatcggtgaacctccggtatcaacgctggaaggtgacgctaacgcagatg'},
    {'class': '-', 'id': '1019', 'sequencia': 'tctcgtggatggacgttcaacattgaggaaggcataacgctactacctgatgtttac'},
    {'class': '-', 'id': '  19', 'sequencia': 'tattggcttgctcaagcatgaactcaaggctgatacggcgagacttgcgagccttgt'},
    {'class': '-', 'id': '1320', 'sequencia': 'tagagggtgtactccaagaagaggaagatgaggctagacgtctctgcatggagtatg'},
    {'class': '-', 'id': '  91', 'sequencia': 'cagcggcagcacgtttccacgcggtgagagcctcaggattcatgtcgatgtcttccg'},
    {'class': '-', 'id': ' 217', 'sequencia': 'ttacgttggcgaccgctaggactttcttgttgattttccatgcggtgttttgcgcaa'},
    {'class': '-', 'id': ' 957', 'sequencia': 'acgctaacgcagatgcagcgaacgctcggcgtattctcaacaagattaaccgacaga'},
    {'class': '-', 'id': ' 260', 'sequencia': 'ggtgttttgcgcaatgttaatcgctttgtacacctcaggcatgtaaacgtcttcgta'},
    {'class': '-', 'id': ' 557', 'sequencia': 'aaccattccggttgactcaatgagcatctcgatgcagcgtactcctacatgaataga'},
    {'class': '-', 'id': '1355', 'sequencia': 'agacgtctctgcatggagtatgagatggactacggtgggtacaatatgctggatgga'},
    {'class': '-', 'id': ' 244', 'sequencia': 'tgttgattttccatgcggtgttttgcgcaatgttaatcgctttgtacacctcaggca'},
    {'class': '-', 'id': ' 464', 'sequencia': 'tgcacgggttgcgatagcctcagcgtattcaggtgcgagttcgatagtctcagagtc'},
    {'class': '-', 'id': ' 296', 'sequencia': 'aggcatgtaaacgtcttcgtagcgcatcagtgctttcttactgtgagtacgcaccag'},
    {'class': '-', 'id': ' 648', 'sequencia': 'ccgagtagacccttagagagcatgtcagcctcgacaacttgcataaatgctttcttg'},
    {'class': '-', 'id': ' 230', 'sequencia': 'cgctaggactttcttgttgattttccatgcggtgttttgcgcaatgttaatcgcttt'},
    {'class': '-', 'id': '1163', 'sequencia': 'tatgaccgaacgagtcaatcagaccgctttgactctggtattactgtgaacattatt'},
    {'class': '-', 'id': '1321', 'sequencia': 'agagggtgtactccaagaagaggaagatgaggctagacgtctctgcatggagtatga'},
    {'class': '-', 'id': ' 663', 'sequencia': 'gagagcatgtcagcctcgacaacttgcataaatgctttcttgtagacgtgccctacg'},
    {'class': '-', 'id': ' 799', 'sequencia': 'cctcaatggcctctaaacgggtcttgaggggttttttgctgaaaggaggaactatat'},
    {'class': '-', 'id': ' 987', 'sequencia': 'gtattctcaacaagattaaccgacagattcaatctcgtggatggacgttcaacattg'},
    {'class': '-', 'id': '1226', 'sequencia': 'cgcgactacgatgagatgcctgagtgcttccgttactggattgtcaccaaggcttcc'},
    {'class': '-', 'id': ' 794', 'sequencia': 'ctcgtcctcaatggcctctaaacgggtcttgaggggttttttgctgaaaggaggaac'},
    {'class': '-', 'id': '1442', 'sequencia': 'taacattaataaataaggaggctctaatggcactcattagccaatcaatcaagaact'}
]

# df = pd.DataFrame(dados)
df = pd.DataFrame.from_dict(dados)
# print(sorted(df['sequencia']))

# seq_enconded = []
# for index, row in df.iterrows():
#     # print(row['id'], row['sequencia'])
#     seq_enconded.append(one_hot_encode(row['sequencia']))

# df['seq_enc'] = seq_enconded


# Lista de sequências de DNA
labels = df['id']

# Codificando as sequências
encoded_sequences = np.array([one_hot_encode(seq) for seq in df['sequencia']])

# Preparando os dados para entrada no modelo
X = np.array(encoded_sequences)
y = np.array(labels)

print('X: %s; y: %d' % (len(X),  len(y)))
# print('X: %s; y: %d' % (X[0], y[0]))
# Definindo o modelo
# model = Sequential([
#     Flatten(input_shape=(4, 4)),  # Ajuste o input_shape conforme necessário
#     Dense(32, activation='relu'),
#     Dense(1, activation='sigmoid')  # Para classificação binária
# ])

# model = tf.keras.models.Sequential([
#   tf.keras.layers.Flatten(input_shape=(4, 28)),
#   tf.keras.layers.Dense(128, activation='relu'),
#   tf.keras.layers.Dropout(0.2),
#   tf.keras.layers.Dense(10, activation='softmax')
# ])



# # Compilando o modelo
# model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# # Treinando o modelo
# model.fit(X, y, epochs=10, batch_size=1)





# print(X)
# print(y)


# print(df)

# for index, row in df.iterrows():
#     print(row)


# for value in df.items():
#     print(value)

# columns = ['Student ID', 'Course ID', 'Marks']
# data = [(103, 201, 67), (103, 203, 67), (103, 204, 89)]
# df = pd.DataFrame(data, columns=columns)
# df.at['Total marks', 'Marks'] = df['Marks'].sum()

# print(df)

# print(df.shape)
# print(df)

# print(df.items)

# df.to_csv('list.csv', index=False)


# print(df)



# sequencia_str = 'ATCGTAGCGATAGATCGCTAGCTATCGT CATGACTGACTGATGCATGACTGACTGACTG'
# sequencia_str2 = ''

# for char in sequencia_str:
#     if char.upper() == 'A' or char.upper() == 'C' or char.upper() == 'G' or char.upper() == 'T': sequencia_str2 += char
#     else: sequencia_str2 += '-'
# # print(sequencia_str2)

# classes = list(range(1, (len(sequencia_str)) +1))

# sequencia = list(sequencia_str2)
# print(sequencia)
# dataset = {}

# # loop through sequences and split into individual nucleotides
# for i, seq in enumerate(sequences):

#     # split into nucleotides, remove tab characters
#     nucleotides = list(seq)
#     nucleotides = [x for x in nucleotides if x != '\t']

#     # append class assignment
#     nucleotides.append(classes[i])

#     # add to dataset
#     dataset[i] = nucleotides

# print(dataset[0])


# for seq in range(5,10):
#     print(seq)


# dframe = pd.DataFrame(dataset)



