#!/gpfs/data/yarmarkovichlab/Frank/JN_BAT_5BDHT/gseapy_env/bin/python3.10

import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import gseapy as gp
from Bio.Align import substitution_matrices
from ast import literal_eval

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# # immunoverse per cancer
# df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/stats/final_all_ts_antigens.txt',sep='\t')
# ultimate_data = []
# for cancer,sub_df in df.groupby(by='cancer'):
#     peps = ','.join(list(sub_df['pep'].unique()))
#     ultimate_data.append(('immunoverse_{}_specific'.format(cancer),peps,'phenotype'))

# # differential condition

# # distance to star molecule
# star = {
#     'PRAME_A0201_1':'SLLQHLIGL',
#     'MAGEA4_A0201_1':'GVYDGREHTV',
#     'NYESO1_A0201_1':'SLLMWITQC',
#     'KRAS_G12V_A1101_1':'VVVGAVGVGK',
#     'KRAS_G12V_A1101_2':'VVGAVGVGK',
#     'KRAS_G12D_C0802_1':'GADGVGKSAL',
#     'KRAS_G12D_C0802_2':'GADGVGKSA',
#     'KRAS_G12D_A1101_1':'VVVGADGVGK',
#     'MART1_A0201_1':'EAAGIGILTV',
#     'HAVCR1_1':'DLSRRDVSL',
#     'sPMEL_1':'KTWDQVPFSV',
#     'PMEL_can1':'KTWGQYWQV',
#     'PMEL_can2':'ITDQVPFSV',
#     'PMEL_ImmTac_can3':'YLEPGPVTA',
#     'ORF2_A2402_1':'ILPKVIYRF'
# }

# normal_path = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/code/hla_ligand_atlas_now_0.05_tesorai.txt'
# normal = pd.read_csv(normal_path,sep='\t')

# def similarity(pep):
#     blosum = substitution_matrices.load("BLOSUM62") 
#     alphabet = np.array(list("ARNDCQEGHILKMFPSTWYV"))
#     A = len(alphabet)
#     M = np.zeros((A, A), dtype=np.int16)
#     aa_to_pos_in_bio = {aa:i for i, aa in enumerate(alphabet)}
#     for i, ai in enumerate(alphabet):
#         for j, aj in enumerate(alphabet):
#             M[i, j] = blosum[ai, aj]

#     # get subset
#     l = len(pep)
#     cond = [True if len(item) == l else False for item in normal['peptide']]
#     tmp = normal.loc[cond,:]
#     tmp = tmp.drop_duplicates(subset='peptide')
#     array = []
#     for item in tmp['peptide']:
#         array.append(list(item))
#     nps = np.array(array)
#     D = np.zeros_like(nps,dtype=np.int16)

#     # pep
#     cp = np.array(list(pep))
#     for j in np.arange(len(cp)):
#         ref = aa_to_pos_in_bio[cp[j]]
#         queries = [aa_to_pos_in_bio[item] for item in nps[:,j]]
#         D[:,j] = M[queries,ref]
#     hotspots = [4,5,6,7,8]
#     hotspots = [item-1 for item in hotspots]
#     D = D[:,hotspots]
#     distances = D.sum(axis=1)
#     indices = np.flip(np.argsort(distances))
#     sorted_nps = nps[indices,:]
#     distances = distances[indices]
#     return sorted_nps,distances

# for title,pep in star.items():
#     sorted_nps, distances = similarity(pep)
#     lis = []
#     for i in np.arange(100):
#         tmp = ''.join(sorted_nps[i].tolist())
#         lis.append(tmp)
#     if pep in lis:
#         lis.remove(pep)
#     lis = ','.join(lis)
#     ultimate_data.append(('{}_{}_top50_similar'.format(title,pep),lis,'distance_to_star'))



# # with open('result.txt','w') as f:
# #     title = '\t'.join(['pos{}'.format(i+1) for i in range(len(pep))])
# #     f.write('{}\n'.format(title))
# #     ref = '\t'.join(list(pep))
# #     f.write('{}\n'.format(ref))
# #     for row,d in zip(sorted_nps,distances):
# #         q = '\t'.join(list(row))
# #         f.write('{}\t{}\n'.format(q,d))

# # inherent mSigDB
# rootdir = '/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas'
# cancers = [
#     'BRCA',
#     'KIRC',
#     'COAD',
#     'STAD',
#     'MESO',
#     'LIHC',
#     'ESCA',
#     'CESC',
#     'BLCA',
#     'RT',
#     'AML',
#     'DLBC',
#     'GBM',
#     'NBL',
#     'PAAD',
#     'HNSC',
#     'OV',
#     'LUSC',
#     'LUAD',
#     'CHOL',
#     'SKCM'
# ] 

# data = []
# for c in cancers:
#     p = os.path.join(rootdir,c,'antigen','0.01','final_enhanced_all.txt')
#     final = pd.read_csv(p,sep='\t')
#     final = final.loc[final['typ']=='self_gene',:]
#     for row in final.itertuples():
#         peptide = row.pep
#         weight = row.relative_abundance
#         if row.ensgs.startswith('ENSG'):
#             data.append((peptide,weight,row.ensgs,c))
#         else:
#             ensgs = literal_eval(row.ensgs)
#             for ensg in ensgs:
#                 data.append((peptide,weight,ensg,c))
# df = pd.DataFrame.from_records(data,columns=['peptide','weight','ensg','cancer'])
# dic = {}
# for ensg,sub_df in df.groupby(by='ensg'):
#     dic[ensg] = sub_df['peptide'].unique().tolist()

# gene_lfc = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/atlas/SKCM/gene_lfc.txt',sep='\t',index_col=0)
# ensg2symbol = gene_lfc['gene_symbol'].to_dict()
# symbol2ensg = {v:k for k,v in ensg2symbol.items()}

# df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/ImmunoSigDB/mSigDB/h.all.v2026.1.Hs.symbols.gmt',sep='\t',header=None)
# for i in range(df.shape[0]):
#     lis = []
#     title = df.iloc[i].iloc[0]
#     for gene in df.iloc[i].iloc[2:]:
#         ensg = symbol2ensg.get(gene,'unknown')
#         if ensg != 'unknown':
#             lis.extend(dic.get(ensg,[]))
#     ultimate_data.append((title,','.join(lis),'mSigDB'))


# final = pd.DataFrame.from_records(ultimate_data,columns=['ontology','peps','category'])
# final.to_csv('db/final.txt',sep='\t',index=None)

# diff pep
df = pd.read_csv('/gpfs/data/yarmarkovichlab/Frank/pan_cancer/NYU_Tesorai_all_searches/tesorai_peptide_fdr_ALL_BALL_P1.tsv',sep='\t')
meta = pd.read_csv('/gpfs/data/yarmarkovichlab/public/ImmunoVerse/ImmunoVerse_Hub/BALL_metadata.txt',sep='\t')
meta = meta.loc[meta['biology'].str.startswith('SupB15'),:]
dic = {}
valid_samples = []
for biology,sub_df in meta.groupby(by='biology'):
    dic[biology] = sub_df['sample'].unique().tolist()
    valid_samples.extend(dic[biology])

df = df.loc[df['filename'].isin(valid_samples),:]
pep = 'GLSLPLPPK'
df = df.loc[df['clean_sequence']==pep,:]








# test gsea
query = pd.read_csv('/gpfs/data/yarmarkovichlab/JH_AML/tesorai_results/patient6_run_pep_fdr.tsv',sep='\t')
# rnk = query.loc[:,['pep','relative_abundance']].set_index(keys='pep')
rnk = query.loc[:,['clean_sequence','score']].drop_duplicates(subset='clean_sequence').set_index(keys='clean_sequence')
final = pd.read_csv('db/final.txt',sep='\t')
gene_sets = {row.ontology:row.peps.split(',') for row in final.itertuples()}
pre_res = gp.prerank(rnk=rnk,
                     gene_sets=gene_sets,
                     threads=4,
                     min_size=1,
                     max_size=1000,
                     permutation_num=100,
                     outdir=None,
                     seed=6,
                     verbose=True)
pre_res.res2d.to_csv('result.txt',sep='\t')
from gseapy import gseaplot2
terms = pre_res.res2d['Term'][0:10]
gseaplot2(terms=terms,
          RESs=[pre_res.results[term]['RES'] for term in terms],
          hits=[pre_res.results[term]['hits'] for term in terms],
          rank_metric=pre_res.ranking,
          ofname='result.pdf')