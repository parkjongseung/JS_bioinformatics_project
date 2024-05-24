import pandas as pd
cnts = pd.read_csv('read-counts.txt', sep='\t', comment='#', index_col=0)

## RPM 구하기
total_reads_CLIP = cnts['binfo1-datapack1/CLIP-35L33G.bam'].sum()
total_reads_RNA = cnts['binfo1-datapack1/RNA-control.bam'].sum()
total_reads_RPF_Lin = cnts['binfo1-datapack1/RPF-siLin28a.bam'].sum()
total_reads_RNA_Lin = cnts['binfo1-datapack1/RNA-siLin28a.bam'].sum()
total_reads_RPF_Luc = cnts['binfo1-datapack1/RPF-siLuc.bam'].sum()
total_reads_RNA_Luc = cnts['binfo1-datapack1/RNA-siLuc.bam'].sum()
cnts['CLIP-RPM'] = (cnts['binfo1-datapack1/CLIP-35L33G.bam']/ total_reads_CLIP) * 1e6
cnts['RNA-RPM'] = (cnts['binfo1-datapack1/RNA-control.bam']/ total_reads_RNA) * 1e6
cnts['RPF-Lin-RPM'] = (cnts['binfo1-datapack1/RPF-siLin28a.bam']/ total_reads_RPF_Lin) * 1e6
cnts['RNA-Lin-RPM'] = (cnts['binfo1-datapack1/RNA-siLin28a.bam']/total_reads_RNA_Lin) * 1e6
cnts['RPF-Luc-RPM'] = (cnts['binfo1-datapack1/RPF-siLuc.bam']/total_reads_RPF_Luc) * 1e6
cnts['RNA-Luc-RPM'] = (cnts['binfo1-datapack1/RNA-siLuc.bam']/total_reads_RNA_Luc) * 1e6

cnts['clip_enrichment'] = cnts['CLIP-RPM'] / cnts['RNA-RPM']
cnts['rden_change'] = (cnts['RPF-Lin-RPM'] / cnts['RNA-Lin-RPM']) / (cnts['RPF-Luc-RPM'] / cnts['RNA-Luc-RPM'])

#low read count(<30)와 low ribosome foot prints(<80) 삭제
idx1 = cnts[cnts['binfo1-datapack1/CLIP-35L33G.bam'] < 30].index
cnts.drop(idx1, inplace=True)
idx2 = cnts[cnts['binfo1-datapack1/RPF-siLin28a.bam'] < 80].index
cnts.drop(idx2, inplace=True)
cnts.head(20)

from matplotlib import pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.scatter(np.log2(cnts['clip_enrichment']),
               np.log2(cnts['rden_change']), s = 1) # dot 크기 조절

#논문처럼 x축, y축 범위 정함
plt.xlim([-6, 4])
plt.ylim([-2, 2])


