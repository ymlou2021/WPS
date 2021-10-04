import pandas as pd
import time

# timing point
t0=time.time()

# read wps file
wps = pd.read_table("EE85756.hg19.wps.mapq30.bedGraph", names=['chr','start','end','value'])
# read gene file
gene_anno = pd.read_table("human_gene_hg38.txt")
# sort gene file by chr and start position
gene_anno_sort = gene_anno.sort_values(by=['Chromosome/scaffold name','Gene start (bp)'])

# timing point
t1=time.time()

# test on chr 21 and chr 1
gene_groupby = gene_anno_sort.groupby("Chromosome/scaffold name")
gene_chr_21 = gene_groupby.get_group("21")
# gene_chr_1 = gene_groupby.get_group("1")

wps_groupby = wps.groupby("chr")
wps_chr=wps_groupby.get_group("chr21")
# wps_chr=wps_groupby.get_group("chr1")

gene_chr=gene_chr_21[0:1]
# gene_chr=gene_chr_1[0:1]
# wps_chr=wps[0:1000]

# define the dictionary storing gene version ID and its avg. WPS
dic_gene_score={}

# loops the gene
for i in gene_chr.index:
    bp_no=0
    gene_all_score=0
    gene_str = gene_chr.loc[i,"Gene start (bp)"]
    gene_end = gene_chr.loc[i,"Gene end (bp)"]
    # change the wps str and end as gene str and end
    wps_gene=pd.DataFrame(
        {
            'start':wps_chr['start'].mask( wps_chr['start'] < gene_str, gene_str),
            'end':wps_chr['end'].mask( wps_chr['end'] > gene_end, gene_end),
            'value':wps_chr['value']
        }
    )
    # t3=time.time()
    for k in wps_gene.index:
        wps_start = wps_gene.loc[k,'start']
        wps_end = wps_gene.loc[k,'end']
        if wps_start <= wps_end:
            value = wps_gene.loc[k,'value']
            wps_score = (wps_end - wps_start + 1) * value
            bp_no = bp_no + wps_end - wps_start + 1
            gene_all_score = gene_all_score + wps_score
    if bp_no == 0:
        gene_score = 'NaN'
    else:
        gene_score = gene_all_score/bp_no
    dic_gene_score[gene_chr.loc[i,"Gene stable ID version"]] = gene_score

print (dic_gene_score)

# timing point
t2=time.time()

print (t2-t1)