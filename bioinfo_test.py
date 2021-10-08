import pandas as pd
import time

# start time
t0=time.time()

# read wps file
wps = pd.read_table("C:\Work\Radiooncology\EE85756.hg38.wps.mapq30.bedgraph", 
                    names=['chr','start','end','value'],
                    dtype={"start": 'int32', "end": 'int32',
                           "value": 'int16'}) #setting datatypes reduces memory
# read gene file
gene_anno = pd.read_table("C:\Work\Radiooncology\human_gene_hg38.txt", 
                          dtype={"Gene start (bp)": 'int32',
                           "Gene end (bp)": 'int32'})
gene_anno = gene_anno.drop(columns=['Gene stable ID version'])
# sort gene file by chr and start position
gene_anno = gene_anno.sort_values(by=['Chromosome/scaffold name',
                                      'Gene start (bp)'])

# time when data loaded
t1=time.time()

# test on chr 21 and chr 1

#gene_chr_21 = gene_anno.groupby("Chromosome/scaffold name").get_group("21")
gene_chr_1 = gene_anno.groupby("Chromosome/scaffold name").get_group("1")

#wps_chr = wps.groupby("chr").get_group("chr21")
wps_chr = wps.groupby("chr").get_group("chr1")

#gene_chr = gene_chr_21
gene_chr = gene_chr_1
# gene_chr=gene_chr_1[0:1]
# wps_chr=wps[0:1000]

# define the dictionary storing gene version ID and its avg. WPS
dic_gene_score={}

# time before looping
t2=time.time()

# loops the gene
for i in gene_chr.index:
    gene_all_score=0
    gene_str = gene_chr.loc[i,"Gene start (bp)"]
    gene_end = gene_chr.loc[i,"Gene end (bp)"]
    # change the wps str and end as gene str and end
    close_to_start = wps_chr[wps_chr["start"] <= gene_str].start.max() # the highest number of the start locations that is still lesser than (or equal to) the gene's start location
    close_to_end = wps_chr[wps_chr["end"] >= gene_end].end.min()
    wps_gene = wps_chr.loc[(wps_chr.start >= close_to_start) &
                           (wps_chr.end <= close_to_end), :]
    wps_gene=pd.DataFrame(
        {'start':wps_gene['start'].mask(wps_gene['start'] < gene_str, gene_str),
         'end':wps_gene['end'].mask(wps_gene['end'] > gene_end, gene_end),
         'value':wps_gene['value']}
    ) # Note that if the gene's start or end location is in an area where there is no WPS entry, then the masking will lead to funny results
    wps_gene = wps_gene.loc[(wps_gene.start < wps_gene.end), :] # this line removes those funny results
# Dataframes are empty if there are no entries in the WPS file along the gene
    if wps_gene.empty:
        gene_score = None # don't use "NaN"! I don't actually know if it would be interpreted the correctly, but even if it is in this case, it'd still be bad practice.
    else:
        wps = (wps_gene.end - wps_gene.start) * wps_gene.value
        gene_score = wps.sum() / (wps_gene.end.max() - wps_gene.start.min())
    dic_gene_score[gene_chr.loc[i, "Gene stable ID"]] = gene_score

print(dic_gene_score)
# finish time
t3=time.time()
print(t3-t2)
print(t3-t1)
print(t3-t0)
