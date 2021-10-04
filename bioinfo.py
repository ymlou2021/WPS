# pandas useful func.
'''
1. data types:
https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.astype.html
df.dtypes
>> show data types of the data frame
df_new=df.astype({'chr': 'string'}) # cast data type by column
>> the dtype in new df is changed

2. sort files
df_sort = df.sort_values(by=['Chromosome/scaffold name','Gene start (bp)']) # sort the chr first and then sort start inside the chr
'''
import pandas as pd
# read wps file
wps = pd.read_table("EE85756.hg19.wps.mapq30.bedGraph", names=['chr','start','end','value'])
# read gene file
gene_anno = pd.read_table("human_gene_hg38.txt")
# sort gene file by chr and start position
gene_anno_sort = gene_anno.sort_values(by=['Chromosome/scaffold name','Gene start (bp)'])

# group by chr
gene_groupby=gene_anno_sort.groupby("Chromosome/scaffold name") 
wps_groupby=wps.groupby("chr")
chr_name=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']

# define the dictionary storing gene version ID and its avg. WPS score
dic_gene_score={}

# loops the chromosome
for j in chr_name:
    gene_chr=gene_groupby.get_group(j)
    a='chr'+j
    wps_chr=wps_groupby.get_group(a)
    # loops through the gene in chrj
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
        for k in wps_gene.index:
            wps_start=wps_gene.loc[k,'start']
            wps_end=wps_gene.loc[k,'end']
            if wps_start <= wps_end:
                value=wps_gene.loc[k,'value']
                wps_score=(wps_end-wps_start+1)*value
                bp_no=bp_no+wps_end-wps_start+1
                gene_all_score=gene_all_score+wps_score
        if bp_no == 0:
            gene_score = 'NaN'
        else:
            gene_score=gene_all_score/bp_no
        dic_gene_score[gene_chr.loc[i,"Gene stable ID version"]]= gene_score

print (dic_gene_score)