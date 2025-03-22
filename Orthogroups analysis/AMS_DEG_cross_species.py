import re
import pandas as pd

# read AM response genes in different species
AMs_Mt_data = pd.read_csv("/home/wuzefeng/MyResearch/Medicago_RNA/analysis_result/Mt_CK_AM.deg.csv",index_col=0)
AMs_Mt = list(AMs_Mt_data[AMs_Mt_data.change!="NOT"].index)

AMs_Gm_data = pd.read_csv("/home/wuzefeng/MyResearch/Gmax_RNA/analysis_result/DEG_AM_CK.csv",index_col=0)
AMs_Gm = list(AMs_Gm_data[AMs_Gm_data.change!="NOT"].index)

AMs_Lt_data = pd.read_csv("/home/wuzefeng/MyResearch/Lotus_RNA/analysis_result/Lt_CK_AM.deg.csv",index_col=0)
AMs_Lt = list(AMs_Lt_data[AMs_Lt_data.change!="NOT"].index)

AMs_Zm_data = pd.read_csv("/home/wuzefeng/MyResearch/Maize_RNA/analysis_result/DEG_AM_CK.csv",index_col=0)
AMs_Zm = list(AMs_Zm_data[AMs_Zm_data.change!="NOT"].index)

AMs_Ri_data = pd.read_csv("/home/wuzefeng/MyResearch/Rice_RNA/analysis_result/DEG_AM_CK.csv",index_col=0)
AMs_Ri = list(AMs_Ri_data[AMs_Ri_data.change!="NOT"].index)

AMs_Sy_data = pd.read_csv("/home/wuzefeng/MyResearch/Sly_RNA/analysis_result/DEG_AM_CK.csv",index_col=0)
AMs_Sy = list(AMs_Sy_data[AMs_Sy_data.change!="NOT"].index)

All_AMs = AMs_Mt + AMs_Gm + AMs_Lt + AMs_Zm + AMs_Ri+AMs_Sy

#outfile = open("/home/wuzefeng/MyResearch/2cross_species_compare/AM_orthogous","w")
outfile2 = open("/home/wuzefeng/MyResearch/2cross_species_compare/AM_orthogous_genes","w")
#outfile.write("\t".join(["Orthogroups","Glycine_max","Glycine_max_AM","Lotus_japonicus","Lotus_japonicus_AM","Medicago_truncatula","Medicago_truncatula_AM","Oryza_sativa","Oryza_sativa_AM","Soly","Soly_AM","Maize","Maize_AM"])+"\n")

with open("/home/wuzefeng/MyResearch/1orthofinder/0proteins/OrthoFinder/Results_Sep27/Orthogroups/Orthogroups.txt") as fh:
    for row in fh:
        data = row.strip().split(" ")
        ortho_group = data[0]
        temp_key = {"Glycine_max":[0,0],"Lotus_japonicus":[0,0],"Medicago_truncatula":[0,0],"Oryza_sativa":[0,0],"Soly":[0,0],"Maize":[0,0]}
        temp_key2 = {"Glycine_max":[],"Lotus_japonicus":[],"Medicago_truncatula":[],"Oryza_sativa":[],"Soly":[],"Maize":[]}
                
        sp_genes = [m.split("|") for m in data[1:]]
        
        # modify lotus gene names
        for sp_gene in sp_genes:
            if sp_gene[0] == "Lotus_japonicus":
                sp_gene[1] = sp_gene[1].split(".")[0]
        #print (sp_genes)
        
        # filling elements into temp_key
        for sp_gene in sp_genes:
            temp_key[sp_gene[0]][0]+=1
            if sp_gene[1] in All_AMs:
                temp_key[sp_gene[0]][1]+=1
                temp_key2[sp_gene[0]].append(sp_gene[1])
        #outfile.write(ortho_group+"\t"+"\t".join([str(m[0])+"\t"+str(m[1]) for m in temp_key.values()])+"\n")
        print(temp_key2.values())
        
        if sum([len(x) for x in temp_key2.values()])>=1:
            outfile2.write(ortho_group+"\t"+"\t".join(["\t".join(m) for m in temp_key2.values()if len(m)!=0])+"\n")
print ("OK")
        
