# -*- coding: utf-8 -*-

inputmuts = open ("OmicsSomaticMutations.csv","r")

#skip header line
skipline = inputmuts.readline() 





#create lists for each column
rowlist = []
chrnumlist = []
poslist = []
reflist= []
altlist = []
variantinfolist = []
hugoidlist = []
geneidlist = []
modelidlist = []

for line in inputmuts:
    
    line = line.split(",")
    rowlist.append(line)
    
    chrnum = line[0]
    chrnumlist.append(chrnum)
    
    pos = line[1]
    poslist.append(pos)
    
    ref = line[2]
    reflist.append(ref)
    
    alt = line[3]
    altlist.append(alt)
    
    variantinfo = line[11]
    variantinfolist.append(variantinfo)
    
    hugoid = line[14]
    hugoidlist.append(hugoid)
    
    geneid = line[15]
    geneidlist.append(geneid)
    
    modelid = line[62]
    modelidlist.append (modelid)
    
    
#identify cell lines with mutations in jak2
#look at mutation info for mutated cell lines
jak2mutlinelist = []
jak2mutposlist = []
jak2mutreflist = []
jak2mutaltlist = []
jak2mutvariantinfolist = []

for index,gid in enumerate(geneidlist): 
    if "ENSG00000096968" in gid:
        jak2mutlinelist.append(modelidlist[index])
        jak2mutposlist.append(poslist[index])
        jak2mutreflist.append(reflist[index])
        jak2mutaltlist.append(altlist[index])
        jak2mutvariantinfolist.append(variantinfolist[index])
        
        
        
        
        
        
#V617F specifically

        
jak2V617Flist = []
for item in rowlist:    
    if "5073770" in item:
        jak2V617Flist.append(item)

jak2V617Fmutlines = []
for b in jak2V617Flist:
    mutline = b[67]
    jak2V617Fmutlines.append(mutline)

print ("Cell Lines with JAK2 V617F Mutation:",*jak2V617Fmutlines, sep = "\n")       

#print (jak2mutlinelist)
#print (jak2mutposlist)
#print (jak2mutreflist)
#print (jak2mutaltlist)
#print (jak2mutvariantinfolist)

#chr9:g.5073770G>T

inputmuts.close()





#info for jak2 mutated cell lines
input_cell_line_info = open ("primary-screen-cell-line-info.csv","r")

#skip header
skiplineagainn = input_cell_line_info.readline()

#create lists for each cell line
cell_line_info = []
for roww in input_cell_line_info:
    crow = roww.split(',')
    cell_line_info.append(crow)

#isolate cell line info for jak2 mutant lines
jak2mutinfolist = []
for ccell_line in cell_line_info:
    for jak2mut in jak2V617Fmutlines:
        if jak2mut == ccell_line[0]:
            jak2mutinfolist.append(ccell_line)
#none of them are leukemia?

input_cell_line_info.close()



import matplotlib.pyplot as plt
import numpy as np


input_drug_viability = open ("sanger-viability.csv","r")

skip_header = input_drug_viability.readline()

drugrowlist = []
for brow in input_drug_viability:
    drugrow = brow.split(',')
    drugrowlist.append(drugrow)
    
jak2mutdruginfolist = []
for drugitem in drugrowlist:
    for cell_line_item in jak2V617Fmutlines:
        if cell_line_item == drugitem[5]:
            jak2mutdruginfolist.append(drugitem)
            
half_viable_list = []
for dritem in jak2mutdruginfolist:
    if float(dritem[3]) < 0.50:
        half_viable_list.append(dritem)
            
half_viable_drug_list = []
for hvdritem in half_viable_list:
    half_viable_drug_list.append(hvdritem[6])
    
half_viable_scores = []
for hv_item in half_viable_list:
    half_viable_scores.append(hv_item[3])
    
half_dose_values = []
for hv_item_i in half_viable_list:
    half_dose_values.append(hv_item_i[2])

#plot a histogram of the dose values for those with viability under 50%    
plt.hist(half_dose_values, bins = 150)
plt.xticks(visible=False)

#make histogram of viability scores for those with viability under 50%
#plt.hist(half_viable_scores)

#remove duplicates
half_viable_drug_list_nodup = []
for hv_drug in half_viable_drug_list:
    if hv_drug not in half_viable_drug_list_nodup:
        half_viable_drug_list_nodup.append(hv_drug)
        
print ("Drugs that will kill 50% of JAK2 V617 mutant cells:",*half_viable_drug_list_nodup, sep = "\n")
#print (half_viable_list)

#LESTAURTINIB is a JAK2 inhibitor and is in clinical trials for PV treatment and AML treatment
#AIRCAR used in treatment of other cancers but no testing for PV
#VORINOSTAT effective in murine PV models


#compare drugs that will kill 50% of JAK2 V617 mutant cells to drugs that will kill 50% of any cancer cell
all_half_viable_drug_list = []
for drugitem_i in drugrowlist:
    if float(drugitem_i[3]) < 0.50:
        if drugitem_i[6] not in all_half_viable_drug_list:
            all_half_viable_drug_list.append (drugitem_i[6])
    
jak2_unique_hv_drugs = []
for drugitem_ii in half_viable_drug_list_nodup:
    if drugitem_ii not in all_half_viable_drug_list:
        jak2_unique_hv_drugs.append (drugitem_ii)
#there are no drugs specific to just JAK2 V617 mutants; this is unsurprising considering these drugs could work on similar cells with similar cancer types but not this specific mutation

#https://www.cell.com/fulltext/S0092-8674(16)30746-2

input_drug_viability.close()





