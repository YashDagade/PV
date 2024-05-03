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
half_viable_drug_ids = []
for hvdritem in half_viable_list:
    half_viable_drug_list.append(hvdritem[6])
    half_viable_drug_ids.append(hvdritem[0])
    
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


#analyze pharmacological profiles of drugs with maximum 50% viability
input_pharm_profiles = open ("Drug_listThu Apr 25 17_24_43 2024.txt",'r')

skipp = input_pharm_profiles.readline()

all_pharm_profiles = []
for pharm in input_pharm_profiles:
    pharm_row = pharm.split('\t')
    all_pharm_profiles.append(pharm_row)
    
#isolate pharm profiles for half viability drugs
half_viable_pharm_profiles = []
for pharmprof in all_pharm_profiles:
    for drugid in half_viable_drug_ids:
        if drugid in pharmprof[0]:
            half_viable_pharm_profiles.append(pharmprof)

half_viable_pharm_profiles_nodups = []
for pharmprof_i in half_viable_pharm_profiles:
    if pharmprof_i not in half_viable_pharm_profiles_nodups:
        half_viable_pharm_profiles_nodups.append(pharmprof_i)

#make a list of targets and target pathways for half viability drugs
half_viable_targets = []
for pharmprof_ii in half_viable_pharm_profiles_nodups:
    target = pharmprof_ii[3]
    half_viable_targets.append(target.strip('"'))
half_viable_targets = list(filter(None, half_viable_targets))

half_viable_targets_split = []
for target_i in half_viable_targets:
    target_i = target_i.split(',')
    half_viable_targets_split.extend(target_i)

half_viable_targets_final = []
for target_ii in half_viable_targets_split:
    target_ii = target_ii.strip()
    half_viable_targets_final.append(target_ii)
   
#count number of occurences for each target
target_counts = {}
for target_item in half_viable_targets_final:
    if target_item in target_counts:
        target_counts[target_item] += 1
    else:
        target_counts[target_item] = 1
        

#make a pie chart for the drug target
target_labels = []
target_counts_pie = []
for c,d in target_counts.items():
    target_labels.append(c)
    target_counts_pie.append(d)

plt.pie(target_counts_pie, labels = target_labels)
plt.show()

for a, b in target_counts.items():
    print(a,b)   

    
half_viable_target_pathways = []
for pharmprof_iii in half_viable_pharm_profiles_nodups:
    target_pathway = pharmprof_iii[4]
    half_viable_target_pathways.append(target_pathway.strip('"'))

#count number of occurences for each pathway
target_pathway_counts = {}
for targetpath in half_viable_target_pathways:
    if targetpath in target_pathway_counts:
        target_pathway_counts[targetpath] += 1
    else:
        target_pathway_counts[targetpath] = 1
        
for key,value in target_pathway_counts.items():
    print(key,value)
    
#make a pie chart for the drug target pathways
pathway_labels = []
pathway_counts = []
for x,y in target_pathway_counts.items():
    pathway_labels.append(x)
    pathway_counts.append(y)

plt.pie(pathway_counts, labels = pathway_labels)
plt.show()

print("Aside from Unclassified and Other, PI3K/MTOR Signaling is the most frequently targeted pathway in drugs that reduce cell viability in JAK2 V617F mutant cells.")

#identify drugs that target PI3K/MTOR Signaling

pi3k_drug_ids = []
for profile in half_viable_pharm_profiles_nodups:
    if profile[4] == "PI3K/MTOR signaling":
        pi3k_drug_ids.append(profile[0])
        
pi3k_drugs = []
for drug_info in half_viable_list:
    for p_drugid in pi3k_drug_ids:
        if p_drugid in drug_info[0]:
            pi3k_drugs.append(drug_info[6])
            
pi3k_drugs_nodups = []
for p_drug in pi3k_drugs:
    if p_drug not in pi3k_drugs_nodups:
        pi3k_drugs_nodups.append(p_drug)
        
print(pi3k_drugs_nodups)

input_pharm_profiles.close()



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












#look at drugs with lower viabilities and compare

third_viable_list = []
for dr_item in jak2mutdruginfolist:
    if float(dr_item[3]) < 0.33:
        third_viable_list.append(dr_item)
            
third_viable_drug_list = []
third_viable_drug_ids = []
for third_dritem in third_viable_list:
    third_viable_drug_list.append(third_dritem[6])
    third_viable_drug_ids.append(third_dritem[0])

third_viable_drug_list_nodup = []
for third_drug in third_viable_drug_list:
    if third_drug not in third_viable_drug_list_nodup:
        third_viable_drug_list_nodup.append(third_drug)
        
        
        
tenth_viable_list = []
for drugitem_iii in jak2mutdruginfolist:
    if float(drugitem_iii[3]) < 0.1:
        tenth_viable_list.append(drugitem_iii)
        
tenth_viable_drug_list = []
tenth_viable_drug_ids = []
for tenth_dritem in tenth_viable_list:
    tenth_viable_drug_list.append(tenth_dritem[6])
    tenth_viable_drug_ids.append(tenth_dritem[0])

tenth_viable_drug_list_nodup = []
for tenth_drug in tenth_viable_drug_list:
    if tenth_drug not in tenth_viable_drug_list_nodup:
        tenth_viable_drug_list_nodup.append(tenth_drug)

tenth_viable_drug_list_nodup_strip = []
for ten_drug in tenth_viable_drug_list_nodup:
    ten_drug = ten_drug.strip('"')
    tenth_viable_drug_list_nodup_strip.append(ten_drug)
    
#isolate pharm profiles for tenth viability drugs
tenth_viable_pharm_profiles = []
for pharmprof_iii in all_pharm_profiles:
    for drugid_i in tenth_viable_drug_ids:
        if drugid_i in pharmprof_iii[0]:
            tenth_viable_pharm_profiles.append(pharmprof_iii)

tenth_viable_pharm_profiles_nodups = []
for pharmprof_iv in tenth_viable_pharm_profiles:
    if pharmprof_iv not in tenth_viable_pharm_profiles_nodups:
        tenth_viable_pharm_profiles_nodups.append(pharmprof_iv)

#make a list of targets and target pathways for tenth viability drugs
tenth_viable_targets = []
for pharmprof_v in tenth_viable_pharm_profiles_nodups:
    target_i = pharmprof_v[3]
    tenth_viable_targets.append(target_i.strip('"'))
tenth_viable_targets = list(filter(None, tenth_viable_targets))

tenth_viable_targets_split = []
for target_ii in tenth_viable_targets:
    target_ii = target_ii.split(',')
    tenth_viable_targets_split.extend(target_ii)

tenth_viable_targets_final = []
for target_iii in tenth_viable_targets_split:
    target_iii = target_iii.strip()
    tenth_viable_targets_final.append(target_iii)
    
tenth_viable_target_pathways = []
for pharmprof_vi in tenth_viable_pharm_profiles_nodups:
    target_pathway_i = pharmprof_vi[4]
    tenth_viable_target_pathways.append(target_pathway_i.strip('"'))

    
#count number of occurences for each pathway
target_pathway_tenth_counts = {}
for targetpath_tenth in tenth_viable_target_pathways:
    if targetpath_tenth in target_pathway_tenth_counts:
        target_pathway_tenth_counts[targetpath_tenth] += 1
    else:
        target_pathway_tenth_counts[targetpath_tenth] = 1
        
for e,f in target_pathway_tenth_counts.items():
    print(e,f)
    
#make a pie chart for the drug target pathways
pathway_tenth_labels = []
pathway_tenth_counts = []
for g,h in target_pathway_tenth_counts.items():
    pathway_tenth_labels.append(g)
    pathway_tenth_counts.append(h)

plt.pie(pathway_tenth_counts, labels = pathway_tenth_labels)
plt.show()


#identify most toxic drugs
drug_viability_scores = []
for drug_row in jak2mutdruginfolist:
    drug_viability = drug_row[2]
    drug_viability_scores.append(drug_viability)
    
enum_viability_scores = list(enumerate(drug_viability_scores))
sorted_viability_scores = sorted(enum_viability_scores, key = lambda x: x[1])

top_thirty_scores = sorted_viability_scores[:31]
top_thirty_index = []
for drug_pair in top_thirty_scores:
    drug_index = drug_pair[0]
    top_thirty_index.append(drug_index)
    
top_thirty_profiles = []
for top_index in top_thirty_index:
    top_thirty_profiles.append(jak2mutdruginfolist[top_index])
    
top_thirty_drugs = []
for top_profile in top_thirty_profiles:
    top_drug = top_profile[6]
    top_thirty_drugs.append(top_drug)
    
top_thirty_drug_ids = []
for top_profile_i in top_thirty_profiles:
    top_drug_id = top_profile_i[0]
    top_thirty_drug_ids.append(top_drug_id)
    
top_thirty_pharm_profiles = []
for top_drug_id_i in top_thirty_drug_ids:
    for prof in all_pharm_profiles:
        if prof[0] == top_drug_id_i:
            top_thirty_pharm_profiles.append(prof)
    


all_top_index = []
for drug_pair_i in sorted_viability_scores:
    drug_index_i = drug_pair_i[0]
    all_top_index.append(drug_index_i)
    
all_top_profiles = []
for top_index_i in all_top_index:
    all_top_profiles.append(jak2mutdruginfolist[top_index_i])

all_top_drugs = []
for top_profile_i in all_top_profiles:
    top_drug_i = top_profile_i[6]
    all_top_drugs.append(top_drug_i)


input_drug_doses = open ("sanger-dose-response.csv","r")

skippp = input_drug_doses.readline()

all_drug_doses = []
for rrow in input_drug_doses:
    rrow = rrow.split(',')
    all_drug_doses.append(rrow)

jakmut_drug_doses = []
for drug_dose in all_drug_doses:
    for mut_line in jak2V617Fmutlines:
        if drug_dose[9] == mut_line:
            jakmut_drug_doses.append(drug_dose)
            
drug_ic50s = []
for mut_drug_dose in jakmut_drug_doses:
    drug_ic50s.append(mut_drug_dose[7])
    
enum_drug_ic50s = list(enumerate(drug_ic50s))
sorted_drug_ic50s = sorted(enum_drug_ic50s, key = lambda x: x[1])

drug_dose_index = []
for index_drug in sorted_drug_ic50s:
    drug_dose_index.append(index_drug[0])
    
sorted_drug_doses = []
for drug_dose_index_item in drug_dose_index:
    sorted_drug_doses.append(jakmut_drug_doses[drug_dose_index_item])
    
sorted_drug_names = []
for drug_dose_item in sorted_drug_doses:
    sorted_drug_names.append(drug_dose_item[10])
    
sorted_drug_ic50s_noindex = []
for index_drug_i in sorted_drug_ic50s:
    sorted_drug_ic50s_noindex.append(index_drug_i[1])

input_drug_doses.close()


input_drug_viability.close()





