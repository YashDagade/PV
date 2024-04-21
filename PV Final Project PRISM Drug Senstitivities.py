# -*- coding: utf-8 -*-

inputmuts = open ("OmicsSomaticMutations.csv","r")

#skip header line
skipline = inputmuts.readline() 

#create lists for each column
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

#for iindex,cell_line in enumerate(jak2mutlinelist):
  #  if "ACH" not in cell_line:
       # jak2mutlinelist.remove(cell_line)
       # jak2mutposlist.remove(jak2mutposlist[index])

#jak2mutlinelist.remove('')
#empty = jak2mutlinelist[34]
#jak2mutlinelist.remove(empty)        

print (jak2mutlinelist)
print (jak2mutposlist)
print (jak2mutreflist)
print (jak2mutaltlist)
print (jak2mutvariantinfolist)


inputmuts.close()

inputdrugsensitivities = open ("primary-screen-replicate-collapsed-logfold-change.csv",'r')

#skip header
skiplineagain = inputdrugsensitivities.readline()

#create lists for each cell line
cell_line_info = []
for row in inputdrugsensitivities:
    
    srow = row.split(',')
    cell_line_info.append(srow)
    
#isolate drug sensitivities for jak2 mutant lines
jak2mutinfolist = []
for cell_line in cell_line_info:
    for jak2mut in jak2mutlinelist:
        if jak2mut in cell_line:
            jak2mutinfolist.append(cell_line)
            


inputdrugsensitivities.close()





