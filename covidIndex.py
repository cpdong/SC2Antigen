import argparse,os;
import pandas as pd;
from os.path import dirname, abspath
scriptDir = dirname(abspath(__file__));# get parent dir of the scripts
# =============================================================================
# initial setting for future update
# =============================================================================
#
parser = argparse.ArgumentParser()
parser.add_argument('-hla', metavar = 'inputstring', dest='inputHLA', help='Give arcasHLA calling result fullname');
parser.add_argument('-len', metavar = 'pep length', type=str, dest='pepLen', help='Give the peptide files');
parser.add_argument('-o', metavar = 'output', dest='outputFile',help='Give output file fullname with path');

args = parser.parse_args();
inputHLA  = args.inputHLA;
#pepLen= args.pepLen;
outputFile = args.outputFile;

path, filename = os.path.split(outputFile)
basename, ext = os.path.splitext(filename)

#================================================================
#pepLen_list= pepLen.split(',');
virus_protein = ["ORF1ab-protein","S-protein","ORF3a-protein","E-protein","M-protein","ORF6-protein",
"ORF7a-protein","ORF7b-protein","ORF8-protein","N-protein","ORF10-protein"]

netMHC_result= scriptDir + '/data/SARS-COV2-pMHC-I_9mer_result.txt.gz'
#================================================================
#####
hla_ref = [l.rstrip() for l in open(scriptDir + '/data/allelenames')]
if inputHLA: #get HLA from user input
    hla_list= inputHLA.split(',')
    if 'HLA' in hla_list[0]:
        hla_list= [x.replace('*','') for x in hla_list]
    else:
        hla_list= ['HLA-' + x.replace('*','') for x in hla_list]
    hla_list= [ x for x in hla_list if x in hla_ref ];
    hla_list = list(set(hla_list)); # combine the identical HLA
    hla_list.sort();
    #hla_String= [x[:5] + '*' + x[5:] for x in hla_list]
else:
    hla_list=None;
#================================================================
#================================================================

df= pd.read_csv(netMHC_result, header=0, sep='\t', compression='gzip', quotechar='"', error_bad_lines=False);

header= ["peptide","protein","pos"] + hla_list;
data_list= df[header].values.tolist();

weakList= [sublist for sublist in data_list if any(x for x in sublist[3:] if float(x) <= 2)];
weakList= [header] + weakList
df_weak = pd.DataFrame(weakList[1:],columns=weakList[0]);
df_weak.to_csv(path + '/' + basename + '.weakBinder.txt', index=False, sep='\t');

strongList= [sublist for sublist in data_list if any(x for x in sublist[3:] if float(x) <= 0.5)];
strongList= [header] + strongList
df_strong = pd.DataFrame(strongList[1:],columns=strongList[0]);
df_strong.to_csv(path + '/' + basename + '.strongBinder.txt', index=False, sep='\t');
#
result=[['protein'] + hla_list + ['bestRank', 'phbr_Score','weakBinder','strongBinder']]
for protein in virus_protein:
    proX= [protein];
    protein_list= [sublist for sublist in data_list if sublist[1] == protein];

    protein_weak_num= len([sublist for sublist in protein_list if any(x for x in sublist[3:] if float(x) <= 2)]);
    protein_strong_num= len([sublist for sublist in protein_list if any(x for x in sublist[3:] if float(x) <= 0.5)]);

    hla_stat=[];
    bestScore={}
    for z in range(3,len(hla_list)+3):
        tmp_min = min([ sublist[z] for sublist in protein_list ]);
        #print(tmp_min)
        tmp_list= [sublist for sublist in protein_list if float(sublist[z]) == tmp_min];
        #print(tmp_list)
        hla_stat.append(tmp_list[0][0] + '_' + str(tmp_min));
        bestScore[hla_list[z-3]]= tmp_min;

    proX= proX + hla_stat;
    bestRank= min(bestScore.values());
    proX.append(bestRank);
    
    recip_list=[]
    for g in ['A','B','C']: # defined for MHC-I genes
        gene = [x for x in hla_list if 'HLA-'+ g in x];
        if len(gene)==1:
            geneScore = [bestScore[gene[0]]]*2;
        elif len(gene)==2:
            geneScore = [bestScore[ge] for ge in gene];
        elif len(gene)==0:
            geneScore = [];
        #print(g, geneScore);
        recip_list= recip_list + geneScore;
    recip_list= [ 1/x for x in recip_list ];
    protein_phbr= round(len(recip_list)/sum(recip_list),4);
    proX.append(protein_phbr);
    proX.append(protein_weak_num);
    proX.append(protein_strong_num);
    
    result.append(proX);

df_result = pd.DataFrame(result[1:],columns=result[0]);
df_result.to_csv(outputFile, index=False, sep='\t');
#
