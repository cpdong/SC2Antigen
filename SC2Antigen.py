import argparse,os,json;
import sys,re;
import pandas as pd;
import numpy as np;
import time,subprocess;
from multiprocessing import Pool;
from os.path import dirname, abspath
dir = dirname(abspath(__file__));# get parent dir of the scripts
#
#print(time.strftime("%Y-%m-%d %H:%M:%S"))
# =============================================================================
# initial setting for future update
# =============================================================================
#
parser = argparse.ArgumentParser()
parser.add_argument('-hla', metavar = 'inputstring', dest='inputHLA', help='Give arcasHLA calling result fullname');
parser.add_argument('-genotype', metavar = 'input', dest='genotype', help='Give the intron calling result files');
parser.add_argument('-l', metavar = 'pep length', type=str, dest='pepLen', help='Give the peptide files');
parser.add_argument('-t', dest='thread', type=int, help='number of multiple thread number,>0');
parser.add_argument('-o', metavar = 'output', dest='output',help='Give output file fullname');

args = parser.parse_args();
genotype  = args.genotype;
inputHLA  = args.inputHLA;
pepLen= args.pepLen;
output = args.output;
thread = args.thread;

allelenames = dir + '/data/allelenames'
hla_pseudo = dir + '/data/MHC_pseudo.dat'
vfile = dir + '/data/sc2_proteins.fasta'
virus_protein = ["ORF1ab-protein","S-protein","ORF3a-protein","E-protein","M-protein","ORF6-protein",
"ORF7a-protein","ORF7b-protein","ORF8-protein","N-protein","ORF10-protein"]

if thread:
    thread=min(os.cpu_count(), thread);
else:
    thread=8; # defualt thread number setting as 4

if pepLen:
    pepLen_list= pepLen.split(',');
else:
    pepLen_list=['9']; # default using 9mer!

if output[-1]=='_':
    outdir =output;
else:
    outdir =output + '_';

# python version check
pyversion= str(sys.version)[0]
#================================================================
hla_ref = [l.rstrip() for l in open(allelenames)]
if genotype: # get hla results from arcasHLA json file
    hla_json = open(genotype);
    data= json.load(hla_json);
    # genes= list(data.keys());
    # for gene in genes:
        # hla_result= data[gene];
    hla_result= data['A'] + data['B'] + data['C']
    hla_list= ['HLA-' + x.replace('*','') for x in hla_result]
    hla_list= [x.split(':')[0] + ':' + x.split(':')[1] for x in hla_list];
    hla_list= [ x for x in hla_list if x in hla_ref ];
    hla_list = list(set(hla_list)); # combine the identical HLA
    hla_list.sort();
elif inputHLA: #get HLA from user input
    hla_list= inputHLA.split(',')
    if 'HLA' in hla_list[0]:
        hla_list= [x.replace('*','') for x in hla_list]
    else:
        hla_list= ['HLA-' + x.replace('*','') for x in hla_list]; # format the HLA input
    hla_list= [x.split(':')[0] + ':' + x.split(':')[1] for x in hla_list]; # 4 digit setting!
    hla_list= [ x for x in hla_list if x in hla_ref ];
    hla_list = list(set(hla_list)); # combine the identical HLA
    hla_list.sort();
else:
    hla_list=None;
#================================================================
#path, filename = os.path.split(inputFile)
#basename, ext = os.path.splitext(filename)

def single_cmd_run(hla_kmer):
    hla= hla_kmer.split('_')[0];
    pep_len= hla_kmer.split('_')[1];
    print(hla,pep_len)
    # pass hla from main app
    hla_string= hla[:5] + '*' + hla[5:]
    p= subprocess.Popen(['netMHCpan','-a', hla,'-l',pep_len,'-f', vfile],shell=False,stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
    if pyversion =='3':
        netMHCpan_out= p.decode('utf-8');
    elif pyversion =='2':
        netMHCpan_out = p;
    result=[];
    netMHCpan_data= netMHCpan_out.splitlines();
    for line in netMHCpan_data:
        if hla_string in line and 'protein' in line and not 'Number' in line:
            line_data= [x.strip() for x in line.split(' ') if x != ''];
            if float(line_data[12])!=0: # some very strong signal with 0.000
                result.append([line_data[2],hla,line_data[10],line_data[0],str(round(float(line_data[12]),3))]);
            else:
                result.append([line_data[2],hla,line_data[10],line_data[0],str(round(1- float(line_data[11]),8))]);
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Finish with ' + hla + '!');
    return result;

if __name__ == '__main__':
    ###
    hla_kmer= [h + '_' + l for h in hla_list for l in pepLen_list ];
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Doing with Multiprcessing jobs')
    pool=Pool(processes=thread); # defualt thread number setting as 4
    process = pool.map(single_cmd_run, hla_kmer); #multiple process

    summaryList = [item for result_list in process for item in result_list];
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'Summary all the piece file results!')
    pool.close();
    pool.join();

    data_list=[]
    for l in pepLen_list:
        kmer_list= [sublist for sublist in summaryList if len(sublist[0])==int(l)]
        peplist = list(set([x[0] for x in kmer_list]))
        peplist.sort();
        res_tmp= []
        for h in hla_list:
            kmerHLA_list= [sublist for sublist in kmer_list if sublist[1]==h];
            kmerHLA_list.sort();
            protein_list= [x[2] for x in kmerHLA_list];
            pos_list = [x[3] for x in kmerHLA_list];
            hla_value= [x[4] for x in kmerHLA_list];
            res_tmp.append(hla_value);
        res_tmp= [peplist,protein_list,pos_list] + res_tmp;
        kmerSum_list=np.transpose(res_tmp).tolist();
        data_list= data_list + kmerSum_list;
    
    header = ['peptide', 'protein','pos'] + hla_list;
    weakList = [sublist for sublist in data_list if any(x for x in sublist[3:] if float(x) < 2)];
    weakList = [header] + weakList;
    df = pd.DataFrame(weakList[1:],columns=weakList[0]);
    df.to_csv(outdir + 'weak_bind_peptide.txt', index=False, sep='\t');

    strongList= [sublist for sublist in data_list if any(x for x in sublist[3:] if float(x) < 0.5)];
    strongList= [header] + strongList;
    df = pd.DataFrame(strongList[1:],columns=strongList[0]);
    df.to_csv(outdir + 'strong_bind_peptide.txt', index=False, sep='\t');

    result=[['protein'] + hla_list + ['bestRank', 'phbr_Score','weakBinder','strongBinder']]
    for protein in virus_protein:
        proX= [protein];
        protein_list= [sublist for sublist in data_list if sublist[1] == protein];
        print(protein_list)
        protein_weak_num= len([sublist for sublist in protein_list if any(x for x in sublist[3:] if float(x) < 2)]);
        protein_strong_num= len([sublist for sublist in protein_list if any(x for x in sublist[3:] if float(x) < 0.5)]);

        hla_stat=[];
        bestScore={}
        for z in range(3,len(hla_list)+3):
            tmp_min = min([ float(sublist[z]) for sublist in protein_list ]);
            print(tmp_min)
            tmp_list= [sublist for sublist in protein_list if float(sublist[z]) == tmp_min];
            print(tmp_list)
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
    df_result.to_csv(outdir + 'Final-output.txt', index=False, sep='\t');
    #
    print(time.strftime("%Y-%m-%d %H:%M:%S"), 'All jobs done!')
#
#
