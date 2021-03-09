import os
import collections
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord,CircularGraphicRecord
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# assign color for gene in nif island
color = {'fabG' : '#FFFFCC', #Catalyzes the NADPH-dependent reduction of beta-ketoacyl-ACP substrates to beta-hydroxyacyl-ACP products, the first reductive step in the elongation cycle of fatty acid biosynthesis.
         'nifA' : '#33CC99', #synthesize nitrogenase
         'nifB' : '#33CC99',
         'nifD' : '#33CC99',
         'nifK' : '#33CC99',
         'nifH' : '#33CC99',
         'nifQ' : '#33CC99',
         'nifE' : '#33CC99',
         'nifN' : '#33CC99',
         'nifT' : '#33CC99',
         'nifU-like protein' : '#33CC99',
         'nifV' : '#33CC99',
         'nifW' : '#33CC99',
         'nifX' : '#33CC99',
         'nifZ' : '#33CC99',
         'fixA' : '#99CC00', # transfer electron to nitrogenase
         'fixB' : '#99CC00',
         'fixC' : '#99CC00',
         'fixX' : '#99CC00',
         'fixL' : '#99CC00',
         'fixJ' : '#99CC00',
         'sufB' : '#FFD700', # FeS cluster assembly protein
         'sufC' : '#FFD700',
         'sufS' : '#FFD700',
         'sufD' : '#FFD700',
         'sufE' : '#FFD700',
         'sufX' : '#FFD700',
         'iscA' : '#FFD700',
         'iscS' : '#FFD700',
         'hcaC' : '#FFFF00', # like ferredoxin
         'fdx' : '#FFFF00', # ferredoxin  transfer electron 
         'per' : '#FFC0CB', # synthesize sugar
         'cysE' : '#FF7F50', # synthesizes L-cysteine
         'paaD' : '#CD853F', # phenylacetate degradation, part of Aromatic compound metabolism
         'bolA': '#66CCFF', # DNA-binding transcriptional regulator BolA
         'grxD' : '#FFD700', # Monothiol glutaredoxin involved in the biogenesis of iron-sulfur clusters.
         'irr' : '#66CCFF', #Acts as a global negative controlling element, employing Fe2+ as a cofactor to bind the operator of the repressed genes.
         'hyaA' : '#99FFFF',
         'hyaB' : '#99FFFF',
         'hyaC' : '#99FFFF',
         'hyaD' : '#99FFFF',
         'hyaF' : '#99FFFF',
         'hypA' : '#99FFFF',
         'hypB' : '#99FFFF',
         'hypBA1': '#99FFFF',
         'hypC' : '#99FFFF',
         'hypD' : '#99FFFF',
         'hypE' : '#99FFFF',
         'hypF' : '#99FFFF',
         'hypV' : '#99FFFF',
         'hupV' : '#99FFFF',
         'glbN' : '#CC0000',
         'glbO' : '#CC0000',
         'hspQ' : '#1E90FF',
         'modA' : '#CC99FF',
         'modB' : '#CC99FF',
         'modC' : '#CC99FF',
         'modD' : '#CC99FF',
         'modE' : '#CC99FF',
         'nodA' : '#F08080',
         'nodB' : '#F08080',
         'nodC' : '#F08080',
         'nodI' : '#F08080',
         'nodJ' : '#F08080',         
          }

# get nifH locus, several genomes have more than one nifH gene, their genoem was named like species_name|1
def get_nifH_id(df,spe_order):
    nifH_id = defaultdict(lambda:list)
    for s in spe_order:
        ns = s.replace('.','')
        if ns not in df.columns.tolist():
            print(s,'not in merged_hmm_info')
        if pd.isna(df.loc['K02588',ns]):
            print(s,'doesn\'t have nifH gene')
        else:
            locus = df.loc['K02588',ns].split(',')
            dup = 0
            for locu in locus:
                if len(locus) > 1:
                    dup = dup +1
                    spedup = s + '|' + str(dup)
                else:
                    spedup = s     
                nifH_id[spedup] = locu.split('|')[-1]
    return nifH_id

#get gene name
def get_locas2ko(df):
    kotbl = [line.strip() for line in open("/home-user/thliao/data/protein_db/kegg/ko_info.tab", 'r')]
    koname = {}
    for line in kotbl:
        ko = line.split('\t')[0].strip('ko:')
        fullname = line.split('\t')[1]
        name = fullname.split(';')[0]
        koname[ko] = name    
    locas2ko = defaultdict()
    for i in df.index.tolist():
        row = df.loc[i]
        for info in row:
            if not pd.isna(info):
                if ',' in info:
                    para = info.split(',')
                    for a in para:
                        spe = a.split('|')[0]
                        locus = a.split('|')[1]          
                        locas2ko[locus]= i
                else:
                    spe = info.split('|')[0]
                    locus = info.split('|')[1]
                    locas2ko[locus]= i
    return koname,locas2ko

def get_correct_name(key, indir):
    filelist = os.listdir(indir)
    if key + '.gbk' in filelist:
        fh = indir + key + '.gbk'
        return fh
    else:    
        key = key.replace('sp', 'sp.') 
        if key + '.gbk' in filelist:
            fh = indir + key + '.gbk'
            return fh
        else:
            key = '_'.join(key.split('_')[:-1]) 
            fh = indir + key + '.gbk'
            return fh
        
def trans_features_to_location(CDS):
    '''transform features to location'''
    location = defaultdict(lambda: defaultdict())
    for feature in CDS:
        locus_tag = feature.qualifiers['locus_tag'][0]
        location[locus_tag]['start'] = int(feature.location.start)        
        location[locus_tag]['end'] = int(feature.location.end)
        location[locus_tag]['strand'] = int(feature.location.strand)
        if locus_tag in locas2ko.keys():
            location[locus_tag]['label'] = koname[locas2ko[locus_tag]]
        else:
            location[locus_tag]['label'] = 'unknown'
    return location

def select_contig(nifH_id):
    locations = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    tqdm.write("reading gbk files")
    for spe,nifhid in nifH_id.items():
            contig = []
            #spe = 'Bradyrhizobium_sp._AT1'
            #nifHid = 'SE92_10630'
            filename = spe.split('|')[0]
            fh = get_correct_name(filename, indir)
            record = list(SeqIO.parse(fh,format='genbank'))
            for n in range(len(record)): 
                cds_num = 0         
                for (index, feature) in enumerate(record[n].features):
                    if feature.type == feature_type:
                        cds_num = cds_num +1
                        if feature.qualifiers['locus_tag'][0] == nifhid: #find the contig
                            contig = record[n]
                            locate = cds_num 
            CDS = []                                   
            for (index_, feature_) in enumerate(contig.features):
                if feature_.type == feature_type:
                    CDS.append(feature_)
            if len(CDS) <= 180:
                location = trans_features_to_location(CDS)
            else:
                if locate < 80:
                    CDS = CDS[:locate+80]
                elif locate > (len(CDS)-80):
                    CDS = CDS[locate-80:]
                else:
                    CDS = CDS[locate-80:locate+80]
                location = trans_features_to_location(CDS)   
            locations[spe] = location
            print(spe, len(locations[spe]))
    return locations

#reverse chromosome and redefine the starting point
def reverse_contig(locations):
    locations_sort = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for spe in locations.keys():
        #spe = 'Bradyrhizobium_sp._S23321'
        nifhid = nifH_id[spe]
        origin = list(locations[spe].values())[0]['start']
        terminal = list(locations[spe].values())[-1]['end']
        if locations[spe][nifhid]['strand'] == 1:
            print(spe,'forward', origin)
            for locus_tag in locations[spe].keys():
                locations_sort[spe][locus_tag]['start'] = locations[spe][locus_tag]['start']- origin
                locations_sort[spe][locus_tag]['end'] = locations[spe][locus_tag]['end'] - origin
                locations_sort[spe][locus_tag]['strand'] = locations[spe][locus_tag]['strand']
                locations_sort[spe][locus_tag]['label'] = locations[spe][locus_tag]['label']
        else:   
            print(spe,'backward',terminal)
            for locus_tag in sorted(locations[spe].keys(),reverse = True):
                locations_sort[spe][locus_tag]['start'] = terminal - locations[spe][locus_tag]['end']
                locations_sort[spe][locus_tag]['end'] = terminal - locations[spe][locus_tag]['start'] 
                locations_sort[spe][locus_tag]['strand'] = -(locations[spe][locus_tag]['strand'])
                locations_sort[spe][locus_tag]['label'] = locations[spe][locus_tag]['label']
    return locations_sort

#cut genome region based on gene numbers 
def cut_contig_by_genenumber(locations_sort,lf,rf):
    cut_n = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for spe in locations_sort.keys(): 
        nifhid = nifH_id[spe]
        count = 0
        for locus_tag in locations_sort[spe].keys():
            count = count +1
            if locus_tag == nifhid:
                locate = count                                         
        for (index, locus_tag) in enumerate(locations_sort[spe]):  
            if locate < lf:
                left = 0
                right = locate+rf
            elif locate > (len(locations_sort[spe].keys())-rf):
                left = locate-lf
                right = len(locations_sort[spe].keys())
            else:
                left = locate-lf
                right = locate+rf      
            for k,v in list(locations_sort[spe].items())[left:right]:
                cut_n[spe][k] = v
    return cut_n

#draw picture 
def plot(start_sites,cut,length,group_name,fig_size):
    fig, ax = plt.subplots(len(start_sites),1,figsize=(fig_size,len(start_sites)*1.5))
    i = 0
    for spe in start_sites:
        fea = []
        if start_sites[spe] == 0:
            origin = list(cut[spe].values())[0]['end']
        else: 
            origin = start_sites[spe]
        for locus_tag in cut[spe].keys():
            s = cut[spe][locus_tag]['end']
            t = cut[spe][locus_tag]['end']
            if t-origin > length or s < origin:
                pass
            else:
                a = cut[spe][locus_tag]['start'] - start_sites[spe]
                z = cut[spe][locus_tag]['end'] - start_sites[spe]
                d = cut[spe][locus_tag]['strand']
                l = cut[spe][locus_tag]['label'].split(',')[0]
                if l == 'unknown':
                    c = "#FFFFCC"
                    gf = GraphicFeature(start=a,
                                        end=z,
                                        strand=d,
                                        color=c)          
                elif l not in color.keys():            
                    c = "#FFFFCC"
                    gf = GraphicFeature(start=a,
                                        end=z,
                                        strand=d,
                                        color=c,
                                        label = l 
                                        )
                
                else:
                    c = color[l]           
                    gf = GraphicFeature(start=a,
                                        end=z,
                                        strand=d,
                                        color=c,
                                        label=l)
                fea.append(gf)
        record = GraphicRecord(sequence_length=length, features=fea)
        for f in record.features:
            f.data['fixed_level'] = 0
        ax[i].set_xlabel(spe.replace('_',' '))
        record.plot(ax = ax[i]) 
        i = i+1
    plt.savefig(f'{group_name}.pdf')
    return('done')

############################################# main #############################################
## species list
spe_order = PB + flnif + flwithsym + SymBasal + sym

## nifH id
annotation = '/home-user/jjtao/Rhizobiales/kegg_hmmsearch/Brady345/after_rename_3Div/1e10/merged_hmm_info.tab'
protein_info = pd.read_csv(annotation,sep='\t',header=0,index_col=0,low_memory=False)
nifH_id = get_nifH_id(protein_info,spe_order)
koname,locas2ko = get_locas2ko(protein_info)

##read gbk and get contig 
indir = '/home-user/jjtao/Rhizobiales/data/Rhizobiales/gbk/' 
feature_type = 'CDS'
locations = select_contig(nifH_id)

locations_sort = reverse_contig(locations)
lf = 65 #gene numbers in left flank of nifH
rf = 50 #gene numbers in right flank of nifH
cut_n = cut_contig_by_genenumber(locations_sort,lf,rf)

###################################### refine start site ######################################
for s in cut_n.keys():
    labels = []
    found = False
    for l in cut_n[s].keys():
        label = cut_n[s][l]['label']
        if 'fabG' in label:
            found = True
            print(s,l,cut_n[s][l]['start'])
            break
    if not found:
        print(s,' not find')

################## draw picture for PB group ##################
PB = {'Bradyrhizobium_AUGA_SZCCT0283|1':-1000,
                #'Bradyrhizobium_AUGA_SZCCT0283|2':0, the same to |1
                #'Bradyrhizobium_sp_STM_3843':,
                #'Bradyrhizobium_oligotrophicum_S58|1':62758, the same to |2
                #'Bradyrhizobium_oligotrophicum_S58|2':19999,
                #'Bradyrhizobium_denitrificans_SZCCT0094':69070,
                'Bradyrhizobium_sp_BTAi1|1':20998,
                #'Bradyrhizobium_sp_BTAi1|2':76529,#the same island with|2
                'Bradyrhizobium_sp_ORS_278|1':20438,
                #'Bradyrhizobium_sp_ORS_285|1':65693,
                #'Bradyrhizobium_sp_ORS_285|2':20969,
                #'Bradyrhizobium_sp_ORS_375':,
                #'Bradyrhizobium_sp_STM_3809': 
                }
plot(PB,cut_n,55000,'PB',56)

################## draw picture for flnif group ##################
flnif = {
      'Bradyrhizobium_sp._S23321|1':31195,
      #'Bradyrhizobium_sp._S23321|2':62979,        
      #'Bradyrhizobium_japonicum_22':28505,
      'Bradyrhizobium_sp_DOA9_chromosome':35549,
      'Bradyrhizobium_guangxiense_CCBAU_53363_CP022219':29168,      
      'Bradyrhizobium_iriomotense_SZCCT0007':-3225,
      'Bradyrhizobium_sacchari_p9-20_LWIG01000010':-3225,
      'Bradyrhizobium_sp._AT1':46651,
      #'Bradyrhizobium_sp._BM-T':35390,
      #'Bradyrhizobium_sp_DOA9_chromosome|2':39546
      }
plot(flnif,cut_n,50000,'flnif',56)

################# draw picture for flwithsym group ##################
flwithsym = {'Bradyrhizobium_mercantei_SEMIA_6399':18709,
         'Bradyrhizobium_yuanmingense_P10_130':0,
         'Bradyrhizobium_liaoningense_CCBAU_83689':0
         }
plot(flwithsym,cut_n,73000,'flwithsym',60)

################## draw picture for symbiotic group ##################
lf = 80 #gene numbers in left flank of nifH
rf = 80 #gene numbers in right flank of nifH
cut_n = cut_contig_by_genenumber(locations_sort,lf,rf)
sym = {
             'Bradyrhizobium_sp_DOA1':71688,            #symbasal
             'Bradyrhizobium_stylosanthis_BR_446':70911,#symbasal
             'Bradyrhizobium_sp_NAS96.2':0,             #sym that colse to fl
             'Bradyrhizobium_sp_CCGE-LA001':53248,          #sym that colse to fl
             'Bradyrhizobium_sacchari_p9-20_LWIG01000056':25891,
              #'Bradyrhizobium_elkanii_USDA_76': 24194,
             'Bradyrhizobium_genosp_SA-4_str_CB756':40371,
             'Bradyrhizobium_sp_DOA9_plasmid': 50865,
             'Bradyrhizobium_guangxiense_CCBAU_53363_CP022220': 52231,
             'Bradyrhizobium_diazoefficiens_USDA_110': 38842,
            }
plot(sym,cut_n,73000,'sym',60)


################## draw picture for outgroup ##################
outgroup = ['Azorhizobium_doebereinerae_UFLA1-100':0,
              'Azorhizobium_caulinodans_ORS_571':0,
              'Rhodopseudomonas_palustris_DX-1':0]
plot(outgroup,cut_n,80000,'outgroup',60)

################## draw picture for flwithsym&sym ##################
fl_sym_compare = {'Bradyrhizobium_sp_NAS96.2':0,
                  'Bradyrhizobium_mercantei_SEMIA_6399':18709,
                  'Bradyrhizobium_yuanmingense_P10_130':0,
                  'Bradyrhizobium_sp_CCGE-LA001':53248
                }
plot(fl_sym_compare,cut_n,73000,'fl_sym_compare',60)
