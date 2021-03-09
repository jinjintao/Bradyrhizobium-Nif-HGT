import pandas as pd
from collections import defaultdict
import sys
from os.path import *
import os
from tqdm import tqdm

#reformat
def reformat(fp):
    df = pd.read_csv(fp,sep=',',header=0)
    df.columns = ['db_id','extant','ancestor','change','genename','function','transition','appeartimes']
    df = df.sort_values(by=['appeartimes', 'genename','transition'], ascending=[False,True,True])
    trans = list(set(df.transition))
    trans.sort()
    columns = ['genename','function','database','db_ID','gain','loss'] + trans
    result_df = pd.DataFrame(columns=columns)
    tqdm.write('formating csv files')
    for i in tqdm(df.index):
        #row = df.loc[row,]
        db_id = df.loc[i,'db_id']
        extant = int(df.loc[i,'extant'])
        ancestor = int(df.loc[i,'ancestor'])
        change = int(df.loc[i,'change'])
        if change > 0:
            gain_or_loss = ','.join(([str(extant),str(ancestor),str(change),'(G)']))
        else:
            gain_or_loss = ','.join(([str(extant),str(ancestor),str(change),('(L)')]))
        genename = df.loc[i,'genename']
        function = df.loc[i,'function']
        tran = df.loc[i,'transition']
        #appeartime = df.loc[i,'appeartimes']
        result_df.loc[db_id,'genename'] = genename
        result_df.loc[db_id,'function'] = function
        result_df.loc[db_id,'db_ID'] = db_id
        result_df.loc[db_id,tran] = gain_or_loss
    #add genes not undergo change
    res = '/'.join(fp.split('/')[:-1])  + '/res/'
    files = os.listdir(res)
    for tran in files:
        kf = open(res+tran)
        KD = {}
        iter_kf = iter(kf)
        for line in iter_kf:
            line_l = line.split('\t')
            key = line_l[0].strip()
            value = ','.join(line_l[1:]).strip()
            KD[key] = value
        for _, arow in result_df.iterrows():
            if pd.isna(result_df.loc[_,tran]):
                ko = result_df.loc[_,'db_ID'].strip()
                added = KD[ko]
                result_df.loc[_,tran] = added
    for i in result_df.index:
        gain = 0
        loss = 0 
        for tran in trans:
            if not pd.isna(result_df.loc[i,tran]):
                if 'G' in result_df.loc[i,tran]:
                    gain += 1
                elif 'L' in result_df.loc[i,tran]:
                    loss += 1
        result_df.loc[i,'gain'] = gain
        result_df.loc[i,'loss'] = loss
    #add database
    if result_df.iloc[0,3].startswith('K'):
        result_df['database']  = 'KEGG'
    elif result_df.iloc[0,3].startswith('C'): 
        result_df['database']  = 'COG'
    elif result_df.iloc[0,3].startswith('T'): 
        result_df['database']  = 'TIGRFAM'
    return result_df

############################### main ###############################
#three lifestyle
#for COG
fpc1 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/FLnonnif-FLnif/gain_and_loss/COG-BDI-FR-CWP/FLnonnif-FLnif_COG-BDI-FR-CWP.csv'
fpc2 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/FLnonnif-Sym/gain_and_loss/COG-BDI-FR-CWP/FLnonnif-Sym_COG-BDI-FR-CWP.csv'
fpc3 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/Sym-FLnif/gain_and_loss/cog/COG-BDI-FR-CWP/Sym-FLnif_COG-BDI-FR-CWP.csv'
fpc4 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/Sym-FLnonnif/gain_and_loss/COG-BDI-FR-CWP/Sym-FLnonnif_COG_BDI-FR-CWP.csv'
reformated_df1 = reformat(fpc1)
reformated_df2 = reformat(fpc2)
reformated_df3 = reformat(fpc3)
reformated_df4 = reformat(fpc4)
writer = pd.ExcelWriter('3COG-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FLnonnif-FLnif')
reformated_df2.to_excel(writer,sheet_name='FLnonnif-Sym')
reformated_df3.to_excel(writer,sheet_name='Sym-FLnif')
reformated_df4.to_excel(writer,sheet_name='Sym-FLnonnif')
writer.save()

#for TIGRFAM
fpt1 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/FLnonnif-FLnif/gain_and_loss/TIGRFAM-BDI-FR-CWP/FLnonnif-FLnif_TIGRFAM-BDI-FR-CWP.csv'
fpt2 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/FLnonnif-Sym/gain_and_loss/TIGRFAM-BDI-FR-CWP/FLnonnif-Sym_TIGRFAM-BDI-FR-CWP.csv'
fpt3 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/Sym-FLnif/gain_and_loss/TIGRFAM/TIGRFAM-BDI-FR-CWP/Sym-FLnif_TIGRFAM_BDI-FR-CWP.csv'
fpt4 = '/home-user/jjtao/Rhizobiales/badirate/3Lifestyle/Sym-FLnonnif/gain_and_loss/TIGRFAM-BDI-FR-CWP/Sym-FLnonnif_TIGRFAM_BDI-FR-CWP.csv'
reformated_df1 = reformat(fpt1)
reformated_df2 = reformat(fpt2)
reformated_df3 = reformat(fpt3)
reformated_df4 = reformat(fpt4)
writer = pd.ExcelWriter('3TIGRFAM-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FLnonnif-FLnif')
reformated_df2.to_excel(writer,sheet_name='FLnonnif-Sym')
reformated_df3.to_excel(writer,sheet_name='Sym-FLnif')
reformated_df4.to_excel(writer,sheet_name='Sym-FLnonnif')
writer.save()

#for KEGG
fpk1 = '/home-user/jjtao/Rhizobiales/kegg_hmmsearch/Brady345/before_rename/e20/badirate/FLnonnif-FLnif/FLnonnif-FLnif_kegg-BDI-FR-CWP/FLnonnif-FLnif_hmmkegg-BDI-FR-CWP.csv'
fpk2 = '/home-user/jjtao/Rhizobiales/kegg_hmmsearch/Brady345/before_rename/e20/badirate/FLnonnif-Sym/FLnonnif-Sym_kegg-BDI-FR-CWP/FLnonnif-Sym_hmmkegg-BDI-FR-CWP.csv'
fpk3 = '/home-user/jjtao/Rhizobiales/kegg_hmmsearch/Brady345/before_rename/e20/badirate/Sym-FLnif/Sym-FLnif_hmmkegg-BDI-FR-CWP/Sym-FLnif_hmmkegg-BDI-FR-CWP.csv'
fpk4 = '/home-user/jjtao/Rhizobiales/kegg_hmmsearch/Brady345/before_rename/e20/badirate/Sym-FLnonnif/Sym-FLnonnif_hmmkegg-BDI-FR-CWP/Sym-FLnonnif_hmmkegg-BDI-FR-CWP.csv'
reformated_df1 = reformat(fpk1)
reformated_df1.to_csv(fpk1.replace('.csv','-reformat.csv'))
reformated_df2 = reformat(fpk2)
reformated_df2.to_csv(fpk2.replace('.csv','-reformat.csv'))
reformated_df3 = reformat(fpk3)
reformated_df3.to_csv(fpk3.replace('.csv','-reformat.csv'))
reformated_df4 = reformat(fpk4)
reformated_df4.to_csv(fpk4.replace('.csv','-reformat.csv'))
writer = pd.ExcelWriter('3KEGG-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FLnonnif-FLnif')
reformated_df2.to_excel(writer,sheet_name='FLnonnif-Sym')
reformated_df3.to_excel(writer,sheet_name='Sym-FLnif')
reformated_df4.to_excel(writer,sheet_name='Sym-FLnonnif')
writer.save()


#two lifestyle
#for COG
fpc1 = '/home-user/jjtao/Rhizobiales/badirate/2Lifestyle/FL-Sym/gainloss/COG-BDI-FR-CWP/FL-Sym-COG-BDI-FR-CWP.csv'
fpc2 = '/home-user/jjtao/Rhizobiales/badirate/2Lifestyle/Sym-FL/gainlosss/COG-BDI-FR-CWP/Sym-FL-COG-BDI-FR-CWP.csv'
reformated_df1 = reformat(fpc1)
reformated_df2 = reformat(fpc2)
writer = pd.ExcelWriter('2COG-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FL-Sym')
reformated_df2.to_excel(writer,sheet_name='Sym-FL')
writer.save()

#for TIGRFAM
fpt1 = '/home-user/jjtao/Rhizobiales/badirate/2Lifestyle/FL-Sym/gainloss/TIGRFAM-BDI-FR-CWP/FL-Sym-TIGRFAM-BDI-FR-CWP.csv'
fpt2 = '/home-user/jjtao/Rhizobiales/badirate/2Lifestyle/Sym-FL/gainlosss/TIGRFAM-BDI-FR-CWP/Sym-FL-TIGRFAM-BDI-FR-CWP.csv'
reformated_df1 = reformat(fpt1)
reformated_df2 = reformat(fpt2)
writer = pd.ExcelWriter('2TIGRFAM-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FL-Sym')
reformated_df2.to_excel(writer,sheet_name='Sym-FL')
writer.save()

#for KEGG
fpk1 = '/home-user/jjtao/Rhizobiales/badirate/2Lifestyle/FL-Sym/gainloss/hmmKEGG-BDI-FR-CWP/FL-Sym-hmmKEGG-BDI-FR-CWP.csv'
fpk2 = '~/Rhizobiales/badirate/2Lifestyle/Sym-FL/gainlosss/hmmKEGG-BDI-FR-CWP/Sym-FL-hmmKEGG-BDI-FR-CWP.csv'
reformated_df1 = reformat(fpt1)
reformated_df2 = reformat(fpt2)
writer = pd.ExcelWriter('2KEGG-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FL-Sym')
reformated_df2.to_excel(writer,sheet_name='Sym-FL')
writer.save()


#for ref_seq
def reformat(fp):
    df = pd.read_csv(fp,sep=',',header=0)
    df.columns = ['genename','extant','ancestor','change','function','transition','appeartimes']
    df = df.sort_values(by=['appeartimes', 'genename','transition'], ascending=[False,True,True])
    trans = list(set(df.transition))
    trans.sort()
    columns = ['genename','function','database','gain','loss'] + trans
    result_df = pd.DataFrame(columns=columns)
    tqdm.write('formating csv files')
    for i in tqdm(df.index):
        #row = df.loc[row,]
        extant = int(df.loc[i,'extant'])
        ancestor = int(df.loc[i,'ancestor'])
        change = int(df.loc[i,'change'])
        if change > 0:
            gain_or_loss = ','.join(([str(extant),str(ancestor),str(change),'(G)']))
        else:
            gain_or_loss = ','.join(([str(extant),str(ancestor),str(change),('(L)')]))
        genename = df.loc[i,'genename']
        function = df.loc[i,'function']
        tran = df.loc[i,'transition']
        #appeartime = df.loc[i,'appeartimes']
        result_df.loc[genename,'genename'] = genename
        result_df.loc[genename,'function'] = function
        result_df.loc[genename,tran] = gain_or_loss
    #add genes not undergo change
    res = '/'.join(fp.split('/')[:-1])  + '/res/'
    files = os.listdir(res)
    for tran in files:
        kf = open(res+tran)
        KD = {}
        iter_kf = iter(kf)
        for line in iter_kf:
            line_l = line.split('\t')
            key = line_l[0].strip()
            value = ','.join(line_l[1:]).strip()
            KD[key] = value
        for _, arow in result_df.iterrows():
            if pd.isna(result_df.loc[_,tran]):
                ko = result_df.loc[_,'genename'].strip()
                added = KD[ko]
                result_df.loc[_,tran] = added
    for i in result_df.index:
        gain = 0
        loss = 0 
        for tran in trans:
            if not pd.isna(result_df.loc[i,tran]):
                if 'G' in result_df.loc[i,tran]:
                    gain += 1
                elif 'L' in result_df.loc[i,tran]:
                    loss += 1
        result_df.loc[i,'gain'] = gain
        result_df.loc[i,'loss'] = loss
    return result_df

#three lifestyle
fpk1 = '/home-user/jjtao/Rhizobiales/refseq/badirate/gain_and_loss/FLnonnif-FLnif/ref-FLnonnif-FLnif.csv'
fpk2 = '/home-user/jjtao/Rhizobiales/refseq/badirate/gain_and_loss/FLnonnif-Sym/ref-FLnonnif-Sym.csv'
fpk3 = '/home-user/jjtao/Rhizobiales/refseq/badirate/gain_and_loss/Sym-FLnif/ref-Sym-FLnif.csv'
fpk4 = '/home-user/jjtao/Rhizobiales/refseq/badirate/gain_and_loss/Sym-FLnonnif/ref-Sym-FLnonnif.csv'
reformated_df1 = reformat(fpk1)
reformated_df1.to_csv(fpk1.replace('ref-','reformat-ref-'))
reformated_df2 = reformat(fpk2)
reformated_df2.to_csv(fpk2.replace('ref-','reformat-ref-'))
reformated_df3 = reformat(fpk3)
reformated_df3.to_csv(fpk3.replace('ref-','reformat-ref-'))
reformated_df4 = reformat(fpk4)
reformated_df4.to_csv(fpk4.replace('ref-','reformat-ref-'))

writer = pd.ExcelWriter('3ref-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FLnonnif-FLnif')
reformated_df2.to_excel(writer,sheet_name='FLnonnif-Sym')
reformated_df3.to_excel(writer,sheet_name='Sym-FLnif')
reformated_df4.to_excel(writer,sheet_name='Sym-FLnonnif')
writer.save()

#two lifestyle

fpk1 = '/home-user/jjtao/Rhizobiales/refseq/badirate/gain_and_loss2/FL-Sym/FL-Sym_ref-BDI-FR-CWP.csv'
fpk2 = '/home-user/jjtao/Rhizobiales/refseq/badirate/gain_and_loss2/Sym-FL/Sym-FL_ref-BDI-FR-CWP.csv'
reformated_df1 = reformat(fpk1)
reformated_df1.to_csv(fpk1.replace('ref-','reformat-ref-'))
reformated_df2 = reformat(fpk2)
reformated_df2.to_csv(fpk2.replace('ref-','reformat-ref-'))
writer = pd.ExcelWriter('2ref-BDI-FR-CWP.xlsx')
reformated_df1.to_excel(writer,sheet_name='FL-Sym')
reformated_df2.to_excel(writer,sheet_name='Sym-FL')
writer.save()
