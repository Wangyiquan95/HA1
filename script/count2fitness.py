import os
import glob
import pandas as pd
import numpy as np

def norm_amp(count_DF,col_name):
    count_DF[col_name]=count_DF[col_name]+1
    count_DF[col_name+'_norm'] = ((count_DF[col_name]) / count_DF.groupby('Amplicon')[col_name].transform('sum'))
    return count_DF
def cal_fit(count_DF,column_name,input):
    count_DF[column_name+'_enri'] = np.log10((count_DF[column_name])/(count_DF[input]))
    silent_DF=count_DF[count_DF.Mutation.str.contains('silent')]
    silent_mean=silent_DF[column_name+'_enri'].mean()
    silent_median = silent_DF[column_name + '_enri'].median()
    count_DF[column_name + '_enri']=(count_DF[column_name + '_enri']-silent_median)
    count_DF=count_DF.rename(columns={column_name + '_enri':column_name + '_fitness'})
    return count_DF


def cal_mean(count_DF,mean,column1,column2,column3):
    count_DF[mean] = (count_DF[column1]+count_DF[column2]+count_DF[column3])/3
    return count_DF

def grab_files_as_Dict(path,suffix):
    #path = 'result/*/*'
    #pattern = '_count.tsv'
    files = glob.glob(path+suffix)
    #filenames = [os.path.basename(i).split(suffix, 1)[0] for i in files]
    file_dict={}
    for file in files:
        filename=os.path.basename(file).split(suffix, 1)[0]
        file_dict[filename]=file
    return file_dict

def main():
    path = 'result/*/*'
    suffix = '_count.tsv'
    file_dict=grab_files_as_Dict(path,suffix)
    DF=pd.DataFrame(columns=['Mutation'])
    for filename,file in file_dict.items():
        df= pd.read_csv(file, header=0, sep= '\t')
        df.columns=['Mutation',filename]
        DF=pd.merge(DF,df,how='outer',on='Mutation')
    # split amplicon and mutation
    DF[['Amplicon','Mutation']] = DF.Mutation.str.split("|", n=1, expand=True)
    DF=DF.fillna(0)
    print(DF)
    amplicon_ls = ['Amp1','Amp2','Amp3']
    contacted_df=pd.DataFrame()
    for amplicon in amplicon_ls:
        amp_DF = DF[DF.Amplicon==amplicon]
        print('writting',amplicon)
        #DF.to_csv('result/HA1_all_count.csv')
        norm_ls=['Plasmid_DMS']
        for col in norm_ls: amp_DF = norm_amp(amp_DF, col)
        #only one mutation & filter low input
        onemut_df = amp_DF[(~amp_DF["Mutation"].str.contains('-'))& (amp_DF['Plasmid_DMS'] >=10)]
        #filter out WT high error mutants
        # onemut_df = onemut_df[onemut_df['HA_DMS_Plasmid_norm'] > onemut_df['HA_WT_plasmid_norm']*8]
        fit_ls=['R_1_DMS','R_2_DMS','R_3_DMS','WT_1_DMS','WT_2_DMS','WT_3_DMS','Y_1_DMS','Y_2_DMS','Y_3_DMS']
        for c in fit_ls:
            onemut_df = norm_amp(onemut_df, c)
            onemut_df = cal_fit(onemut_df, c + '_norm', 'Plasmid_DMS_norm')
        cal_mean(onemut_df, 'R_Fitness', 'R_1_DMS_norm_fitness', 'R_2_DMS_norm_fitness','R_3_DMS_norm_fitness')
        cal_mean(onemut_df, 'Y_Fitness', 'Y_1_DMS_norm_fitness', 'Y_2_DMS_norm_fitness', 'Y_3_DMS_norm_fitness')
        cal_mean(onemut_df, 'WT_Fitness', 'WT_1_DMS_norm_fitness', 'WT_2_DMS_norm_fitness', 'WT_3_DMS_norm_fitness')
        #classify mutation type(silent;missense;nonsense)
        onemut_df['mutation_type'] = 'missense'
        onemut_df['mutation_type'][onemut_df.Mutation.str.contains('silent')]='silent'
        onemut_df['mutation_type'][(~onemut_df.Mutation.str.contains('silent'))&onemut_df.Mutation.str.contains('_')] = 'nonsense'
        contacted_df=pd.concat([onemut_df,contacted_df])
    # DF.to_csv('result/HA1_fit_raw.csv')
    pos_df = contacted_df.Mutation.str.extract('(\d+)')
    one_mut_sortdf = contacted_df.join(pos_df, lsuffix='_caller', rsuffix='_other') \
        .set_index(0) \
        .fillna(0) \
        .sort_index()
    one_mut_sortdf.to_csv('result/HA1_fit.csv')


if __name__ == "__main__":
    main()