# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:16:30 2024

@author: sanchofernandez
"""

#%% GENERAL SET-UP 

#set a save directory 
savedir_images="J:\\Laura\\Jillytest\\analysis\\images\\take2" #change this to be where analyzed images live 
savedir_tables="J:\\Laura\\Jillytest\\analysis\\tables\\take2" #change this to be where analysis tables live 

import os 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt

#empty dataframe
all_data=pd.DataFrame()

#read in all csv files from output tables folder
for filename in os.listdir(savedir_tables):
    f=pd.read_csv(os.path.join(savedir_tables,filename))
    f['Filename']=filename
    all_data=pd.concat([all_data,f])

#extract relevant information from file name- line, treatment, and antibody used 
fname_split=all_data['Filename'].str.split("_",expand=True)
all_data['line']=fname_split[0] 
all_data['treatment']=fname_split[1]
all_data['antibody']=fname_split[2]
all_data['antibody']=all_data['antibody'].str.strip('bodipy')

#calculate bodipy count per micron2 of astrocyte area
all_data['bodipycountpermicron']=all_data['bodipy_count'].div(all_data['area'], axis=0)

#%% PLOTTING PARAMETERS 

import seaborn as sns 

#set some plotting preferences 
my_palette= ['lightgray', 'violet', 'darkviolet'] #can change to whatever colors 
order=['notx','lxr','soat']  #can change to whatever order you like 
sns.set_theme("talk") #my personal preference 
sns.set_theme(style="ticks") #i like ticks 



  #%% PLOT S100B 

#working directory   
os.chdir("J:\\Laura\\Jillytest\\analysis\\take2")

#useful to split if you have different antibodies 
all_data_s100b= all_data[all_data['antibody'].isin(['S100'])]

#replace NaNs with 0s 
all_data_s100b=all_data_s100b.fillna(0)
all_data_s100b.reset_index(inplace=True, drop=True) 



#-------------BODIPY COUNT PER ASTRO---------------------------
#plot bodipy count per astrocyte 

sns.boxplot(data=all_data_s100b,x='line',y='bodipy_count',
     hue='treatment', hue_order=order, whis=[0,100],palette=my_palette)

#plot the points (all images per line per condition)
sns.stripplot(data=all_data_s100b,x='line', y='bodipy_count', dodge=True,
        hue= 'treatment', hue_order=order, size=6, color= ".3")
plt.title('Number of Lipid droplets per astrocyte')          
sns.despine() 
plt.savefig('s100B_bodipycount_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()
plt.figure().clf()


  #----------BODIPY MEAN AREA PER ASTRO-----------------
#plot mean area of bodipy puncta per astrocyte   
  
sns.boxplot(data=all_data_s100b,x='line',y='mean_bodipy_area',
     hue='treatment', hue_order=order, whis=[0,100],palette=my_palette)

 #plot the points 
sns.stripplot(data=all_data_s100b,x='line', y='mean_bodipy_area', dodge=True,
        hue= 'treatment', hue_order=order, size=6, color= ".3")
plt.title('Lipid droplet mean area per astrocyte')          
sns.despine() 
plt.savefig('s100B_bodipyarea_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()
plt.figure().clf()

  #----------BODIPY PERCENT AREA PER ASTRO-----------------
#plot percent bodipy area per astrocyte

sns.boxplot(data=all_data_s100b,x='line',y='percent_bodipy_area',
     hue='treatment', hue_order=order, whis=[0,100],palette=my_palette)

 #plot the points 
sns.stripplot(data=all_data_s100b,x='line', y='percent_bodipy_area', dodge=True,
        hue= 'treatment', hue_order=order, size=6, color= ".3")
plt.title('Lipid droplet percent area per astrocyte')          
sns.despine() 
plt.savefig('s100B_bodipypercentarea_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()
plt.figure().clf()

#%% Plot by group (CT vs AD)

group_mapping={'UCI107': 'CT',
               'UCI34':'CT',
               'UCI48': 'CT',
               'UCI70':'CT',
               'UCI39': 'AD',
               'UCI66': 'AD',
               'UCI7': 'AD'}

all_data_s100b=all_data_s100b.assign(group=all_data_s100b.line.map(group_mapping))

#%% s100b staining by group (CT vs AD)

order2=['CT','AD']
x_order=['notx', 'lxr', 'soat']

palette2=['lightgray','orange']

#bodipy count 
sns.boxplot(data=all_data_s100b, x='treatment', y='bodipy_count',order=x_order,
            hue='group', hue_order=order2, whis=[0,100],
            palette=palette2)

ax=sns.stripplot(data=all_data_s100b,x='treatment', y='bodipy_count', dodge=True,order=x_order,
                    hue='group', hue_order=order2,size=6, color= ".3",
                    )
h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])

plt.title('Lipid droplets per astrocyte')
plt.savefig('ByGroup_s100b_bodipycount_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()

#bodipy mean area 
sns.boxplot(data=all_data_s100b, x='treatment', y='mean_bodipy_area',order=x_order,
            hue='group', hue_order=order2, whis=[0,100],
            palette=palette2)

ax=sns.stripplot(data=all_data_s100b,x='treatment', y='mean_bodipy_area', dodge=True,order=x_order,
                    hue='group', hue_order=order2,size=6, color= ".3",
                    )
h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])


plt.title('Lipid droplet mean area per astrocyte')
plt.savefig('ByGroup_s100b_bodipymeanarea_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()

#bodipy percent area 
 
sns.boxplot(data=all_data_s100b,x='treatment', y='percent_bodipy_area',order=x_order,
                     hue='group', hue_order=order2,
                     whis=[0,100], 
                     palette=palette2,
                     )
#plot the points 
ax=sns.stripplot(data=all_data_s100b,x='treatment', y='percent_bodipy_area', dodge=True,order=x_order,
                    hue='group', hue_order=order2,size=6, color= ".3",
                    )

h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])

plt.title('Lipid droplet percent area per astrocyte')
plt.savefig('ByGroup_s100b_percent_bodipy_area_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()

#bodipy count per micron of astrocyte
sns.boxplot(data=all_data_s100b,x='treatment', y='bodipycountpermicron',order=x_order,
                     hue='group', hue_order=order2,
                     whis=[0,100], 
                     palette=palette2,
                     )
    #plot the points 
ax= sns.stripplot(data=all_data_s100b,x='treatment', y='bodipycountpermicron', dodge=True,order=x_order,
                    hue='group', hue_order=order2,size=6, color= ".3",
                    )

h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])

plt.title('Lipid droplets per micron astrocyte')
plt.savefig('ByGroup_s100b_bodipycountpermicron_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()


#%%
# no treatment only 
#bodipy count 
sns.boxplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'], x='treatment', y='bodipy_count',
            hue='group', hue_order=order2, whis=[0,100],
            palette=palette2)

ax=sns.stripplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'], x='treatment', y='bodipy_count', dodge=True,
                    hue='group', hue_order=order2,size=6, color= ".3",
                    )

h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])

plt.title('Lipid droplets per astrocyte')
plt.savefig('ByGroup_NoTx_s100b_bodipycount_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()

#bodipy mean area 
sns.boxplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'], x='treatment', y='mean_bodipy_area',
            hue='group', hue_order=order2, whis=[0,100],
            palette=palette2)

ax=sns.stripplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'], x='treatment', y='mean_bodipy_area', dodge=True,
                    hue='group', hue_order=order2,size=6, color= ".3",
                    )

h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])

plt.title('Lipid droplet mean area per astrocyte')
plt.savefig('ByGroup_NoTx_s100b_bodipymeanarea_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()

#bodipy percent area 
 
sns.boxplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'],x='treatment', y='percent_bodipy_area',
                     hue='group', hue_order=order2,
                     whis=[0,100], 
                     palette=palette2,
                     )
    #plot the points 
ax=sns.stripplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'],x='treatment', y='percent_bodipy_area', dodge=True,
              hue='group', hue_order=order2,size=6, color= ".3",     
              )

h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])

plt.title('Lipid droplet percent area per astrocyte')
plt.savefig('ByGroup_NoTx_s100b_percent_bodipy_area_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()

#bodipy count per micron astrocyte
sns.boxplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'],x='treatment', y='bodipycountpermicron',
                     hue='group', hue_order=order2,
                     whis=[0,100], 
                     palette=palette2,
                     )
    #plot the points 
ax=sns.stripplot(data=all_data_s100b.loc[all_data_s100b['treatment']=='notx'],x='treatment', y='bodipycountpermicron', dodge=True,
              hue='group', hue_order=order2,size=6, color= ".3",     
              )

h,l= ax.get_legend_handles_labels()
ax.legend(h[0:2], l[0:2])

plt.title('Lipid droplets per micron astrocyte')
plt.savefig('ByGroup_NoTx_s100b_bodipycountpermicron_boxplot.pdf',format='pdf', bbox_inches="tight")
plt.show()

#%% STATISTICS 

import statsmodels.api as sm
from statsmodels.formula.api import ols 
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison



#-----------------ANOVAs--------------------------------
#bodipy count 
model = ols ('bodipy_count ~ C(treatment) + C(group) +\
             C(treatment):C(group)', 
             data=all_data_s100b).fit()
result=sm.stats.anova_lm(model, type=2)

print("Bodipy count\n", result, "\n")    
result.to_csv(os.path.join('s100b_bodipy_count_ANOVA.csv'))

#bodipy mean area 
model2 = ols ('mean_bodipy_area ~ C(treatment) + C(group) +\
             C(treatment):C(group)', 
             data=all_data_s100b).fit()
result2=sm.stats.anova_lm(model2, type=2)

print("Mean bodipy area\n", result2,"\n")    
result2.to_csv(os.path.join('s100b_mean_bodipy_area_ANOVA.csv'))

#percent_bodipy_area
model3 = ols ('percent_bodipy_area ~ C(treatment) + C(group) +\
             C(treatment):C(group)', 
             data=all_data_s100b).fit()
result3=sm.stats.anova_lm(model3, type=2)

print("Percent bodipy area\n", result3, "\n")    
result3.to_csv(os.path.join('s100b_percent_bodipy_area_ANOVA.csv'))

#bodipy count per micron astrocyte
model4 = ols ('bodipycountpermicron ~ C(treatment) + C(group) +\
             C(treatment):C(group)', 
             data=all_data_s100b).fit()
result4=sm.stats.anova_lm(model4, type=2)

print("Bodipy count per micron\n",result4,"\n")    
result3.to_csv(os.path.join('s100b_bodipycountpermicron_ANOVA.csv'))

#--------------------post-hoc tukey's tests-----------------------------------
all_data_s100b['combo']=all_data_s100b['treatment'] + all_data_s100b['group']

#bodipy count 
bodipy_count_posthoc_results=pairwise_tukeyhsd(endog=all_data_s100b['bodipy_count'],groups=all_data_s100b['combo'],alpha=0.05)
print("Bodipy_count\n", bodipy_count_posthoc_results, "\n")

df_bodipycount_results=pd.DataFrame(data=bodipy_count_posthoc_results._results_table.data[1:], 
                                    columns=bodipy_count_posthoc_results._results_table.data[2:])
df_bodipycount_results.to_csv(os.path.join('s100b_bodipy_count_tukey.csv'))


#bodipy mean area
bodipy_meanarea_posthoc_results=pairwise_tukeyhsd(endog=all_data_s100b['mean_bodipy_area'],groups=all_data_s100b['combo'],alpha=0.05)
print("Bodipy mean area\n", bodipy_meanarea_posthoc_results, "\n")

df_bodipy_meanarea_posthoc_results=pd.DataFrame(data=bodipy_meanarea_posthoc_results._results_table.data[1:], 
                                    columns=bodipy_meanarea_posthoc_results._results_table.data[2:])
df_bodipy_meanarea_posthoc_results.to_csv(os.path.join('s100b_bodipy_mean_area_tukey.csv'))


#bodipy percent area
bodipy_percentarea_posthoc_results=pairwise_tukeyhsd(endog=all_data_s100b['percent_bodipy_area'],groups=all_data_s100b['combo'],alpha=0.05)
print("Bodipy percent area\n", bodipy_percentarea_posthoc_results, "\n")

df_bodipy_percentarea_posthoc_results=pd.DataFrame(data=bodipy_percentarea_posthoc_results._results_table.data[1:], 
                                    columns=bodipy_percentarea_posthoc_results._results_table.data[2:])
df_bodipy_percentarea_posthoc_results.to_csv(os.path.join('s100b_bodipy_percent_area_tukey.csv'))


#bodipy count per astro 
bodipy_countperastro_posthoc_results=pairwise_tukeyhsd(endog=all_data_s100b['bodipycountpermicron'],groups=all_data_s100b['combo'],alpha=0.05)
print("Bodipy count per micron astrocyte area\n", bodipy_countperastro_posthoc_results, "\n")

df_bodipy_countperastro_posthoc_results=pd.DataFrame(data=bodipy_countperastro_posthoc_results._results_table.data[1:], 
                                    columns=bodipy_countperastro_posthoc_results._results_table.data[2:])
df_bodipy_countperastro_posthoc_results.to_csv(os.path.join('s100b_bodipy_count_permicron_tukey.csv'))







