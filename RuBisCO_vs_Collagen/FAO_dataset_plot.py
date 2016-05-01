## I used this paper to calculate the mass of different animals:
### chrome-extension://ecnphlgnajanjnkcmbpancdjoidceilk/http://faostat.fao.org/Portals/_Faostat/documents/pdf/ManureManagement.pdf
## in combination with:
## chrome-extension://ecnphlgnajanjnkcmbpancdjoidceilk/http://www.fao.org/climatechange/41521-0373071b6020a176718f15891d3387559.pdf
import pandas as pd
import matplotlib.pyplot as plt

region_mappings_subcat =  {'Africa': 'Africa + (Total)',
                    'Eastern Europe':' --Eastern Europe + (Total)',
                    'Western Europe': ' --Western Europe + (Total)',
                    'Oceania':'Oceania + (Total)',
                    'Northern America': ' --Northern America + (Total)',
                    'Latin America': 'Americas + (Total)',
                    'Indian Subcontinent': ' --Southern Asia + (Total)',
                    'Asia': 'Asia + (Total)',
                    #'Middle east': ''
                    }

region_mappings_numbers =  {'Africa': 'Africa',
                    'Eastern Europe':'Eastern Europe',
                    'Western Europe': 'Western Europe',
                    'Oceania':'Oceania',
                    'Northern America': 'Northern America',
                    'Latin America': 'Americas',
                    'Indian Subcontinent': 'Southern Asia',
                    'Asia': 'Asia',
                    #'Middle east': ''
                    }
region_mappings_numbers_reverse = dict ( (v,k) for k, v in region_mappings_numbers.items() )
region_mappings_subcat_reverse= dict ( (v,k) for k, v in region_mappings_subcat.items() )

animal_categories_mass = ['Asses','Buffaloes','Camelids, other','Camels','Cattle - dairy','Cattle - non-dairy','Chickens - Broilers','Chickens - Layers','Ducks','Goats','Horses','Mules','Swine - market','Swine - breeding','Sheep','Turkeys']
animal_categories_final = ['Asses','Buffaloes','Camelids, other','Camels','Cattle','Chickens','Ducks','Goats','Horses','Mules','Pigs','Sheep','Turkeys']
data = pd.read_csv('/home/yinonbaron/git/RuBisCO_vs_Collagen/Production_Livestock_E_All_Data_(Norm).csv')
weight = pd.read_csv('/home/yinonbaron/git/RuBisCO_vs_Collagen/ipcc_animal_weight.csv',index_col=0)
subcategories = pd.read_csv('/home/yinonbaron/git/RuBisCO_vs_Collagen/dairy_egg_global_data.csv',index_col=0)
filtered_subcat = subcategories.loc[[r for r in region_mappings_subcat.values() if len(r) >0]]
chicken_cat = filtered_subcat.loc[filtered_subcat['element'] == 'Laying (1000 Head)'][filtered_subcat.columns[5:]].reset_index()
chicken_melt = pd.melt(chicken_cat,id_vars=['countries'],value_vars=[x for x in chicken_cat.columns[1:]])
chicken_melt = chicken_melt.loc[chicken_melt['variable']<>'Unnamed: 58']
chicken_melt['variable'] = chicken_melt['variable'].astype('int')
chicken_melt = chicken_melt.replace({'countries':region_mappings_subcat_reverse})
chicken_melt.columns = ['Country','Year','Chickens - Layers']
#chicken_melt = chicken_melt.replace({'Country':{'Indian Subcontinent':'Southern Asia'}})
dairy_cat = filtered_subcat.loc[filtered_subcat['element'] == 'Milk Animals (Head)'][filtered_subcat.columns[5:]].reset_index()
dairy_melt = pd.melt(dairy_cat,id_vars=['countries'],value_vars=[x for x in dairy_cat.columns[1:]])
dairy_melt = dairy_melt.loc[dairy_melt['variable']<>'Unnamed: 58']
dairy_melt['variable'] = dairy_melt['variable'].astype('int')
dairy_melt = dairy_melt.replace({'countries':region_mappings_subcat_reverse})
dairy_melt.columns = ['Country','Year','Cattle - dairy']
#dairy_melt = dairy_melt.replace({'Country':{'Indian Subcontinent':'Southern Asia'}})
world_data = data.loc[data['Country'] == 'World']

region_data = data.loc[data['Country'].isin(region_mappings_numbers.values())][['Country','Item','Year','Value']]
region_data = region_data.loc[region_data['Item'].isin(animal_categories_final)]
region_data_piv = pd.pivot_table(region_data,values='Value',index=['Country','Year'],columns='Item')
region_data_piv = region_data_piv.reset_index()
#region_data_piv = region_data_piv.replace({'Country':region_mappings_numbers_reverse})
region_data_piv = region_data_piv.replace({'Country':{'Southern Asia':'Indian Subcontinent'}})
region_data_piv = region_data_piv.replace({'Country':{'Americas':'Latin America'}})

merged_data = pd.merge(region_data_piv,chicken_melt,how='left',on=['Country','Year'])
merged_data = pd.merge(merged_data,dairy_melt,how='left',on=['Country','Year'])

## replace Asia with Asia- Southern Asia
asia = merged_data.loc[merged_data.Country == 'Asia'].reset_index()
india = merged_data.loc[merged_data.Country == 'Indian Subcontinent'].reset_index()
merged_data.loc[merged_data.Country == 'Asia'][asia.columns[3:]] = asia[asia.columns[3:]]-india[india.columns[3:]]

## replace Latin America with Americas - Northern America
americas = merged_data.loc[merged_data.Country == 'Latin America'].reset_index()
north_america = merged_data.loc[merged_data.Country == 'Northern America'].reset_index()
merged_data.loc[merged_data.Country == 'Latin America'][asia.columns[3:]] = americas[americas.columns[3:]]-north_america[north_america.columns[3:]]

## replace Cattle with Cattle - non-dairy
merged_data['Cattle'] = merged_data['Cattle']- merged_data['Cattle - dairy']

##replace Chickens with Chickens - Chickens - Broilers
merged_data['Chickens'] = merged_data['Chickens']- merged_data['Chickens - Layers']

mass_data = merged_data[merged_data.columns[:-2]].copy()
for animal in animal_categories_final:
    for region in weight.index:
        if animal == 'Cattle':
            mass_data.loc[mass_data.Country == region,animal] = merged_data.loc[merged_data.Country == region]['Cattle']*weight.loc[region]['Cattle - non-dairy'] + merged_data.loc[merged_data.Country == region]['Cattle - dairy']*weight.loc[region]['Cattle - dairy']
        elif animal == 'Chickens':
            mass_data.loc[mass_data.Country == region,animal] = merged_data.loc[merged_data.Country == region]['Chickens']*weight.loc[region]['Chicken - Broilers'] + merged_data.loc[merged_data.Country == region]['Chickens - Layers']*weight.loc[region]['Chicken - Layers']
        elif animal == 'Pigs':
            mass_data.loc[mass_data.Country == region,animal] = merged_data.loc[merged_data.Country == region]['Pigs']*(0.9*weight.loc[region]['Swine - market'] + 0.1*weight.loc[region]['Swine - breeding'])
        elif animal == 'Camelids, other':
            mass_data.loc[mass_data.Country == region,animal] = merged_data.loc[merged_data.Country == region][animal]*weight.loc[region]['Llamas']
        else:
            mass_data.loc[mass_data.Country == region,animal] = merged_data.loc[merged_data.Country == region][animal]*weight.loc[region][animal]


mass_pivot = mass_data.groupby(['Year']).sum()
mass_pivot = mass_pivot.loc[mass_pivot.index<2012]
simple_pivot = mass_pivot[['Cattle','Chickens']].copy()
simple_pivot['Cattle'] = mass_pivot['Cattle'] + mass_pivot['Buffaloes']
simple_pivot['Chickens'] = mass_pivot.sum(axis=1) - simple_pivot['Cattle']
simple_pivot.columns = ['Cattle','Other Livestock']


import wbpy
api = wbpy.IndicatorAPI()
total_population = "SP.POP.TOTL"
dataset = api.get_dataset(total_population, date="1961:2012")
e = dataset.as_dict()
w = pd.Series(e['1W'],name='Humans').reset_index().astype('int')
w = w.set_index(w['index'])['Humans']
w = pd.DataFrame(w)
av_hum_weight = 50
hum_pop = 70018668738
collagen_frac = 0.05
simple_pivot = simple_pivot.merge(w*av_hum_weight,left_index=True,right_index=True)
mass_table = simple_pivot
simple_pivot  = simple_pivot *collagen_frac*1000
mega = pd.DataFrame(pd.Series(data = [hum_pop],index= simple_pivot.index),columns=['Wild Megafauna'])
simple_pivot = simple_pivot.merge(mega*collagen_frac*1000,left_index=True,right_index=True)
simple_pivot =simple_pivot[['Wild Megafauna','Humans','Cattle','Other Livestock']]
simple_pivot=simple_pivot/10**12

## ploting
simple_pivot.plot(kind='area',cmap=plt.get_cmap('viridis'),fontsize=30)
plt.hold(True)
plt.hlines(5*10**13/10**12,1961,2012,linewidth=100.0,colors='grey',alpha=0.7)
plt.hlines(5*10**13/10**12,1961,2012,linewidth=10.0,color='#393933')


plt.text(1971,7.75*10**13/10**12,'RuBisCO mass',weight='bold',size='x-large',color='#393933')
plt.xlabel('Year',fontsize=30,weight='bold')
plt.xticks(range(1961,2012,10))
plt.ylim([0,100])
plt.ylabel('Mass (Mt)',fontsize=30,weight='bold')
