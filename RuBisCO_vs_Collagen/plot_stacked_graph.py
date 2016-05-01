# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 18:03:11 2016

@author: yinonbaron
"""

import pandas as pd
import pygal
data = pd.read_csv('/home/yinonbaron/git/RuBisCO_vs_Collagen/biomass_table.csv')

line_chart = pygal.StackedLine(fill=True)
line_chart.title = 'RuBisCO vs. Collagen Throughout History'
line_chart.x_labels = data['Time'].loc[5::-1].astype('int')

cattle_fraction =          663149950000.0/988565250000.0;
sheep_fraction =            77871734000.0/988565250000.0;
pig_fraction =              95835247000.0/988565250000.0;
horses_fraction =           53797867000.0/988565250000.0;
buffalo_camel_fraction =    71653831000.0/988565250000.0;
sheep_fraction =            77871734000.0/988565250000.0;
chicken_fraction =          26082000000.0/988565250000.0;



line_chart.add('Wild Megafauna',data['Wild_mega'].loc[5::-1]*6)
line_chart.add('Humans',data['Humans'].loc[5::-1]*6)
line_chart.add('Cattle',data['Livestock'].loc[5::-1]*cattle_fraction*6)
line_chart.add('Sheep',data['Livestock'].loc[5::-1]*sheep_fraction*6)
line_chart.add('Pigs',data['Livestock'].loc[5::-1]*pig_fraction*6)
line_chart.add('Chicken',data['Livestock'].loc[5::-1]*chicken_fraction*6)
line_chart.add('Other Livestock',data['Livestock'].loc[5::-1]*(1-cattle_fraction-sheep_fraction-pig_fraction-chicken_fraction)*6)
line_chart.render_to_png('/home/yinonbaron/git/RuBisCO_vs_Collagen/plot.png')

fig =figure()
hold(True)
total_sum = data['Wild_mega'].loc[5::-1]*6 + data['Humans'].loc[5::-1]*6 + data['Livestock'].loc[5::-1]*6
plot(data['Time'].loc[5::-1],total_sum,'o')
ax1 = fig.add_subplot(111)
ax1.fill_between(data['Time'].loc[5::-1], 0, data['Wild_mega'].loc[5::-1]*6, facecolor="#ff0000", alpha=.7)
ax1.fill_between(data['Time'].loc[5::-1], data['Wild_mega'].loc[5::-1]*6, data['Wild_mega'].loc[5::-1]*6+data['Humans'].loc[5::-1]*6, facecolor="#3333cc", alpha=.7)
'''plot(data['Time'].loc[5::-1],data['Wild_mega'].loc[5::-1]*6)
plot(data['Time'].loc[5::-1],data['Humans'].loc[5::-1]*6)
plot(data['Time'].loc[5::-1],data['Livestock'].loc[5::-1]*cattle_fraction*6)
plot(data['Time'].loc[5::-1],data['Livestock'].loc[5::-1]*sheep_fraction*6)
plot(data['Time'].loc[5::-1],data['Livestock'].loc[5::-1]*pig_fraction*6)
plot(data['Time'].loc[5::-1],data['Livestock'].loc[5::-1]*chicken_fraction*6)
plot(data['Time'].loc[5::-1],data['Livestock'].loc[5::-1]*(1-cattle_fraction-sheep_fraction-pig_fraction-chicken_fraction)*6)
'''
