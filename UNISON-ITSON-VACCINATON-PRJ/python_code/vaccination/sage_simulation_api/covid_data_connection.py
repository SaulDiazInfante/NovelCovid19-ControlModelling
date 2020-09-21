import requests
import numpy as np
import pandas as pd
import googlemaps
import geopandas as gpd
import seaborn
from shapely.geometry import Point, Polygon
url = 'https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/'
url = url + '455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv'
myfile = requests.get(url)
open('data/conposcovidloc.csv', 'wb').write(myfile.content)

covid_19 = pd.read_csv('covid_19_data_path.csv')
cv0 = pd.get_dummies(covid_19, columns=['OUTCOME1', 'CLIENT_GENDER', 'Age_Group'])
# this creates dummy variables for the categorical variables names
cv0 = cv0.drop(cv0.columns[0], axis=1)
cv0 = cv0.replace(0, np.nan)
cv1 = cv0.groupby(['Reporting_PHU',
                   'Reporting_PHU_Latitude',
                   'Reporting_PHU_Longitude',
                   'Reporting_PHU_City']).count().reset_index()
Counties = gpd.read_file('geo_data_path.shp')
Counties = Counties.sort_values('OFF_NAME')
for i, v in enumerate(city_order):
    city = cv1['Reporting_PHU_City'][v]
    county_city.append(city)

Counties['CITY'] = county_city
# longitude must always come before latitude
geometry = [Point(xy) for xy in zip(cv1['Reporting_PHU_Longitude'],
                                    cv1['Reporting_PHU_Latitude'])]
geo_cv = gpd.GeoDataFrame(cv1, geometry=geometry)
