"""The module that splits the data in to training and testing parts."""

import pandas as pd
import os

# Function to remove rows with all zeros (excluding specified columns by index)
def remove_all_zero_rows_by_index(dataframe):

    # firstly filter out the rows with no data or non-numerical values
    float_dataframe = dataframe.select_dtypes(include=['float64'])
    
    # Create a mask to filter out rows where all the specified columns (except those ignored) are zero
    mask = (float_dataframe != 0).any(axis=1)
    filtered_dataframe = dataframe[mask]
    filtered_dataframe = filtered_dataframe.dropna(how='all')
    
    return filtered_dataframe

# Load the first dataset (circular polysolar)
circular = pd.read_csv('weather_data/circular_polysolar_kent_cellwise_irradiance.csv')
circular.drop(['Unnamed: 0', 'hour'], inplace=True, axis=1)

# Load the second dataset (ninja pv polysolar)
ninja = pd.read_csv('weather_data/ninja_pv_polysolar_kent.csv', skiprows=3)
ninja.drop(['time', 'local_time', 'temperature'], inplace=True, axis=1)

# Remove rows where all values are zero in both datasets
circular_non_zero = remove_all_zero_rows_by_index(circular)
ninja_non_zero = remove_all_zero_rows_by_index(ninja)

# Split and get both datasets
circular_train = circular_non_zero.sample(frac=0.8, random_state=40)
circular_test = circular_non_zero.drop(circular_train.index)

ninja_train = ninja_non_zero.sample(frac=0.8, random_state=40)
ninja_test = ninja_non_zero.drop(circular_train.index)


