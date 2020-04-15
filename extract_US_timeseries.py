import os
import pandas as pd

US_CONFIRMED_CSV_PATH = '/home/jmc678/covid_data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv'
US_DEATHS_CSV_PATH = '/home/jmc678/covid_data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv'

def extract_csv(df, file_prefix, main_columns):
    main_df = df[main_columns]
    date_set = set(df.columns) - set(main_columns)
    new_files = []
    for date in date_set:
        file_name = file_prefix + '_' + date.replace('/','_') + '.csv'
        if os.path.exists(file_name):
            continue
        date_df = main_df.copy()
        date_df['count'] = pd.Series(df[date], index = date_df.index)
        date_df['date'] = pd.Series([date] * date_df.shape[0], index = date_df.index)
        date_df.to_csv(file_name,index_label=False, index=False)

        new_files.append(file_name)

    return new_files

def extract_sql_command(file_name, main_columns, table_name):
    return "copy {table_name}({column_list}) from '{file_name}' delimiter ',' csv header;".format(
            table_name=table_name, 
            column_list = ','.join(main_columns + ['count','date']),
            file_name=file_name)


main_columns_deaths = [u'UID', u'iso2', u'iso3', u'code3', u'FIPS', u'Admin2',
       u'Province_State', u'Country_Region', u'Lat', u'Long_', u'Combined_Key',
       u'Population']
deaths_df = pd.read_csv(US_DEATHS_CSV_PATH)
file_prefix = '/nfs01/covid_extracted_csvs/time_series_deaths_US'
new_files = extract_csv(deaths_df, file_prefix, main_columns_deaths)
for file_name in new_files:
    print(extract_sql_command(file_name, main_columns_deaths, 'time_series_covid19_deaths'))

main_columns_confirmed = [u'UID', u'iso2', u'iso3', u'code3', u'FIPS', u'Admin2',
       u'Province_State', u'Country_Region', u'Lat', u'Long_', u'Combined_Key']
confirmed_df = pd.read_csv(US_CONFIRMED_CSV_PATH)
file_prefix = '/nfs01/covid_extracted_csvs/time_series_confirmed_US'
new_files = extract_csv(confirmed_df, file_prefix, main_columns_confirmed)
for file_name in new_files:
    print(extract_sql_command(file_name, main_columns_confirmed, 'time_series_covid19_confirmed'))
