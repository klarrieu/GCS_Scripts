import os
from file_functions import *
from classify_landforms import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import arcpy
arcpy.env.overwriteOutput = True
from Tkinter import *
import logging
init_logger(__file__)

'''
Data structure:

Input:
list of .csv tables, each containing same number of rows, and columns named dist_down, Z_s, W_s, Z_s_W_s
index for division between reaches

Store in a dict (like JSON objects) with key = filename (or discharge, if given), value = pandas dataframe from corresponding .csv


Analyses between flows vs independent analyses for each flow

between all flows:
correlation between C(Z,W)s

between some flows:
hierarchical landform abundances/sequencing

independent:
% of Zs/Ws/Czw above/below certain values

outputs:


'''


def analysis_1(flows, reaches = None, zs = 'Z_s', ws = 'W_s', cov = 'Z_s_W_s'):
    '''Computes mean covariance, percent of covariance above/below 0, Pearson's correlation coefficient between Z_s and W_s

    Args:
        flows (list): a list of dataframes for each flow/stage
        reaches (dict, optional):  names of reaches as keys and dataframe index for beginning of each reach as values
        zs (str): name of the column in each flow dataframe for standardized bed elevation
        ws (str): name of the column in each flow dataframe for standardized width
        cov (str): name of the column in each flow dataframe for covariance series between Z_s and W_s

    Returns:
        list: a list of titled dataframes with flow/stages as rows and reaches as columns.
              One titled dataframe is returned for each of the following variables:
                  mean Cov(Z_s,W_s)
                  % of Cov(Z_s,W_s)>0
                  Pearsons's correlation between Z_s and W_s
    '''
    
    if reaches == None:
        mean_covs = []
        percent_above_0s = []
        correlations = []
        
        for flow in flows:
            mean_cov = np.mean(flow[cov])
            mean_covs.append(mean_cov)
            count_above_0 = flow[flow[cov]>0].count()[cov]
            percent_above_0 = count_above_0 * 1.0 / len(flow) * 100
            percent_above_0s.append(percent_above_0)
            correlation = np.corrcoef(flow[zs],flow[ws])[0][1]
            correlations.append(correlation)

        df_list = []
        variables = {'mean Czw':mean_covs,
                     r'% Czw>0':percent_above_0s,
                     'Corr(Zs,Ws)':correlations
                     }

        for title,data in zip(variables.keys(),variables.values()):
            df = DF(data)
            df.title = title
            df.index.name = 'flow'
            df.columns = [title]
            df_list.append(df)
        
        return df_list
            
    else:
        reach_df_lists = []
        for i in range(len(reaches)):
            reach_flows = [split_reaches(flow, reaches.values())[i] for flow in flows]
            reach_df_list = analysis_1(reach_flows, zs = zs, ws = ws, cov = cov)
            reach_df_lists.append(reach_df_list)
            
        df_list = []
        #for each variable
        for i in range(len(reach_df_lists[0])):
            #concatenate dataframes from each reach (use pandas dataframes to concat, then convert back to titled dfs)
            title = reach_df_lists[0][i].title
            df = pd.concat([pd.DataFrame(reach_df_list[i]) for reach_df_list in reach_df_lists], axis = 1)
            df = DF(df)
            df.title = title
            df.columns = reaches.keys()
            df.columns.name = 'reach'
            df_list.append(df)
        return df_list


def runs_test(series):
    '''Does WW runs test for values above/below median of series

    Args:
        series (list): a list of values for which to perform the Wald-Wolfowitz runs test for values below/above median

    Returns:
        A dataframe containing the following:
            number of runs
            number of expected runs (if random)
            expected standard deviation of number of runs (if random)
            Z: number of standard deviations difference between actual and expected number of run (standard deviation of num. of runs if random)
    '''

    m = np.median(series)
    #omit values from series equal to the median
    test_series = [x for x in series if x != m]
    run_lengths = []
    count = 0
    for i, vals in enumerate(zip(test_series, test_series[1:])):
        x1,x2 = vals
        count += 1
        #if transition between value above median to value equal or below median, end the run
        if (x1 > m and x2 < m) or (x1 < m and x2 > m):
            run_lengths.append(count)
            count = 0
        #if on the last value, but no transition, then last value is part of current run
        if i == len(test_series)-2:
            count += 1
            run_lengths.append(count)
            count = 0

    #total number of values (excluding median values)
    n = len(test_series)
    #num of values above median
    n_plus = sum([1 for x in test_series if x>m])
    #num of values below median
    n_minus = n-n_plus
    #expected number of runs if random
    exp_runs = 2*n_plus*n_minus*1.0/n + 1
    #actual number of runs
    num_runs = len(run_lengths)
    #standard deviation of expected num of runs if random
    exp_run_std = np.sqrt((exp_runs-1)*(exp_runs-2)*1.0/(n-1))
    #number of standard deviations (of exected run length) that the actual run count differs from WW expected run count
    z_diff_expected = (num_runs - exp_runs)*1.0/exp_run_std
    
    data = {'num. of runs':[num_runs],
            'expected runs':[exp_runs],
            'expected run StDev':[exp_run_std],
            'Z':[z_diff_expected]
            }
    
    return pd.DataFrame.from_dict(data)


def analysis_2(df, zs = 'Z_s', ws = 'W_s', cov = 'Z_s_W_s'):
    '''Does WW runs test on Z_s and W_s series in input dataframe'''
    series_list = [zs,ws]
    data = [runs_test(df[series]) for series in series_list]
    data = pd.concat(data)
    data.index = series_list

    return data


def analysis_3(flows, reaches = None, zs = 'Z_s', ws = 'W_s', cov = 'Z_s_W_s'):
    '''Computes mean(Ws), mean(Zs), % Abs(Ws)>0.5, % Abs(Zs)>0.5, % Abs(Ws)>1, % Abs(Zs)>1

    Args:
        flows (list): a list of dataframes for each flow/stage
        reaches (dict, optional):  names of reaches as keys and dataframe index for beginning of each reach as values
        zs (str): name of the column in each flow dataframe for standardized bed elevation
        ws (str): name of the column in each flow dataframe for standardized width
        cov (str): name of the column in each flow dataframe for covariance series between Z_s and W_s

    Returns:
        list: a list of titled dataframes with flow/stages as rows and reaches as columns.
              One titled dataframe is returned for each of the following variables:
                  mean Ws
                  mean Zs
                  % Abs(Ws)>0.5
                  % Abs(Zs)>0.5
                  % Abs(Ws)>1
                  % Abs(Zs)>1
    '''

    if reaches == None:
        mean_wss = []
        mean_zss = []
        abs_ws_percent_above_0pt5s = []
        abs_zs_percent_above_0pt5s = []
        abs_ws_percent_above_1s = []
        abs_zs_percent_above_1s = []

        for flow in flows:
            mean_ws = np.mean(flow[ws])
            mean_wss.append(mean_ws)
            mean_zs = np.mean(flow[zs])
            mean_zss.append(mean_zs)
            abs_ws_percent_above_0pt5 = sum([1 for x in flow[ws] if x>0.5])*1.0/len(flow)*100
            abs_ws_percent_above_0pt5s.append(abs_ws_percent_above_0pt5)
            abs_zs_percent_above_0pt5 = sum([1 for x in flow[zs] if x>0.5])*1.0/len(flow)*100
            abs_zs_percent_above_0pt5s.append(abs_zs_percent_above_0pt5)
            abs_ws_percent_above_1 = sum([1 for x in flow[ws] if x>1])*1.0/len(flow)*100
            abs_ws_percent_above_1s.append(abs_ws_percent_above_1)
            abs_zs_percent_above_1 = sum([1 for x in flow[zs] if x>1])*1.0/len(flow)*100
            abs_zs_percent_above_1s.append(abs_zs_percent_above_1)

        df_list = []
        variables = {'mean Ws':mean_wss,
                     'mean Zs':mean_zss,
                     '% Abs(Ws)>0.5':abs_ws_percent_above_0pt5s,
                     '% Abs(Zs)>0.5':abs_zs_percent_above_0pt5s,
                     '% Abs(Ws)>1':abs_ws_percent_above_1s,
                     '% Abs(Ws)>1':abs_zs_percent_above_1s
                     }
        
        for title,data in zip(variables.keys(),variables.values()):
            df = DF(data)
            df.title = title
            df.index.name = 'flow'
            df.columns = [title]
            df_list.append(df)

        return df_list
        
    else:
        reach_df_lists = []
        for i in range(len(reaches)):
            reach_flows = [split_reaches(flow, reaches.values())[i] for flow in flows]
            reach_df_list = analysis_3(reach_flows, zs = zs, ws = ws, cov = cov)
            reach_df_lists.append(reach_df_list)
            
        df_list = []
        #for each variable
        for i in range(len(reach_df_lists[0])):
            #concatenate dataframes from each reach (use pandas dataframes to concat, then convert back to titled dfs)
            title = reach_df_lists[0][i].title
            df = pd.concat([pd.DataFrame(reach_df_list[i]) for reach_df_list in reach_df_lists], axis = 1)
            df = DF(df)
            df.title = title
            df.columns = reaches.keys()
            df.columns.name = 'reach'
            df_list.append(df)
        return df_list


def analysis_4(df, zs = 'Z_s', ws = 'W_s', cov = 'Z_s_W_s', code = 'code'):
    '''Computes percentage of channel each landform class occupies'''
    n = len(df)
    percent_oversized = sum([1 for x in df[code] if x == -2])*1.0/n*100
    percent_const_pool = sum([1 for x in df[code] if x == -1])*1.0/n*100
    percent_normal = sum([1 for x in df[code] if x == 0])*1.0/n*100
    percent_wide_bar = sum([1 for x in df[code] if x == 1])*1.0/n*100
    percent_nozzle = sum([1 for x in df[code] if x == 2])*1.0/n*100

    data = {r'% oversized':[percent_oversized],
            r'% constricted pool':[percent_const_pool],
            r'% normal channel':[percent_normal],
            r'% wide bar':[percent_wide_bar],
            r'% nozzle':[percent_nozzle]
            }

    return pd.DataFrame.from_dict(data)


def analysis_5(df, zs = 'Z_s', ws = 'W_s', cov = 'Z_s_W_s', code = 'code'):
    '''Computes percentage of time each landform type was followed by each other landform type'''
    transitions = []
    for c1,c2 in zip(df[code],df[code][1:]):
        if c1 != c2:
            transitions.append([c1, c2])

    code_dict = {-2:'oversized',
                 -1:'const_pool',
                 0:'normal',
                 1:'wide_bar',
                 2:'nozzle'
                 }
    a = []
    for row_code in code_dict.keys():
        counts = []
        for col_code in code_dict.keys():
            count = sum([ 1 for x in transitions if x == [row_code, col_code] ])
            counts.append(count)
        row = [x*1.0/sum(counts)*100 for x in counts]
        a.append(row)
            
    df = pd.DataFrame(a, index = code_dict.values(), columns = code_dict.values())
    df.index.name = 'starting unit'
    df.columns.name = '% of time unit follows starting unit'
    return df
    

def analysis_6(flows, zs = 'Z_s', ws = 'W_s', cov = 'Z_s_W_s', code = 'code'):
    '''Identifies abundance and percentage of each hierarchical landform nesting permutation'''
    
    #combine code from each flow to create hierarchical landform series
    #must have wetted rectangle at every station so indices line up!
    #check that each flow dataframe is same length
    if not all(len(flow) == len(flows[0]) for flow in flows):
        raise Exception('Dataframes need to be same length for all flows.')

    hierarchy_codes = zip(*[flow[code] for flow in flows])

    #get unique hierarchical nesting permutations
    unique_nests = list(set(hierarchy_codes))

    #count number of times each unique nesting permutation occurs
    unique_nest_counts = list(np.zeros(len(unique_nests), dtype = int))
    for hier_code in hierarchy_codes:
        i = unique_nests.index(hier_code)
        unique_nest_counts[i] += 1

    #rank the permutations by count
    ranked_hierarchies = sorted(zip(unique_nest_counts,unique_nests))[::-1]

    #rank permutations
    
    return ranked_hierarchies
    

def complete_analysis(tables):
    '''Executes various analyses on tables (one table per discharge) containing columns for dist_down, Z_s, W_s, and Z_s_W_s'''


if __name__ == '__main__':

    tables = []
    flows = []
    for table in tables:
        clean_in_table(table)
        standardize(table)
        std_covar_series(table)
        df = classify_landforms(table)
        df = classify_landform_polygons(table, table.replace('_joined_table.csv', '.shp') )
        flows.append(df)

    df_list = analysis_1(flows)


    #make the GUI window
    root = Tk()
    root.wm_title('GCS Analysis')

    L1 = Label(root, text = 'GCS Data Tables')
    L1.grid(sticky = E, row = 0, column = 1)
    E1 = Entry(root, bd = 5)
    E1.grid(row = 0, column = 2)
    b1 = Button(root, text = 'Browse', command = lambda: browse(root, E1, select = 'files', ftypes = [('Comma-delimited text','.csv'),
                                                                                                      ('All files','*')]
                                                                )
                )
    b1.grid(sticky = W, row = 0, column = 3)
