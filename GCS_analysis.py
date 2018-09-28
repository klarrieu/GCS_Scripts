from file_functions import *
from classify_landforms_GUI import *
import pandas as pd
import numpy as np
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
correlation between C(Z_S,W_s)s

between some flows:
hierarchical landform abundances/sequencing

independent:
% of Zs/Ws/Czw above/below certain values

outputs:


'''


def flow_cov_corrs(data):
    '''
    Computes binary correlation between C(Z_s,W_s)'s for each pair of flows
    '''
    output = []

    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    # create dataframes of correlations between flows for each reach
    for reach_name in reach_names:
        # the (titled) dataframe to output for each reach
        reach_df = DF(index=flow_names, columns=flow_names, title='Corr(Czw_i,Czw_j) %s' % reach_name)

        # loop for rows, loop for columns
        for flow_rkey in flow_names:
            for flow_ckey in flow_names:
                # covariance series corresponding to flow row
                czw_r = data[flow_rkey][reach_name]['Z_s_W_s'].tolist()
                # covariance series corresponding to flow column
                czw_c = data[flow_ckey][reach_name]['Z_s_W_s'].tolist()
                # add correlation between the series to dataframe
                reach_df.loc[flow_rkey, flow_ckey] = np.corrcoef(czw_r, czw_c)[0][1]

        output.append(reach_df)

    return output


def zw_corrs(data):
    '''
    Computes the correlation between Z_s and W_s for each flow and reach
    '''

    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Corr(Z_s,W_s)')

    for flow in flow_names:
        for reach in reach_names:
            zs = data[flow][reach]['Z_s'].tolist()
            ws = data[flow][reach]['W_s'].tolist()
            corr = np.corrcoef(zs, ws)[0][1]
            output.loc[flow, reach] = corr

    return output


def czw_means(data):
    '''
    Computes mean C(Z_s,W_s) values for each flow and reach
    '''

    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Mean C(Z_s,W_s)')

    for flow in flow_names:
        for reach in reach_names:
            czw = data[flow][reach]['Z_s_W_s'].tolist()
            output.loc[flow, reach] = np.mean(czw)

    return output


def czw_pos_percents(data):
    '''
    Computes percent of C(Z_s,W_s) values above 0 for each flow and reach
    '''

    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Percent of C(Z_s,W_s) > 0')

    for flow in flow_names:
        for reach in reach_names:
            czw = data[flow][reach]['Z_s_W_s'].tolist()
            n = len(czw)
            num_above_0 = len([x for x in czw if x > 0])
            percent = num_above_0 * 100.0 / n
            output.loc[flow, reach] = percent

    return output


def zs_percents(data):
    '''
    Computes percent of abs(Z_s) values above 1 for each flow and reach
    '''
    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Percent of |Z_s| > 1')

    for flow in flow_names:
        for reach in reach_names:
            zs = data[flow][reach]['Z_s'].tolist()
            n = len(zs)
            num_abs_above_1 = len([x for x in zs if abs(x) > 1])
            percent = num_abs_above_1 * 100.0 / n
            output.loc[flow, reach] = percent

    return output


def ws_percents(data):
    '''
    Computes percent of abs(W_s) values above 1 for each flow and reach
    '''
    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Percent of |W_s| > 1')

    for flow in flow_names:
        for reach in reach_names:
            ws = data[flow][reach]['W_s'].tolist()
            n = len(ws)
            num_abs_above_1 = len([x for x in ws if abs(x) > 1])
            percent = num_abs_above_1 * 100.0 / n
            output.loc[flow, reach] = percent

    return output


def ww_runs_z(data):
    '''
    Computes runs and expected runs for values above/below median Z_s
    '''
    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = []

    # df for each reach, row for each flow, columns for runs, expected runs, expected run stdev, and z (number of expected stdevs difference between runs and expected runs)
    for reach in reach_names:
        # the (titled) dataframe to output for each reach
        reach_df = DF(index=flow_names, columns=['Runs', 'Expected Runs', 'Expected Run StDev', 'Z'],
                      title='Z_s WW runs test %s' % reach)

        for flow in flow_names:
            zs = data[flow][reach]['Z_s'].tolist()
            run_df = runs_test(zs)
            run_df.index = [flow]
            reach_df.append(run_df)

        output.append(reach_df)

    return output


def ww_runs_w(data):
    '''
    Computes runs and expected runs for values above/below median W_s
    '''
    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = []

    # df for each reach, row for each flow, columns for runs, expected runs, expected run stdev, and z (number of expected stdevs difference between runs and expected runs)
    for reach in reach_names:
        # the (titled) dataframe to output for each reach
        reach_df = DF(index=flow_names, columns=['Runs', 'Expected Runs', 'Expected Run StDev', 'Z'],
                      title='W_s WW runs test %s' % reach)

        for flow in flow_names:
            ws = data[flow][reach]['W_s'].tolist()
            run_df = runs_test(ws)
            run_df.index = [flow]
            reach_df.append(run_df)

        output.append(reach_df)

    return output


def analysis_1(flows, reaches=None, zs='Z_s', ws='W_s', cov='Z_s_W_s'):
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
            count_above_0 = flow[flow[cov] > 0].count()[cov]
            percent_above_0 = count_above_0 * 1.0 / len(flow) * 100
            percent_above_0s.append(percent_above_0)
            correlation = np.corrcoef(flow[zs], flow[ws])[0][1]
            correlations.append(correlation)

        df_list = []
        variables = {'mean Czw': mean_covs,
                     r'% Czw>0': percent_above_0s,
                     'Corr(Zs,Ws)': correlations
                     }

        for title, data in zip(variables.keys(), variables.values()):
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
            reach_df_list = analysis_1(reach_flows, zs=zs, ws=ws, cov=cov)
            reach_df_lists.append(reach_df_list)

        df_list = []
        # for each variable
        for i in range(len(reach_df_lists[0])):
            # concatenate dataframes from each reach (use pandas dataframes to concat, then convert back to titled dfs)
            title = reach_df_lists[0][i].title
            df = pd.concat([pd.DataFrame(reach_df_list[i]) for reach_df_list in reach_df_lists], axis=1)
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
    # omit values from series equal to the median
    test_series = [x for x in series if x != m]
    run_lengths = []
    count = 0
    for i, vals in enumerate(zip(test_series, test_series[1:])):
        x1, x2 = vals
        count += 1
        # if transition between value above median to value equal or below median, end the run
        if (x1 > m and x2 < m) or (x1 < m and x2 > m):
            run_lengths.append(count)
            count = 0
        # if on the last value, but no transition, then last value is part of current run
        if i == len(test_series) - 2:
            count += 1
            run_lengths.append(count)
            count = 0

    # total number of values (excluding median values)
    n = len(test_series)
    # num of values above median
    n_plus = sum([1 for x in test_series if x > m])
    # num of values below median
    n_minus = n - n_plus
    # expected number of runs if random
    exp_runs = 2 * n_plus * n_minus * 1.0 / n + 1
    # actual number of runs
    num_runs = len(run_lengths)
    # standard deviation of expected num of runs if random
    exp_run_std = np.sqrt((exp_runs - 1) * (exp_runs - 2) * 1.0 / (n - 1))
    # number of standard deviations (of exected run length) that the actual run count differs from WW expected run count
    z_diff_expected = (num_runs - exp_runs) * 1.0 / exp_run_std

    data = {'Runs': [num_runs],
            'Expected Runs': [exp_runs],
            'Expected Run StDev': [exp_run_std],
            'Z': [z_diff_expected]
            }

    return pd.DataFrame.from_dict(data)


def analysis_2(df, zs='Z_s', ws='W_s', cov='Z_s_W_s'):
    '''Does WW runs test on Z_s and W_s series in input dataframe'''
    series_list = [zs, ws]
    data = [runs_test(df[series]) for series in series_list]
    data = pd.concat(data)
    data.index = series_list

    return data


def analysis_3(flows, reaches=None, zs='Z_s', ws='W_s', cov='Z_s_W_s'):
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
            abs_ws_percent_above_0pt5 = sum([1 for x in flow[ws] if x > 0.5]) * 1.0 / len(flow) * 100
            abs_ws_percent_above_0pt5s.append(abs_ws_percent_above_0pt5)
            abs_zs_percent_above_0pt5 = sum([1 for x in flow[zs] if x > 0.5]) * 1.0 / len(flow) * 100
            abs_zs_percent_above_0pt5s.append(abs_zs_percent_above_0pt5)
            abs_ws_percent_above_1 = sum([1 for x in flow[ws] if x > 1]) * 1.0 / len(flow) * 100
            abs_ws_percent_above_1s.append(abs_ws_percent_above_1)
            abs_zs_percent_above_1 = sum([1 for x in flow[zs] if x > 1]) * 1.0 / len(flow) * 100
            abs_zs_percent_above_1s.append(abs_zs_percent_above_1)

        df_list = []
        variables = {'mean Ws': mean_wss,
                     'mean Zs': mean_zss,
                     '% Abs(Ws)>0.5': abs_ws_percent_above_0pt5s,
                     '% Abs(Zs)>0.5': abs_zs_percent_above_0pt5s,
                     '% Abs(Ws)>1': abs_ws_percent_above_1s,
                     '% Abs(Ws)>1': abs_zs_percent_above_1s
                     }

        for title, data in zip(variables.keys(), variables.values()):
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
            reach_df_list = analysis_3(reach_flows, zs=zs, ws=ws, cov=cov)
            reach_df_lists.append(reach_df_list)

        df_list = []
        # for each variable
        for i in range(len(reach_df_lists[0])):
            # concatenate dataframes from each reach (use pandas dataframes to concat, then convert back to titled dfs)
            title = reach_df_lists[0][i].title
            df = pd.concat([pd.DataFrame(reach_df_list[i]) for reach_df_list in reach_df_lists], axis=1)
            df = DF(df)
            df.title = title
            df.columns = reaches.keys()
            df.columns.name = 'reach'
            df_list.append(df)
        return df_list


def analysis_4(df, zs='Z_s', ws='W_s', cov='Z_s_W_s', code='code'):
    '''Computes percentage of channel each landform class occupies'''
    n = len(df)
    percent_oversized = sum([1 for x in df[code] if x == -2]) * 1.0 / n * 100
    percent_const_pool = sum([1 for x in df[code] if x == -1]) * 1.0 / n * 100
    percent_normal = sum([1 for x in df[code] if x == 0]) * 1.0 / n * 100
    percent_wide_bar = sum([1 for x in df[code] if x == 1]) * 1.0 / n * 100
    percent_nozzle = sum([1 for x in df[code] if x == 2]) * 1.0 / n * 100

    data = {r'% oversized': [percent_oversized],
            r'% constricted pool': [percent_const_pool],
            r'% normal channel': [percent_normal],
            r'% wide bar': [percent_wide_bar],
            r'% nozzle': [percent_nozzle]
            }

    return pd.DataFrame.from_dict(data)


def analysis_5(df, zs='Z_s', ws='W_s', cov='Z_s_W_s', code='code'):
    '''Computes percentage of time each landform type was followed by each other landform type'''
    transitions = []
    for c1, c2 in zip(df[code], df[code][1:]):
        if c1 != c2:
            transitions.append([c1, c2])

    code_dict = {-2: 'oversized',
                 -1: 'const_pool',
                 0: 'normal',
                 1: 'wide_bar',
                 2: 'nozzle'
                 }
    a = []
    for row_code in code_dict.keys():
        counts = []
        for col_code in code_dict.keys():
            count = sum([1 for x in transitions if x == [row_code, col_code]])
            counts.append(count)
        row = [x * 1.0 / sum(counts) * 100 for x in counts]
        a.append(row)

    df = pd.DataFrame(a, index=code_dict.values(), columns=code_dict.values())
    df.index.name = 'starting unit'
    df.columns.name = '% of time unit follows starting unit'
    return df


def analysis_6(flows, zs='Z_s', ws='W_s', cov='Z_s_W_s', code='code'):
    '''Identifies abundance and percentage of each hierarchical landform nesting permutation'''

    # combine code from each flow to create hierarchical landform series
    # must have wetted rectangle at every station so indices line up!
    # check that each flow dataframe is same length
    if not all(len(flow) == len(flows[0]) for flow in flows):
        raise Exception('Dataframes need to be same length for all flows.')

    hierarchy_codes = zip(*[flow[code] for flow in flows])

    # get unique hierarchical nesting permutations
    unique_nests = list(set(hierarchy_codes))

    # count number of times each unique nesting permutation occurs
    unique_nest_counts = list(np.zeros(len(unique_nests), dtype=int))
    for hier_code in hierarchy_codes:
        i = unique_nests.index(hier_code)
        unique_nest_counts[i] += 1

    # rank the permutations by count
    ranked_hierarchies = sorted(zip(unique_nest_counts, unique_nests))[::-1]

    # rank permutations

    return ranked_hierarchies


def clean_in_data(tables, reach_breaks=None):
    '''
    Structures data for all further analyses
    
    Args:
        tables: a list of filenames, each a .csv table for a particular discharge with the same number of rows, containing columns for dist_down, Z (detrended), W
        reach_breaks: a list of indices where new reaches begin. If reach_breaks == None, it will just look at the entire dataset
    
    Returns:
        a dict of dicts containing dataframes:
            top level keys: original filenames/discharges
            top level values: a dict
            second level keys: section of river, i.e. "All", "Reach 1", "Reach 2", ...
            second level values: pandas dataframes for the discharge specified by key 1 and reach specified by key 2
            
            In each dataframe, columns are added for standardized Z_s, W_s, covariance series Z_s_W_s, and code (GCS landform classification)
    '''

    # make sure tables all have same number of rows!

    data = {}

    for table in tables:

        # ensures each table has columns named dist_down, Z, and W for distance downstream, detrended elevation, and channel width
        clean_in_table(table)
        # creates columns for standardized Z_s and W_s
        standardize(table)
        # creates column Z_s_W_s for covariance series between standardized detrended elevation and standardized channel width
        std_covar_series(table)
        # adds column "code" for landforms based on covariance value
        landforms(table)

        # split up into second level dict based on reach
        if reach_breaks != None:
            df = pd.read_csv(table)
            table_dict = {'All': df}

            reach_indices = split_list(df.index.tolist(), reach_breaks)
            i = 1  # reach number
            for reach in reach_indices:
                reach_df = df.loc[reach]
                table_dict['Reach %i' % i] = reach_df
                i += 1

            data[table] = table_dict



        # if no reach breaks, just have one second level key: "All"
        else:
            df = pd.read_csv(table)
            table_dict = {'All': df}
            data[table] = table_dict

    return data


def complete_analysis(tables, reach_breaks=None):
    '''Executes various analyses'''

    logging.info('Cleaning input data...')
    data = clean_in_data(tables, reach_breaks=reach_breaks)
    logging.info('OK')

    logging.info('Computing C(Z_s,W_s) correlation between flows...')
    czw_corrs = flow_cov_corrs(data)
    logging.info('OK')

    logging.info('Computing Corr(Z_s,W_s)...')
    corrs = zw_corrs(data)
    logging.info('OK')

    logging.info('Computing mean C(Z_s,W_s)...')
    means = czw_means(data)
    logging.info('OK')

    logging.info('Computing percentage of C(Z_s,W_s) values > 0')
    percents = czw_pos_percents(data)
    logging.info('OK')

    logging.info('Computing percentage of abs(Z_s) values > 1')
    zs_percntz = zs_percents(data)
    logging.info('OK')

    logging.info('Computing percentage of abs(W_s) values > 1')
    ws_percntz = ws_percents(data)
    logging.info('OK')

    logging.info('Conducting Wald-Wolfowitz runs tests...')
    ww_test_z = ww_runs_z(data)
    ww_test_w = ww_runs_w(data)
    logging.info('OK')

    logging.info('Writing outputs to files...')

    # save tables to excel files/sheets
    writer = pd.ExcelWriter('GCS_output_data.xlsx', engine='xlsxwriter')

    for df_list in [czw_corrs, ww_test_z, ww_test_w]:
        for df in df_list:
            df.to_excel(writer, sheet_name=df.title)

    for df in [corrs, means, percents, zs_percntz, ws_percntz]:
        df.to_excel(writer, sheet_name=df.title)

    writer.save()

    # save images of plots

    logging.info('OK')

    return data


if __name__ == '__main__':
    '''
    tables = []
    flows = []
    for table in tables:
        df = classify_landform_polygons(table, table.replace('_joined_table.csv', '.shp') )
    '''

    # make the GUI window
    root = Tk()
    root.wm_title('GCS Analysis')

    L1 = Label(root, text='GCS Data Tables: ')
    L1.grid(sticky=E, row=0, column=1)
    E1 = Entry(root, bd=5)
    E1.grid(row=0, column=2)
    b1 = Button(root, text='Browse',
                command=lambda: browse(root, E1, select='files', ftypes=[('Comma-delimited text', '.csv'),
                                                                         ('All files', '*')]
                                       )
                )
    b1.grid(sticky=W, row=0, column=3)

    L2 = Label(root, text='Reach Breaks: ')
    L2.grid(sticky=E, row=1, column=1)
    E2 = Entry(root, bd=5)
    E2.grid(row=1, column=2)

    b = Button(root, text='    Run    ', command=lambda: complete_analysis(tables=list(root.tk.splitlist(E1.get())),
                                                                           reach_breaks=map(int, E2.get().split(
                                                                               ',')) if E2.get() != '' else None
                                                                           )
               )
    b.grid(sticky=W, row=2, column=2)
    root.mainloop()
