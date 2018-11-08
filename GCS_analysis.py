from file_functions import *
from classify_landforms_GUI import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Tkinter import *
import logging

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


def runs_test(series):
    '''
    Does WW runs test for values above/below median of series

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

    data = {'Runs': num_runs,
            'Expected Runs': exp_runs,
            'Expected Run StDev': exp_run_std,
            'Z': z_diff_expected
            }

    return data


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
            runs_data = runs_test(zs)
            reach_df.loc[flow, runs_data.keys()] = runs_data.values()

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
            runs_data = runs_test(ws)
            reach_df.loc[flow, runs_data.keys()] = runs_data.values()

        output.append(reach_df)

    return output


def landform_occupation(data):
    '''Computes percentage of channel each landform class occupies (not nested)'''
    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())
    # code number and corresponding MU
    code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}

    output = []

    # table for each reach, row for each flow, columns for landform class, percentage as values
    for reach in reach_names:

        reach_df = DF(index=flow_names, columns=['O', 'CP', 'NC', 'WB', 'NZ'],
                      title='%s %% Landform Composition' % reach)

        for flow in flow_names:
            code_series = data[flow][reach]['code'].tolist()
            n = len(code_series)
            for code in code_dict.keys():
                mu = code_dict[code]
                count = sum([1 for x in code_series if x == code])
                percentage = round(count*100.0/n, 2)
                reach_df.loc[flow, mu] = percentage

        output.append(reach_df)

    return output


def landform_following(data, nc=True):
    '''Computes percentage of times each landform type is followed by each other landform type type'''
    flow_names = sorted(data.keys())
    # code number and corresponding MU
    if nc:
        code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}
    else:
        code_dict = {-2: 'O', -1: 'CP', 1: 'WB', 2: 'NZ'}

    output = []

    # do not split into reaches, do for entire river segment
    # one table per flow, starting unit as row, following unit as col, percent of times starting unit is followed by following unit as values
    for flow in flow_names:

        code_series = data[flow]['All']['code'].tolist()

        # code pairs for each transition
        transitions = []
        if nc:
            flow_df = DF(index=code_dict.values(), columns=code_dict.values(),
                         title='%s %% times unit follows starting unit' % flow)
            flow_df.index.name = 'starting unit'

            for c1, c2 in zip(code_series, code_series[1:]):
                if c1 != c2:
                    transitions.append([c1, c2])

        else:
            # exclude normal channel condition
            flow_df = DF(index=code_dict.values(), columns=code_dict.values(),
                         title='%s NoNC %% times unit follows starting unit' % flow)
            flow_df.index.name = 'starting unit'
            no_nc = filter(lambda c: c != 0, code_series)
            for c1, c2 in zip(no_nc, no_nc[1:]):
                if c1 != c2:
                    transitions.append([c1, c2])

        for start_code in code_dict.keys():
            # n = total number of transitions for the start_code
            n = sum([1 for x in transitions if x[0] == start_code])

            for follow_code in code_dict.keys():
                # number of times start_code is followed by follow_code
                count = sum([1 for x in transitions if x == [start_code, follow_code]])

                start_mu = code_dict[start_code]
                follow_mu = code_dict[follow_code]
                percentage = round(count*100.0/n, 2)
                flow_df.loc[start_mu, follow_mu] = percentage

        output.append(flow_df)

    return output


def landform_nesting(data):
    '''Identifies abundance of each hierarchically nested landform'''
    flow_names = sorted(data.keys())
    # code number and corresponding MU
    code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}

    output = DF(columns=sorted(flow_names, reverse=True) + ['Count', '%% of River'],
                title='Nested landform abundance')

    # output table containing flood, bankfull, basefull landforms, count, and percentage as columns, row for each sequence in descending order of abundance
    flood_codes = data[flow_names[2]]['All']['code'].tolist()
    bankfull_codes = data[flow_names[1]]['All']['code'].tolist()
    baseflow_codes = data[flow_names[0]]['All']['code'].tolist()

    nested_landforms = list(zip(flood_codes, bankfull_codes, baseflow_codes))
    unique_nests = list(set(nested_landforms))

    # initialize list of lists to count abundance
    unique_nest_counts = list(np.zeros(len(unique_nests), dtype=int))

    for nest in nested_landforms:
        i = unique_nests.index(nest)
        unique_nest_counts[i] += 1

    nest_abundances = list(zip(unique_nests, unique_nest_counts))
    nest_abundances.sort(key=lambda x: x[1], reverse=True)

    for i, nest in enumerate(nest_abundances):
        n = len(nested_landforms)
        p = round(nest[1]*100.0/n, 2)
        output.loc[i] = [code_dict[nst] for nst in nest[0]] + [nest[1], p]

    return output


def GCS_plots(data):
    '''Returns list of plots for Z_s, W_s at each flow, and C(Z_s,W_s) at each flow'''

    flow_names = sorted(data.keys())

    output = []

    # Z_s plot for each flow

    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        zs = data[flow]['All']['Z_s'].tolist()

        fig = plt.figure()
        plt.title(r'$Z_s$' + ' %s' % flow)
        plt.xlabel('Distance downstream (m)')
        plt.ylabel('$Z_s$')
        plt.grid()
        plt.plot(dist, zs)
        fig.savefig('Zs %s.png' % flow)
        output.append(fig)
        plt.close(fig)

    # W_s plot for each flow

    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        zs = data[flow]['All']['W_s'].tolist()

        fig = plt.figure()
        plt.title(r'$W_s$' + ' %s' % flow)
        plt.xlabel('Distance downstream (m)')
        plt.ylabel(r'$W_s$')
        plt.grid()
        plt.plot(dist, zs)
        fig.savefig('Ws %s.png' % flow)
        output.append(fig)
        plt.close(fig)

    # GCS at each flow, same plot (lines colored by flow)

    fig = plt.figure()
    plt.title(r'$C(Z_s, W_s)$')
    plt.xlabel('Distance downstream (m)')
    plt.ylabel(r'$Z_s \cdot W_s$')
    plt.grid()

    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        czw = data[flow]['All']['Z_s_W_s'].tolist()
        plt.plot(dist, czw, label=flow)

    plt.legend()
    fig.savefig('Czw.png')
    output.append(fig)
    plt.close(fig)

    # GCS plot for each flow (points colored by landform type)

    color_code = {-2: 'black', -1: 'blue', 0: 'grey', 1: 'orange', 2: 'red'}
    code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}

    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        czw = data[flow]['All']['Z_s_W_s'].tolist()
        code = data[flow]['All']['code'].tolist()

        fig = plt.figure()
        plt.title(r'$C(Z_s, W_s)$' + ' %s' % flow)
        plt.xlabel('Distance downstream (m)')
        plt.ylabel(r'$Z_s \cdot W_s$')
        plt.grid()
        plt.scatter(dist, czw, c=[color_code[x] for x in code], s=4, edgecolors='none')
        legend_elements = [Line2D([0], [0], marker='o', color='w', label=mu, markerfacecolor=color, markeredgecolor='none', markersize=10)
                           for color, mu in [[color_code[x], code_dict[x]] for x in code_dict.keys()]]
        plt.legend(handles=legend_elements, ncol=5)
        fig.savefig('Czw landforms %s.png' % flow)
        output.append(fig)
        plt.close(fig)

    # Fourier transform of GCS at each flow

    fig = plt.figure()
    plt.title(r'$\hat{f}(Czw)$')
    plt.xlabel('Frequency' + r'$(m^{-1})$')
    plt.grid()
    for flow in flow_names:
        x = data[flow]['All']['dist_down'].tolist()
        gcs = data[flow]['All']['Z_s_W_s'].tolist()
        xf, yf = ft(x, gcs)
        plt.semilogx(xf, yf, label=flow)
    plt.legend()
    fig.savefig('GCSFourier.png')
    output.append(fig)
    plt.close(fig)

    # GCS autocorrelation at each flow
    '''
    fig = plt.figure()
    plt.title('Autocorrelation')
    plt.xlabel('Lag')
    plt.grid()
    for flow in flow_names:
        # x = data[flow]['All']['dist_down'].tolist()
        gcs = data[flow]['All']['Z_s_W_s'].tolist()
        plt.acorr(gcs, usevlines=True, normed=True, label=flow)
    plt.legend()
    plt.savefig('GCSacorr.png')
    output.append(fig)
    plt.close()
    '''
    # GCS cross-correlation between flows
    # 2D Fourier transform for Ws and Zs at each flow
    # power spectral density

    return output


def clean_in_data(tables, reach_breaks=None):
    '''
    Structures data for all further analyses
    
    Args:
        tables: a list of filenames, each a .csv table for a particular discharge with the same number of rows, containing columns for dist_down, Z (detrended), W
        reach_breaks: a list of dist_down values where new reaches begin. If reach_breaks == None, it will just look at the entire dataset
    
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

        flow_name = os.path.basename(table)
        try:
            flow_name = os.path.basename(table).split('_')[0]
            logging.info('Using %s as flow name for %s'%(flow_name,table))
        except:
            pass

        # ensures each table has columns named dist_down, Z, and W for distance downstream, detrended elevation, and channel width
        clean_in_table(table)
        # creates columns for standardized Z_s and W_s
        standardize(table)
        # creates column Z_s_W_s for covariance series between standardized detrended elevation and standardized channel width
        std_covar_series(table)
        # adds column "code" for landforms based on covariance value
        landforms(table)

        # split up into second level dict based on reach
        if reach_breaks is not None:
            df = pd.read_csv(table)
            table_dict = {'All': df}

            # dataframe indices corresponding to reach break distance values
            reach_break_indices = [df.index[df['dist_down'] == reach_break].tolist()[0] for reach_break in reach_breaks]
            reach_indices = split_list(df.index.tolist(), reach_break_indices)
            i = 1  # reach number
            for reach in reach_indices:
                reach_df = df.loc[reach]
                table_dict['Reach %i' % i] = reach_df
                i += 1

            data[flow_name] = table_dict

        # if no reach breaks list provided, look for 'Reach' attribute
        else:
            df = pd.read_csv(table)
            try:
                # try using 'Reach' attribute
                reach_names = list(set(df['Reach'].tolist()))
                reaches = [df[df['Reach'] == reach_name] for reach_name in reach_names]
                reach_names = ['Reach '+str(reach) for reach in reach_names]
                # if only one reach == df, use this one reach for table dict ('All' would be redundant)
                if len(reaches) == 1 and reaches[0] == df:
                    table_dict = dict(zip(reach_names, reaches))
                else:
                    table_dict = dict(zip(reach_names, reaches))
                    table_dict['All'] = df
                data[flow_name] = table_dict
            except:
                # if no reach breaks can be identified, just use the whole enchilada
                logging.info('Couldn\'t identify reach breaks.')
                table_dict = {'All': df}
                data[flow_name] = table_dict

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

    logging.info('Determining landform abundance...')
    landform_percents = landform_occupation(data)
    logging.info('OK')

    logging.info('Determining landform sequencing...')
    follows = landform_following(data)
    follows_no_nc = landform_following(data, nc=False)
    logging.info('OK')

    logging.info('Determining nested landform abundance...')
    nests = landform_nesting(data)
    logging.info('OK')

    logging.info('Writing outputs to files...')

    # save tables to excel files/sheets
    writer = pd.ExcelWriter('GCS_output_data.xlsx', engine='xlsxwriter')

    for df_list in [czw_corrs, ww_test_z, ww_test_w, landform_percents, follows, follows_no_nc]:
        for df in df_list:
            # 31 character limit for excel sheet name
            if len(df.title) > 31:
                logging.info('%s title too long for excel sheet name. Shortening to 31 characters.' % df.title)
            df.to_excel(writer, sheet_name=df.title[:31])

    for df in [corrs, means, percents, zs_percntz, ws_percntz, nests]:
        if len(df.title) > 31:
            logging.info('%s title too long for excel sheet name. Shortening to 31 characters.' % df.title)
        df.to_excel(writer, sheet_name=df.title[:31])

    writer.save()
    logging.info('OK')

    # save images of plots

    logging.info('Making plots...')
    plots = GCS_plots(data)
    logging.info('OK')

    return data


if __name__ == '__main__':
    '''
    tables = []
    flows = []
    for table in tables:
        df = classify_landform_polygons(table, table.replace('_joined_table.csv', '.shp') )
    '''

    # initialize logger
    init_logger(__file__)
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

    b = Button(root, text='    Run    ',
               command=lambda: complete_analysis(tables=list(root.tk.splitlist(E1.get())))
               )
    b.grid(sticky=W, row=2, column=2)

    root.mainloop()
