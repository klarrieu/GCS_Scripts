from file_functions import *
from classify_landforms_GUI import *
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.tri as tri
import scipy.signal as sig
from Tkinter import *
import logging

# Set the default color cycle
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['deepskyblue', 'royalblue', 'navy', 'black'])


def flow_cov_corrs(data, field_1, field_2):
    '''
    Computes binary correlation between C(field_1,field_2)'s for each pair of flows
    '''
    output = []

    flow_names = sorted(data.keys())
    # sort flow names by numerical value instead of alphabetically
    flow_nums = [float(flow.replace('pt', '.').replace('cms', '')) for flow in flow_names]
    flow_names = zip(*sorted(zip(flow_nums, flow_names)))[1]
    reach_names = sorted(data.values()[0].keys())

    # create dataframes of correlations between flows for each reach
    for reach_name in reach_names:
        # the (titled) dataframe to output for each reach
        field_abbrev = field_1[0].lower() + field_2[0].lower()
        reach_df = DF(index=flow_names, columns=flow_names, title='Corr(C%s_i,C%s_j) %s' % (field_abbrev, field_abbrev, reach_name))

        # loop for rows, loop for columns
        for flow_rkey in flow_names:
            for flow_ckey in flow_names:
                # covariance series corresponding to flow row
                gcs_r = data[flow_rkey][reach_name]['%s_%s' % (field_1, field_2)].tolist()
                # covariance series corresponding to flow column
                gcs_c = data[flow_ckey][reach_name]['%s_%s' % (field_1, field_2)].tolist()
                # add correlation between the series to dataframe
                reach_df.loc[flow_rkey, flow_ckey] = np.corrcoef(gcs_r, gcs_c)[0][1]

        output.append(reach_df)

    return output


def agg_corrs(data, field_1, field_2):
    '''
    Computes the (aggregate) correlation between field_1 and field_2 for each flow and reach
    '''

    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Corr(%s,%s)' % (field_1, field_2))

    for flow in flow_names:
        for reach in reach_names:
            f1 = data[flow][reach][field_1].tolist()
            f2 = data[flow][reach][field_2].tolist()
            corr = np.corrcoef(f1, f2)[0][1]
            output.loc[flow, reach] = corr

    return output


def gcs_means(data, field_1, field_2):
    '''
    Computes mean C(field_1,field_2) values for each flow and reach
    '''

    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Mean C(%s,%s)' % (field_1, field_2))

    for flow in flow_names:
        for reach in reach_names:
            gcs = data[flow][reach]['%s_%s' % (field_1, field_2)].tolist()
            output.loc[flow, reach] = np.mean(gcs)

    return output


def gcs_pos_percents(data, field_1, field_2):
    '''
    Computes percent of C(field_1,field_2) values above 0 for each flow and reach
    '''

    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Percent of C(%s,%s) > 0' % (field_1, field_2))

    for flow in flow_names:
        for reach in reach_names:
            gcs = data[flow][reach]['%s_%s' % (field_1, field_2)].tolist()
            n = len(gcs)
            num_above_0 = len([x for x in gcs if x > 0])
            percent = num_above_0 * 100.0 / n
            output.loc[flow, reach] = percent

    return output


def percents(data, field):
    '''
    Computes percent of abs(field) values above 1 for each flow and reach
    '''
    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = DF(index=flow_names, columns=reach_names, title='Percent of |%s| > 1' % field)

    for flow in flow_names:
        for reach in reach_names:
            series = data[flow][reach][field].tolist()
            n = len(series)
            num_abs_above_1 = len([x for x in series if abs(x) > 1])
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


def ww_runs(data, field):
    '''
    Computes runs and expected runs for values above/below median field value
    '''
    flow_names = sorted(data.keys())
    reach_names = sorted(data.values()[0].keys())

    output = []

    # df for each reach, row for each flow, columns for runs, expected runs, expected run stdev, and z (number of expected stdevs difference between runs and expected runs)
    for reach in reach_names:
        # the (titled) dataframe to output for each reach
        reach_df = DF(index=flow_names, columns=['Runs', 'Expected Runs', 'Expected Run StDev', 'Z'],
                      title='%s WW runs test %s' % (field, reach))

        for flow in flow_names:
            series = data[flow][reach][field].tolist()
            runs_data = runs_test(series)
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


def series_plots(data, field, odir=''):
    '''Returns list of plots for the given field at each flow'''

    flow_names = sorted(data.keys())
    output = []
    figsize=(12, 6)

    # plot series for each flow
    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        series = data[flow]['All'][field].tolist()
        fig = plt.figure(figsize=figsize)
        plt.title(r'$%s$' % field + ' %s' % flow.replace('pt', '.').replace('cms', ' cms'))
        plt.xlabel('Distance downstream (m)')
        plt.ylabel(r'$%s$' % field)
        plt.grid()
        plt.plot(dist, series)
        fig.savefig(odir + '%s %s.png' % (field, flow.replace('pt', '.').replace('cms', ' cms')), bbox_inches='tight', pad_inches=0.1)
        output.append(fig)
        plt.close(fig)

    # one plot for all flows
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.set_title(r'$%s$' % field)
    ax.set_xlabel('Distance downstream (m)')
    ax.set_ylabel(r'$%s (\sigma)$' % field)
    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        series = data[flow]['All'][field].tolist()
        ax.plot(dist, series, label=flow.replace('pt', '.').replace('cms', ' cms'))
        ax.set_xlim(dist[0], dist[-1])
    ax.legend()
    ax.grid()
    fig.savefig(odir + '%s_all_discharges.png' % field, bbox_inches='tight', pad_inches=0.1)
    output.append(fig)
    plt.close(fig)


def GCS_plots(data, field_1, field_2, odir=''):
    '''Returns list of plots for field_1, field_2 at each flow, and C(field_1,field_2) at each flow'''
    field_abbrev = field_1[0].lower() + field_2[0].lower()
    flow_names = sorted(data.keys())
    qs = [float(flow.replace('pt', '.').replace('cms', '')) for flow in flow_names]
    output = []
    figsize = (12, 6)

    # GCS for all flows
    fig = plt.figure(figsize=figsize)
    plt.ylabel(r'$C(%s, %s)$' % (field_1, field_2))
    plt.xlabel('Distance downstream (m)')
    plt.grid()
    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        gcs = data[flow]['All']['%s_%s' % (field_1, field_2)].tolist()
        plt.plot(dist, gcs, label=flow.replace('pt', '.').replace('cms', ' cms'))
    plt.legend()
    plt.axhline(0, linestyle='--', color='grey')
    plt.xlim(dist[0], dist[-1])
    fig.savefig(odir + 'C%s.png' % field_abbrev, bbox_inches='tight', pad_inches=0.1)
    output.append(fig)
    plt.close(fig)

    # field_1, field_2, GCS in one plot for each flow
    for flow in flow_names:
        dist = data[flow]['All']['dist_down'].tolist()
        s1 = data[flow]['All'][field_1].tolist()
        s2 = data[flow]['All'][field_2].tolist()
        gcs = data[flow]['All']['%s_%s' % (field_1, field_2)].tolist()
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_title(flow.replace('pt', '.').replace('cms', ' cms'))
        ax.plot(dist, s1, label=r'$%s$' % field_1)
        ax.plot(dist, s2, label=r'$%s$' % field_2)
        ax.plot(dist, gcs, label=r'$C(%s, %s)$' % (field_1, field_2))
        ax.legend()
        ax.grid()
        ax.set_xlabel('Distance downstream (m)')
        fig.savefig(odir + 'C%s_%s.png' % (field_abbrev, flow), bbox_inches='tight', pad_inches=0.1)
        output.append(fig)
        plt.close(fig)

    # GCS plot for each flow (points colored by landform type)
    if 'W_s' in [field_1, field_2] and 'Z_s' in [field_1, field_2]:
        color_code = {-2: 'black', -1: 'blue', 0: 'grey', 1: 'orange', 2: 'red'}
        code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}
        for flow in flow_names:
            dist = data[flow]['All']['dist_down'].tolist()
            czw = data[flow]['All']['W_s_Z_s'].tolist()
            code = data[flow]['All']['code'].tolist()
            fig = plt.figure(figsize=figsize)
            plt.title(r'$C(Z_s, W_s)$' + ' %s' % flow.replace('pt', '.').replace('cms', ' cms'))
            plt.xlabel('Distance downstream (m)')
            plt.ylabel(r'$Z_s \cdot W_s$')
            plt.grid()
            plt.scatter(dist, czw, c=[color_code[x] for x in code], s=4, edgecolors='none')
            legend_elements = [Line2D([0], [0], marker='o', color='w', label=mu, markerfacecolor=color, markeredgecolor='none', markersize=10)
                               for color, mu in [[color_code[x], code_dict[x]] for x in code_dict.keys()]]
            plt.legend(handles=legend_elements, ncol=5)
            fig.savefig(odir + 'Czw landforms %s.png' % flow, bbox_inches='tight', pad_inches=0.1)
            output.append(fig)
            plt.close(fig)


    # Bivariate Correlation vs discharge
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    corr_list = []
    for flow in flow_names:
        s1 = data[flow]['All'][field_1]
        s2 = data[flow]['All'][field_2]
        corr = np.corrcoef(s1, s2)[0][1]
        corr_list.append(corr)
    qs, corr_list = zip(*sorted(zip(qs, corr_list)))

    n = data[data.keys()[0]]['All'].shape[0]
    lower_cband, upper_cband = r_confidence_interval(0, n)

    ax.semilogx(qs, corr_list, '-o', color='black', markersize=2)
    ax.axhline(lower_cband, linestyle='--', color='grey')
    ax.axhline(upper_cband, linestyle='--', color='grey')
    qticks = list(qs[:3]) + [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    tick_labels = [float(str(q).rstrip('0').rstrip('.')) for q in qticks]
    ax.set_xticks(qticks)
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(min(qs), max(qs))
    plt.xticks(rotation=-45)
    ax.set(xlabel=r'Discharge $(m^3/s)$', ylabel='Correlation')
    ax.grid()
    fig.savefig(odir + 'bivar_corrs_%s.png' % field_abbrev, bbox_inches='tight', pad_inches=0.1)

    # Autocorrelations
    fig, ax = plt.subplots(3, 3, sharex=True, sharey=True, figsize=figsize)
    for i, flow in enumerate(flow_names[:3]):
        s1 = data[flow]['All'][field_1]
        s2 = data[flow]['All'][field_2]
        gcs = data[flow]['All']['%s_%s' % (field_1, field_2)]
        dist_down = data[flow]['All']['dist_down']
        spacing = abs(dist_down[1] - dist_down[0])
        maxlags = int(len(dist_down)/2)
        lags, lower_white, upper_white = white_noise_acf_ci(s1, maxlags=maxlags)

        lags, ar1_acorrs, lower_red, upper_red = ar1_acorr(s1, maxlags=maxlags)
        ax[0][i].plot(*cox_acorr(s1, maxlags=maxlags))
        ax[0][i].plot(lags, ar1_acorrs, color='red')
        ax[0][i].plot(lags, lower_red, '--', color='salmon')
        ax[0][i].plot(lags, upper_red, '--', color='salmon')
        ax[0][i].plot(lags, lower_white, '--', color='grey')
        ax[0][i].plot(lags, upper_white, '--', color='grey')
        ax[0][i].set_ylabel(r'$%s$' % field_1 + ' Autocorrelation')

        lags, ar1_acorrs, lower_red, upper_red = ar1_acorr(s2, maxlags=maxlags)
        ax[1][i].plot(*cox_acorr(s2, maxlags=maxlags))
        ax[1][i].plot(lags, ar1_acorrs, color='red')
        ax[1][i].plot(lags, lower_red, '--', color='salmon')
        ax[1][i].plot(lags, upper_red, '--', color='salmon')
        ax[1][i].plot(lags, lower_white, '--', color='grey')
        ax[1][i].plot(lags, upper_white, '--', color='grey')
        ax[1][i].set_ylabel(r'$%s$' % field_2 + ' Autocorrelation')

        lags, ar1_acorrs, lower_red, upper_red = ar1_acorr(gcs, maxlags=maxlags)
        ax[2][i].plot(*cox_acorr(gcs, maxlags=maxlags))
        ax[2][i].plot(lags, ar1_acorrs, color='red')
        ax[2][i].plot(lags, lower_red, '--', color='salmon')
        ax[2][i].plot(lags, upper_red, '--', color='salmon')
        ax[2][i].plot(lags, lower_white, '--', color='grey')
        ax[2][i].plot(lags, upper_white, '--', color='grey')
        ax[2][i].set_ylabel('GCS Autocorrelation')

        ax[0][i].set_title(flow.replace('pt', '.').replace('cms', ' cms'))
        ax[2][i].set_xlabel('Lag (m)')
        for j in range(3):
            ax[j][i].grid()
            ax[j][i].set_xlim(0, maxlags)
            ticks = map(int, ax[j][i].get_xticks() * spacing)
            ax[j][i].set_xticklabels(ticks)
    fig.savefig(odir + 'Acorrs_%s.png' % field_abbrev, bbox_inches='tight', pad_inches=0.1)

    # Autocorrelation Heat Map
    x = []
    y = []
    z = []
    lower_whites = []
    upper_whites = []
    for i, flow in enumerate(flow_names):
        s1 = data[flow]['All'][field_1]
        s2 = data[flow]['All'][field_2]
        gcs = data[flow]['All']['%s_%s' % (field_1, field_2)]
        dist_down = data[flow]['All']['dist_down']
        spacing = abs(dist_down[1] - dist_down[0])
        maxlags = int(len(dist_down) / 2)
        lags, lower_white, upper_white = white_noise_acf_ci(gcs, maxlags=maxlags)
        lags, acorrs = cox_acorr(gcs, maxlags=maxlags)
        # only include positive autocorrelations
        for j, acorr_val in enumerate(acorrs):
            if acorr_val > 0:
                x.append(np.log10(qs[i]))
                y.append(lags[j])
                z.append(acorrs[j])
                lower_whites.append(lower_white[j])
                upper_whites.append(upper_white[j])

    triang = tri.Triangulation(x, y)
    # mask where acorr (z) is below the white noise threshold at that lag (y)
    mask = []
    for triangle in triang.triangles:
        z_vals = [z[vertex] for vertex in triangle]
        lws = [lower_whites[vertex] for vertex in triangle]
        uws = [upper_whites[vertex] for vertex in triangle]
        cond_1 = all(z_val > lw for z_val, lw in zip(z_vals, lws))
        cond_2 = all(z_val < uw for z_val, uw in zip(z_vals, uws))
        if cond_1 and cond_2:
            mask_val = 1
        else:
            mask_val = 0
        mask.append(mask_val)
    triang.set_mask(mask)

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    tpc = ax.tripcolor(triang, z, cmap='jet')
    ax.set(xlabel='Discharge (cms)', ylabel='Lag (m)')
    qlabels = list(qs[:3]) + [2, 3, 4, 5, 6, 8, 10]  # show less discharge ticks to make plot less cluttered
    xticks = [np.log10(q) for q in qlabels]
    ax.set_xticks(xticks)
    xticklabels = [q for q in qlabels]
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=-45)
    yticklabels = map(int, ax.get_yticks() * spacing)
    ax.set_yticklabels(yticklabels)
    fig.colorbar(tpc)
    fig.savefig(odir + 'Acorr_heatmap_%s.png' % field_abbrev, bbox_inches='tight', pad_inches=0.1)


    # power spectral density
    x = []
    y = []
    z = []
    for i, flow in enumerate(flow_names):
        s1 = data[flow]['All'][field_1]
        s2 = data[flow]['All'][field_2]
        gcs = data[flow]['All']['%s_%s' % (field_1, field_2)]
        dist_down = data[flow]['All']['dist_down']
        spacing = abs(dist_down[1] - dist_down[0])
        frequencies, psd = sig.periodogram(gcs, 1.0/3, window=sig.get_window('hamming', len(gcs)))
        for j, psd_val in enumerate(psd):
            x.append(np.log10(qs[i]))
            y.append(frequencies[j])
            z.append(psd[j])
    std_z = [z_val*1.0/np.std(z) for z_val in z]

    triang = tri.Triangulation(x, y)

    mask = []
    for triangle in triang.triangles:
        z_vals = [z[vertex] for vertex in triangle]
        lws = [lower_whites[vertex] for vertex in triangle]
        uws = [upper_whites[vertex] for vertex in triangle]
        cond_1 = all(z_val > lw for z_val, lw in zip(z_vals, lws))
        cond_2 = all(z_val < uw for z_val, uw in zip(z_vals, uws))
        if cond_1 and cond_2:
            mask_val = 1
        else:
            mask_val = 0
        mask.append(mask_val)
    triang.set_mask(mask)

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    tpc = ax.tripcolor(triang, std_z, cmap='jet')
    ax.set(xlabel='Discharge (cms)', ylabel='Spatial Frequency (cycles/m)')
    qlabels = list(qs[:3]) + [2, 3, 4, 5, 6, 8, 10]  # show less discharge ticks to make plot less cluttered
    xticks = [np.log10(q) for q in qlabels]
    ax.set_xticks(xticks)
    xticklabels = [q for q in qlabels]
    ax.set_xticklabels(xticklabels)
    plt.xticks(rotation=-45)
    ax.set_ylim(min(y), 0.01)
    fig.colorbar(tpc)
    fig.savefig(odir + 'PSD_heatmap_%s.png' % field_abbrev, bbox_inches='tight', pad_inches=0.1)

    # Landform stratified velocity

    # quadrant histogram/ MU histogram

    return output


def clean_in_data(tables, fields, reach_breaks=None):
    '''
    Structures data for all further analyses
    
    Args:
        tables: a list of filenames, each a .csv table for a particular discharge with the same number of rows, containing columns for dist_down, Z, W, V
        fields: list of fields to use for standardizing/calculating GCS
        reach_breaks: a list of dist_down values where new reaches begin. If reach_breaks == None, it will just look at the entire dataset
    
    Returns:
        a dict of dicts containing dataframes:
            top level keys: original filenames/discharges
            top level values: a dict
            second level keys: section of river, i.e. "All", "Reach 1", "Reach 2", ...
            second level values: pandas dataframes for the discharge specified by key 1 and reach specified by key 2
            
            In each dataframe, columns are added for standardized Z_s, W_s, V_s, GCS between all pairs, and code (GCS landform classification)
    '''

    std_fields = [field + '_s' for field in fields]
    std_pairs = itertools.combinations(std_fields, 2)
    data = {}

    for table in tables:

        flow_name = os.path.basename(table)
        try:
            flow_name = os.path.basename(table).split('_')[0]
            logging.info('Using %s as flow name for %s' % (flow_name, table))
        except:
            pass

        # ensures each table has columns named dist_down, W, Z, V for distance downstream, channel width, detrended elevation, velocity
        clean_in_table(table)
        # creates columns for standardized fields
        standardize(table, fields=fields)
        # creates column for GCS series between all standardized field pairs
        for std_pair in std_pairs:
            std_covar_series(table, std_pair[0], std_pair[1])
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
    '''Executes various analyses

    Args:
        tables: list of file names for cross-section data tables
        reach_breaks: list of distances downstream to set reach breaks
    '''

    odir = os.path.dirname(__file__)+'\\results\\'
    if os.path.isdir(odir) == False:
        os.mkdir(odir)

    fields = ['W', 'Z', 'V']
    std_fields = [field + '_s' for field in fields]
    std_pairs = list(itertools.combinations(std_fields, 2))

    df_lol = []

    logging.info('Cleaning input data...')
    data = clean_in_data(tables, fields=fields, reach_breaks=reach_breaks)
    logging.info('OK')

    logging.info('Computing GCS correlation between flows...')
    for std_pair in std_pairs:
        gcs_corrs = flow_cov_corrs(data, std_pair[0], std_pair[1])
        df_lol.append(gcs_corrs)
    logging.info('OK')

    logging.info('Computing correlations between field pairs...')
    corrs_list = []
    for std_pair in std_pairs:
        corrs = agg_corrs(data, std_pair[0], std_pair[1])
        corrs_list.append(corrs)
    df_lol.append(corrs_list)
    logging.info('OK')

    logging.info('Computing GCS means...')
    means_list = []
    for std_pair in std_pairs:
        means = gcs_means(data, std_pair[0], std_pair[1])
        means_list.append(means)
    df_lol.append(means_list)
    logging.info('OK')

    logging.info('Computing percentage of GCS values > 0')
    pos_percents_list = []
    for std_pair in std_pairs:
        pos_percents = gcs_pos_percents(data, std_pair[0], std_pair[1])
        pos_percents_list.append(pos_percents)
    df_lol.append(pos_percents_list)
    logging.info('OK')

    logging.info('Computing percentage of abs(std_field) values > 1')
    percntz_list = []
    for std_field in std_fields:
        percntz = percents(data, field=std_field)
        percntz_list.append(percntz)
    df_lol.append(percntz_list)
    logging.info('OK')

    logging.info('Conducting Wald-Wolfowitz runs tests...')
    for std_field in std_fields:
        ww_test = ww_runs(data, field=std_field)
        df_lol.append(ww_test)
    logging.info('OK')

    logging.info('Determining landform abundance...')
    landform_percents = landform_occupation(data)
    df_lol.append(landform_percents)
    logging.info('OK')

    logging.info('Determining landform sequencing...')
    follows = landform_following(data)
    df_lol.append(follows)
    follows_no_nc = landform_following(data, nc=False)
    df_lol.append(follows_no_nc)
    logging.info('OK')

    logging.info('Determining nested landform abundance...')
    if len(data.keys()) == 3:
        nests = landform_nesting(data)
    else:
        nests = landform_nesting({col: data[col] for col in data.keys() if col in data.keys()[:3]})
    df_lol.append([nests])
    logging.info('OK')

    logging.info('Writing outputs to files...')

    # save tables to excel files/sheets
    writer = pd.ExcelWriter(odir + 'GCS_output_data.xlsx', engine='xlsxwriter')

    for df_list in df_lol:
        for df in df_list:
            # 31 character limit for excel sheet name
            if len(df.title) > 31:
                logging.info('%s title too long for excel sheet name. Shortening to 31 characters.' % df.title)
            df.to_excel(writer, sheet_name=df.title[:31])

    writer.save()
    logging.info('OK')

    # save images of plots

    logging.info('Making plots...')
    for std_field in std_fields:
        series_plots(data, std_field, odir=odir)
    for std_pair in std_pairs:
        GCS_plots(data, std_pair[0], std_pair[1], odir=odir)
    logging.info('OK')
    logging.info('Finished. Results in: %s' % odir)

    return data


if __name__ == '__main__':

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
