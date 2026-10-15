import pandas as pd
from shutil import copy
import os
import TraceGUI as tg
import winsound
import TraceIO as tio
import matplotlib.pyplot as plt
import numpy as np
import ForceSpectroscopy as fs


# list of measurements

def measurements():
    measurements = []
    measurements.append(['190829', '032', '2', 'HMfAmut', 'Mutants'])
    measurements.append(['171107', '018', '29', 'HMfAwt', 'Mutants'])  # comparison-candidate --> winner
    measurements.append(['200224', '007', ' 75', 'HMfBmut', 'Mutants'])
    measurements.append(['200224', '023', ' 28', 'HMfBmut', 'Mutants'])
    measurements.append(['171201', '011', ' 5', 'HMfBwt', 'Mutants'])  # comparison candidate --> winner
    measurements.append(['171107', '005', ' 17', 'DNA', 'Mutants'])  # DNA ONLY - winner

    # DNA ONLY
    measurements.append(['171201', '005', ' 8', 'DNA', 'DNA'])  # DNA ONLY
    measurements.append(['171201', '004', ' 18', 'DNA', 'DNA'])  # DNA ONLY
    measurements.append(['171201', '005', ' 18', 'DNA', 'DNA'])  # DNA ONLY
    measurements.append(['171201', '004', ' 7', 'DNA', 'DNA'])  # DNA ONLY
    measurements.append(['171201', '005', ' 7', 'DNA', 'DNA'])  # DNA ONLY
    measurements.append(['171201', '004', ' 18', 'DNA', 'DNA'])  # DNA ONLY
    measurements.append(['171201', '005', ' 18', 'DNA', 'DNA'])  # DNA ONLY
    measurements.append(['171107', '005', ' 17', 'DNA', 'DNA'])  # DNA ONLY - winner

    # HMFA

    # fast

    measurements.append(['171107', '011', '4', 'HMfA', 'HMfA'])
    measurements.append(['171107', '011', '7', 'HMfA', 'HMfA'])
    measurements.append(['171107', '011', '13', 'HMfA', 'HMfA'])
    measurements.append(['171107', '011', '14', 'HMfA', 'HMfA'])
    measurements.append(['171107', '011', '19', 'HMfA', 'HMfA'])
    measurements.append(['171107', '012', '4', 'HMfA', 'HMfA'])
    measurements.append(['171107', '012', '13', 'HMfA', 'HMfA'])
    measurements.append(['171107', '012', '14', 'HMfA', 'HMfA'])
    measurements.append(['171107', '012', '19', 'HMfA', 'HMfA'])
    measurements.append(['171107', '013', '5', 'HMfA', 'HMfA'])
    measurements.append(['171107', '013', '18', 'HMfA', 'HMfA'])
    measurements.append(['171107', '015', '12', 'HMfA', 'HMfA'])
    measurements.append(['171107', '017', '19', 'HMfA', 'HMfA'])
    measurements.append(['171107', '027', '3', 'HMfA', 'HMfA'])
    measurements.append(['171107', '029', '1', 'HMfA', 'HMfA'])
    measurements.append(['171107', '039', '46', 'HMfA', 'HMfA'])
    measurements.append(['171107', '043', '13', 'HMfA', 'HMfA'])
    measurements.append(['171107', '045', '3', 'HMfA', 'HMfA'])
    measurements.append(['171107', '047', '23', 'HMfA', 'HMfA'])
    measurements.append(['171107', '051', '16', 'HMfA', 'HMfA'])
    measurements.append(['171107', '051', '28', 'HMfA', 'HMfA'])
    measurements.append(['171107', '051', '31', 'HMfA', 'HMfA'])
    measurements.append(['171107', '051', '40', 'HMfA', 'HMfA'])
    measurements.append(['171107', '054', '27', 'HMfA', 'HMfA'])
    measurements.append(['171107', '054', '31', 'HMfA', 'HMfA'])
    measurements.append(['171107', '054', '37', 'HMfA', 'HMfA'])
    measurements.append(['171107', '054', '43', 'HMfA', 'HMfA'])
    measurements.append(['171107', '055', '7', 'HMfA', 'HMfA'])
    measurements.append(['171107', '055', '15', 'HMfA', 'HMfA'])
    measurements.append(['171107', '055', '22', 'HMfA', 'HMfA'])
    measurements.append(['171107', '056', '20', 'HMfA', 'HMfA'])
    measurements.append(['171107', '056', '22', 'HMfA', 'HMfA'])
    measurements.append(['171107', '059', '10', 'HMfA', 'HMfA'])
    measurements.append(['171107', '059', '22', 'HMfA', 'HMfA'])
    measurements.append(['171107', '060', '10', 'HMfA', 'HMfA'])
    measurements.append(['171107', '060', '22', 'HMfA', 'HMfA'])
    measurements.append(['171107', '060', '24', 'HMfA', 'HMfA'])
    measurements.append(['171107', '063', '32', 'HMfA', 'HMfA'])
    measurements.append(['171107', '063', '35', 'HMfA', 'HMfA'])
    measurements.append(['171107', '064', '17', 'HMfA', 'HMfA'])

    # moderate

    measurements.append(['171107', '007', '2', 'HMfA', 'HMfA'])
    measurements.append(['171107', '007', '15', 'HMfA', 'HMfA'])  # tether broken
    measurements.append(['171107', '007', '17', 'HMfA', 'HMfA'])  # compare with DNA --> winner
    measurements.append(['171107', '007', '21', 'HMfA', 'HMfA'])
    measurements.append(['171107', '008', '11', 'HMfA', 'HMfA'])  # weird
    measurements.append(['171107', '009', '4', 'HMfA', 'HMfA'])
    measurements.append(['171107', '009', '5', 'HMfA', 'HMfA'])
    measurements.append(['171107', '009', '7', 'HMfA', 'HMfA'])
    measurements.append(['171107', '009', '14', 'HMfA', 'HMfA'])
    measurements.append(['171107', '009', '19', 'HMfA', 'HMfA'])
    # measurements.append(['171107', '009', '21', 'HMfA', 'HMfA'])  # bad tether
    measurements.append(['171107', '010', '4', 'HMfA', 'HMfA'])
    measurements.append(['171107', '010', '7', 'HMfA', 'HMfA'])
    measurements.append(['171107', '010', '14', 'HMfA', 'HMfA'])
    measurements.append(['171107', '014', '5', 'HMfA', 'HMfA'])
    measurements.append(['171107', '014', '18', 'HMfA', 'HMfA'])  # excellent, comparison-candidate
    measurements.append(['171107', '016', '12', 'HMfA', 'HMfA'])  # comparison-candidate
    measurements.append(['171107', '018', '29', 'HMfA', 'HMfA'])  # comparison-candidate --> winner
    measurements.append(['171107', '018', '40', 'HMfA', 'HMfA'])  # comparison-candidate
    measurements.append(['171107', '028', '3', 'HMfA', 'HMfA'])  # weird
    measurements.append(['171107', '030', '1', 'HMfA', 'HMfA'])
    measurements.append(['171107', '040', '22', 'HMfA', 'HMfA'])
    measurements.append(['171107', '042', '6', 'HMfA', 'HMfA'])
    measurements.append(['171107', '042', '13', 'HMfA', 'HMfA'])  # a bit strange
    measurements.append(['171107', '044', '13', 'HMfA', 'HMfA'])  # prob a double tether
    measurements.append(['171107', '048', '0', 'HMfA', 'HMfA'])
    measurements.append(['171107', '048', '23', 'HMfA', 'HMfA'])
    measurements.append(['171107', '052', '16', 'HMfA', 'HMfA'])  # prob a double tether
    measurements.append(['171107', '052', '30', 'HMfA', 'HMfA'])
    measurements.append(['171107', '052', '31', 'HMfA', 'HMfA'])
    measurements.append(['171107', '052', '40', 'HMfA', 'HMfA'])  # something sticks
    measurements.append(['171107', '053', '16', 'HMfA', 'HMfA'])  # ruptures
    measurements.append(
        ['171107', '053', '30', 'HMfA', 'HMfA'])  # overstretching - does not fit in our program unfortunately
    measurements.append(['171107', '057', '11', 'HMfA', 'HMfA'])
    measurements.append(
        ['171107', '058', '7', 'HMfA', 'HMfA'])  # overstretching - does not fit in our program unfortunately
    measurements.append(
        ['171107', '058', '11', 'HMfA', 'HMfA'])  # overstretching - does not fit in our program unfortunately
    measurements.append(['171107', '058', '14', 'HMfA', 'HMfA'])
    measurements.append(['171107', '061', '10', 'HMfA', 'HMfA'])
    measurements.append(['171107', '061', '32', 'HMfA', 'HMfA'])  # drift in fittable region
    measurements.append(['171107', '062', '10', 'HMfA', 'HMfA'])
    measurements.append(['171107', '062', '22', 'HMfA', 'HMfA'])
    measurements.append(['171107', '062', '24', 'HMfA', 'HMfA'])  # partial overstretching
    measurements.append(['171107', '065', '30', 'HMfA', 'HMfA'])
    measurements.append(['171107', '065', '47', 'HMfA', 'HMfA'])
    measurements.append(['171107', '066', '30', 'HMfA', 'HMfA'])

    # HMfB

    # fast

    # experiments June (different pulling rate, HMfB concentration --> hysteresis)

    measurements.append(['170622', '030', ' 41', 'HMfB', 'HMfB'])
    measurements.append(['170622', '030', ' 95', 'HMfB', 'HMfB'])
    measurements.append(['170622', '030', ' 96', 'HMfB', 'HMfB'])
    measurements.append(['170622', '037', ' 9', 'HMfB', 'HMfB'])
    measurements.append(['170622', '054', ' 26', 'HMfB', 'HMfB'])
    measurements.append(['170622', '054', ' 47', 'HMfB', 'HMfB'])
    measurements.append(['170622', '054', ' 60', 'HMfB', 'HMfB'])
    measurements.append(['170622', '055', ' 24', 'HMfB', 'HMfB'])
    measurements.append(['170622', '075', ' 3', 'HMfB', 'HMfB'])
    measurements.append(['170622', '077', ' 12', 'HMfB', 'HMfB'])
    measurements.append(['170622', '082', ' 7', 'HMfB', 'HMfB'])

    measurements.append(['170628', '037', ' 0', 'HMfB', 'HMfB'])
    measurements.append(['170628', '037', ' 6', 'HMfB', 'HMfB'])
    measurements.append(['170628', '037', ' 9', 'HMfB', 'HMfB'])
    measurements.append(['170628', '075', ' 3', 'HMfB', 'HMfB'])
    measurements.append(['170628', '077', ' 12', 'HMfB', 'HMfB'])
    measurements.append(['170628', '078', ' 11', 'HMfB', 'HMfB'])
    measurements.append(['170628', '082', ' 7', 'HMfB', 'HMfB'])

    # experiments December

    measurements.append(['171201', '006', ' 7', 'HMfB', 'HMfB'])
    measurements.append(['171201', '006', ' 8', 'HMfB', 'HMfB'])
    measurements.append(['171201', '006', ' 18', 'HMfB', 'HMfB'])
    measurements.append(['171201', '006', ' 22', 'HMfB', 'HMfB'])
    measurements.append(['171201', '010', ' 16', 'HMfB', 'HMfB'])
    measurements.append(['171201', '016', ' 14', 'HMfB', 'HMfB'])
    measurements.append(['171201', '016', ' 21', 'HMfB', 'HMfB'])
    measurements.append(['171201', '016', ' 22', 'HMfB', 'HMfB'])
    measurements.append(['171201', '017', ' 10', 'HMfB', 'HMfB'])
    measurements.append(['171201', '017', ' 13', 'HMfB', 'HMfB'])
    measurements.append(['171201', '017', ' 14', 'HMfB', 'HMfB'])
    measurements.append(['171201', '017', ' 22', 'HMfB', 'HMfB'])
    measurements.append(['171201', '020', ' 13', 'HMfB', 'HMfB'])
    measurements.append(['171201', '020', ' 14', 'HMfB', 'HMfB'])
    measurements.append(['171201', '020', ' 23', 'HMfB', 'HMfB'])
    measurements.append(['171201', '022', ' 18', 'HMfB', 'HMfB'])
    measurements.append(['171201', '025', ' 18', 'HMfB', 'HMfB'])
    measurements.append(['171201', '026', ' 7', 'HMfB', 'HMfB'])
    measurements.append(['171201', '029', ' 2', 'HMfB', 'HMfB'])
    measurements.append(['171201', '029', ' 29', 'HMfB', 'HMfB'])
    measurements.append(['171201', '034', ' 7', 'HMfB', 'HMfB'])
    measurements.append(['171201', '034', ' 11', 'HMfB', 'HMfB'])
    measurements.append(['171201', '037', ' 6', 'HMfB', 'HMfB'])
    measurements.append(['171201', '037', ' 7', 'HMfB', 'HMfB'])
    measurements.append(['171201', '037', ' 9', 'HMfB', 'HMfB'])
    measurements.append(['171201', '038', ' 3', 'HMfB', 'HMfB'])
    measurements.append(['171201', '038', ' 11', 'HMfB', 'HMfB'])
    measurements.append(['171201', '042', ' 4', 'HMfB', 'HMfB'])
    measurements.append(['171201', '045', ' 4', 'HMfB', 'HMfB'])
    measurements.append(['171201', '045', ' 18', 'HMfB', 'HMfB'])

    # moderate

    measurements.append(['171201', '007', ' 8', 'HMfB', 'HMfB'])
    measurements.append(['171201', '007', ' 18', 'HMfB', 'HMfB'])  # excellent
    measurements.append(['171201', '007', ' 22', 'HMfB', 'HMfB'])  # tether breaks
    measurements.append(['171201', '008', ' 7', 'HMfB', 'HMfB'])  # tether breaks
    measurements.append(['171201', '008', ' 18', 'HMfB', 'HMfB'])  # tether breaks
    measurements.append(['171201', '011', ' 5', 'HMfB', 'HMfB'])  # comparison candidate --> winner
    measurements.append(['171201', '011', ' 16', 'HMfB', 'HMfB'])  # fitting is a problem
    measurements.append(['171201', '012', ' 5', 'HMfB', 'HMfB'])  # excellent
    measurements.append(['171201', '012', ' 16', 'HMfB', 'HMfB'])  # tether breaks
    measurements.append(['171201', '013', ' 24', 'HMfB', 'HMfB'])  # fast -> might be outlier
    measurements.append(['171201', '014', ' 24', 'HMfB', 'HMfB'])
    measurements.append(['171201', '018', ' 10', 'HMfB', 'HMfB'])  # comparison candidate (some sticking at low force)
    measurements.append(['171201', '018', ' 13', 'HMfB', 'HMfB'])
    measurements.append(['171201', '018', ' 14', 'HMfB', 'HMfB'])
    measurements.append(['171201', '018', ' 22', 'HMfB', 'HMfB'])
    measurements.append(['171201', '018', ' 23', 'HMfB', 'HMfB'])  # something is amiss
    measurements.append(['171201', '019', ' 13', 'HMfB', 'HMfB'])  # very low stiffness fiber -> outlier
    measurements.append(
        ['171201', '019', ' 22', 'HMfB', 'HMfB'])  # overstretching - does not fit in our program unfortunately
    measurements.append(['171201', '023', ' 18', 'HMfB', 'HMfB'])
    measurements.append(['171201', '024', ' 18', 'HMfB', 'HMfB'])
    measurements.append(['171201', '027', ' 8', 'HMfB', 'HMfB'])
    measurements.append(['171201', '027', ' 29', 'HMfB', 'HMfB'])
    measurements.append(['171201', '028', ' 29', 'HMfB', 'HMfB'])
    measurements.append(
        ['171201', '028', ' 40', 'HMfB', 'HMfB'])  # overstretching - does not fit in our program unfortunately
    measurements.append(['171201', '035', ' 7', 'HMfB', 'HMfB'])
    measurements.append(['171201', '035', ' 11', 'HMfB', 'HMfB'])  # sticking at low force
    measurements.append(['171201', '036', ' 6', 'HMfB', 'HMfB'])  # comparison-candidate
    measurements.append(['171201', '036', ' 7', 'HMfB', 'HMfB'])  # comparison-candidate
    measurements.append(['171201', '036', ' 9', 'HMfB', 'HMfB'])
    measurements.append(['171201', '036', ' 11', 'HMfB', 'HMfB'])  # comparison-candidate -> actually tether breaking
    measurements.append(['171201', '039', ' 3', 'HMfB', 'HMfB'])
    measurements.append(['171201', '039', ' 11', 'HMfB', 'HMfB'])
    measurements.append(['171201', '043', ' 4', 'HMfB', 'HMfB'])
    measurements.append(['171201', '044', ' 4', 'HMfB', 'HMfB'])
    measurements.append(['171201', '044', ' 18', 'HMfB', 'HMfB'])  # double tether

    return measurements


def move_selected_files(start_dir=r'D:\Data\HMf', end_dir=r'D:\Data\HMf\paper_resubmission'):
    df = pd.DataFrame(measurements(), columns=['dir', 'filenr', 'trace', 'type', 'dest'])
    for _, m in df.iterrows():
        try:
            src = rf"{start_dir}\{m['dir']}\data_{m['filenr']}.dat"
            dst = rf"{end_dir}\{m['dest']}\{m['type']}_{m['dir']}_{m['filenr']}.dat"
            print(f'{src} > {dst}')
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            copy(src, dst)
        except FileNotFoundError:
            pass


def write_selected_h5(dir=r'D:\Data\HMf\paper_resubmission', dest=None):
    df = pd.DataFrame(measurements(), columns=['dir', 'filenr', 'trace', 'type', 'dest'])
    if dest is not None:
        df = df[df['dest'] == dest]

    for _, m in df.iterrows():
        try:
            filename = rf"{dir}\{m['dest']}\{m['type']}_{m['dir']}_{m['filenr']}.h5"
            print(f'{filename} : {m["trace"].strip()}')
            with tio.Traces(filename) as dataset:
                dataset.write_par(m['trace'].strip(), ' selected', 1, type='local')

        except FileNotFoundError:
            pass


def aggregate_params(dir=r'D:\Data\HMf\paper_resubmission', dest=None):
    df = pd.DataFrame(measurements(), columns=['dir', 'filenr', 'trace', 'type', 'dest'])
    if dest is not None:
        df = df[df['dest'] == dest]

    filenames = []
    for _, m in df.iterrows():
        try:
            filenames.append(rf"{dir}\{m['dest']}\{m['type']}_{m['dir']}_{m['filenr']}.h5")
        except FileNotFoundError:
            pass
    filenames = list(set(filenames))
    all_params = pd.DataFrame()
    for filename in filenames:
        print(filename)
        with tio.Traces(filename) as dataset:
            params = dataset.parameters
        params = params[params[' selected'] == 1].reset_index()
        params.insert(0, 'filename', filename)
        all_params = pd.concat([all_params, params], axis=0, sort=False)
    all_params.reset_index(inplace=True, drop=True)
    print(all_params)
    file_out = fr'{dir}\Parameters_{dest}.csv'
    all_params.to_csv(file_out, index=False)


def analyse_params(filenames):
    data_to_csv = pd.DataFrame()
    fileout = filenames[0].replace('Parameters_HMfA.csv', 'Fits_X.jpg')
    colors = ['red', 'blue']
    plots = []
    plots.append(['Ndimer (-)', 0, 150])
    plots.append(['k (pN_nm)', 0, 1])
    plots.append(['G1 (kT)', 0, 10])
    plots.append(['G2 (kT)', 0, 10])
    plots.append(['angle (deg)', 0, 25])
    import string
    alphabet = string.ascii_uppercase[3:]
    results = pd.DataFrame(index=['HMfA', 'HMfB'], columns=[t[0] for t in plots])
    for plot, ref in zip(plots, alphabet):
        par = plot[0]
        for filename, color in zip(filenames, colors):
            prot_name = filename.split('Parameters_')[-1].split('.')[0]
            df = pd.read_csv(filename)
            df = df[df['G1 (kT)'] > 0.5]
            df = df[df['k (pN_nm)'] < 1.5]

            bins = np.linspace(plot[1], plot[2], 20)
            plt.hist(df[par].values, bins=bins, color=color, alpha=0.5, lw=1, edgecolor=color, histtype='barstacked',
                     rwidth=0.7, label=prot_name)

            res = tio.round_significance(np.mean(df[par].values), np.std(df[par].values), delimiter=' +/- ')
            results.loc[prot_name][par] = res

            his, bins = np.histogram(df[par].values, bins=bins)
            bins = 0.5 * (bins[1:] + bins[:-1])
            # df = pd.DataFrame({par: bins, 'N': his})
            # csv_filename = filenames[0].replace('Parameters_HMfA.csv', f"Figure3{ref}_{filename.split('_')[-1]}")
            # print(csv_filename)
            # df.to_csv(csv_filename , index=False)
            label =filename.split('_')[-1].split('.')[0]
            if filename.split('.')[0][-1] == 'A':
                data_to_csv = pd.concat([data_to_csv, pd.DataFrame(bins, columns=[par])], axis=1)
            data_to_csv = pd.concat([data_to_csv, pd.DataFrame(his, columns=[label])], axis=1)


        if ref == 'D':
            legend = ['HMfA', 'HMfB']
        else:
            legend = None
        tio.format_plot(par.replace('_', r'/'), 'N', scale_page=0.3, aspect=0.5, ref=ref, yrange=[0, 35],
                        save=fileout.replace('X', par).split(' (')[0] + '.jpg', legend=legend)

    pd.DataFrame(results).T.to_csv(fileout.replace('Fits_X.jpg', 'Aggregate_results.csv'))
    print(data_to_csv)
    data_to_csv.to_csv(fileout.replace('Fits_X.jpg', 'Figure3D_H.csv'), index=False)

def plot_curves(folder, measurements, figure=None):
    data_to_csv = pd.DataFrame()
    results = None
    legend = ['_', '_']
    for m in measurements:
        filename = rf'{folder}\{m[0]}'
        print(filename)

        data = tg.CurrentData()
        data.check(filename, m[1])

        if m[2] == 'DNA':
            model, params, names = tio.prepare_fit('zWLC')
        else:
            model, params, names = tio.prepare_fit('zHMf2')

        if results is None:
            results = pd.DataFrame(index=[name[2] for name in measurements],
                                   columns=[names[key] for key in names if key not in ['x', 'y']])

        for key in names:
            if key not in ['x', 'y']:
                stored_value = data.pars.loc[names[key]].value
                if str(stored_value) != 'nan':
                    params[key].value = data.pars.loc[names[key]].value
                if data.pars.loc[names[key]].type == 'fit':
                    params[key].vary = True
                if str(data.pars.loc[names[key]].error) == 'nan':
                    res = str(data.pars.loc[names[key]].value)
                else:
                    res = tio.round_significance(data.pars.loc[names[key]].value, data.pars.loc[names[key]].error,
                                                  delimiter=' +/- ')
                results.loc[m[2]][names[key]] = res

        t = data.traces.loc['Time (s)'].values
        z = data.traces.loc['Z (um)'].values - tg.calculate_drift(data)
        f = data.traces.loc['F (pN)'].values
        selection = data.traces.loc['Selection'].values == 1.0
        selection = np.diff(f, prepend=0) > 0

        z -= data.pars.loc['Zbead (um)'].value
        plt.scatter(z, f, s=20, facecolors='none', edgecolors='lightgrey', linewidths=0.5)
        legend.append('_')
        plt.scatter(z[selection], f[selection], s=20, facecolors='none', edgecolors=m[-1], linewidths=0.5)
        data_to_csv = pd.concat([data_to_csv, pd.DataFrame(f, columns=['F (pN)'])], axis=1)
        data_to_csv = pd.concat([data_to_csv, pd.DataFrame(z, columns=['z (um)'])], axis=1)

        legend.append(m[-2])
        model_f = np.linspace(5e-2, 50, 1000)
        data_to_csv = pd.concat([data_to_csv, pd.DataFrame(model_f, columns=['F (pN)'])], axis=1)
        zfit = model.eval(x=model_f, params=params)
        zfit -= data.pars.loc['Zbead (um)'].value

        if figure != 'figure3A':
            plt.plot(zfit, model_f, color='black')
            data_to_csv = pd.concat([data_to_csv, pd.DataFrame(zfit, columns=['z (um)'])], axis=1)

        if m[2] == 'DNA':
            plt.plot(zfit, model_f, color='black', linestyle='dashed')
            plt.scatter(z * -100, f, s=20, facecolors='none', edgecolors='lightgrey', linewidths=0.5)
            legend.append('release')
            data_to_csv = pd.concat([data_to_csv, pd.DataFrame(zfit, columns=['z (um)'])], axis=1)

        if figure == 'figure3A' and m[2] != 'DNA':
            legend.append('_')
            params['g1'].value = -1
            zfit = model.eval(x=model_f, params=params)
            zfit -= data.pars.loc['Zbead (um)'].value
            plt.plot(zfit, model_f, color='black', linestyle='dashed')
            data_to_csv = pd.concat([data_to_csv, pd.DataFrame(zfit, columns=['z (um)'])], axis=1)

            legend.append('_')
            params['g2'].value = -1
            # params['angle'].value = 70
            zfit = model.eval(x=model_f, params=params)
            zfit -= data.pars.loc['Zbead (um)'].value
            plt.plot(zfit, model_f, color='black', linestyle='dashed')
            data_to_csv = pd.concat([data_to_csv, pd.DataFrame(zfit, columns=['z (um)'])], axis=1)

    if figure in ['figure3B', 'figure4A', 'figure4B']:
        plt.scatter(z / 0, f, s=20, facecolors='none', edgecolors='lightgrey', linewidths=0.5)
        legend.append('release')

    if figure == 'figure3A':
        legend = ['_' for i, l in enumerate(legend)]
        legend[4] = 'HMfB'
        legend[6] = 'DNA'
        legend[7] = 'release'

        print(legend)

    # filename_plot = rf'{folder}\{figure}.jpg'
    # tio.format_plot(r'z ($\mu$m)', 'F (pN)', ref=figure[-1],
    #                 title='model', scale_page=0.42, aspect=1.4, ylog=True, xrange=[-0.25, 1.35],
    #                 yrange=[0.05, 70], legend=legend, save=filename_plot)
    # results.T.to_csv(filename_plot.replace('jpg', 'csv'))

    print(figure)
    print(data_to_csv)
    filename_data = rf'{folder}\{figure}_data.csv'
    data_to_csv.to_csv(filename_data, index=False)


# aggregate analysis
if True:
    filenames = []
    filenames.append(r'D:\Data\HMf\paper_resubmission\Parameters_HMfA.csv')
    filenames.append(r'D:\Data\HMf\paper_resubmission\Parameters_HMfB.csv')
    analyse_params(filenames)

# F-z plots
if False:
    folder = r"D:\Data\HMf\paper_resubmission\Mutants"
    measurements = []
    measurements.append(['HMfAwt_171107_018.h5', '29', 'HMfA wt', 'lightcoral'])
    measurements.append(['HMfAmut_190829_032.h5', '2', 'mutant', 'red'])
    plot_curves(folder, measurements, figure='figure4A')

    measurements = []
    measurements.append(['HMfBwt_171201_011.h5', '5', 'HMfB wt', 'lightblue'])
    measurements.append(['HMfBmut_200224_023.h5', '28', 'mutant', 'mediumblue'])
    plot_curves(folder, measurements, figure='figure4B')

    # measurements = []
    # measurements.append(['HMfAwt_171107_018.h5', '29', 'HMfA', 'red'])
    # measurements.append(['HMfBwt_171201_011.h5', '5', 'HMfB', 'mediumblue'])
    # measurements.append(['DNA_171107_005.h5', '17', 'DNA', 'orange'])
    # plot_curves(folder, measurements, figure='figure3A')

    # measurements = []
    # measurements.append(['HMfAwt_171107_018.h5', '29', 'HMfA', 'red'])
    # measurements.append(['HMfBwt_171201_011.h5', '5', 'HMfB', 'mediumblue'])
    # plot_curves(folder, measurements, figure='figure3B')

winsound.MessageBeep()
