# -*- coding: utf-8 -*-
import sys
import glob
import time

import peptide_fragmentor
import pymzml

import pandas as pd
import numpy as np

import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html

from dash.dependencies import Input, Output, State

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
# allow output to not yet existing elements
app.config['suppress_callback_exceptions'] = True


mzml_files = [path for path in glob.glob('../data/*.mzML')]
mzml_files += [path for path in glob.glob('../data/*.gz')]
ident_files = [path for path in glob.glob('../data/*.csv')]

# global global_data
# gobal_data = {
#     'run': None,
#     'df' : None,
# }


app.layout = html.Div(children=[
    html.H1(
        children='pymzML GUI',
        style={
            'textAlign': 'center',
        }
    ),
    html.Div(
        children='pymzML GUI, yo!',
        style={
            'textAlign': 'center',
        }
    ),
    html.Label('mzML file:'),
    dcc.Dropdown(
        id='input-mzml',
        options=[{'label': i, 'value': i} for i in mzml_files],
    ),
    html.Label('Ident file:'),
    dcc.Dropdown(
        id='input-ident',
        options=[{'label': i, 'value': i} for i in ident_files],
    ),
    html.Div(id='table', children=[dash_table.DataTable(id='table-obj')]),
    # dash_table.DataTable(id='table'),
    html.Label('Peptide Annotation:'),
    dcc.Input(id='peptide-input', placeholder='ELVISLIVES', type='text', debounce=True),
    dcc.Graph(
        id='spec-graph',
        figure={
            'data': [
                {'x': [1, 2, 3], 'y': [4, 1, 2], 'type': 'line', 'name': 'SF'},
            ],
            'layout': {
                'title': 'Spectrum 1'
            }
        }
    ),
    html.Div(id='spec-slider', children=[dcc.Slider(id='spec-slider-ele')]),
    html.Div(id='intermediate-value'),
    html.Div(id='intermediate-value2')
])


def match_df_to_mz_i_list(fragment_df, mz_i_list, dalton_cutoff=0.02, rel_cutoff=None, test=False):
    if not isinstance(mz_i_list, np.ndarray):
        mz_i_list = np.ndarray(mz_i_list)
    mz = mz_i_list[:,0]
    matches = []
    mismatches = []
    # print('\nMATCHES')
    for i, row in fragment_df.iterrows():
        # print(row['mass'])
        if rel_cutoff:
            cutoff = row['mass'] * rel_cutoff
        else:
            cutoff = dalton_cutoff
        res =  (row['mass'], 'max_intensity', row['name'])
        if (abs(mz - row['mass']) < cutoff).any():
            matches.append(res)
        else:
            mismatches.append(res)
        # print(cutoff)
        # print(mz[abs(mz - row['mass']) < cutoff])
    if test is True:
        matches.append(
            (131.2, 'max_intensity', 'TEST')
        )
    # print('\n\nMATCHED\n\n')
    # print(matches, mismatches)
    return matches, mismatches


@app.callback(
    Output('spec-graph', 'figure'),
    [
        Input('spec-slider-ele', 'value'),
        Input('peptide-input', 'value'),
        Input('table-obj', "derived_virtual_data"),
        Input('table-obj', "derived_virtual_selected_rows")
    ],
    # [
    #     State('table-obj', "derived_virtual_data"),
    #     State('table-obj', "derived_virtual_selected_rows")
    # ]
)
def update_figure(spec_id, seq, df_data, df_row):
    fac = pymzml.plot.Factory()
    fac.new_plot(MS_precision=20e-6)
    data = {'data': [], 'layout': []}
    print('UPDATE FIGURE')
    # print(f'INPUT {spec_id}, {seq}, {df_data}, {df_row}')
    if spec_id is not None:
        print('-'*100)
        print('Plot data')
        start = time.time()
        spec = run[spec_id]
        print('Took {0} seconds to get spectrum'.format(time.time() - start))
        spec_peaks = spec.peaks('centroided')
        fac.add(data=spec_peaks, name='data spec_id: {0}'.format(spec_id))
        data_tmp = fac.save()
        data['data'] += data_tmp['data']
        print('-'*100)

    if df_row is not None and len(df_row) > 0:
        # fac = pymzml.plot.Factory()
        print('-'*100)
        print(df_row)
        print(df_data)
        table_seq = df_data[df_row[0]]['Sequence']
        table_mods = df_data[df_row[0]]['Modifications']
        table_seq_mods = '#'.join([table_seq, table_mods])
        fragger = peptide_fragmentor.PeptideFragment0r(table_seq_mods)
        df = fragger.df
        if len(df) > 0:
            df = df[df['series'].isin(['y', 'b'])][['seq', 'series', 'mass', 'name']]
            matches, mismatches = match_df_to_mz_i_list(df, spec_peaks, rel_cutoff=20e-6)
            if len(matches) > 0:
                fac.add(
                    matches,
                    style='label.triangle.MS_precision',
                    name='fragment matches',
                    color=(0, 255, 0)
                )
            if len(mismatches) > 0:
                fac.add(
                    mismatches,
                    style='label.triangle.MS_precision',
                    name='fragment mismatches',
                    color=(255, 0, 0)
                )

            data_tmp = fac.save()
            data['data'] += data_tmp['data']
        print('-'*100)
    if seq is not None:
        print('-'*100)
        print('input anno')
        print('-'*100)
        try:
            # fac = pymzml.plot.Factory()
            fragger = peptide_fragmentor.PeptideFragment0r(seq)
            df = fragger.df
            if len(df) > 0:
                df = df[df['series'].isin(['y', 'b'])][['seq', 'series', 'mass', 'name']]
                matches, mismatches = match_df_to_mz_i_list(df, spec_peaks, rel_cutoff=20e-6)
                if len(matches) > 0:
                    fac.add(
                        matches,
                        style='label.triangle.MS_precision',
                        name='fragment matches',
                        color=(0, 255, 0)
                    )
                if len(mismatches) > 0:
                    fac.add(
                        mismatches,
                        style='label.triangle.MS_precision',
                        name='fragment mismatches',
                        color=(255, 0, 0)
                    )

                data_tmp = fac.save()
                data['data'] += data_tmp['data']
        except NameError as ne:
            print(ne)
    return data


@app.callback(
    Output('table', 'children'),
    [Input('input-ident', 'value'), Input('spec-slider-ele', 'value')]
)
def create_psm_table(path_to_file, spec_id):
    global df
    # print('CREATE_PSM_TABLE')
    # print(path_to_file, spec_id)
    if path_to_file is not None:
        try:
            df
        except NameError:
            # print('READ TABLE')
            df = pd.read_csv(path_to_file)
        if spec_id:
            # print('GIVE SPEC ID')
            # print(spec_id)
            df_tmp = df[df['Spectrum ID'] == spec_id][['Spectrum ID', 'Sequence', 'Modifications', 'q-value']]
            table = dash_table.DataTable(
                id='table-obj',
                columns=[{"name": i, "id": i} for i in df_tmp.columns],
                data=df_tmp.to_dict("rows"),
                row_selectable="multi",
            )
        else:
            table = dash_table.DataTable(id='table-obj',)
    else:
        table = dash_table.DataTable(id='table-obj',)
    print(table)
    return table

@app.callback(
    Output('intermediate-value', 'children'),
    [
        Input('table-obj', "derived_virtual_data"),
        Input('table-obj', "derived_virtual_selected_rows")
    ],
)
def selected_rows(data, rows):
    print('SELECTED_ROWS')
    print(data, rows)
    return rows

@app.callback(
    Output('spec-slider', 'children'),
    [Input('input-mzml', 'value')])
def set_runner(value):
    global run
    global fac
    if value:
        run = pymzml.run.Reader(value)
        spec_count = run.get_spectrum_count()
        fac = pymzml.plot.Factory()
        fac.new_plot(MS_precision=20e-6)
    else:
        spec_count = 0
    # return {'max': spec_count, 'marks': {i: str(i) for i in range(1, spec_count)}}
    return dcc.Slider(id='spec-slider-ele', min=1, max=spec_count)


# @app.callback(
#     Output('intermediate-value2', 'children'),
#     [Input('input-mzml', 'value'), Input('input-ident', 'value')]
# )
# def init_data(mzml_path, ident_path):
#     global_data = {}
#     if mzml_path is not None:
#         run = pymzml.run.Reader(mzml_path)
#         global_data['run'] = mzml_path
#     if ident_path is not None:
#         df = pd.read_csv(ident_path)
#         global_data['df'] = ident_path
#     print(global_data)
#     return global_data


if __name__ == '__main__':
    app.run_server(debug=True)
