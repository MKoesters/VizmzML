# -*- coding: utf-8 -*-
import sys
import glob

import peptide_fragmentor
import pymzml

import numpy as np

import dash
import dash_core_components as dcc
import dash_html_components as html

from dash.dependencies import Input, Output

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


mzml_files = [path for path in glob.glob('../data/*.mzML')]
ident_files = [path for path in glob.glob('../data/*.csv')]


app.layout = html.Div(children=[
    html.H1(
        children='pymzML GUI',
        style={
            'textAlign': 'center',
        }
    ),
    html.Div(
        children='pymzML GUI.',
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
    # dcc.Slider(
    #     id='spec-slider',
    #     min=1,
    #     max=10,
    #     # marks={i: str(i) for i in range(1, 11)},
    #     value=1,
    # ),
    html.Div(id='spec-slider', children=[dcc.Slider(id='spec-slider-ele')]),
    html.Div(id='intermediate-value', style={'display': 'none'}),
    html.Div(id='intermediate-value2', style={'display': 'none'})
])


def match_df_to_mz_i_list(fragment_df, mz_i_list, dalton_cutoff=0.02, rel_cutoff=None, test=False):
    if not isinstance(mz_i_list, np.ndarray):
        mz_i_list = np.ndarray(mz_i_list)
    mz = mz_i_list[:,0]
    matches = []
    mismatches = []
    print('\nMATCHES')
    for i, row in fragment_df.iterrows():
        print(row['mass'])
        if rel_cutoff:
            cutoff = row['mass'] * rel_cutoff
        else:
            cutoff = dalton_cutoff
        res =  (row['mass'], 'max_intensity', row['name'])
        if (abs(mz - row['mass']) < cutoff).any():
            # print(mz[], row['mass'])
            matches.append(res)
        else:
            mismatches.append(res)
        print(cutoff)
        print(mz[abs(mz - row['mass']) < cutoff])
    if test is True:
        matches.append(
            (131.2, 'max_intensity', 'TEST')
        )
    print('\n\nMATCHED\n\n')
    print(matches, mismatches)
    return matches, mismatches


@app.callback(
    Output('spec-graph', 'figure'),
    [Input('spec-slider-ele', 'value'), Input('peptide-input', 'value')])
def update_figure(spec_id, seq):
    fac = pymzml.plot.Factory()
    fac.new_plot(MS_precision=20e-6)
    print('app prec :', fac.MS_precisions)
    data = {'data': [], 'layout': []}
    print('\n\n')
    print('1: ', len(data['data']))
    try:
        spec = run[spec_id]
        spec_peaks = spec.peaks('raw')
        fac.add(data=spec_peaks, name='data spec_id: {0}'.format(spec_id))
        data_tmp = fac.save()
        data['data'] += data_tmp['data']
    except NameError as e:
        pass
    print('2: ', len(data['data']))
    if seq is not None:
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
                    )
                if len(mismatches) > 0:
                    fac.add(
                        mismatches,
                        style='label.triangle.MS_precision',
                        name='fragment mismatches',
                    )

                data_tmp = fac.save()
                data['data'] += data_tmp['data']
        except NameError as ne:
            print(ne)
    print('3: ', len(data['data']))
    return data


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

if __name__ == '__main__':
    app.run_server(debug=True)
