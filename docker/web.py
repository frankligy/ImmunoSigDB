import os,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import gseapy as gp
from Bio.Align import substitution_matrices
from ast import literal_eval
from dash import Dash, Input, Output, callback, dash_table, html, dcc, State
import subprocess
import dash

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

@callback(Output('result_table','data'),
          Output('gsea_img','src'),
          State('peptide_input','value'),
          State('db_type','value'),
          Input('submit_button','n_clicks'))
def run_gsea(peptides,db,n_clicks):

    if n_clicks > 0:
        # textarea is a string with /n
        peptides = peptides.split('\n')
        data = []
        n = len(peptides)
        for i,peptide in enumerate(peptides):
            data.append((peptide,n-i))
        rnk = pd.DataFrame.from_records(data,columns=['peptide','score']).set_index(keys='peptide')

        final = pd.read_csv('/static/final.txt',sep='\t')
        mapping = {
            'phenotype':'phenotype',
            'distance':'distance_to_star',
            'mSigDB':'mSigDB'
        }
        if db != 'All':
            final = final.loc[final['category']==mapping[db],:]
        gene_sets = {row.ontology:row.peps.split(',') for row in final.itertuples()}

        pre_res = gp.prerank(rnk=rnk,
                             gene_sets=gene_sets,
                             threads=4,
                             min_size=1,
                             max_size=1000,
                             permutation_num=100,
                             outdir=None,
                             seed=6,
                             verbose=True)
        pre_res.res2d.to_csv('/static/result.txt',sep='\t')
        from gseapy import gseaplot2
        terms = pre_res.res2d['Term'][0:5]
        gseaplot2(terms=terms,
                RESs=[pre_res.results[term]['RES'] for term in terms],
                hits=[pre_res.results[term]['hits'] for term in terms],
                rank_metric=pre_res.ranking,
                ofname='/static/result.png')
        
        return pre_res.res2d.to_dict('records'), '/static/result.png'

    else:
        return dash.no_update, dash.no_update

@callback(Output('downloader','data'),
          Input('download_button','n_clicks'))
def download(n_clicks):

    result_txt = '/static/result.txt'
    result_png = '/static/result.png'
    zip_path = '/static/ImmunoSigDB_results.zip'


    if n_clicks > 0:
        if os.path.exists(zip_path):
            os.remove(zip_path)
        subprocess.run(['zip','-j',zip_path,result_txt,result_png])

        return dcc.send_file(zip_path)
    
    else:
        return dash.no_update
    





app = Dash(__name__,assets_folder='/assets')
server = app.server
app.title = 'ImmunoSigDB'

app.layout = html.Div(
    [
        html.Div(
        [
            html.Div(html.Img(src='/static/logo.png',className='actual_logo'), className='logo'),
            html.Div(html.P(['Developed by ', html.A('Guangyuan(Frank) Li',href='https://github.com/frankligy',target='_blank')]),className='develper'),
            html.Div(dcc.Textarea(id='peptide_input',className='peptide_textarea'),className='textarea'),
            html.Div(dcc.Dropdown(id='db_type',className='peptide_dropdown',options=['All','phenotype','distance','mSigDB','differential'],value='All'),className='dropdown'),
            html.Div(html.Button('Run',id='submit_button',n_clicks=0,className='actual_button'),className='run_button')
        ],className='left_div'),

        html.Div(
        [
            html.Div(dash_table.DataTable(id='result_table',
                                          page_size=10,page_current=0,page_action='native',
                                          style_table={'height': '100%','overflowY': 'auto','overflowX': 'auto'}),className='table'),
            html.Div(html.Img(id='gsea_img',className='actual_img'),className='image'),
            html.Div([html.Button('Download',id='download_button',n_clicks=0,className='actual_download_button'),dcc.Download(id='downloader')],className='download_button')
        ],className='right_div')
    ],className='main_div'

)



if __name__ == '__main__':
    app.run(debug=False,host='0.0.0.0',port=int(os.environ.get('PORT',5000)))

