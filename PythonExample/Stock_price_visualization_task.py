# Perform imports here:
import plotly.graph_objs as go
import numpy as np
import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from datetime import datetime as dt
from datetime import datetime
import pandas_datareader.data as web
import os
import json

########################################################################################################
# Specify Your IEX API Token here
os.environ["IEX_API_KEY"] = ""
########################################################################################################

df = pd.read_csv("./data/NASDAQcompanylist.csv")
option_dict = [{'label': i, 'value': i} for i in df['Name'].unique()]

# Launch the application:
app = dash.Dash()

# Create a Dash layout that contains input components
# and at least one output. Assign IDs to each component:
app.layout = html.Div([
    html.Div("Stock Ticker Dashboard",
             style={"fontSize":"30", "marginLeft":"10px"}),
    html.Label("Select Stock symbles", style={'width': '30%', 'display': 'inline-block'}),
    html.Label("Select Start and End date", style={'width':' 40%', 'display': 'inline-block'}),
    html.Div(),
    html.Div(
        dcc.Dropdown(id='my_company',
                     options=option_dict,
                     multi=True,
                     value=[]),
        style={'width': '30%', 'display': 'inline-block'}),
    html.Div(dcc.DatePickerRange(id='my_date_picker',
                                 min_date_allowed=datetime(2015, 1, 1),
                                 max_date_allowed=datetime.today(),
                                 start_date=datetime(2018, 1, 1),
                                 end_date=datetime.today()),
             style={'width': '40%', 'display': 'inline-block'}),
    html.Div(html.Button(id='submit_button', n_clicks=0, children='Submit'),
             style={'width': '30%', 'display': 'inline-block'}),
    html.Div(),
    dcc.Graph(id='result_plot',
              figure={'data': go.Scatter(),
                      'layout': go.Layout(hovermode="closest", title="my_plot")}),
    html.Markdown(id="my_explanation", children=""),
])


# Create a Dash callback:
@app.callback(Output('result_plot', 'figure'),
              [Input('submit_button', 'n_clicks')],
              [State('my_company', 'value'),
               State('my_date_picker', 'start_date'),
               State('my_date_picker', 'end_date')])
def update_figure(n_clicks, companies, start_date, end_date):

    traces =[]
    for i in companies:
        df = web.DataReader(i, 'iex', start_date, end_date, access_key=os.environ["IEX_API_KEY"])
        traces.append({'x':df.index, 'y': df.close, 'name':i})

    fig = {
        'data': traces,
        'layout': {'title': ', '.join(companies) + ' Closing Prices'}
    }
    return fig

@app.callback(Output('my_explanation', 'children'),
              [Input('result_plot', 'clickData')])
def update_figure2(selected_data):
    my_json_comment = json.dump(selected_data)

    comment = ""
    for i in range(my_json_comment["points"]):
        comment += "the index for {} is in {} \n".format(my_json_comment["points"][i]['traces'],
                                                      my_json_comment["points"][i]['y'])

    return comment

# Add the server clause:
app.run_server()
