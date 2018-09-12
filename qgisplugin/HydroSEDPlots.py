from wmf import wmf 
import numpy as np
import plotly.graph_objs as go
import pandas as pd
from plotly import tools
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

class PlotRainfall():

    def __init__(self, pathRain):
        self.path_rainBin, self.path_rainHdr = wmf.__Add_hdr_bin_2route__(pathRain) 
        self.rainData = wmf.read_rain_struct(self.path_rainHdr)
    
    def Plot_Rainfall(self, pathFigure):
        '''Hace el plot de la serie de tiempo de lluvia e incluye en esta toda la info necesaria'''
        #obtiene los records
        Rs = self.rainData.copy()
        Records = ['Record: %d'%i for i in Rs[' Record'].values]
        #Datos de lluvia 
        trace1 = go.Scatter(
            x = Rs.index.to_pydatetime(),
            y = Rs[' Lluvia'].values,
            text = Records,
            name = 'Lluvia [mm]',
            line = {'width':3},
            fill='tozeroy'
        )
        #Datos de lluvia acumulada
        trace2 = go.Scatter(
            x = Rs.index.to_pydatetime(),
            y = Rs[' Lluvia'].values.cumsum(),
            name = 'Acumulado [mm]',
            line = {'width':3},
            yaxis = 'y2'
        )
        #Datos y formato de la figura
        data = [trace1, trace2]
        layout = dict(
            width=1020,
            height=400,
            showlegend = False,
            margin=dict(
                l=50,
                r=50,
                b=50,
                t=50,
                pad=4
            ),
            yaxis=dict(
                title='Precipitacion',
                titlefont=dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                ),
                tickangle=-90,
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                ),
            ),
            
            xaxis = dict(
                tickfont = dict(
                    size = 16
                ), 
            ),
            
            yaxis2=dict(
                title='Acumulada',
                titlefont=dict(
                    color='rgb(255, 153, 51)',
                    size = 15
                ),
                tickangle=90,
                tickfont=dict(
                    color='rgb(255, 153, 51)',
                    size = 16,            
                ),
                overlaying='y',
                side='right',
                range = [0,Rs[' Lluvia'].values.sum()]
            ))
        #Figura 
        fig = dict(data=data, layout=layout)
        #Guarda el html
        plot(fig,filename=pathFigure, auto_open = False)
        
    def Plot_Histogram(self, pathFigure):
        '''Hace un plot del histograma de la distribucion de la lluvia en la cuenca.'''
        #Hace una cipia de la informacion
        Data = self.rainData.copy()
        Data = Data[' Lluvia'].values
        Data = Data[Data>0]
        step = (np.percentile(Data,95) - np.percentile(Data,5))/7.
        #Genera los datos de la figura
        trace1 = go.Histogram(
            x = Data,
            name = 'Lluvia [mm]',
            xbins = dict(
                start = np.percentile(Data,5),
                end = np.percentile(Data,95),
                size = step)
        )
        #Establece la configuracion de la misma
        layout = dict(
            width=400,
            height=400,
            showlegend = False,
            margin=dict(
                l=50,
                r=50,
                b=70,
                t=50,
                pad=4
            ),
            yaxis=dict(
                title='PDF',
                titlefont=dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                    ),
                tickangle=45,
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                    ),),
            xaxis = dict(
                title = 'Lluvia [mm]',
                titlefont =dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                    ),
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                    )
                )
            )
        #Monta la figura 
        data = [trace1]
        fig = dict(data = data, layout = layout)
        plot(fig,filename=pathFigure, auto_open = False)
        
    def Plot_MediaMensual(self, pathFigure):
        '''Hace el plot de la media mensual multi-anual de la serie de lluvia de la cuenca'''
        Rs = self.rainData.copy()
        Rain = Rs[' Lluvia'].resample('M').sum()
        #Obtiene media y desviqacion
        Mean = []
        Desv = []
        for i in range(1,13):
            Mean.append(Rain[Rain.index.month == i].mean())
            Desv.append(np.std(Rain[Rain.index.month == i]))
        Mean = np.array(Mean)
        Desv = np.array(Desv)
        #Hace la figura
        trace1 = go.Scatter(
            x = np.arange(1,13,1),
            y = Mean,
            name = 'Media',
            line = {'width':5,
                'color':('rgb(0, 51, 102)')},    
        )
        
        trace2 = go.Scatter(
            x = np.arange(1,13,1),
            y = Mean + Desv,
            name = 'm + s',
            line = dict(color = ('rgb(255, 153, 51)'),
                width = 3),
        )
        
        trace3 = go.Scatter(
            x = np.arange(1,13,1),
            y = Mean - Desv,
            name = 'm - s',
            fill='tonexty',
            line = dict( width=3,
                 color = ('rgb(255, 153, 51)')),
        )
        data = [trace2, trace3, trace1]
        layout = dict(
            width=800,
            height=400,
            showlegend = False,
            margin=dict(
                l=50,
                r=50,
                b=50,
                t=50,
                pad=4
            ),
            xaxis=dict(
                dtick=1,
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16),
                ticktext = ['Ene','Feb','Mar','Abr','May','Jun','Jul','Ago','Sep','Oct','Nov','Dic'],
                tickvals = range(1,13)
            ),
            yaxis=dict(
                title='Precipitacion',
                titlefont=dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                ),
                tickangle=-90,
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                ),
            ),
        )
        fig = dict(data=data, layout=layout)
        plot(fig,filename=pathFigure, auto_open = False)

class PlotCaudal():
    def __init__(self):
        self.a = 0
        
    def Plot_Caudal(self, pathFigure,DataQ,id_est,colorline):
        '''Hace el plot de la serie de tiempo de lluvia e incluye en esta toda la info necesaria'''
        #obtiene los records
        Q = DataQ[id_est].copy()
        #Datos de lluvia 
        trace1 = go.Scatter(
            x = Q.index.to_pydatetime(),
            y = Q.values,
            name = 'Caudal Observado',
            line = dict(color = (colorline),width=3),
            fill='tozeroy'
        )
        #Datos y formato de la figura
        data = [trace1]
        layout = dict(
            width=1020,
            height=400,
            showlegend = False,
            margin=dict(
                l=50,
                r=50,
                b=50,
                t=50,
                pad=4
            ),
            yaxis=dict(
                title='Caudal Observado [m3/s]',
                titlefont=dict(
                    color='black',
                    size = 15
                ),
                tickangle=-90,
                tickfont=dict(
                    color='black',
                    size = 16,            
                ),
            ),
            
            xaxis = dict(
                tickfont = dict(
                    size = 16
                ), 
            ),
            )
        #Figura 
        fig = dict(data=data, layout=layout)
        #Guarda el html
        plot(fig,filename=pathFigure, auto_open = False)
        
    def Plot_Caudal_Simu(self,pathFigure,dfQobs,dfQsim,rain_mean):
    
        trace_high = go.Scatter(
                        x=dfQsim.index,
                        y=dfQsim.values,
                        name = "Q Simulado",
                        line = dict(color = 'red'),
                        opacity = 0.8)

        trace_low = go.Scatter(
                        x=dfQsim.index,
                        y=dfQobs.values,
                        name = "Q Observado",
                        line = dict(color = 'blue'),
                        opacity = 0.8)
        trace_inv = go.Scatter(
                        x=dfQsim.index,
                        y=rain_mean.values,
                        name = "Q Observado",
                        line = dict(color = 'blue'),
                        opacity = 0.8,
                        yaxis='y2',
                        fill='tozeroy',
                        fillcolor='rgba(0, 102, 153,0.2)')

        data = [trace_high,trace_low,trace_inv]

        layout = dict(showlegend = False,
            xaxis = dict(
                title='Fechas'),
            yaxis=dict(
                title="Caudal [m3/s]"),
            yaxis2=dict(autorange="reversed",
                title="Caudal [m3/s]",
                overlaying ='y',
                side='right')
            
        )

        fig = dict(data=data, layout=layout)
        #Guarda el html
        plot(fig,filename=pathFigure, auto_open = False)

    def Plot_CDC_caudal(self,pathFigure,dfQobs,dfQsim):
        
        Qs=np.sort(np.array(dfQsim.values))
        Qo=np.sort(np.array(dfQobs.values))
        porcen_s=[]
        porcen_o=[]
        for i in range(len(Qo)):
            porcen_s.append((len(Qs[Qs>Qs[i]]))/float(len(Qs))*100)
            porcen_o.append((len(Qo[Qo>Qo[i]]))/float(len(Qo))*100)
            
        trace_high = go.Scatter(
                x=porcen_s,
                y=Qs,
                name = "Q simulado",
                line = dict(color = 'red'),
                opacity = 0.8)

        trace_low = go.Scatter(
                        x=porcen_o,
                        y=Qo,
                        name = "Q observado",
                        line = dict(color = 'blue'),
                        opacity = 0.8)

        data = [trace_high,trace_low]

        layout = dict(showlegend = False,
            width=500,
            height=400,
            xaxis = dict(
                title='Porcentaje de Excedencia'),
            yaxis=dict(
                title='$Caudal [m^3/s]$')
    
        )
        fig = dict(data=data, layout=layout)
        #Guarda el html
        plot(fig,filename=pathFigure, auto_open = False)    
        
    def Plot_Series_Sedimentos(self,pathFigure,dfArenas,dfLimos,dfArcillas,dfSed_total,dfDataQobs):   
             
        trace_high = go.Scatter(
                        x=dfArenas.index,
                        y=dfArenas.values,
                        name = "Arenas",
                        line = dict(color = 'rgb(255, 255, 0)',dash = 'dashdot'),
                        opacity = 0.8)

        trace_1 = go.Scatter(
                        x=dfArenas.index,
                        y=dfLimos.values,
                        name = "Limos",
                        line = dict(color = 'rgb(179, 143, 0)',dash = 'dash'),
                        opacity = 0.8)
                        
        trace_2 = go.Scatter(
                        x=dfArenas.index,
                        y=dfArcillas.values,
                        name = "Arcillas",
                        line = dict(color = 'rgb(85, 128, 0)',dash = 'dot'),
                        opacity = 0.8)
                        
        trace_3 = go.Scatter(
                        x=dfArenas.index,
                        y=dfSed_total.values,
                        name = "Qsed total",
                        line = dict(color = 'red',width=3),
                        opacity = 0.8)
                        
        trace_4 = go.Scatter(
                        x=dfArenas.index,
                        y=dfDataQobs.values,
                        name = "Qsed Observado",
                        line = dict(color = ('rgb(255, 71, 26)'),width=3),
                        opacity = 0.8)        
                                                   

        data = [trace_high,trace_1,trace_2,trace_3,trace_4]

        layout = dict(showlegend = False,
            xaxis = dict(
                title='Fechas'),
            yaxis=dict(
                title="Caudal Solido [m3/s]"),
            
        )

        fig = dict(data=data, layout=layout)
        #Guarda el html
        plot(fig,filename=pathFigure, auto_open = False)
            
   
class PlotGeomorphology():
    
    def __init__(self):
        self.a = 0
    
    def VarHistogram(self, var, pathFigure, ncells):
        '''Hace el plot de un histograma de una variable seleccionada'''
        #Seleccion de datos
        if ncells > 10000:
            pos = np.random.choice(len(var), 10000)
        else:
            pos = range(10000)
        #Obtiene la pdf y cdf de la variable    
        var = var[pos]
        p1 = np.percentile(var, 0.1)
        p99 = np.percentile(var, 99.9)
        h,b = np.histogram(var, bins=np.linspace(p1,p99,10))
        pdf = h.astype(float)/h.sum()
        cdf = pdf.cumsum()
        b = (b[:-1]+b[1:])/2.
        #Create a trace
        trace1 = go.Scatter(
            x = b,
            y = pdf,
            name = 'PDF',
            line = {'width':3},
            fill='tozeroy'
        )
        trace2 = go.Scatter(
            x = b,
            y = cdf,
            name = 'CDF',
            line = {'width':3},
            yaxis = 'y2'
        )
        #Datos del plot
        data = [trace1, trace2]
        #Configuracion del plot
        layout = dict(
            width=400,
            height=400,
            showlegend = False,
            margin=dict(
                l=50,
                r=50,
                b=50,
                t=50,
                pad=4
            ),
            xaxis = dict(
                titlefont =dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                    ),
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                    )
                ),
            yaxis=dict(
                title='PDF',
                titlefont=dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                ),
                tickangle=-90,
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                ),
            ),
            
            yaxis2=dict(
                title='CDF',
                titlefont=dict(
                    color='rgb(255, 153, 51)',
                    size = 15
                ),
                tickangle=90,
                tickfont=dict(
                    color='rgb(255, 153, 51)',
                    size = 16,            
                ),
                overlaying='y',
                side='right',
            ))
        #Figura y guardado
        fig = dict(data=data, layout=layout)
        plot(fig,filename=pathFigure, auto_open = False)
    
    def TiempoConcentracion(self, DicTiempos, pathFigure):
        '''Plot de los tiempos de concentracion de la cuenca'''
        #Los datos
        Names = [i[0][:7] for i in DicTiempos.items()]
        Texto = [i[0] for i in DicTiempos.items()]
        Times = ['%.3f'%i[1] for i in DicTiempos.items()]
        #El plot
        data = [go.Bar(
                    x = Names,
                    y = Times,
                    text = Texto
            )]
        #Configuracion del plot
        layout = dict(
            width=600,
            height=400,
            showlegend = False,
            margin=dict(
                l=50,
                r=50,
                b=50,
                t=50,
                pad=4
            ),          
            xaxis = dict(
            title = 'Variable',
            titlefont =dict(
                color='rgb(0, 102, 153)',
                size = 15
                ),
            tickfont=dict(
                color='rgb(0, 102, 153)',
                size = 16,            
                )
            ),
            yaxis=dict(
                title='Tiempo viaje [hrs]',
                titlefont=dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                ),
                tickangle=0,
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                ),
            ),)
        #Figura y guardado
        fig = dict(data=data, layout=layout)
        plot(fig,filename=pathFigure, auto_open = False)
    
    def StreamProfileHipsometric(self, pathFigure, HipsoData, StreamData):
        '''Hace el plot del perfil longitudinal y de la curva hipsometrica'''
        #plot# Create a trace
        trace1 = go.Scatter(
            x = HipsoData[0], 
            y = HipsoData[1],
            name = 'Curva hipso',
            line = {'width':3},
        )
        x = StreamData[0]
        y = StreamData[1]
        trace2 = go.Scatter(
            x = x[30:],
            y = y[30:],
            name = 'Perfil',
            line = {'width':3},
            xaxis = 'x2'
        )
        
        data = [trace1, trace2]
        
        layout = dict(
            width=600,
            height=400,
            showlegend = False,
            margin=dict(
                l=50,
                r=50,
                b=50,
                t=50,
                pad=4
            ),
            
            xaxis = dict(
            title = 'Porcentaje del area [%]',
            titlefont =dict(
                color='rgb(0, 102, 153)',
                size = 15
                ),
            tickfont=dict(
                color='rgb(0, 102, 153)',
                size = 16,            
                )
            ),
            
            yaxis=dict(
                title='Elevacion [m.s.n.m]',
                titlefont=dict(
                    color='rgb(0, 102, 153)',
                    size = 15
                ),
                tickangle=-90,
                tickfont=dict(
                    color='rgb(0, 102, 153)',
                    size = 16,            
                ),
            ),
            
            xaxis2=dict(
                title='Distancia a la salida [km]',
                titlefont=dict(
                    color='rgb(255, 153, 51)',
                    size = 15
                ),
                tickfont=dict(
                    color='rgb(255, 153, 51)',
                    size = 16,            
                ),
                overlaying='x',
                side='top',
            ))
        fig = dict(data=data, layout=layout)
        plot(fig,filename=pathFigure, auto_open = False)
    
class PlotStorage():
    
    def __init__(self):
        self.foo = 1
    
    def Plot_Storages(self, StoragePath, PathFigure):
        '''Hace un plot del storage en el periodo de simulacion'''
        #Lee los datos 
        Data = pd.read_csv(StoragePath, skiprows=4, index_col=6, parse_dates=True)
        #Hace la figura
        fig = tools.make_subplots(rows=5, cols=1)
        for c,key in enumerate(Data.columns.values[1:].tolist()):
            trace1 = go.Scatter(
                x = Data.index.to_pydatetime(),
                y = Data[key].values,
                name = key,
                line = {'width':3},
                fill='tozeroy',
            )
            fig.append_trace(trace1, c+1, 1)
        fig['layout'].update(height=600, width=600,
            showlegend = False,
            yaxis=dict(title='Estado [mm]',),
            margin=dict(
                l=50,
                r=50,
                b=50,
                t=50,
                pad=4
            ))
        plot(fig,filename=PathFigure, auto_open = False)
    
