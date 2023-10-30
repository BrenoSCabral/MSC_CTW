import xarray as xr
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sys import path
import time
from scipy.fft import fft, fftfreq
path.append(
    '/Users/breno/Documents/Mestrado/tese/scripts'
)
from filtro import filtra_dados
fig_folder = '/Users/breno/Documents/Mestrado/tese/figs/teste/'


path_dado = '/Users/breno/dados_mestrado/dados/'


# le dados nc:
def le_dado_nc(filename, path_dado = path_dado):
    '''
        filename: (string) nome do arquivo csv

        retorna um Dataset com os dados
    '''
    return (xr.open_dataset(path_dado + filename))

def le_dado_csv(filename, path_dado = path_dado, separador=';', pula_linhas = None):
    '''
        OBS: esse pode mudar de acordo com a estrutura do csv
        filename: (string) nome do arquivo csv

        retorna um Dataset com os dados
    '''
    df_file = pd.read_csv(path_dado + filename,
                          sep=separador, header=None, skiprows=pula_linhas)
    df_file.columns = ['time', 'sea_level','drop']
    df_file.index = pd.to_datetime(df_file['time'], format="%d/%m/%Y %H:%M")
    df_file = df_file.drop('drop', axis=1)
    df_file = df_file.drop('time', axis=1)

    return(xr.Dataset.from_dataframe(df_file))

def corta_dado(dataini, datafim, dado):
    '''
        dataini: (string) data incial do corte
        datafim: (string) data final do corte
        dado: (Dataset) xarray do daddo que vc deseja cortar

        retorna um DataArray com os dados cortados e outro com
            as datas utilizadas
    '''
    variaveis = list(dado.data_vars)
    variavel = ' '
    while variavel not in variaveis:
        print('-', '\n- '.join(variaveis))
        variavel = input('\nQual das variáveis acima você deseja acessar? ')
        if variavel not in variaveis:
            print("\n",
                  "Opção inválida. Por favor, escolha uma das variáveis listadas.",
                  "\n")
            time.sleep(1)
    dado_cortado = dado[variavel].sel(time=slice(dataini, datafim))
    datas_cortadas = dado.time.sel(time=slice(dataini, datafim))

    return(dado_cortado, datas_cortadas)

def reamostra_dado(dado, dado_tempo):
    '''
        Pega um array numpy de dado e outro de tempo e reamostra para 1 dia.
        dado: (numpy array) array com dados (de ssh)
        dado_tempo: (numpy array) array de datas
    '''

    da = xr.DataArray(data=dado,
                      dims=['time'],
                      coords=dict(time=dado_tempo))
    '''
    metodos de interpolacao:
    https://xarray.pydata.org/en/v0.14.1/generated/xarray.DataArray.resample.html
    '''
    da_diario = da.resample(time='1D').interpolate('cubic')
    return da_diario

def pega_fft(ssh, nome_serie, fig_folder):
    N = 365*24 # numero de medicoes
    T = 1.0 / 24 # frequencia das medicoes (nesse caso = 1 medicao a cada 24h)
    yf = fft(ssh) # faz a transformada de fourier
    xf = fftfreq(N, T) # freequencia de amostragem da transformada
    # por algum motivo, o exemplo usa esse xfa, mas eu nao achei que valeria a pena
    # xfa = xf[:N//2]
    # tfa = 1/xfa
    tf = 1/xf

    # fig, [ax1, ax2] = plt.subplots(nrows=2, ncols=1)
    # ax1.semilogx(tfa, 2.0/N * np.abs(yf[0:N//2]))
    # ax1.set_ylim([-0.1, 40])
    # ax1.grid()
    fig = plt.figure(figsize=(20, 7))
    plt.semilogx(tf, 2.0/N * np.abs(yf))
    plt.ylim([-0.1, 40])
    plt.ylabel('Energia')
    plt.xlabel('Dias')
    plt.title('Decomposição das frequências da série de '+nome_serie)
    plt.grid()
    plt.savefig(f'{fig_folder}{nome_serie}_frequencias.png')

def pega_fft_diario(ssh, nome_serie, fig_folder):
    N = 365 # numero de medicoes
    T = 1.0 # frequencia das medicoes (nesse caso = 1 medicao a cada 24h)
    yf = fft(ssh) # faz a transformada de fourier
    xf = fftfreq(N, T) # freequencia de amostragem da transformada
    # por algum motivo, o exemplo usa esse xfa, mas eu nao achei que valeria a pena
    # xfa = xf[:N//2]
    # tfa = 1/xfa
    tf = 1/xf

    # fig, [ax1, ax2] = plt.subplots(nrows=2, ncols=1)
    # ax1.semilogx(tfa, 2.0/N * np.abs(yf[0:N//2]))
    # ax1.set_ylim([-0.1, 40])
    # ax1.grid()
    fig = plt.figure(figsize=(20, 7))
    plt.semilogx(tf, 2.0/N * np.abs(yf))
    plt.ylim([-0.1, 40])
    plt.ylabel('Energia')
    plt.xlabel('Dias')
    plt.title('Decomposição das frequências da série de '+nome_serie)
    plt.grid()
    plt.savefig(f'{fig_folder}{nome_serie}_frequencias.png')

def roda_analise(filename, formato, nome, metodo, fig_folder=fig_folder,
                 dataini = None, datafim = None):
    formatos = ['csv', 'nc']
    while formato not in formatos:
        formato = input('Vc deseja abrir um csv ou um nc? ')
        formato = formato.lower()
        if formato not in formatos:
            print('\nfavor escolher um formato válido.\n')
            time.sleep(1)
    if formato == 'csv':
        dado = le_dado_csv(filename)
    elif formato == 'nc':
        dado = le_dado_nc(filename)
    
    if dataini == None:
        dataini = dado.time[0].values
    if datafim == None:
        datafim = dado.time[-1].values

    dado_ssh, dado_tempo = corta_dado(dataini, datafim, dado)

    dado_filtrado = filtra_dados(dado_ssh, dado_tempo, nome, fig_folder, metodo)

    dado_reamostrado = reamostra_dado(dado_filtrado, dado_tempo)

    return dado_reamostrado



################# FIM DO SCRIPT #################
################# DAQUI PRA BAIXO EH SO RASCUNHO #################

# ##### AINDA VOU APAGAR

# df = xr.DataArray(data=mdp_filtrado,
#                   dims=['time'],
#                   coords=dict(time=tempo_mdp)
# )
# ddf = df.resample(time='1D').interpolate('cubic')
# # plt.plot(smdp.sea_level)

# ### INTERPOLACAO:
# '''
# mdp1d = mdp.resample(time='1D').interpolate('linear')
# metodos de interpolacao = https://xarray.pydata.org/en/v0.14.1/generated/xarray.DataArray.resample.html

# '''


# # le dados imbituba 2012
# # # imb = pd.read_csv(path_dado + 'imbituba.csv')
# # imb = xr.open_dataset(path_dado + 'imbituba.nc')

# # nivel_mar_imb = (imb.sea_level.sel(
# #     time=slice('2005-01-01', '2005-12-31')).values[0])/1000

# # tempo_imb = imb.time.sel(
# #     time=slice('2005-01-01', '2005-12-31')).values

# # filtra_dados(nivel_mar_imb, tempo_imb)

# # le dados ilha fiscal 2012
# # ilf = pd.read_csv(path_dado + 'ilha_fiscal.csv')
# ilf = xr.open_dataset(path_dado + 'ilha_fiscal.nc')

# nivel_mar_ilf = (ilf.sea_level.sel(
#     time=slice('2000/01/01', '2013/12/31')).values[0])/1000

# tempo_ilf = ilf.time.sel(
#     time=slice('2000/01/01', '2013/12/31')).values


# ilf_filtrado = filtra_dados(nivel_mar_ilf, tempo_ilf, 'ilha_fiscal')


# #############
# # LE DADO CSV ILHA FISCAL
# # - TEM QUE LIMPAR O DADO
# #############
# #
# # LATITUDE:22∫ 53'.8 S
# # LONGITUDE:043∫ 10'.0 W
# #

# dfilf = pd.read_csv(path_dado + 'Ilha Fiscal 2014.txt', sep=';', header=None)
# dfilf.columns = ['time', 'sea_level','drop']
# dfilf.index = pd.to_datetime(dfilf['time'], format="%d/%m/%Y %H:%M")
# dfilf = dfilf.drop('drop', axis=1)
# dfilf = dfilf.drop('time', axis=1)

# ilf = xr.Dataset.from_dataframe(dfilf)

# ssh_ilf = ilf.sea_level.values/100 # passando pra cm
# time_ilf = ilf.time.values
# ilf_filtrado = filtra_dados(ssh_ilf, time_ilf, 'ilha_fiscal')

# ilfilt = pd.DataFrame({'ssh':ilf_filtrado})
# ilfilt.index = time_ilf



# ax = plt.axes()
# ax.plot(ilfilt['05/17/2014':'08/03/2014'])
# # Set ylim
# ax.set_ylim(-0.6, 0.6)
# # Set xlim
# ax.grid()

# plt.savefig(fig_folder+'compara_pedro.png')


# df = xr.DataArray(data=ilf_filtrado,
#                   dims=['time'],
#                   coords=dict(time=time_ilf)
# )

# ddf = df.resample(time='1D').interpolate('cubic')

# ssh_diario = ddf.sel(time=slice('05/17/2014','08/03/2014'))
# time_plot = pd.date_range('2014-05-17','2014-08-03')
# # time_plot = ilfilt['05/17/2014':'08/03/2014'].index

# ax = plt.axes()
# ax.plot(time_plot, ssh_diario)
# # Set ylim
# ax.set_ylim(-0.6, 0.6)
# # Set xlim
# ax.grid()
# plt.show()

# # plt.savefig(fig_folder+'compara_pedro.png')