import xarray as xr
import numpy as np
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import itertools
import cartopy
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

path_reanalise = '/Users/breno/dados_mestrado/dados/reanalise/GLOR4/'
figs_folder = '/Users/breno/Documents/Mestrado/tese/figs/teste/'
# long1 = (302.47 + 180) % 360 - 180

##################################################################
#                           FUNCTIONS                            #
##################################################################

# BRAN: TEMPO = 365 - LAT = 951 - LON = 1401
# HYCOM: TEMPO: 365 - LAT = 1751 - LON = 1751
# GLOR4: TEMPO = 365, LAT = 401 - LON = 561

# ----> dataset.sel({time_dim: time, lat_dim: lat, lon_dim: lon}) <---

def read_reanalisys(year, path_reanalise):
    '''
        Opens the reanalisys for the year of interest.
        :path_reanalise: str with the path of the reanalisys
        :year: year of interest
    '''
    reanalisys = xr.open_mfdataset(f'{path_reanalise}{year}/*.nc')

    return reanalisys

def cut_reanalisys(reanalisys, ti, tf, lat, lon, dimensions):
    '''
        Makes a selection from the realanilsys domain in time and space
        :reanalisys: Dataarray
        :ti: str - initial time in the format yyyy-mm-dd
        :tf: str - final time in the format yyyy-mm-dd
        :lat: float - latitude of the poin of interest
        :lon: float - longitude of the point of interest
    '''
    variable = list(reanalisys.data_vars)[0]

    if dimensions['lat'] == 'lat' and dimensions['lon'] == 'lon':
        ssh_raw = reanalisys[variable].sel(lat=lat,
                                    lon=lon,
                                    method='nearest')
    elif dimensions['lat'] == 'latitude' and dimensions['lon'] == 'longitude':
        ssh_raw = reanalisys[variable].sel(latitude=lat,
                                    longitude=lon,
                                    method='nearest')
    elif dimensions['lat'] == 'yt_ocean' and dimensions['lon'] == 'xt_ocean':
        ssh_raw = reanalisys[variable].sel(yt_ocean=lat,
                                    xt_ocean=lon,
                                    method='nearest')        
    else:
        raise('Check latitude/longitude dimensions name.')

    try:
        ssh = ssh_raw.sel(time=slice(ti, tf))
    except:
        ssh = ssh_raw.sel(Time=slice(ti, tf))

    return ssh

def get_reanalisys_dims(reanalisys):
    for i in list(reanalisys.dims):
        if i[:2] == 'la':
            lat_dim = i
        elif i[:2] == 'lo':
            lon_dim = i
        elif i[:2] == 'ti':
            time_dim = i
    return({'time':time_dim, 'lat':lat_dim, 'lon':lon_dim})

def get_lat_lon(reanalisys, dimensions):
    if dimensions['lat'] == 'lat' and dimensions['lon'] == 'lon':
        return(reanalisys.lat.values, reanalisys.lon.values)
    elif dimensions['lat'] == 'latitude' and dimensions['lon'] == 'longitude':
        return(reanalisys.latitude.values, reanalisys.longitude.values)
    else:
        raise('Check latitude/longitude dimensions name.')

def plot_domain(lats, lons, nome_modelo, figs_folder):  
    proj = ccrs.PlateCarree()

    fig, ax = plt.subplots(subplot_kw=dict(projection=proj), figsize=(12,12))

    lon_max, lon_min, lat_max, lat_min = lons[-1], lons[0], lats[-1], lats[0]
    ax.set_extent([lon_max, lon_min, lat_max, lat_min], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND, facecolor='0.3')
    ax.add_feature(cfeature.LAKES, alpha=0.9)  
    ax.add_feature(cfeature.BORDERS, zorder=10)
    ax.add_feature(cfeature.COASTLINE, zorder=10)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',  name='admin_1_states_provinces_lines',
            scale='50m', facecolor='none')
    ax.add_feature(states_provinces, edgecolor='black', zorder=10)   

    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right=True
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator(lons)
    gl.ylocator = mticker.FixedLocator(lats)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'red'}
    ax.set_title('Domínio do ' + nome_modelo)

    plt.savefig(figs_folder + f'{nome_modelo}_dominio.png')
    plt.close()

    # plt.show()

def plot_point(lat, lon, nome_modelo, nome_ponto, figs_folder, lons, lats):
    plot_point_broad(lat, lon, nome_ponto, figs_folder, lons, lats)
    proj = ccrs.PlateCarree()

    fig, ax = plt.subplots(subplot_kw=dict(projection=proj), figsize=(12,12))

    lon_max, lon_min, lat_max, lat_min = lon + 0.5 , lon - 0.5, lat + 0.5, lat - 0.5
    ax.set_extent([lon_max, lon_min, lat_max, lat_min], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND, facecolor='0.3')
    ax.add_feature(cfeature.LAKES, alpha=0.9)  
    ax.add_feature(cfeature.BORDERS, zorder=10)
    ax.add_feature(cfeature.COASTLINE, zorder=10)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',  name='admin_1_states_provinces_lines',
            scale='50m', facecolor='none')
    ax.add_feature(states_provinces, edgecolor='black', zorder=10) 

    plt.plot(lon, lat,
            color='red', linewidth=2, marker='o',
            transform=ccrs.PlateCarree()
            )  
    plt.text(lon - 0.005, lat - 0.005, nome_ponto,
          horizontalalignment='right', color = 'red', weight = 'bold',
          transform=ccrs.PlateCarree())  

# # to_remove
#     plt.plot(-43.2, -22.88,
#             color='red', linewidth=2, marker='o',
#             transform=ccrs.PlateCarree()
#             )
#     plt.text(-43.2 - 0.005, -22.88 - 0.005, 'ponto mais proximo hycom (s/ dados)',
#           horizontalalignment='right', color = 'red', weight = 'bold',
#           transform=ccrs.PlateCarree())
#     plt.plot(-43.12, -22.96,
#             color='red', linewidth=2, marker='o',
#             transform=ccrs.PlateCarree()
#             )
#     plt.text(-43.12 - 0.005, -22.96 - 0.005, 'ponto sem dados',
#           horizontalalignment='right', color = 'red', weight = 'bold',
#           transform=ccrs.PlateCarree())
#     plt.plot(-43.12, -23.04,
#             color='green', linewidth=2, marker='o',
#             transform=ccrs.PlateCarree()
#             )
#     plt.text(-43.12 - 0.005, -23.04 - 0.005, 'ponto mais proximo com dados',
#           horizontalalignment='right', color = 'green', weight = 'bold',
#           transform=ccrs.PlateCarree())


    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5,
                      color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right=True
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator(lons)
    gl.ylocator = mticker.FixedLocator(lats)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'red'}
    ax.set_title(f'Pontos do {nome_modelo} vs ponto {nome_ponto}')

    plt.tight_layout()
    plt.savefig(figs_folder + f'{nome_modelo}vs{nome_ponto}.png')

def plot_point_broad(lat, lon, nome_ponto, figs_folder, lons, lats):
    proj = ccrs.PlateCarree()

    fig, ax = plt.subplots(subplot_kw=dict(projection=proj), figsize=(12,12))

    lon_max, lon_min, lat_max, lat_min = lon + 20 , lon - 20, lat + 15, lat - 15
    ax.set_extent([lon_max, lon_min, lat_max, lat_min], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND, facecolor='0.3')
    ax.add_feature(cfeature.LAKES, alpha=0.9)  
    ax.add_feature(cfeature.BORDERS, zorder=10)
    ax.add_feature(cfeature.COASTLINE, zorder=10)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',  name='admin_1_states_provinces_lines',
            scale='50m', facecolor='none')
    ax.add_feature(states_provinces, edgecolor='black', zorder=10) 


    plt.plot(lon, lat,
            color='red', linewidth=2, marker='o',
            transform=ccrs.PlateCarree()
            )  
    plt.text(lon - 0.05, lat - 0.05, nome_ponto,
          horizontalalignment='right', color = 'red', weight = 'bold',
          transform=ccrs.PlateCarree())    

    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5,
                      color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylabels_right=True
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator(lons)
    gl.ylocator = mticker.FixedLocator(lats)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'red'}
    ax.set_title(f'Ponto {nome_ponto}')

    plt.tight_layout()
    plt.savefig(figs_folder + f'{nome_ponto}_dominio.png')

def roda_tudo(year, lat, lon, model_name, point_name, ti=None, tf=None,
              path_reanalise = path_reanalise, figs_folder = figs_folder):
    
    reanalisys = read_reanalisys(str(year), path_reanalise)

    dimensions = get_reanalisys_dims(reanalisys)

    lats, lons = get_lat_lon(reanalisys, dimensions)
    plot_domain(lats, lons, model_name, figs_folder)
    plot_point(lat, lon, model_name, point_name, figs_folder, lons, lats)
    plt.close('all')

    if ti == None:
        ti = reanalisys.time.values[0]
    if tf == None:
        tf = reanalisys.time.values[-1]

    reanalisys_point = cut_reanalisys(reanalisys, ti, tf, lat, lon, dimensions)
    return (reanalisys_point)


################# FIM DO SCRIPT #################
################# DAQUI PRA BAIXO EH SO RASCUNHO #################

# # precisa melhorar isso aqui
# mercator = xr.open_mfdataset(path_reanalise + 'GLOR12/2012/*.nc')
# nivel_mar_mercator_mdp = (mercator.zos.sel(latitude=-38.035,
#                                            longitude=-57.52999999999997, method='nearest'))

# nivel_mar_mercator_mdp = nivel_mar_mercator_mdp.sel(time=slice('2012-01-01', '2012-12-31'))

# ssh_mdp = nivel_mar_mercator_mdp

# lats = mercator.latitude.values
# lons = mercator.longitude.values
# latlon = list(itertools.product(lats, lons))

# ###############
# #HYCOM####
# ##############

# hycom = xr.open_mfdataset(path_reanalise+'/HYCOM/2014/*.nc')
# lats = hycom.lat.values
# lons = hycom.lon.values
# if_lat = -22.897
# if_lon = -43.16500000000002

# ifh = hycom.surf_el.sel(lat=-23.04, lon=-43.12, method='nearest')
# ifhp = ifh.sel(time=slice('05/17/2014','08/03/2014'))
# ifhpd = ifhp.values
# tifhp = ifhp.time.values


# #### plots

# import cartopy
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.pyplot as plt
# import matplotlib.ticker as mticker

# from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# figs_folder = '/Users/breno/Documents/Mestrado/tese/figs/rep1/'

# def plot_domain(lats, lons, nome_modelo, figs_folder = figs_folder):  
#     proj = ccrs.PlateCarree()

#     fig, ax = plt.subplots(subplot_kw=dict(projection=proj), figsize=(12,12))

#     lon_max, lon_min, lat_max, lat_min = lons[-1], lons[0], lats[-1], lats[0]
#     ax.set_extent([lon_max, lon_min, lat_max, lat_min], crs=ccrs.PlateCarree())

#     ax.add_feature(cfeature.LAND, facecolor='0.3')
#     ax.add_feature(cfeature.LAKES, alpha=0.9)  
#     ax.add_feature(cfeature.BORDERS, zorder=10)
#     ax.add_feature(cfeature.COASTLINE, zorder=10)

#     states_provinces = cfeature.NaturalEarthFeature(
#             category='cultural',  name='admin_1_states_provinces_lines',
#             scale='50m', facecolor='none')
#     ax.add_feature(states_provinces, edgecolor='black', zorder=10)   

#     gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
#     gl.xlabels_top = False
#     gl.ylabels_left = False
#     gl.ylabels_right=True
#     gl.xlines = True
#     gl.xlocator = mticker.FixedLocator(lons)
#     gl.ylocator = mticker.FixedLocator(lats)
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlabel_style = {'color': 'red'}
#     ax.set_title('Domínio do ' + nome_modelo)

#     plt.savefig(figs_folder + f'{nome_modelo}_dominio.png')
#     plt.close()

#     # plt.show()

# # long1 = (302.47 + 180) % 360 - 180

# def plot_point(lat, lon, nome_modelo, nome_ponto, figs_folder = figs_folder):
#     plot_point_broad(lat, lon, nome_modelo, nome_ponto, figs_folder = figs_folder)
#     proj = ccrs.PlateCarree()

#     fig, ax = plt.subplots(subplot_kw=dict(projection=proj), figsize=(12,12))

#     lon_max, lon_min, lat_max, lat_min = lon + 0.5 , lon - 0.5, lat + 0.5, lat - 0.5
#     ax.set_extent([lon_max, lon_min, lat_max, lat_min], crs=ccrs.PlateCarree())

#     ax.add_feature(cfeature.LAND, facecolor='0.3')
#     ax.add_feature(cfeature.LAKES, alpha=0.9)  
#     ax.add_feature(cfeature.BORDERS, zorder=10)
#     ax.add_feature(cfeature.COASTLINE, zorder=10)

#     states_provinces = cfeature.NaturalEarthFeature(
#             category='cultural',  name='admin_1_states_provinces_lines',
#             scale='50m', facecolor='none')
#     ax.add_feature(states_provinces, edgecolor='black', zorder=10) 

#     plt.plot(lon, lat,
#             color='red', linewidth=2, marker='o',
#             transform=ccrs.PlateCarree()
#             )  
#     plt.text(lon - 0.005, lat - 0.005, nome_ponto,
#           horizontalalignment='right', color = 'red', weight = 'bold',
#           transform=ccrs.PlateCarree())  

# # # to_remove
# #     plt.plot(-43.2, -22.88,
# #             color='red', linewidth=2, marker='o',
# #             transform=ccrs.PlateCarree()
# #             )
# #     plt.text(-43.2 - 0.005, -22.88 - 0.005, 'ponto mais proximo hycom (s/ dados)',
# #           horizontalalignment='right', color = 'red', weight = 'bold',
# #           transform=ccrs.PlateCarree())
# #     plt.plot(-43.12, -22.96,
# #             color='red', linewidth=2, marker='o',
# #             transform=ccrs.PlateCarree()
# #             )
# #     plt.text(-43.12 - 0.005, -22.96 - 0.005, 'ponto sem dados',
# #           horizontalalignment='right', color = 'red', weight = 'bold',
# #           transform=ccrs.PlateCarree())
# #     plt.plot(-43.12, -23.04,
# #             color='green', linewidth=2, marker='o',
# #             transform=ccrs.PlateCarree()
# #             )
# #     plt.text(-43.12 - 0.005, -23.04 - 0.005, 'ponto mais proximo com dados',
# #           horizontalalignment='right', color = 'green', weight = 'bold',
# #           transform=ccrs.PlateCarree())


#     gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5,
#                       color='black', alpha=0.5, linestyle='--', draw_labels=True)
#     gl.xlabels_top = False
#     gl.ylabels_left = False
#     gl.ylabels_right=True
#     gl.xlines = True
#     gl.xlocator = mticker.FixedLocator(lons)
#     gl.ylocator = mticker.FixedLocator(lats)
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlabel_style = {'color': 'red'}
#     ax.set_title(f'Pontos do {nome_modelo} vs ponto {nome_ponto}')

#     plt.tight_layout()
#     plt.savefig(figs_folder + f'{nome_modelo}vs{nome_ponto}.png')

# def plot_point_broad(lat, lon, nome_modelo, nome_ponto, figs_folder = figs_folder):
#     proj = ccrs.PlateCarree()

#     fig, ax = plt.subplots(subplot_kw=dict(projection=proj), figsize=(12,12))

#     lon_max, lon_min, lat_max, lat_min = lon + 20 , lon - 20, lat + 15, lat - 15
#     ax.set_extent([lon_max, lon_min, lat_max, lat_min], crs=ccrs.PlateCarree())

#     ax.add_feature(cfeature.LAND, facecolor='0.3')
#     ax.add_feature(cfeature.LAKES, alpha=0.9)  
#     ax.add_feature(cfeature.BORDERS, zorder=10)
#     ax.add_feature(cfeature.COASTLINE, zorder=10)

#     states_provinces = cfeature.NaturalEarthFeature(
#             category='cultural',  name='admin_1_states_provinces_lines',
#             scale='50m', facecolor='none')
#     ax.add_feature(states_provinces, edgecolor='black', zorder=10) 


#     plt.plot(lon, lat,
#             color='red', linewidth=2, marker='o',
#             transform=ccrs.PlateCarree()
#             )  
#     plt.text(lon - 0.05, lat - 0.05, nome_ponto,
#           horizontalalignment='right', color = 'red', weight = 'bold',
#           transform=ccrs.PlateCarree())    

#     gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5,
#                       color='black', alpha=0.5, linestyle='--', draw_labels=True)
#     gl.xlabels_top = False
#     gl.ylabels_left = False
#     gl.ylabels_right=True
#     gl.xlines = True
#     gl.xlocator = mticker.FixedLocator(np.arange(
#        lon_min, lon_max, 1
#     ))
#     gl.ylocator = mticker.FixedLocator(np.arange(
#        lat_min, lat_max, 1
#     ))
#     gl.xformatter = LONGITUDE_FORMATTER
#     gl.yformatter = LATITUDE_FORMATTER
#     gl.xlabel_style = {'color': 'red'}
#     ax.set_title(f'Ponto {nome_ponto}')

#     plt.tight_layout()
#     plt.savefig(figs_folder + f'{nome_ponto}_dominio.png')

# plot_point(-38.035, -57.52999999999997, 'mercator', 'plata')


# # IFHPD1 eh o dado modelado, mas deslocado em 1 dia
# fig= plt.figure(figsize=(16,4))
# trans = fig.transFigure
# plt.title('Dado Medido (filtrado) x Modelado')
# plt.plot(tifhp, ssh_diario, label='Dado medido reamostrado')
# plt.plot(tifhp, ifhpd1, label='Modelado')
# plt.legend()
# plt.text(0.60, 0.85, 'CORR: 80%', transform = trans, color = 'red',
#          weight = 'bold')
# plt.grid()

# plt.tight_layout()
# plt.savefig(figs_folder + 'comparacao_dado_modelo_deslocado.png')


# plt.close('all')

# #### correlacao = np.corrcoef(arr1, arr2)