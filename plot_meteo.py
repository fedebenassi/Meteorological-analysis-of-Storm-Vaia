#################################
# Librerie necessarie:
# - numpy e matplotlib (che immagino tu abbia già)
# - netCDF4
# - cartopy versione 0.20
# - wrf-python

# Riferimenti:
# https://wrf-python.readthedocs.io/en/latest/plot.html
# https://wrf-python.readthedocs.io/en/latest/basic_usage.html

## IMPORTAZIONE LIBRERIE ##

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.colors import Normalize
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel)
import numpy as np
from matplotlib.colors import ListedColormap

## CREAZIONE MAPPA PER PLOT ##
# Praticamente scarica da un sito la mappa con risoluzione intermedia.
# Dato che ci serve per tutti i plot che facciamo, la scarico una volta sola
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_0_boundary_lines_land")

## PROIEZIONI GIUSTE DEI DOMINI ##
# Dominio 1 #
ncfile = Dataset("wrfout_d01_2018-10-26_060000")
gpt = getvar(ncfile, "z")
cart_proj_d01 = get_cartopy(interp_gpt)

# Dominio 3  (il 2 per ora non lo usiamo)#
ncfile = Dataset("wrfout_d03_2018-10-26_060000")
gpt = getvar(ncfile, "z")
cart_proj_d02 = get_cartopy(interp_gpt)

###############################################################################

## GEOPOTENZIALE VS TEMPERATURA a 500 hPa# #
# Ciclo for sui giorni disponibili, alle 18:00.

# Definizione della figura: numero di righe, colonne, grandezza e proiezione usata
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 17), subplot_kw = {'projection' : cart_proj_d01})

# Definizione della lista dei giorni
days_list = ["26", "27", "28", "29", "30"]


for i, ax in zip(days_list, axs.ravel()):

# Importazione dei dati delle 18
    namefile = "wrfout_d01_2018-10-"+i+"_180000"
    ncfile = Dataset(namefile)

# Qui prende rispettivamente il geopotenziale, la pressione e la temperatura
# in Celsius
    gpt = getvar(ncfile, "z")
    p = getvar(ncfile, "pressure")
    t = getvar(ncfile, "tc")


# I dati che abbiamo sono tridimensionali (uno strato per ogni height).
# Queste funzioni interpolano i dati ai livelli desiderati, in questo caso 500 hPa
    interp_gpt = interplevel(gpt, p, 500)
    interp_t = interplevel(t, p, 500)

# Prende latitudine e longitudine dai dati che gli abbiamo dato
    lats, lons = latlon_coords(interp_gpt)

# Inizia a disegnare il plot. Qui mette gli stati e la costa
    ax.add_feature(states, linewidth=.6, edgecolor="black", zorder = 2, alpha = 0.8)
    ax.coastlines('50m', linewidth=1, zorder = 2, alpha = 0.8)

# Disegna le linee di geopotenziale
    cs = ax.contour(to_np(lons), to_np(lats), to_np(interp_gpt), 10, colors="black",
            transform=crs.PlateCarree())
# Aggiunge i numeri alle linee
    ax.clabel(cs, inline=1, fontsize=20, colors = 'black', zorder = 3)
# Colora i livelli di temperatura (vmin e vmax sono importanti e devono essere uguali in questa
# funzione e nel Normalize della colorbar sotto)
    ax.contourf(to_np(lons), to_np(lats), to_np(interp_t), 10,
             transform=crs.PlateCarree(),
             cmap=get_cmap("jet"), alpha = 0.6,
               vmin = -35, vmax = -10)

# Mette i limiti agli assi
    ax.set_xlim(cartopy_xlim(interp_gpt))
    ax.set_ylim(cartopy_ylim(interp_gpt))

# Mette le linee di griglia
    ax.gridlines(color="black", linestyle="dotted")

# Titolo con orario
    ax.set_title(i + " ott 2018, 18:00",
                fontsize = 30)

# Questo toglie il sesto plot (dato che sono cinque). Se effettivamente
# togliamo il 26 ottobre, si può togliere questa cosa
fig.delaxes(axs[1,2])

# Queste cose servono per la colorbar orizzontale. Se togliamo il 26
# ottobre c'è da lavorarci di nuovo un po' sopra

# Mappa colori scelta
cmap = get_cmap("jet")

# Normalizzazione della mappa colori tra i due limiti che vogliamo (-35 e -10)
# Devono essere uguali in vmin e vmax del contourf!
cmappable = ScalarMappable(norm=Normalize(-35,-10), cmap=cmap)

# Stupidaggini tecniche
p0 = axs[0,0].get_position().get_points().flatten()
p1 = axs[0,2].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.02])

# Settaggio della colorbar
cbar = fig.colorbar(cmappable, cax=ax_cbar, orientation='horizontal', shrink=.5, alpha = 0.7)
cbar.ax.tick_params(labelsize=20)
cbar.set_label("Temperatura [°C]", fontsize=25)

# Aggiunstamento dei plot
fig.subplots_adjust(top = 0.92, bottom = 0.15, hspace = 0.1, wspace = 0.01, left = 0.01, right = 0.99)

# Salva la figura
fig.savefig('gph+t.pdf', bbox_inches='tight')

###############################################################################

## PRESSIONE A LIVELLO DEL MARE ##

# Settaggio figura
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 17), subplot_kw = {'projection' : cart_proj_d01})

# Lista giorni
days_list = ["26", "27", "28", "29", "30"]

for i, ax in zip(days_list, axs.ravel()):

# Importazione dei file
    namefile = "wrfout_d01_2018-10-"+i+"_180000"
    ncfile = Dataset(namefile)

# Prende la pressione a livello del mare
    slp = getvar(ncfile, "slp")

# Questo applica uno smoothing alla pressione (era già fatto nell'esempio)
# sul sito
    smooth_slp = smooth2d(slp, 3, cenweight=4)

# Latitudine e longitudine
    lats, lons = latlon_coords(slp)

# Set del plot con le coste e gli stati
    ax.add_feature(states, linewidth=.6, edgecolor="black", zorder = 2)
    ax.coastlines('50m', linewidth=1, zorder = 2)

# Disegna le linee di pressione
    cs = ax.contour(to_np(lons), to_np(lats), to_np(smooth_slp), 10, colors="black",
            transform=crs.PlateCarree())
# Mette il valore dei livelli per ogni livello disegnato
    ax.clabel(cs, inline=1, fontsize=20, colors = 'black', zorder = 3)
# Colora i livelli (vmin e vmax sono importanti e devono essere uguali in questa
# funzione e nel Normalize della colorbar sotto)
    ax.contourf(to_np(lons), to_np(lats), to_np(smooth_slp), 10,
             transform=crs.PlateCarree(), alpha = 0.9,
             cmap=get_cmap("Reds"), vmin = 980, vmax = 1040)

# Limiti della mappa
    ax.set_xlim(cartopy_xlim(smooth_slp))
    ax.set_ylim(cartopy_ylim(smooth_slp))

# Griglia
    ax.gridlines(color="black", linestyle="dotted")

# Titolo del plot con giorno e ora
    ax.set_title(i + " ott 2018, 18:00",
                fontsize = 30)

# Elimina il sesto plot
fig.delaxes(axs[1,2])

# Mappa colori scelta
cmap = get_cmap("Reds")

# Settaggio dei valori massimo e minimo della colorbar. Devono essere uguali a vmin
# e vmax in contourf!
cmappable = ScalarMappable(norm=Normalize(980,1040), cmap=cmap)

# Stupidaggiini tecniche
p0 = axs[0,0].get_position().get_points().flatten()
p1 = axs[0,2].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.02])

# Colorbar
cbar = fig.colorbar(cmappable, cax=ax_cbar, orientation='horizontal', shrink=.5, alpha = 0.6)
cbar.ax.tick_params(labelsize=20)
cbar.set_label("Pressione [hPa]", fontsize=25)

# Ridimensionamento dei subplot
fig.subplots_adjust(top = 0.92, bottom = 0.15, hspace = 0.1, wspace = 0.01, left = 0.01, right = 0.99)

# Salva il grafico
fig.savefig('slp.pdf', bbox_inches='tight')

################################################################################

## UMIDITA E VENTO A 700 hPa ##

# Solito settaggio della figura con numero righe, colonne, size e proiezione
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 17), subplot_kw = {'projection' : cart_proj_d01})

# Lista dei giorni da mettere
days_list = ["26", "27", "28", "29", "30"]


for i, ax in zip(days_list, axs.ravel()):

# Importazione del file
    namefile = "wrfout_d01_2018-10-"+i+"_180000"
    ncfile = Dataset(namefile)

# Qui andiamo a prendere umidità relativa, pressione, le componenti orizzontali
# del vento e il geopotenziale
    rh = getvar(ncfile, "rh")
    p = getvar(ncfile, "pressure")
    ua = getvar(ncfile, "ua")
    va = getvar(ncfile, "va")
    gpt = getvar(ncfile, "z")

# Le variabili vengono interpolate a 700 hPa
    rh_700 = interplevel(rh, p, 700)
    ua_700 = interplevel(ua, p, 700)
    va_700 = interplevel(va, p, 700)
    wspd_700 = interplevel(wspd, p, 700)
    gpt_700 = interplevel(gpt, p, 700)

# Latitudine e longitudine
    lats, lons = latlon_coords(rh_700)

# Plotta stati e coste
    ax.add_feature(states, linewidth=.5, edgecolor="black", zorder = 2)
    ax.coastlines('50m', linewidth=1, zorder = 2)

# Creazione delle barbs del vento. Dato che plottarle tutte renderebbe il plot
# molto pesante, ne plotto una ogni 10.
    ax.barbs(to_np(lons[::10,::10]), to_np(lats[::10,::10]),
          to_np(ua_700[::10, ::10]), to_np(va_700[::10, ::10]), color = 'red',
          transform=crs.PlateCarree(), length=6, zorder = 2)

# Linee di contorno per il geopotenziale.
    cs = ax.contour(to_np(lons), to_np(lats), to_np(gpt_700), 10,
             transform=crs.PlateCarree(),
             colors="black", vmin = 0, vmax = 100, zorder = 2)
# Mette il valore dei livelli per ogni livello di geopotenziale
    ax.clabel(cs, inline=1, fontsize=20, colors = 'black', zorder = 3)


# Colora aree per umidità. Mettere vmin e vmax uguali in Normalize della
# cbar!!
    ax.contourf(to_np(lons), to_np(lats), to_np(rh_700), 10,
             transform=crs.PlateCarree(), alpha = 0.7,
             cmap=get_cmap("Blues"), vmin = 0, vmax = 100)

# Xlim e ylim per la mappa
    ax.set_xlim(cartopy_xlim(rh_700))
    ax.set_ylim(cartopy_ylim(th_700))

# Griglia
    ax.gridlines(color="black", linestyle="dotted")

# Titolo
    ax.set_title(i + " ott 2018, 18:00",
                fontsize = 30)

# Elimina il sesto plot
fig.delaxes(axs[1,2])

# Colormap
cmap = get_cmap("Blues")

# Mettere i valori in Normalize uguali a vmin e vmax di contourf!
cmappable = ScalarMappable(norm=Normalize(0,100), cmap=cmap)

# Stupidaggini tecniche
p0 = axs[0,0].get_position().get_points().flatten()
p1 = axs[0,2].get_position().get_points().flatten()
ax_cbar = fig.add_axes([p0[0], 0.1, p1[2]-p0[0], 0.02])

# colorbar
cbar = fig.colorbar(cmappable, cax=ax_cbar, orientation='horizontal', shrink=.5, alpha = 0.6)
cbar.ax.tick_params(labelsize=20)
cbar.set_label("Umidità relativa [%]", fontsize=25)

# Aggiustamento dei plot
fig.subplots_adjust(top = 0.92, bottom = 0.15, hspace = 0.1, wspace = 0.01, left = 0.01, right = 0.99)

# Salva la figura
fig.savefig('rh+wind.pdf', bbox_inches='tight')

###############################################################################
