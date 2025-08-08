import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sourcePy import source
import os
import datetime as dt


# Set up the geographic plot
def basicGeoPlot(extent=None,proj=ccrs.PlateCarree(),s=(5,5),states=True,borders=True,coastline=True):
    """Create a basic geographic plot using matplotlib. The precision of the plotted borders adjust automatically to the size and extent of the plot.

    Parameters
    ----------
    extent : list (float), default: None
        The extent of the figure (minLon,maxLon,minLat,maxLat).
    proj : projection, default: ccrs.PlateCarree()
        The figure projection.
    s : list (float), default: (5,5)
        The size of the figure (x,y).
    states : bool, default: True
        Whether or not to include US state borders in the figure.
    borders : bool, default: True
        Whether or not to include country borders in the figure.
    coastline : bool, default: True
        Whether or not to include coastlines in the figure.

    Returns
    -------
    fig
        The requested matplotlib figure.
    """
    
    fig = plt.figure(figsize=s)
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    if states:
        ax.add_feature(cfeature.STATES, edgecolor='k', lw=.5)
    if borders:
        ax.add_feature(cfeature.BORDERS, edgecolor='k', lw= .5)
    if coastline:
        ax.add_feature(cfeature.COASTLINE, edgecolor='k', lw= .5)
    
    return fig


# Create a plot of trajectories or plot them on an existing figure
def plotTrajs(trajs,fProj=ccrs.LambertConformal(),dProj=ccrs.PlateCarree(),startPoints=False,endPoints=False,
             alpha=0.2,c=None,newFig=False,extent=None,s=(5,5)):
    """Plot trajectory paths from a list of trajectories, either on a new geographic plot or the one currently being worked on.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    fProj : projection, default: ccrs.LambertConformal()
        The figure projection.
    dProj : projection, default: ccrs.PlateCarree()
        The data projection, ccrs.PlateCarree() if data is on a regular lat/lon grid.
    startPoints : bool, default: False
        Whether or not to plot the start points of the trajectories.
    endPoints : bool, default: False
        Whether or not to plot the end points of the trajectories.
    alpha : float, default: 0.2
        The alpha for the trajectory paths.
    c : string, default: None
        The color of the trajectories, if desired. Otherwise each trajectory will be assigned a color automatically by matplotlib.
    newFig : bool, default: False
        Whether to create a new figure.
    extent : list (float), default: None
        The extent of the figure (minLon,maxLon,minLat,maxLat), if new.
    s : list (float), default: (5,5)
        The size of the figure (x,y), if new.    

    Returns
    -------
    fig
        The requested matplotlib figure, if new.
    """
    
    # Add handling here for forward trajs eventually?
    
    # Create a figure if desired
    if newFig:
        fig = basicGeoPlot(proj=fProj,extent=extent,s=s)
    
    for tr in trajs:
        if c:
            plt.plot(tr.lons,tr.lats,transform=dProj,alpha=alpha,zorder=0,c=c)
        else:
            plt.plot(tr.lons,tr.lats,transform=dProj,alpha=alpha,zorder=0)
            
        if startPoints:
            plt.scatter(tr.lons[0],tr.lats[0],c='lightgray',transform=dProj,s=8)
        if endPoints:
            plt.scatter(tr.lons[-1],tr.lats[-1],c='k',transform=dProj,s=5,marker='x')
    
    if newFig:
        return fig
    else:
        return
    
    
def plotUniqueRecs(trajs,fProj=ccrs.LambertConformal(),dProj=ccrs.PlateCarree(),
                   size=None,labels=False,newFig=False,extent=None,fSize=None,alpha=1):
    """Plot the unique receptor locations, either on a new geographic plot or the one currently being worked on.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    fProj : projection, default: ccrs.LambertConformal()
        The figure projection.
    dProj : projection, default: ccrs.PlateCarree()
        The data projection, ccrs.PlateCarree() if data is on a regular lat/lon grid.
    size : float, default: None
        The size of the markers for the receptor locations.
    labels : bool, default: False
        Whether to plot the IDs of the receptor sites.
    newFig : bool, default: False
        Whether to create a new figure.
    extent : list (float), default: None
        The extent of the figure (minLon,maxLon,minLat,maxLat), if new.
    fSize : list (float), default: (5,5)
        The size of the figure (x,y), if new.
    alpha : float, default: 1.0
        The alpha for the receptor markers.

    Returns
    -------
    fig
        The requested matplotlib figure, if new.
    """
    
    IDs,lons,lats = source.getUniqueRecs(trajs)
    
    if newFig:
        fig = basicGeoPlot(proj=fProj,extent=extent,s=fSize)
    
    plt.scatter(lons,lats,c='silver',transform=dProj,s=size,alpha=alpha,edgecolors='k',zorder=2)
    if labels:
        for i in range(len(lons)):
            plt.text(lons[i],lats[i],IDs[i],c='k',transform=dProj)
    
    if newFig:
        return fig
    else:
        return

    
    
def sourcePlot(gridLons,gridLats,data,extent=None,cmap="gist_heat_r",fProj=ccrs.LambertConformal(),dProj=ccrs.PlateCarree(),
               trajs=[],pTrajs=False,pRecs=False,recLabels=False,vmin=0,vmax=None,alpha=0.5,newFig=True,
               fSize=(5,5),title=None,eColor='lightgray',eWidth=0.05,norm=None,extend=None):
    """Plot a 2D field using pcolormesh, either on a new geographic plot or the one currently being worked on. Intended to reduce most basic plotting procedures to a single method. 

    Parameters
    ----------
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    data : array-like (float)
        The 2D field to plot.
    extent : list (float), default: None
        The extent of the figure (minLon,maxLon,minLat,maxLat), if new.
    cmap : str, default: "gist_heat_r"
        The colormap to be used for the 2D data.
    fProj : projection, default: ccrs.LambertConformal()
        The figure projection.
    dProj : projection, default: ccrs.PlateCarree()
        The data projection, ccrs.PlateCarree() if data is on a regular lat/lon grid.
    trajs : list (trajectory), default: []
        The list of trajectories to be used.
    pTrajs : bool, default: False
        Whether to plot the trajectories.
    pRecs : bool, default: False
        Whether to plot the unique receptor sites.
    recLabels : bool, default: False
        Whether to plot the IDs of the receptor sites.
    vmin : float, default: 0
        The minimum value to plot.
    vmax : float, default: None
        The maximum value to plot.
    alpha : float, default: 0.5
        The alpha for the 2D data.
    newFig : bool, default: False
        Whether to create a new figure.
    fSize : list (float), default: (5,5)
        The size of the figure (x,y), if new.
    title : str, default: None
        The title of the figure, if desired.
    eColor : str, default: 'lightgray'
        The edgecolor for the pcolormesh plot.
    eWidth : float, default: 0.05
        The width of the edges in the pcolormesh plot.
    norm : normalization object, default: None
        The normalization to be passed to pcolormesh(), if desired.
    extend : str, default: None
        Passed to the extend parameter in colorbar().

    Returns
    -------
    fig
        The requested matplotlib figure, if new.
    """
    
    if not len(trajs) and pRecs:
        print("Trajectories required for receptor and/or trajectory plotting.")
        pTrajs=False
        pRecs=False
    
    # Create a figure if desired
    if newFig:
        fig = basicGeoPlot(proj=fProj,extent=extent,s=fSize)
    
    # Plot the provided data
    plt.pcolormesh(gridLons,gridLats,data,cmap=cmap,vmin=vmin,vmax=vmax,alpha=alpha,
                   edgecolor=eColor,lw=eWidth,transform=dProj,norm=norm)
    plt.colorbar(extend=extend)
    
    # Plot receptors
    if pRecs:
        plotUniqueRecs(trajs,dProj=ccrs.PlateCarree(),labels=recLabels,size=8)
    
    # Plot trajs
    if pTrajs:
        plotTrajs(trajs,fProj=fProj,dProj=dProj,startPoints=False,endPoints=False)
          
    if title:
        plt.title(title)
        
    if newFig:
        return fig
    else:
        return


        
    


