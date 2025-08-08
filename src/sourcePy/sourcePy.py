#####################################################################
#source.py ##########################################################
#####################################################################

##### TODO
# - Add handling (and make default) for timedelta runtimes in trajRequests and measurementToRequest

import numpy as np
import os
import datetime as dt
import copy
from scipy.ndimage import gaussian_filter


class trajectory:
    """A class that contains all the information needed to perform trajectory source analyses by pairing HYSPLIT trajectory output files with measurement data. Not all attributes are necessary for the objects to function but are included as options to assist the user in their analysis, particularly with subsetting and parsing lists of trajectory objects.

    Here and in much of the literature, trajectory "endpoints" refer to all points along the trajectories and not the single final point along the trajectory (in either direction). The "arrival time" is the latest time in a trajectory along the time axis. In other words, for a backwards trajectory the arrival time is also the starting time and the opposite is true for a forwards trajectory.
    

    Attributes
    ----------
    fName : str
        The full path to the associated HYSPLIT trajectory file.
    dir : str
        The direction of the trajectory, either 'forwards' or 'backwards'.
    datetimes : list (datetime)
        A list of datetime objects corresponding to the endpoint times.
    data : list
        The data from the trajectory file.
    lats : list (float)
        A list of endpoint latitudes.
    lons : list (float)
        A list of endpoint longitudes, all negative.
    hs : list (float)
        A list of endpoint heights AGL in meters.
    ps : list (float)
        A list of endpoint pressures in hPa.
    measurements : list (float)
        A list of measurement objects, intended to be those associated with the trajectory.


    Methods (Public)
    ----------------
    readNewTraj(fName):
        Reads in and populates a trajectory using the data in the file fName.
    """
    def __init__(self,fName=None,direct=None,datetimes=None,lats=None,lons=None,
                 hs=None,ps=None,measurements=None,data=None):
        """Constructor"""
        
        self.fName = fName
        self.data = data
        self.dir = direct
        self.datetimes = datetimes
        
        # Locations
        self.lats = lats
        self.lons = lons
        self.hs = hs
        self.ps = ps
        
        self.measurements = measurements
        
    def readNewTraj(self,fName=None):
        """Reads in and processes the data from a trajectory
        
        Parameters
        ----------
        fName : str (optional)
            The full path to the associated HYSPLIT trajectory file. Only needed if the fName attribute has not already been specified.
        """
        
        if fName:
            self.fName = fName
        
        # Read in the direction (forwards or backwards)
        info = np.loadtxt(self.fName,skiprows=2,max_rows=1,dtype=str)
        self.dir = str(info[1])
        
        # Read in the data
        skips=5
        while skips < 30:
            try:
                self.data = np.loadtxt(self.fName,skiprows=skips)
                skips = 30
            except ValueError:
                skips += 1
        
        # Dates and times 
        self.datetimes = []
        for i in range(self.data.shape[0]):
            
            # Breaks after the year 2050 or before 1950!
            # NEEDS AN OVERHAUL
            if self.data[i,2].astype(int) > 50:
                self.datetimes.append(dt.datetime(self.data[i,2].astype(int)+1900,self.data[i,3].astype(int),
                                                  self.data[i,4].astype(int),self.data[i,5].astype(int),
                                                  self.data[i,6].astype(int)))
            else:
                self.datetimes.append(dt.datetime(self.data[i,2].astype(int)+2000,self.data[i,3].astype(int),
                                                  self.data[i,4].astype(int),self.data[i,5].astype(int),
                                                  self.data[i,6].astype(int)))
        
        # Locations
        self.lats = self.data[:,9]
        self.lons = self.data[:,10]
        self.hs = self.data[:,11]
        self.ps = self.data[:,12]
        
        self.measurements = []
        
        # Fix positive longitudes (all longitudes negative)
        newLons = self.lons
        newLons[newLons>=0] -= 360
        self.lons = newLons
        

class measurement:
    """A class to store measurements of pollutants and everything about those measurements.

    Each instance contains a single sample from a single site.

    Users are encouraged to put their data into the standard measurment csv format. If the data is not in the specified format, users will have to read it into measurement objects on their own.
    
    Measurements sites are often referred to as 'receptor sites' or simply 'receptors'.

    Attributes
    ----------
    lat : float
        The latitude of the receptor.
    lon : float
        The longitude of the receptor, negative only.
    h : float
        The height of the receptor AGL in meters.
    datetime : datetime
        The datetime object reflecting the start of the sampling period.
    duration : timedelta
        The timedelta of the length of the sampling period.
    name : str
        A short description of the sample.
    value : float
        The actual measured value, often a concentration or count rate.
    ID : str
        The ID for the receptor site or measurement platform.

    (Planned) Methods
    -----------------
    processMeasurements(fName,sr=1,dateForm='%m/%d/%Y %H:%M'):
        Reads in measurements from a standard measurement file and returns a list of measurment objects.
    """
    def __init__(self,lat,lon,h,datetime,name,value,ID,endtime=None,duration=1,releaseTime=None,releaseLength=None,releaseID=None):
        self.lat = lat
        self.lon = lon
        self.h = h
        self.datetime = datetime
        self.duration = duration
        
        # Get the length into a timeDelta
        if type(duration) in (int, float):
            self.duration = dt.timedelta(hours=duration)
            
        # If an end time is given, overwrite the given length
        if endtime:
            self.duration = endtime-datetime
            
        self.name = name
        self.value = value
        self.ID = ID
        self.releaseTime=releaseTime
        self.releaseLength=releaseLength
        self.releaseID=releaseID
        

class trajRequest:
    """A class that can be passed to the hysplit module that outlines the basic properties of a desired trajectory.

    May be combined with the plumeRequest class in the future.

    Attributes
    ----------
    lat : float
        The latitude of the trajectory starting location.
    lon : float
        The longitude of the trajectory starting location.
    hs : list (float)
        A list of starting heights, meters AGL.
    start : datetime
        The trajectory start time.
    runtime : int or float
        The run time for the trajectory, can be positive or negative.
    end : datetime
        The trajectory end time.
    ID : string (optional)
        An ID string to help further organize requests, if desired.
    measNum : int or float (optional)
        A number to help further organize requests, if desired.
    """
    def __init__(self,lat,lon,hs,start,runtime=-60,end=None,ID="defaultID",measNum=1):
        self.lat = lat
        self.lon = lon
        self.hs = hs
        self.start = start
        self.runtime=runtime

        # If an endtime is given, overwrite the default or given runtime
        if end:
            self.end = end
            self.runtime = int(-(self.start-self.end).total_seconds()/(60*60))
            
        # If an endtime is not given, define it
        if not end:
            self.end = self.start + dt.timedelta(hours=self.runtime)

        self.ID = ID
        self.measNum=measNum
        

class plumeRequest:
    """A class that can be passed to the hysplit module that outlines the basic properties of a desired concentration plume.

    May be combined with the trajRequest class in the future.

    Attributes
    ----------
    lat : float
        The latitude of the concentration starting location.
    lon : float
        The longitude of the concentration starting location.
    hs : list (float)
        A list of starting heights, meters AGL.
    start : datetime
        The concentration start time.
    runtime : int or float
        The run time for the concentration, can be positive or negative.
    end : datetime
        The concentration end time.
    ID : string (optional)
        An ID string to help further organize requests, if desired.
    measNum : int or float (optional)
        A number to help further organize requests, if desired.
    """
    def __init__(self,lat,lon,hs,start,end=None,ID="defaultID",runtime=-60,measNum=1):
        self.lat = lat
        self.lon = lon
        self.hs = hs
        self.start = start
        
        # If an endtime is given, overwrite the default or given runtime
        if end:
            self.end = end
            runtime = int(-(self.start-self.end).total_seconds()/(60*60))

        # If an endtime is not given, define it
        if not end:
            self.end = self.start + dt.timedelta(hours=runtime)

        self.ID = ID
        self.runtime=runtime
        self.measNum=measNum
        

#class plume:
#    def __init__(self,fName):
#        self.fName = fName
#        
#        # Read in the plume
        
        
# Functions        
##########################################################################

#### General Utility ####

# Get all the trajectory names from a directory path
def get_fNames(path):
    """Return a list of the file names (full paths) in a given directory.

    Parameters
    ----------
    path : str
        A string representing the full path to the directory of interest. Should end in '/'.

    Returns
    -------
    output : list (str)
        A list of strings representing the full paths to the files in the directory of interest.
    """
    
    files = os.listdir(path)
    output = []
    for i in range(len(files)):
        output.append(path+files[i])
    return output


# Process measurements from a standard csv file
#def processMeasurements(fName,sr=1,dateForm='%m/%d/%Y %H:%M'):
#
#    headers = np.loadtxt(fName,dtype=str,max_rows=1)
#    data = np.loadtxt(fName,skiprows=sr,dtype=str)
#
#    starts,durs,lats,lons,vals,IDs,hs,descs = None,None,None,None,None,None,None,None
#
#    print(headers)
#    print(data[:,:5])
#
#    # Figure out what columns are included
#    for i in range(len(headers)):
#        header = headers[i]
#        if header in ("start","START","Start"):
#            starts = data[i]
#        elif header in ("dur","duration"):
#            durs = data[i]
#        elif header in ("lat","lats"):
#            lats = data[i]
#        elif header in ("lon","lons"):
#            lons = data[i]
#        elif header in ("val","value"):
#            vals = data[i]
#        elif header in ("Site ID","ID"):
#            IDs = data[i]
#        elif header in ("height","hght","h"):
#            hs = data[i]
#        elif header in ("description","desc"):
#            descs = data[i]
#
#    # Turn into measurement objects
#    for i in range(len(starts)):
#
#        # datetime objects for the times
#        start = dt.datetime.strptime(starts[i], dateForm)
#        print(starts[i],start)
#        #etc
#
#    return
        
        
def measurementToRequests(meas,hs=None,runtime=-60,when="start",dWhen=dt.timedelta(hours=1)):
    """Converts a measurement object to trajectory requests. Users can specify the runtime, heights, and how often trajectories should be launched.

    Future changes: Will eventually include handling for plumes along with a combination of the trajRequest and plumeRequest classes. From there it may be moved to a class method.

    Parameters
    ----------
    meas : measurement
        A single measurement object to be used to generate the requests.
    hs : list (floats), default: None
        A list of heights to start trajectories at. If none are provided, the height of the measurement is used.
    runtime : float, default: -60
        The runtime for the requests in hours. Can be positive or negative.
    when : int or str, default: 'start'
        The strategy for how often and when trajectories should be launched. Valid inputs are described below.
    dWhen : timedelta, int, or float, default: dt.timedelta(hours=1)
        The timedelta (or number of hours) used if when is 'delta'.

    Returns
    -------
    requests : list
        A list of trajRequest objects.

    Notes
    -----
    when strategies:
        'start': Launch a single trajectory at the start of the measurement period.
        'middle': Launch a single trajectory in thhe exact middle of the measurement period.
        'end': Launch a single trajectory at the end of the measurement period.
        'delta': Starting at the start of the measurement period, launch a trajectory every dWhen. Inclusive of the end of the period if appropriate.
        some int: 
            1: Same as 'middle'.
            2: Same as 'start' and 'end' combined.
            3+: Using the start and end times, divide the remaining times (when-2) evenly within the measurement period.
                For example, a value of three would use the start, middle, and end times.
    """
    whenErr = "Invalid when, use an int > 0 or one of the following: (start,middle,end,delta)"
    
    if hs:
        theseHs = hs
    else:
        theseHs = [meas.h,]
        
    requests = []
        
    # Produce the requests for the proper times
    # Use the start time only
    if when == "start":
        requests.append(trajRequest(meas.lat,meas.lon,hs,start=meas.datetime,ID=meas.ID,runtime=runtime))
    
    # Use the end time only
    elif when == "end":
        end = meas.datetime + meas.duration
        requests.append(trajRequest(meas.lat,meas.lon,hs,start=end,ID=meas.ID,runtime=runtime))
    
    # Starting at the start time, use every dWhen 
    elif when == "delta":
        # Make sure dWhen is a timedelta, assume hours if not and correct
        if type(dWhen) is not dt.timedelta:
            dWhen = dt.timedelta(hours=dWhen)
        
        start = meas.datetime
        end = start + meas.duration
        time = start
        while time <= end:
            requests.append(trajRequest(meas.lat,meas.lon,hs,start=time,ID=meas.ID,runtime=runtime))
            time += dWhen

    elif (type(when) in (int,float)) or when == "middle":
        start = meas.datetime
        end = start + meas.duration
        
        # Evenly distribute the times between the start and end times
        # If when is 1, use the middle time between start and end
        if when == 1 or when == "middle":
            dTime = dt.timedelta(hours=(meas.duration.total_seconds()/60/60)/2)
            requests.append(trajRequest(meas.lat,meas.lon,hs,start=start+dTime,ID=meas.ID,runtime=runtime))
        # If when is 2, use the start and end times
        elif when == 2:
            requests.append(trajRequest(meas.lat,meas.lon,hs,start=start,ID=meas.ID,runtime=runtime))
            requests.append(trajRequest(meas.lat,meas.lon,hs,start=end,ID=meas.ID,runtime=runtime))
        # If when is 3 or more, use the start and end times and distribute evenly between them
        elif when >= 3:
            dTime = dt.timedelta(hours=(meas.duration.total_seconds()/60/60)/(when-1))
            requests.append(trajRequest(meas.lat,meas.lon,hs,start=start,ID=meas.ID,runtime=runtime))
            for i in range(when-2):
                requests.append(trajRequest(meas.lat,meas.lon,hs,start=start+(dTime*(i+1)),ID=meas.ID,runtime=runtime))
            requests.append(trajRequest(meas.lat,meas.lon,hs,start=end,ID=meas.ID,runtime=runtime))
        else:
            print(whenErr)
    else:
        print(whenErr)
            
    return requests


def round_to_hour(dateIn):
    """Given a datetime object, return a datetime object rounded to the nearest hour.

    Parameters
    ----------
    dateIn : datetime
        The datetime object to be rounded.

    Returns
    -------
    dateOut : datetime
        The dateIn object rounded to the nearest hour.
    """
    dt_start_of_hour = dateIn.replace(minute=0, second=0, microsecond=0)
    dt_half_hour = dateIn.replace(minute=30, second=0, microsecond=0)

    if dateIn >= dt_half_hour:
        dateOut = dt_start_of_hour + dt.timedelta(hours=1)
    else:
        dateOut = dt_start_of_hour

    return dateOut


# Create a regular lat/lon grid to do pscf on
# The grid is assumed to be the center points of the bins
def makeGrid(minLon,maxLon,minLat,maxLat,dLon,dLat):
    """Given a range of min/max lats/lons and lat/lon resolution, return 2D meshes of lons and lats. These locations are treated as the centers of gridboxes.

    Parameters
    ----------
    minLon,maxLon,minLat,maxLat : floats
        The boundaries (inclusive) of the grid
    dLon,dLat : floats
        The width of the gridboxes in the lon/lat directions

    Returns
    -------
    gridLons,gridLats : 2D arrays (floats)
        2D arrays of lons and lats for the desired grid.
    """
    
    lons = np.arange(minLon,maxLon,dLon)
    lats = np.arange(minLat,maxLat,dLat)
    gridLons,gridLats = np.meshgrid(lons,lats)
    
    return gridLons,gridLats


# From a list of trajs, get a list of unique receptor locations and IDs
def getUniqueRecs(trajs):
    """Given a list of trajectories return unique measurement receptor sites and their lat/lon locations.

    Assumes that receptor locations are static.

    In the future, should also accept lists of measurement objects.

    Parameters
    ----------
    trajs : list (trajectory)
        A list of trajectories of interest.

    Returns
    -------
    uniqueIDs : list (str)
        A list of ID strings for each unique receptor.
    lons : list (float)
        A list of longitudes for each unique receptor.
    lats : list (float)
        A list of latitudes for each unique receptor.
    """
    
    # Only use trajs with an assigned measurement to take the ID from
    subTrajs = [t for t in trajs if len(t.measurements)]
    
    uniqueIDs = []
    for t in subTrajs:
        uniqueIDs.append(t.measurements[0].ID)
          
    uniqueIDs,lats,lons = np.unique(uniqueIDs),[],[]
    
    for ID in uniqueIDs:
        theseTrajs = [t for t in subTrajs if t.measurements[0].ID in [ID,]]
        if len(theseTrajs):
            lons.append(theseTrajs[0].measurements[0].lon)
            lats.append(theseTrajs[0].measurements[0].lat)
    
    return uniqueIDs,lons,lats


def listVals(trajs,im=0):
    """A simple method to return measurement values from a list of trajectories.

    Parameters
    ----------
    trajs : list (trajectory)
        A list of trajectories, all with attached measurements.
    im : int, default: 0
        The index for the measurements of interest, useful if one trajectory has multiple measurements.

    Returns
    -------
    vals : list (float)
        A list of the requested measurement values, one for each trajectory in trajs in the same order.
    """
    
    vals = []
    for traj in trajs:
        vals.append(traj.measurements[im].value)
    return vals


#def normalizeByReceptor(trajs,im=0,norm="mean"):
#    """
#    """
#    
#    # Collect the unique IDs
#    uIDs = getUniqueRecs(trajs)[0]
#    
#    # Create a copy of the trajectories to work on
#    newTrajs = copy.deepcopy(trajs)
#    
#    # For each ID, get the values, calculate the norm, and divide by the norm
#    for ID in uIDs:
#        subTrajs = [t for t in newTrajs if t.measurements[im].ID in [ID,]]
#        vals = listVals(subTrajs,im=im)
#        
#        if norm == "mean":
#            n = np.mean(vals)
#        if norm == "median":
#            n = np.median(vals)
#        if norm == "max":
#            n = np.amax(vals)
#            
#        if np.amax(vals) != 0 and n != 0:
#            for traj in subTrajs:
#                traj.measurements[0].value /= n
#    
#    return newTrajs
    

#### PSCF ####

# Perform unweighted PSCF
def calcPSCF(trajs,gridLons,gridLats,dLon,dLat,critVals,im=0):
    """Primary calculation method for the (unweighted) potential source contribution function (PSCF).
    
    Description of the process and references below.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    dLon : float
        The width of the gridboxes in longitude
    dLat : float
        The width of the gridboxes in latitude
    critVals : list (float)
        The critical values to be used. Length should be the same as trajs.
    im : int, default: 0
        The index for the measurements of interest, useful if one trajectory has multiple measurements.

    Returns
    -------
    PSCF : array-like (float)
        The PSCF values. Shape is (nx, ny).
    nCounts : array-like (float)
        The counts of trajectory endpoints in each gridbox. Same shape as PSCF.
    mCounts : array-like (float)
        The counts of trajectory endpoints associated with measurement values above the critical value in each gridbox. Same shape as PSCF.

    Notes
    -----
    PSCF at any one gridbox is simply the number of significant trajectory endpoints (those with an associated measurement above the critical value) within a gridbox (m) divided by the total number of endpoints within that gridbox (n).

    PSCF = m/n

    Gridboxes with no end points at all (n=0) are assigned a nan value.

    References
    ----------
    Cheng, M.-D., Hopke, P. K., Barrie, L. A., Rippe, A., Olson, M. H., & Landsberger, S. (1993). Qualitative determination of source regions of aerosol in Canadian high Arctic. Environmental Science & Technology, 27(10), 2063–2071. https://doi.org/10.1021/es00047a011
    """
    
    # Count ns and ms
    nCounts = np.zeros_like(gridLats)
    mCounts = np.zeros_like(gridLats)
    
    # Work through each trajectory
    for i in range(len(trajs)):
    
        # Get the lats, lons, and concentration for this trajectory
        theseLons,theseLats,thisConc,critVal = trajs[i].lons,trajs[i].lats,trajs[i].measurements[im].value,critVals[i]
        
        # Create temporary grids of 0s for n and m
        ns = np.zeros_like(nCounts)
        ms = np.zeros_like(nCounts)
        
        for iLoc in range(len(theseLons)):
            thisLat,thisLon = theseLats[iLoc],theseLons[iLoc]

            # distance differences
            #dists = np.absolute((gridLons-thisLon)/dlon)+np.absolute((gridLats-thisLat)/dlat)
            lonDists = np.absolute(gridLons-thisLon)
            latDists = np.absolute(gridLats-thisLat)
            dists = np.absolute(gridLons-thisLon)+np.absolute(gridLats-thisLat)

            # minimize to find the closest bin
            iBin = np.where(dists==np.amin(dists))

            # Add to n at this bin if the bin is within dlon and dlat
            if np.amin(lonDists) <= dLon and  np.amin(latDists) <= dLat:
                ns[iBin[0],iBin[1]]+=1

                # Is this a significant concentration trajectory? If so add to ms
                if thisConc > critVal:

                    ms[iBin[0],iBin[1]]+=1

            # Add to the total n and m arrays
            nCounts += ns
            mCounts += ms
        
    PSCF = mCounts/nCounts
    
    return PSCF, nCounts, mCounts


def weighting_CL2001(nCounts,divisor = 10):
    """Produce weights for a PSCF field using the n field and the method outlined in Cheng and Lin (2001). Using a divisor of 10 follows the Cheng and Lin (2001) method exactly, although other numbers can be used if desired. The output of this method is a 2D set of weights that can be multiplied against the PSCF field.

    Parameters
    ----------
    nCounts : array-like (float)
        The counts of trajectory endpoints in each gridbox.
    divisor : float, default: 10
        The value used to weight the field by.

    Returns
    -------
    weights : array-like (float)
        The output weights in the same shape as nCounts.
    """
    # Weight n by the divisor
    weights = nCounts/divisor
    weights[weights > 1] = 1
    return weights


def weighting_heavy(nCounts,divisor = 3):
    """Produce weights for a PSCF field using the n field and the a more flexible method. The divisor used in Cheng and Lin (2001) is weighted here by the maximum n value in nCounts. The output of this method is a 2D set of weights that can be multiplied against the PSCF field.

    This method has not been tested in any published study to date, use with caution!

    Parameters
    ----------
    nCounts : array-like (float)
        The counts of trajectory endpoints in each gridbox.
    divisor : float, default: 10
        The value used to weight the field by.

    Returns
    -------
    weights : array-like (float)
        The output weights in the same shape as nCounts.
    """
    weights = nCounts/(np.amax(nCounts)/divisor)
    weights[weights > 1] = 1
    return weights


#### CWT ####

def calcCWT(trajs,gridLons,gridLats,dlon,dlat,im=0,fixLowMeas=False):
    """Primary calculation for concentration-weighted trajectories (CWT), aka the concentration field (CF).

    Description of the process and references below.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    dLon : float
        The width of the gridboxes in longitude.
    dLat : float
        The width of the gridboxes in latitude.
    im : int, default: 0
        The index for the measurements of interest, useful if one trajectory has multiple measurements.
    fixLowMeas : bool, default: False
        Whether or not to implement an adjustment to keep all measurements above 1 by adding .01 and multiplying by 100. This ensures the log of the measurement values are all above 0.

    Returns
    -------
    CWT : array-like (float)
        The CWT values. Shape is (nx, ny). Unit is concentration.
    resTimes : array-like (float)
        The residence times for all trajectories in each gridbox. Same shape as CWT. Unit is trajectory time step in hours.
    conResTimes : array-like (float)
        The residence times wieghted by the log of the associated measurement value in each gridbox. Same shape as CWT. Unit is concentration * trajectory time step in hours.

    Notes
    -----
    For now, assumes that trajectories have a time step of 1 hour. This means the units are wrong, but the analysis itself should be self-consistent.

    CWT is a natural evolution of PSCF, and can be more useful in certain cases. Trajectories are weighted by the log of their measurement value, rather than simply with a 1 or 0 depending on significance as in PSCF. For the same trajectories the CWT residence time field is simply the PSCF n field multiplied by the time between trajectory endpoints in hours.

    CWT = sum((residence times)*ln(concentration)) / sum(residence times)

    Because we use the natural log of the measurements, care should be taken to ensure all the measurements above 0 are also above 1 through either unit transformation or the fixLowMeas option. Setting fixLowMeas to True modifies each measurement value by adding 0.01 and multiplying by 100 before taking the natural log.

    PSCF weights can also be applied to CWT if desired using the residence times as input rather than n.

    References
    ----------
    Fleming, Z. L., Monks, P. S., & Manning, A. J. (2012). Review: Untangling the influence of air-mass history in interpreting observed atmospheric composition. Atmospheric Research, 104–105, 1–39. https://doi.org/10.1016/j.atmosres.2011.09.009
    """
    
    # Track residence times and concentration-weighted residence times
    resTimes = np.zeros_like(gridLats,dtype=float)
    conResTimes = np.zeros_like(gridLats,dtype=float)
    
    # Work through each trajectory
    for i in range(len(trajs)):
    
        # Get the lats, lons, and concentration for this trajectory
        theseLons,theseLats,thisConc = trajs[i].lons,trajs[i].lats,trajs[i].measurements[im].value
        
        # Create temporary grids of 0s for residence time
        rTs = np.zeros_like(gridLats,dtype=float)
        
        # Get the residence times for this trajectory
        # Assume dT to be the same between all points/trajs
        # Assume dT to be unitless
        # Only do this if concentration is > 0
        
        for iLoc in range(len(theseLons)):
            thisLat,thisLon = theseLats[iLoc],theseLons[iLoc]

            # distance differences
            lonDists = np.absolute(gridLons-thisLon)
            latDists = np.absolute(gridLats-thisLat)
            dists = np.absolute(gridLons-thisLon)+np.absolute(gridLats-thisLat)

            # minimize to find the closest bin
            iBin = np.where(dists==np.amin(dists))

            # Add to residence time at this bin if the bin is within dlon and dlat
            if np.amin(lonDists) <= dlon and  np.amin(latDists) <= dlat:
                rTs[iBin[0],iBin[1]]+=1
        
            # Sum the residence times
            resTimes += rTs
                
            if fixLowMeas:
                # Make sure all concs are above 1
                w = np.log((thisConc+.01)*100)
            else:
                w = np.log(thisConc)
            
            # Sum the concentration weighted residence times
            conResTimes += rTs*w
        
    # Filter out bins with low residence times
    CWT = conResTimes/resTimes
        
    return CWT,resTimes,conResTimes


#### SQTBA ####

# Q is the transition probability density function
# tPrime should be in hours
def calcQ(gridLons,gridLats,lonPrime,latPrime,tPrime,dSpeed):
    """Calculates the transition probability denity function, or the Q field, for a single trajectory endpoint.
    
    The Q field takes into account the trajectory uncertainty in time. It creates a normally distributed PDF of location around the given endpoint, with a standard deviation weighted by the amount of time since the start of the launch and the disperion speed.

    Parameters
    ----------
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    lonPrime : float
        The longitude of the endpoint.
    latPrime : float
        The latitude of the endpoint.
    tPrime : float
        The time elapsed (in either direction) since the start of the trajectory in hours.
    dSpeed : float
        The dispersion speed in degrees/hour.

    Returns
    -------
    Q : array-like (float)
        The 2D Q field. Shape is (nx, ny). Unitless.
    """
    # Use one dispersion sigma for x and y
    sigma = tPrime*dSpeed
    
    # Break calculation up into parts
    A = (1/(2*np.pi*sigma*sigma))
    B = ((gridLons-lonPrime)/sigma)**2
    C = ((gridLats-latPrime)/sigma)**2
    D = np.exp(-.5*(B+C))

    Q = A*D
    
    return Q
    

def calcSQTBA(trajs,gridLons,gridLats,dSpeed=0.025,t0=None):
    """Calculates the simplified quantitative transport bias analysis (SQTBA) field for a given set of trajectories.

    Description of the process and references below.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    dSpeed : float, default: 0.025
        The dispersion speed in degrees/hour.
    t0 : float, default: None
        The time scaling factor (usually the length of the longest trajectory in time).

    Returns
    -------
    SQTBA : array-like (float)
        The SQTBA values. Shape is (nx, ny). Unit is concentration.

    References
    ----------
    Zhou, C., Zhou, H., Holsen, T. M., Hopke, P. K., Edgerton, E. S., & Schwab, J. J. (2019). Ambient Ammonia Concentrations Across New York State. Atmospheres, 124(14), 8287–8302. https://doi.org/10.1029/2019jd030380
    """
    
    # For now, assume the lat/lon grid can be used as the x/y grid
    # Assume 1 degree = 100km
    # Thus, 5.4 km/hr = 0.054 deg/hr
    tBarSum = np.zeros_like(gridLons,dtype=float)
    tTild = np.zeros_like(gridLons,dtype=float)
    
    # Integrate over all trajs
    for iTraj in range(len(trajs)):
        #print(str(iTraj+1)+"/"+str(len(trajs)))
        traj = trajs[iTraj]
        Q = np.zeros_like(gridLons,dtype=float)
        
        # Calc the 2D probability of natural transport field (Q) for each endpoint, and sum them:
        for i in range(len(traj.lons)-1):
            dt = np.absolute((traj.datetimes[0]-traj.datetimes[i+1]).total_seconds())/60/60
            thisQ = calcQ(gridLons,gridLats,traj.lons[i+1],traj.lats[i+1],dt,dSpeed)
            Q += thisQ
        
        # Calculate the potential transport field, tbar:
        tBar = Q/len(traj.lons)
        
        # Add tBar to the sum:
        tBarSum += tBar
        
        # Calculate the weighted potential transfer field, ttild:
        tTild += tBar*(traj.measurements[0].value)
    
    # if no t0 is provided, use the time length of the last trajectory
    if not t0:
        t0 = len(traj.lons)
    weight = (10*len(traj.lons))/(2*np.pi*(dSpeed*t0)**2)

    SQTBA = (tTild/tBarSum)*(1-np.exp(-tBarSum/weight))
    
    return SQTBA
    

#### RCF ####

def calcCI(trajs,gridLons,gridLats,dLon,dLat,confidence=95,nBoot=1000):
    """Given a set of trajectories, use bootstrapping to calculate a confidence interval of CWT values for each gridbox. Part of the RCF calculation.

    Bootstrapping is done by resampling the trajectories "nBoot" times and recording the CWT values at each gridbox. From there, percentiles are taken "confidence" apart, and these are returned as the upper and lower bounds.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    dLon : float
        The width of the gridboxes in longitude.
    dLat : float
        The width of the gridboxes in latitude.
    confidence : int, default: 95
        The % confidence to calculate the interval for.
    nBoot : int, default: 1000
        The number of times to bootstrap.

    Returns
    -------
    upper : array-like (float)
        The upper CWT values in the interval. Shape is (nx, ny). Unit is concentration.
    lower : array-like (float)
        The lower CWT values. Shape is (nx, ny). Unit is concentration.
    """
    
    # To make this work with bootstrapping, we need to rethink how we do CWT
    # Pre-calculate the residence times and the  for each trajectory individually
    indRes,indConRes = indResTimes(gridLons,gridLats,trajs,dLon,dLat)
    
    # Then, we can re-sample the trajectories quickly
    # Calculate nBoot CWT from random samples of the trajectories
    CWTs = np.zeros(shape=(nBoot,indRes.shape[1],indRes.shape[2]))
    #print(CWTs.shape)
    for i in range(nBoot):
        sample = np.random.choice(range(indRes.shape[0]),size=indRes.shape[0])
        sampleRes = indRes[sample,]
        sampleConRes = indConRes[sample,]
        sampleCWT = np.sum(sampleConRes,axis=0)/np.sum(sampleRes,axis=0)
        CWTs[i] = sampleCWT
        
    dPerc = (100-confidence)/2
    upper = np.percentile(CWTs,100-dPerc,axis=0)
    lower = np.percentile(CWTs,dPerc,axis=0) 
    
    return upper,lower


# Get the residence and conc-weighted residence times for each individual trajectory
def indResTimes(trajs,gridLons,gridLats,dLon,dLat,im=0):
    """Lazily calculates the residence times for each trajectory individually. Part of the RCF calculation.

    Instead of only calculating the residence times, use the calcCWT function for each trajectory and only keep the residence times.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    dLon : float
        The width of the gridboxes in longitude.
    dLat : float
        The width of the gridboxes in latitude.
    im : int, default: 0
        The index for the measurements of interest, useful if one trajectory has multiple measurements.

    Returns
    -------
    resTimes : array-like (float)
        The residence times for all trajectories in each gridbox. Shape is (nTrajs, nx, ny). Unit is trajectory time step in hours.
    conResTimes : array-like (float)
        The residence times wieghted by the log of the associated measurement value in each gridbox. Shape is (nTrajs, nx, ny). Unit is concentration * trajectory time step in hours.
    """
    
    # Make arrays of size (nTrajs,nx,ny)
    resTimes = np.zeros((len(trajs),gridLons.shape[0],gridLons.shape[1]))
    conResTimes = np.zeros_like(resTimes,dtype=float)
    
    for iTraj in range(len(trajs)):
        theseTrajs = [trajs[iTraj]]
        
        CWT,res,conRes = calcCWT(theseTrajs,gridLons,gridLats,dLon,dLat,fixLowMeas=True)
        
        resTimes[iTraj] = res
        conResTimes[iTraj] = conRes
    
    return resTimes,conResTimes


def smoothCWT(CWT,upper,lower,sigma=.5):
    """Given a CWT field and its upper and lower confidence interval bounds, use a gaussian filter to smooth the CWT field. Keep values between the upper and lower bounds. Part of the RCF calculation.

    Parameters
    ----------
    CWT : array-like (float)
        The CWT field to be smoothed.
    upper : array-like (float)
        The upper bounds of the CWT confidence interval.
    lower : array-like (float)
        The lower bounds of the CWT confidence interval.

    Returns
    -------
    smoothed : array-like (float)
        The smoothed CWT field.
    tooHigh : array-like (int)
        The index of where the smoothed field is too high.
    tooLow : array-like (int)
        The index of where the smoothed field is too low.
    corrected : array-like (float)
        The smoothed and bounded CWT field.
    """
    
    smoothed = gaussian_filter(CWT,sigma=sigma)
    uDiff = upper-smoothed # Negative values are out of bounds
    lDiff = smoothed-lower # Same here
    
    tooHigh,tooLow = np.where(uDiff<-.0001),np.where(lDiff<-.0001)
    
    corrected = copy.deepcopy(smoothed)
    corrected[tooHigh] += uDiff[tooHigh]
    corrected[tooLow] -= lDiff[tooLow]
    corrected[np.isnan(corrected)]=0
    
    return smoothed,tooHigh,tooLow,corrected


def redistribute(trajs,thisCWT,gridLons,gridLats,dlon,dlat,im=0):
    """Redistribute the measured values along each trajectory. Part of the RCF calculation.

    Future work: add the ability to control segment length.

    A full explanation is below.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    thisCWT : array-like (float)
        The original 2D CWT field, usually smoothed/bounded.
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    dLon : float
        The width of the gridboxes in longitude.
    dLat : float
        The width of the gridboxes in latitude.
    im : int, default: 0
        The index for the measurements of interest, useful if one trajectory has multiple measurements.

    Returns
    -------
    allSegments : list (trajectory)
        The full list of segments from all trajs, with measurements reflecting the redistributed CWT values.

    Notes
    -----
    This redistribution step re-weights measurements (and thus the CWT values) along trajectories based on the ratio of the average CWT along a segment and the average CWT along the trajectory that segment came from. For now each segment is really just an endpoint, so the average CWT along the segment is simply the CWT of the gridbox the endpoint is in.

    Each trajectory is broken down into a list of trajectories with a single endpoint which serve as the segments, each with a unique measurement object. Once the average CWT along the trajectory is found, each segment has its measurement value redistributed in the following way:

    (segment meas value) = (trajectory meas value) * (mean segment CWT) / (mean trajectory CWT)

    As a result, trajectories with consistent CWT values along the path remain relatively unchanged, and trajectories where the CWT varies will have extra weight added to high-CWT regions and weight removed from low-CWT regions.
    """
    
    # Break down each trajectory into a single trajectory for each point
    # These will be the segments for now, makes the math easy
    # Segments only need to know their lat, lon, and measurement
    
    allSegments = []
    for traj in trajs:
        segments = []
        Xis = []
        Xk = 0
        for i in range(len(traj.lats)):
            # Each segment gets a new measurement object, because we will be modifying them
            t = trajectory(lats=[traj.lats[i],],lons=[traj.lons[i],],
                                  measurements=[copy.deepcopy(traj.measurements[im]),])
            
            # Now we need the CWT values for each segment
            # Find the correct bin
            lonDists = np.absolute(gridLons-t.lons[0])
            latDists = np.absolute(gridLats-t.lats[0])
            dists = lonDists+latDists

            # minimize to find the closest bin
            iBin = np.where(dists==np.amin(dists))

            # Add to Xk at this bin if the bin is within dlon and dlat
            if np.amin(lonDists) <= dlon and  np.amin(latDists) <= dlat:
                
                # Keep this segment
                segments.append(t)
                
                # Calculate Xi
                Xi = np.exp(thisCWT[iBin[0],iBin[1]])/len(traj.lats)
                
                if len(Xi)>1:
                    Xi = Xi[0]
                
                Xis.append(Xi)
                Xk+=Xi
                
        # Now that we have a final value for Xk, we can find the new concentrations for the segments
        for iSeg in range(len(segments)):
            
            newConc = traj.measurements[im].value*(Xis[iSeg]/Xk)
            segments[iSeg].measurements[im].value = newConc
            
            allSegments.append(segments[iSeg])    
    
    return allSegments


# Perform Redistribution Concentration Field analysis (RCF)
def calcRCF(trajs,gridLons,gridLats,dLon,dLat,im=0,sigma=0.5,nBoot=1000,con1=95,con2=99,maxIters=10):
    """Calculate the redistribution concentration field analysis as in Stohl 1996.

    Full documentation is below.

    Parameters
    ----------
    trajs : list (trajectory)
        The list of trajectories to be used.
    gridLons : array-like (float)
        The 2D grid of longitudes.
    gridLats : array-like (float)
        The 2D grid of latitudes.
    dLon : float
        The width of the gridboxes in longitude.
    dLat : float
        The width of the gridboxes in latitude.
    im : int, default: 0
        The index for the measurements of interest, useful if one trajectory has multiple measurements.
    sigma : float, default: 0.5
        The sigma value used for the gaussian filter.
    nBoot : int, default: 1000
        The number of times to bootstrap.
    con1 : int, default: 95
        The % confidence to calculate the interval for the first time.
    con2 : int, default: 99
        The % confidence to calculate the interval for after the first time.
    maxIters : int, default: 10
        The maximum number of iterations.

    Returns
    -------
    corrected : array-like (float)
        The smoothed, bounded, and redistributed CWT field after some number of iterations.
        
    Notes
    -----
    This method attempts to recreate the RCF procedure as outlined in the Stohl 1996 paper. A few liberties had to be taken with the smoothing procedure, as nowhere in the literature (that I could find) is it laid out in detail.

    The first step is to calculate our initial CWT field, as normal. Then, we obtain a 95% confidence interval for this intial field, and smooth using a gaussian filter, keeping the field within the confidence interval bounds.

    The second step involves iteratively redistributing and then smoothing the CWT field (redistribution is described in more detail in the eponymous method) until the % change from the previous step to the current one is less than 0.5%, indicating a solution has been converged to.

    Due to the nature of RCF (and some code inefficiency), this can take a long time to calculate. Additionally, some fields may not converge and will continue to iterate indefinetely, so a hard limit (default 10) on the number of iterations is provided.

    References
    ----------
    Stohl, A. (1996). Trajectory statistics—A new method to establish source–receptor relationships of air pollutants and its application to the transport of particulate sulfate in Europe. Atmospheric Environment, 30(4), 579–587. https://doi.org/10.1016/1352-2310(95)00314-2
    """
    
    # Calculate the first guess CWT
    print("First Guess CWT...")
    CWT,resTimes,conResTimes = calcCWT(trajs,gridLons,gridLats,dLon,dLat,fixLowMeas=True)
    
    # Calculate CI
    upper,lower = calcCI(gridLons,gridLats,trajs,dLon,dLat,confidence=con1,nBoot=nBoot)
    
    # Smooth the original field within the confidence interval
    smooth,tooHigh,tooLow,corrected = smoothCWT(CWT,upper,lower,sigma=sigma)
    
    oldCorrected = corrected
    
    # While the average difference between the previous CWT field and the next is >0.5%:
    i = 0
    diffAvg = 100
    while diffAvg > 0.5 and i <= 5:
        print("Redistribution",i+1)
        segs = redistribute(trajs,oldCorrected,gridLons,gridLats,dLon,dLat)
        RCF,newResTimes,newConResTimes = calcCWT(segs,gridLons,gridLats,dLon,dLat,fixLowMeas=True)
        print("CI/Smoothing...")
        upper,lower = calcCI(gridLons,gridLats,segs,dLon,dLat,confidence=con2,nBoot=nBoot)
        smooth,tooHigh,tooLow,corrected = smoothCWT(RCF,upper,lower,sigma=sigma)
        
        diff = corrected-oldCorrected
        
        diffAvg = np.abs(diff)
        diffAvg = (np.abs(diff)/oldCorrected)*100
        diffAvg[np.isnan(diffAvg)] = 0
        diffAvg[np.isinf(diffAvg)] = 0
        diffAvg = np.mean(diffAvg)
        print(diffAvg)
        i+=1
        
        oldCorrected = corrected
    
    return corrected




################ OLD CAPTEX SPECIFIC STUFF ###################

# Read in the ground measurement data from a file
def groundMeasurements(fName,sr=2,dtype=float):
    
    return np.loadtxt(fName,skiprows=sr,dtype=dtype)


# Process the ground measurement data for the second captex tracer release
# Output is an array of measurement objects
def processCAPTEXgnd(fName,sr=2):
    
    data = groundMeasurements(fName,sr=sr)
    dataStr = groundMeasurements(fName,sr=sr,dtype=str)
    
    measurements = []
    for i in range(data.shape[0]):
        
        # Get the correct datetime object
        hour = int(dataStr[i,3][:2])
        minute = int(dataStr[i,3][3:])
        datetime = dt.datetime(data[i,0].astype(int),data[i,1].astype(int),data[i,2].astype(int),hour,minute)
        
        # Get the correct duration (assuming no minutes)
        duration = int(dataStr[i,4][:2])
        
        # Get the release time
        
        
        # Put it in a measurement object
        measurements.append(measurement(data[i,5],data[i,6],10,datetime,"C7F14 (PMCH) concentration (pg/m3)",
                                        data[i,7],str(data[i,8].astype(int)),duration=duration))
        
    return measurements


def getCAPTEXReleases():
    # Return a list of datetimes for the 7 releases we have
    return [dt.datetime(1983,9,18,17,0),
            dt.datetime(1983,9,25,17,0),
            dt.datetime(1983,10,2,19,0),
            dt.datetime(1983,10,14,16,0),
            dt.datetime(1983,10,26,4,0),
            dt.datetime(1983,10,28,15,30),
            dt.datetime(1983,10,29,6,0)]



#####################################################################
#plotSource.py ######################################################
#####################################################################

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


#####################################################################
#hysplit.py #########################################################
#####################################################################

import os
import numpy as np
import datetime as dt
import subprocess as sp
from sourcePy import source


class met:
    """A class that contains information about a meteorology file in the NOAA ARL format.
    
    Currently not in use, but will be integrated into the sourcePy framework soon.

    Attributes
    ----------
    fName : str
        The full path to the ARL file.
    start : datetime
        The start datetime of the file.
    end : datetime
        The end datetime of the file.
    name : str
        
    """
    def __init__(self,fName,start=None,end=None,name=None):
        self.fName = fName
        self.start = start
        self.end = end
        self.name = name
        

def metDates(metFile,exec_path):
    """Gets all of the datetimes in a given ARL met file.

    Uses the HYSPLIT program filedates.exe

    Attributes
    ----------
    metFile : str
        The full path to the ARL file.
    exec_path : str
        The full path to the HYSPLIT program filedates.exe, usually .../HYSPLIT/exec/filedates.exe.

    Returns
    -------
    dates : list (datetimes)
        A list of the datetimes in the ARL met file.
    """

    # Ensure the metFile path is in the correct format (single / at end)
    # Get the dates, piping the output back to here rather than in the shell
    out = sp.run(exec_path+" "+metFile,shell=True,capture_output=True,text=True).stdout.splitlines()

    # Convert the output strings to datetime objects
    dates = []
    if out:
        for date in out:
            try:
                dates.append(dt.datetime.strptime(date,"%y %m %d %H %M"))
            except:
                None

    return dates

    
# Given a directory path of met (ARL) files, return a list of file names and start/end date pairings
def metRange(met_path,exec_path):
    """Process the met files in a directory and return the names and start/end times of each file.

    Parameters
    ----------
    met_path : str
        The path to the directory containing the met files. Should end with a '/'.
    exec_path : str
        The full path to the HYSPLIT program filedates.exe, usually .../HYSPLIT/exec/filedates.exe.

    Returns
    -------
    mets : list
        A list of met file information, sorted by start date. Each entry contains a tuple with the following structure:
        (full file path (str), start (datetime), end (datetime), time step (timedelta)).
    """
    
    files = source.get_fNames(met_path)
    out = []
    for fName in files:
        allDates = metDates(fName,exec_path)
        dTime = dt.timedelta(hours=0)
        if len(allDates):
            if len(allDates) > 1:
                dTime = (allDates[-1]-allDates[0])/(len(allDates)-1)
            dates = [allDates[0],allDates[-1]]
            out.append((fName,allDates[0],allDates[-1],dTime)) # File path, start time, end time, timedelta

    # Sort the list of mets by start date
    mets = np.asarray(sorted(out,key=lambda x: (x[1])))
    
    return mets


# Given a trajectory/plume request and a directory path of met files, return the necessary met files
# If the start or end time is between files, use both.
def necessaryMets(req,mets,exec_path):
    """Given a trajectory or plume request and a directory of met files, return the met files that HYSPLIT would need to complete those requests.

    Assumes that all met files have the same time step and have no missing times.

    Loops through times from the start of the request to the end by half the met file time step. If the start or end of a met file is within a time step of the current time, include it as necessary.

    Parameters
    ----------
    req : request
        The trajectory or plume request object.
    mets : list
        A list of met file information, sorted by start date. Each entry contains a tuple with the following structure:
        (full file path (str), start (datetime), end (datetime), time step (timedelta)).
    exec_path : str
        The full path to the HYSPLIT program filedates.exe, usually .../HYSPLIT/exec/filedates.exe.

    Returns
    -------
    necessaryFiles : list (str)
        A list of strings containing the full paths to the met files necessary for HYSPLIT to run the requested operation.
    """

    # Determine the full range of coverage needed
    dTime = mets[0][-1]
    
    time = np.amin((req.start,req.end))
    needEnd = np.amax((req.start,req.end))

    
    # For each met file, determine if it is necessary
    necessaryFiles = []
    while time <= needEnd:
        for met in mets:
            # If time is within the file, include it
            if time >= met[1]-(dTime) and time <= met[2]+(dTime):
                necessaryFiles.append(met[0])
        time += dTime/2

    necessaryFiles = np.unique(necessaryFiles)
        
    return necessaryFiles


# Given a request, the necessary met files, a storage path, and a working path, write a control file
def generateTrajCONTROL(req,metFiles,storage_path,working_path,vMotion=0,modelTop=10000.0):
    """Given a request and some further information about the HYSPLIT run, create and write a trajectory CONTROL file in the HYSPLIT working directory.

    Parameters
    ----------
    req : request
        The trajectory request object.
    metFiles : list (str)
        A list of the full paths to the met files necessary for the trajectory.
    storage_path : str
        The full path to the directory where HYSPLIT output files should be written.
    working_path : str
        The full path to the HYSPLIT working directory.
    vMotion : int, default: 0
        The vertical motion setting for HYSPLIT to use. 0 uses the vertical wind from the ARL file.
    modelTop : float, default: 10000.0
        The top of the model in meters AGL.
    """

    # Create the control file
    ctrl = req.start.strftime("%y %m %d %H %M\n")
    ctrl += "1\n"
    ctrl += str(req.lat)+" "+str(req.lon)+" "+str(req.hs[0]) + "\n"
    ctrl += str(req.runtime) + "\n"
    ctrl += str(vMotion) + "\n"
    ctrl += str(modelTop) + "\n"
    ctrl += str(len(metFiles)) + "\n"
    for met in metFiles:
        splits = met.split('/')
        for s in splits[:-1]:
            ctrl += s+"/"
        ctrl += "\n"+splits[-1]+"\n"
    ctrl += storage_path+"\n"
    ctrl += req.ID+"_"+req.start.strftime("%Y%m%d_%H%M_")+str(req.hs[0])+".traj"
    #print(ctrl)

    with open(working_path+"CONTROL", "w") as file:
        file.write(ctrl)
    return


# Given a list of string inputs, write those inputs to the SETUP.CFG file
# Given nothing, write nothing.
# Input should be a list of strings like ["TOUT=60","NVER=0"]
# https://www.ready.noaa.gov/hysplitusersguide/S410.html
def generateSETUP(working_path,setup=[]):
    """Generate a HYSPLIT SETUP file.

    Parameters
    ----------
    working_path : str
        The full path to the HYSPLIT working directory.
    setup : list (str)
        The setup parameters to use. Something like ["TOUT=60","NVER=0"]. A full list of possible parameters can be found here: https://www.ready.noaa.gov/hysplitusersguide/S410.html
    """

    out = "&SETUP\n"
    for setting in setup:
        out+=setting+"\n"
    out += "/\n"

    with open(working_path+"SETUP.CFG", "w") as file:
        file.write(out)
    return
    

# Given a working directory and a path to the hysplit trajectory executable, run hysplit
def runHYSPLITStandard(working_path,exec_path,out = False):
    """Once the CONTROL and SETUP files have been generated, this method runs HYSPLIT. It can be pointed towards any executable, but is generall going to be pointed towards hyts_std.exe for trajectories and hycs_std.exe for plumes.

    Parameters
    ----------
    working_path : str
        The full path to the HYSPLIT working directory.
    exec_path : str
        The full path to the HYSPLIT program, usually hyts_std.exe or hycs_std.exe.
    out : bool
        Whether to return the text from the output of the shell.

    Returns
    -------
    outText : str
        The output text from the shell, if requested.
    """

    # Move to the working directory
    os.chdir(working_path)
    #print(sp.run("dir",shell=True,capture_output=True,text=True).stdout)

    # Run HYSPLIT
    outText = sp.run(exec_path,shell=True,capture_output=True,text=True).stdout
    #print(sp.run("dir",shell=True,capture_output=True,text=True).stdout)

    if out:
        return outText
    return


# Given a request, the necessary met files, a storage path, and a working path, write a control file
def generateConcCONTROL(req,mets,storage_path,working_path,vMotion=0,modelTop=10000.0,polID="TEST",eRate=1.0,eHours=1.0,center=(0.0,0.0),
                        releaseStart=None,samplingStart=None,samplingStop=None):

    # Create the control file
    ctrl = req.start.strftime("%y %m %d %H %M\n")
    ctrl += "1\n"
    ctrl += str(req.lat)+" "+str(req.lon)+" "+str(req.hs[0]) + "\n"
    ctrl += str(req.runtime) + "\n"
    ctrl += str(vMotion) + "\n"
    ctrl += str(modelTop) + "\n"
    ctrl += str(len(mets)) + "\n"
    for met in mets:
        splits = met.split('/')
        for s in splits[:-1]:
            ctrl += s+"/"
        ctrl += "\n"+splits[-1]+"\n"

    # From here, the conc control file diverges from the traj control file
    ctrl += "1\n" # number of pollutants
    # Repeated for each pollutant
    ctrl += str(polID)+"\n"
    ctrl += str(eRate)+"\n"
    ctrl += str(eHours)+"\n"
    if releaseStart:
        ctrl += releaseStart.strftime("%y %m %d %H %M\n")
    else:
        ctrl += req.start.strftime("%y %m %d %H %M\n")

    
    ctrl += "1\n" # number of simultanious concentration grids
    # Repeated for each concentration grid
    ctrl += "35.953131 -84.300788\n" # lat/lon location of grid center
    ctrl += "0.001 0.001\n" # lat/lon grid spacing (degrees)
    ctrl += "1.0 1.0\n" # lat/lon grid span (degrees)
    ctrl += storage_path+"\n"
    ctrl += req.ID+"_"+req.start.strftime("%Y%m%d_%H%M_")+str(req.hs[0])+".conc\n"
    ctrl += "1\n" # Vertical concentration levels
    ctrl += "100\n" # level heights (m)
    if samplingStart:
        ctrl += samplingStart.strftime("%y %m %d %H %M\n")
    else:
        ctrl += req.start.strftime("%y %m %d %H %M\n")
    if samplingStop:
        ctrl += samplingStop.strftime("%y %m %d %H %M\n")
    else:
        ctrl += "00 00 00 00 00\n"
    ctrl += "00 1 00\n" # sampling interval: type, hour, minute

    ctrl += "1\n" # number of pollutants depositing (same as line 10, number of pollutants)
    ctrl += "0.0 0.0 0.0\n" # particle diameter, density, and shape
    ctrl += "0.0 0.0 0.0 0.0 0.0\n" # Other particle stuff
    ctrl += "0.0 0.0 0.0\n" # Even more particle stuff
    ctrl += "0.0\n" # half-life
    ctrl += "0.0\n" # resuspension
    
    

    with open(working_path+"CONTROL", "w") as file:
        file.write(ctrl)
    return