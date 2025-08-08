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
    outName : str, default: None
        The name of the output file. If none, one is generated automatically.
    """

    # Create the control file
    ctrl = req.start.strftime("%y %m %d %H %M\n")
    ctrl += "1\n"
    ctrl += str(req.lat)+" "+str(req.lon)+" "+str(req.hs[0]) + "\n" # Location
    ctrl += str(req.runtime) + "\n" # Run time
    ctrl += str(vMotion) + "\n" # Vertical motion option
    ctrl += str(modelTop) + "\n" # Model top (m)
    ctrl += str(len(metFiles)) + "\n" # Number of met files
    # Met files
    for met in metFiles:
        splits = met.split('/')
        for s in splits[:-1]:
            ctrl += s+"/"
        ctrl += "\n"+splits[-1]+"\n"
    ctrl += storage_path+"\n"
    # Output name
    if outName:
        ctrl += outName+"\n" # Given output file name
    else:
        ctrl += req.ID+"_"+req.start.strftime("%Y%m%d_%H%M_")+str(req.hs[0])+".traj\n" # Automatic output file name

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
                        releaseStart=None,samplingStart=None,samplingStop=None,outName=None):

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
    if outName:
        ctrl += outName+"\n" # Given output file name
    else:
        ctrl += req.ID+"_"+req.start.strftime("%Y%m%d_%H%M_")+str(req.hs[0])+".conc\n" # Automatic output file name
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