from sourcePy import source
import pytest
import datetime as dt
import numpy as np

# Create a test set of measurements
@pytest.fixture
def test_measurements():
    meass = []
    duration = dt.timedelta(hours=3,minutes=0)
    meass.append(source.measurement(0,0,10,dt.datetime(2000,1,1,0,0),'ptch',2,"testSite1",duration=duration))
    meass.append(source.measurement(0,0,10,dt.datetime(2000,1,1,1,0),'ptch',1,"testSite1",duration=duration))
    meass.append(source.measurement(0,0,10,dt.datetime(2000,1,1,2,0),'ptch',0,"testSite1",duration=duration))
    meass.append(source.measurement(1,1,10,dt.datetime(2000,1,1,0,0),'ptch',1,"testSite2",duration=duration))
    meass.append(source.measurement(1,1,10,dt.datetime(2000,1,1,1,0),'ptch',.5,"testSite2",duration=duration))
    meass.append(source.measurement(1,1,10,dt.datetime(2000,1,1,2,0),'ptch',0,"testSite2",duration=duration))
    yield meass

# Create a test set of filled trajectories
@pytest.fixture
def test_valid_trajectories(test_measurements):
    meass = test_measurements
    trajs = []
    for i in range(len(meass)):
        meas = meass[i]
        time = meas.datetime
        j = 0
        datetimes,lats,lons,hs,ps = [],[],[],[],[]
        while time <= meas.datetime+meas.duration:
            datetimes.append(time)
            lats.append(meas.lat+(j*(i-3)*(j-.2)*.2))
            lons.append(meas.lon+(j*(i-3)*(j-.2)*.2))
            hs.append(10)
            ps.append(1000)
            j += .1
            time += dt.timedelta(minutes=10)
        trajs.append(source.trajectory(fName=None,direct=None,datetimes=datetimes,lats=lats,lons=lons,hs=hs,ps=ps,measurements=[meas,],data=None))
    yield trajs

# Create a test set of unfilled trajectories
@pytest.fixture
def test_empty_trajectories(test_measurements):
    trajs = []



##### Little Stuff #####
# Proper unit tests of smaller functions

@pytest.mark.parametrize("date",[(dt.datetime(2000,1,1,0,0)),(dt.datetime(2000,1,1,0,29)),(dt.datetime(1999,12,31,23,30)),(dt.datetime(1999,12,31,23,45))])
def test_round_to_hour(date):
    assert source.round_to_hour(date) == dt.datetime(2000,1,1,0,0)


def test_makeGrid():
    gridLons,gridLats = source.makeGrid(0,6,0,6,5,5)
    assert np.array_equal(gridLons,np.asarray([[0,5],[0,5.0]]))
    assert np.array_equal(gridLats,np.asarray([[0,0],[5,5]]))

    gridLons,gridLats = source.makeGrid(-.5,.5,-2,2,.25,2)
    assert np.array_equal(gridLons,np.asarray([[-0.5,-0.25,0.,0.25],[-0.5,-0.25,0.,0.25]]))
    assert np.array_equal(gridLats,np.asarray([[-2,-2,-2,-2],[0,0,0,0]]))


def test_getUniqueRecs(test_valid_trajectories):
    uniqueIDs,lons,lats = source.getUniqueRecs(test_valid_trajectories)
    assert np.array_equal(uniqueIDs,['testSite1','testSite2'])
    assert np.array_equal(lons,[0,1])
    assert np.array_equal(lats,[0,1])

    uniqueIDs,lons,lats = source.getUniqueRecs([])
    assert np.array_equal(uniqueIDs,[]) and np.array_equal(lons,[]) and np.array_equal(lats,[])


def test_listVals(test_valid_trajectories):
    assert np.array_equal(source.listVals(test_valid_trajectories),[2, 1, 0, 1, 0.5, 0])
    assert np.array_equal(source.listVals([]),[])


##### Big Stuff #####
# More like case testing, some simple stuff
def test_calcPSCF_validTrajs():
    # Answers
    answer_nCounts = np.asarray([[ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,2,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,14,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,41,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,49,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,8,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0]])
    answer_mCounts = np.asarray([[ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,2,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,6,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,11,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,19,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0],
                                 [ 0,0,0,0,0,0,0,0,0,0,0,0]])
    answer_PSCF = np.asarray([[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,1.0,   np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,0.42857143,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,0.26829268,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,0.3877551,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,0.0,   np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                              [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]])
    # Create input
    
