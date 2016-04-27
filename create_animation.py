import numpy as np
import matplotlib.pyplot as plt
import circles
import os
import time
import sys
from matplotlib.collections import LineCollection

def print_fl(x):
    print(x)
    sys.stdout.flush()
    

InFile = sys.argv[1]
OutDir = sys.argv[2]

if not os.path.exists(OutDir):
    os.makedirs(OutDir)

Data = np.load(InFile)

NumSteps = Data["SaveLocation"].shape[0]

StartStep = int(float(sys.argv[3]) / 100.0 * NumSteps)
EndStep = int(float(sys.argv[4]) / 100.0 * NumSteps)

print_fl("Rendering frames %d through %d" % (StartStep, EndStep))

ExtremeBorders = np.zeros((NumSteps, 2, 2))
ExtremeBorders[:,0] = np.min(Data["SaveLocation"], axis=1)
ExtremeBorders[:,1] = np.max(Data["SaveLocation"], axis=1)

Centers = (ExtremeBorders[:,1,:] + ExtremeBorders[:,0,:])/2
Sizes = ExtremeBorders[:,1,:] - ExtremeBorders[:,0,:]
ExtremeSize = np.max(Sizes)

ExtremePositions = np.zeros((2,2)) # minmax x dim
ExtremePositions[0,:] = np.min(ExtremeBorders[:,0,:], axis=0)
ExtremePositions[1,:] = np.max(ExtremeBorders[:,1,:], axis=0)
MigrationReach = np.max(ExtremePositions[1,:] - ExtremePositions[0,:])
MigrationCenter = (ExtremePositions[1,:] + ExtremePositions[0,:])/2

SaveLocation = Data["SaveLocation"]
SaveOrientation = Data["SaveOrientation"]
Radius = Data["Radius"]

NumParticles = SaveLocation.shape[1]

MarginFactorSheet = 1.1
MarginFactorCloseup = 1.0

for t in range(StartStep, EndStep):
    plt.figure(figsize=(18,8))
    axSheet = plt.subplot(121, aspect='equal')
    circles.circles(SaveLocation[t,:,0], SaveLocation[t,:,1], Radius, fc='none')

    plt.xlim(MigrationCenter[0] + np.array([-1, 1]) * MigrationReach/2 * MarginFactorSheet)
    plt.ylim(MigrationCenter[1] + np.array([-1, 1]) * MigrationReach/2 * MarginFactorSheet)

    axCloseUp = plt.subplot(122, aspect='equal')
    circles.circles(SaveLocation[t,:,0], SaveLocation[t,:,1], Radius, fc='none')

    EndPoints = SaveLocation[t,:,:] + np.array([Radius * np.cos(SaveOrientation[t,:]), Radius * np.sin(SaveOrientation[t,:])]).T
    LineSegments = np.zeros((NumParticles, 2, 2))
    LineSegments[:,0,:] = SaveLocation[t,:,:]
    LineSegments[:,1,:] = EndPoints

    lines = LineCollection(LineSegments, color='k')
    axCloseUp.add_collection(lines)

    plt.xlim(Centers[t,0] + np.array([-1, 1]) * ExtremeSize/2 * MarginFactorCloseup)
    plt.ylim(Centers[t,1] + np.array([-1, 1]) * ExtremeSize/2 * MarginFactorCloseup)
    
    plt.savefig(OutDir + ("/frame-%08d.png" % t), dpi=300)
    plt.cla()
    plt.clf()
    plt.close()
    if (t % (0.01) * NumSteps) < 1:
        print_fl("Progress: %d%%" % (t*100/NumSteps))