import pandas as pd
import numpy as np
import time
from sys import exit

def weightStructureFunc(structureType,structureDF,structureWeight,peakPos,tssPos,chrnum):
    #注意这个chrnum是chr1还是1的形式
    if not structureType: 
        return(1)

    if all(isinstance(var, list) for var in [structureDF, structureType, structureWeight]):
        totalweight = 1
        for i in range(len(structureDF)):
            weightone = weightStructureFuncOne(structureType[i],structureDF[i],structureWeight[i],peakPos,tssPos,chrnum)
            totalweight = totalweight * weightone
        return totalweight
    elif all(isinstance(var, str) for var in [structureType]):
        weightone = weightStructureFuncOne(structureType,structureDF,structureWeight,peakPos,tssPos,chrnum)
        return weightone
    else:
        exit(1)

def weightStructureFuncOne(structureType,structureDF,structureWeight,peakPos,tssPos,chrnum):
    # structureDataType in ["boundary","tad","loop","compartment"]
    # boundary is similar to tad, but different.
    # if ture, give a weight to peak-to-gene
    if not structureType: 
        return(1)

    sDataChr = structureDF[structureDF[0]==chrnum]

    if structureType == "boundary":  # cross a boundary
        boundarySite = np.array(sDataChr[1])
        if ((tssPos - boundarySite) * (peakPos - boundarySite) < 0).any():
            return structureWeight
        else:
            return 1
        
    elif structureType == "tad": # within the same tad
        tadstart = np.array(sDataChr[1])
        tadend = np.array(sDataChr[2])
        withinAtad = (peakPos>= tadstart) & (peakPos<= tadend) & (tssPos>=tadstart) & (tssPos<= tadend)
        if sum(withinAtad) >=1:
            return structureWeight
        else:
            return 1
        
    elif structureType == "loop": # connected by a loop
        anchor1_start = np.array(sDataChr[1])
        anchor1_end = np.array(sDataChr[2])
        anchor2_start = np.array(sDataChr[4])
        anchor2_end = np.array(sDataChr[5])
        loopConnect = (
            ((peakPos >= anchor1_start) & (peakPos <= anchor1_end) & 
             (tssPos >= anchor2_start ) & (tssPos <= anchor2_end)) |
            ((peakPos >= anchor2_start ) & (peakPos <= anchor2_end) & 
             (tssPos >= anchor1_start ) & (tssPos <= anchor1_end))
        ).any()
        if loopConnect:
            return structureWeight
        else:
            return 1
        
    elif structureType == "stripe": # connected by a stripe
        anchor1_start = np.array(sDataChr[1])
        anchor1_end = np.array(sDataChr[2])
        anchor2_start = np.array(sDataChr[4])
        anchor2_end = np.array(sDataChr[5])
        stripeConnect = (
            ((peakPos >= anchor1_start) & (peakPos <= anchor1_end) & 
             (tssPos >= anchor2_start ) & (tssPos <= anchor2_end)) |
            ((peakPos >= anchor2_start ) & (peakPos <= anchor2_end) & 
             (tssPos >= anchor1_start ) & (tssPos <= anchor1_end))
        ).any()
        if stripeConnect:
            return structureWeight
        else:
            return 1
        
    elif structureType == "compartmentSameA":  # within a same compartment A
        sDataChr[3] = pd.Categorical(sDataChr[3])
        sDataChr = sDataChr[sDataChr[3]=="compartmentA"]
        structureStart = np.array(sDataChr[1])
        structureEnd = np.array(sDataChr[2])
        withinAstructure = (peakPos>= structureStart) & (peakPos<= structureEnd) & \
                           (tssPos>=structureStart) & (tssPos<= structureEnd)
        if sum(withinAstructure) >=1:
            return structureWeight
        else:
            return 1

    elif structureType == "compartmentSame":  # within a same compartment A or B
        structureStart = np.array(sDataChr[1])
        structureEnd = np.array(sDataChr[2])
        withinAstructure = (peakPos>= structureStart) & (peakPos<= structureEnd) & \
                           (tssPos>=structureStart) & (tssPos<= structureEnd)
        if sum(withinAstructure) >=1:
            return structureWeight
        else:
            return 1

    elif structureType == "compartmentBothA": # both in A but not need to be in the same
        structureStart = np.array(sDataChr[1])
        structureEnd = np.array(sDataChr[2])
        peakCompart = sDataChr[(peakPos>= structureStart) & (peakPos<= structureEnd)][3].values
        tssCompart = sDataChr[(tssPos>= structureStart) & (tssPos<= structureEnd)][3].values
        if peakCompart == "compartmentA" and tssCompart == "compartmentA":
            return structureWeight
        else:
            return 1

    elif structureType == "compartmentBothAvalue": 
        structureStart = np.array(sDataChr[1])
        structureEnd = np.array(sDataChr[2])
        peakCompart = sDataChr[(peakPos>= structureStart) & (peakPos<= structureEnd)][3].values
        tssCompart = sDataChr[(tssPos>= structureStart) & (tssPos<= structureEnd)][3].values
        bothCompartA = (tssCompart>0) and (peakCompart>0)
        if bothCompartA:
            return structureWeight
        else:
            return 1