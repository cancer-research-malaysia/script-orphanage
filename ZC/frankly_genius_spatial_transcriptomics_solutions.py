#!/usr/bin/env python3
# coding: utf-8
"""
@author: Tan Zi Ching       ziching.tan@cancerresearch.my
@last modified: 18/11/2024
@title: Functions for MyBrCa spatial transcriptomics processing/analysis solutions
"""

def adjust(df, d, subplot=False, factor=2): # for plotting 8 subplots
    """
    adjusting figures padding to ensure each fig has the same ratio
    usage: r1, xlim, ylim, pad, addpadx, addpady, dot_size = adjust(df, d)
    """
    TARGET_RATIO = 1.5
    FIXED_X_LEN = 15000
    FIXED_Y_LEN = 10000 
    xlim = None
    ylim = None 
    dx1 = max(df["x"])-min(df["x"])
    dy1 = max(df["y"])-min(df["y"])
    r1 = dx1/dy1 #round(dx1/dy1,1)
    dx2 = dx1
    dy2 = dy1
    dx2 = TARGET_RATIO*dy1
    pad = abs(dx1-dx2)/2
    dy2 = dx1/TARGET_RATIO
    pad = abs(dy1-dy2)/2

    # additional padding for small samples
    addpadx = abs((FIXED_X_LEN-dx2)/2)
    addpady = abs((FIXED_Y_LEN-dy2)/2)
    #factor2 = FIXED_X_LEN/dx2
    factor2 = ((FIXED_X_LEN/dx2) + (FIXED_Y_LEN/dy2))/2
    if subplot:
        dot_size = (.3*factor) * factor2 #3
        if d == "SD1043": dot_size = .2
    else:
        dot_size = (2*factor) * factor2 #3
        if d == "SD1043": dot_size = 1.5

    return r1, xlim, ylim, pad, addpadx, addpady, dot_size

def get_neighbourhood(df, x, y, reach=1):
    """
    returns dataframe of surrounding neighbourhood based on input <reach>
    """
    pad = reach * bin_size
    return df[((df["x"]>=x-pad) & (df["x"]<=x+pad)) & ((df["y"]>=y-pad) & (df["y"]<=y+pad))]

## classify each bin as high/low - based on threshold
## including neighbourhood (immediately sekililing)
def mutations_classify_neighbourhood(thresh=0.8, reach=1):
    """
    classifies and returns dataframe + surrounding neighbourhood based on input <thresh> and <reach>
    """
    all_muts_classified = []
    for muttype in muttypes:
        #print(muttype)
        sub = []
        for i,d in enumerate(datalist):
            #print(d)
            #df = all_muts_combined[muttypes.index(muttype)][all_muts_combined[muttypes.index(muttype)]["sample"]==d]]
            df = all_muts[muttypes.index(muttype)][i].copy()
            #df1 = df[df["normlog%s"%muttype]>thresh]
            df = df[df["normlog%s"%muttype]>0]
            df["mutation_status"] = df["normlog%s"%muttype] > thresh
            df["mutation_status"].replace({True:"high", False:"low"}, inplace=True)
            high = df[df["mutation_status"]=="high"]
            highs = []
            for j in range(len(high)):
                highs.append(get_neighbourhood(df, high.iloc[j].x, high.iloc[j].y, reach=reach))
            all_highs = pd.concat(highs) ## highs + neighbourhood of lows
            #print(len(all_highs))
            all_highs.drop_duplicates(inplace=True)
            #print(len(all_highs))
            all_highs["mutation_status_neighbourhood"] = "high"
            all_lows = df[~df.isin(all_highs)].dropna() ## get the lows
            all_lows["mutation_status_neighbourhood"] = "low"
            newdf = pd.concat([all_highs, all_lows])
            sub.append(newdf)
        all_muts_classified.append(sub)
    
    return all_muts_classified