#!/usr/bin/python
#coding=utf-8
# Created on 
# Author: Zhihua Pei
# Organization: GenePlus
# website: www.zilhua.com
# github: github.com/zilhua
'''
R package:facet -> python
https://github.com/mskcc/facets
'''


import pandas as pd
import numpy as np
import re,sys,time


def readSNP(file=None, err_thresh=0, del_thresh=0):
    dat = pd.read_csv(file, dtype={"Position": np.int32,
                                   "Chromosome": np.str})
    ii = (dat["File1E"] <= err_thresh) & (dat["File1D"] <= del_thresh) \
            & (dat["File2E"] <= err_thresh) & (dat["File2D"] <= \
            del_thresh)
    dat_f = dat[ii].ix[:, 0:2]
    dat_f["NOR.DP"] = dat["File1R"][ii] + dat["File1A"][ii]
    dat_f["NOR.RD"] = dat["File1R"][ii]
    dat_f["TUM.DP"] = dat["File2R"][ii] + dat["File2A"][ii]
    dat_f["TUM.RD"] = dat["File2R"][ii]
    dat_f["Chromosome"] = dat_f["Chromosome"].str.replace("chr", "",
                        flags=re.I)
    dat_f = dat_f[dat_f["Chromosome"].isin([str(i) for i in xrange(23)] + ['X'])]
    dat_f.loc[dat_f["Chromosome"] == "X", "Chromosome"] = "23"
    dat_f["Chromosome"] = dat_f["Chromosome"].astype(int)
    #print dat_f.shape
    preProcSample(dat_f)
    return dat_f


def preProcSample(rcmat, ndepth=35, het_thresh=0.25, snp_window=250,
            cval=25, deltaCN=0, hetscale=True, unmatched=False,
            ndepthmax=5000):
    rcmat = _procSnps(rcmat, ndepth, het_thresh, snp_window, unmatched,
                     ndepthmax)
    _counts2logROR(rcmat, gbuild="hg19", unmatched=unmatched, f=0.2)
    #print rcmat.shape
    #print rcmat


def _procSnps(rcmat, ndepth=35, het_thresh=0.25, snp_window=250,
             unmatched=False, ndepthmax=5000):
    rcmat = rcmat[(rcmat["NOR.DP"] >= ndepth ) & \
            (rcmat["NOR.DP"] <= ndepthmax)]
    out = pd.DataFrame(columns=["chrom", "maploc", "rCountT", "rCountN",
                                "vafT", "vafN", "refN", "het"])
    out["chrom"] = rcmat["Chromosome"]
    out["maploc"] = rcmat["Position"]
    out["rCountT"] = rcmat["TUM.DP"]
    out["rCountN"] = rcmat["NOR.DP"]
    out["vafT"] = 1 - rcmat["TUM.RD"]/rcmat["TUM.DP"]
    out["vafN"] = 1 - rcmat["NOR.RD"]/rcmat["NOR.DP"]
    out["refN"] = rcmat["NOR.RD"]/rcmat["NOR.DP"]
    out = out.sort_values(by=["chrom", "maploc"])
    out = out.reset_index(drop=True)
    if unmatched:
        if het_thresh == 0.25:
            het_thresh = 0.1
        out["het"] = 1 * ((out[["vafN", "refN"]].min(axis=1) > het_thresh) \
                          & (out["rCountT"] > 50))
    else:
        out["het"] = 1 * (out[["vafN", "refN"]].min(axis=1) > het_thresh)
    out["keep"] = _scansnp(out, snp_window=snp_window)
    return out


def _scansnp(data, snp_window=250):
    n = len(data)
    maploc = list(data["maploc"])
    het = list(data["het"])
    keep = np.zeros(n)
    nbhd = snp_window
    i = 0; keep[i] = 1; nsnp = 1.0; nhet = het[i]; isel = i
    for j in xrange(1, n, 1):
        if abs(maploc[j] - maploc[i]) <= nbhd:
            nsnp = nsnp + 1.0
            nhet = nhet + het[j]
            usnp = np.random.uniform()
            if nhet > 0:
                if het[j] == 1:
                    if usnp <= 1.0 / nhet:
                        keep[isel] = 0
                        keep[j] = 1.0
                        isel = j
                    else:
                        keep[j] = 0.0
                else:
                    keep[j] = 0.0
            else:
                if usnp <= 1.0 / nsnp:
                    keep[isel] = 0
                    keep[j] = 1.0
                    isel = j
                else:
                    keep[j] = 0.0
        else:
            i = j
            isel = i
            keep[i] = 1
            nsnp = 1.0
            nhet = het[i]
    return keep


def _counts2logROR(dat, gbuild="hg19", unmatched=False, f=0.2):
    out = dat
    out.loc[:, "gcpct"] = 0
    out = dat[dat["keep"] == 1]
    nchr = max(out["chrom"])
    print "--------->ROR向下"
    print out.shape
    print nchr
    print "--------->ROR"


def _scansnp2(dat, snp_window=250):
    #nsnp, nhet and isel  are initialized
    times = time.time()
    i = 0
    keep = np.zeros(len(dat))
    het = dat["het"]
    nsnp = 1.0
    nhet = het[i]
    isel = i
    keep[i] = 1
    for index, row in dat.iterrows():
        if index == 0:
            continue
        if abs(row["maploc"] - dat.iloc[i]["maploc"]) <= snp_window:
            nsnp = nsnp + 1.0
            nhet = nhet + het[index]
            #generate a random number to decide which one is kept
            usnp = np.random.uniform()
            if nhet > 0:
                if het[index] == 1:
                    if usnp <= 1/nhet:
                        keep[isel] = 0
                        keep[index] = 1.0
                        isel = index
                    else:
                        keep[index] = 0
                else:
                    keep[index] = 0
            else:
                if usnp <= 1/nsnp:
                    keep[isel] = 0
                    keep[index] = 1.0
                    isel = index
                else:
                    keep[index] = 0.0
        else:
            i = index
            isel = i
            keep[i] = 1
            nsnp = 1.0
            nhet = het[i]
        # print "1<<<<<<<<<<<<<<<"
        # print "i = " + str(i)
        # print "keep = " + str(keep[0]) + str(keep[1])
        # print "het = " + str(het[0]) + "\t" + str(het[1])
        # print "nsnp = " + str(nsnp)
        # print "nhet = " + str(nhet)
        # break
    times = time.time() - times
    print "---------------------->1"
    print "time = " + str(times)
    #print ', '.join([str(i) for i in keep])
    return keep


if __name__ == "__main__":
    readSNP(file="stomach.csv")