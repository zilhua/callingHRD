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
import re,sys


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
    preProcSample(dat_f)
    return dat_f


def preProcSample(rcmat, ndepth=35, het_thresh=0.25, snp_window=250,
            cval=25, deltaCN=0, hetscale=True, unmatched=False,
            ndepthmax=5000):
    pmat = procSnps(rcmat, ndepth, het_thresh, snp_window, unmatched, ndepthmax)
    print pmat


def procSnps(rcmat, ndepth=35, het_thresh=0.25, snp_window=250,
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
    if unmatched:
        if het_thresh == 0.25:
            het_thresh = 0.1
        out["het"] = 1 * ((out[["vafN", "refN"]].min(axis=1) > het_thresh) \
                          & (out["rCountT"] > 50))
    else:
        out["het"] = 1 * (out[["vafN", "refN"]].min(axis=1) > het_thresh)
    scansnp(out, snp_window=snp_window)
    print out

def scansnp(dat, snp_window=250):
    pass


if __name__ == "__main__":
    readSNP(file="test.csv")