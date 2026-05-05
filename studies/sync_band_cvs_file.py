#!/usr/bin/env python

import sys
import ROOT
import PlotUtils

ROOT.TH1.AddDirectory(False)

def sync_mnv_band_cvs_in_file(path):
    f = ROOT.TFile.Open(path, "UPDATE")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open file: {}".format(path))

    keys = f.GetListOfKeys().Clone()
    n_synced = 0
    n_skipped = 0

    for key in keys:
        name = key.GetName()
        obj = key.ReadObj()

        try:
            if obj.InheritsFrom("PlotUtils::MnvH1D"):
                cv = ROOT.TH1D(obj)
                cv.SetDirectory(0)

                for bandname in obj.GetVertErrorBandNames():
                    band = obj.GetVertErrorBand(str(bandname))
                    if band:
                        cv.Copy(band)

                obj.Write(name, ROOT.TObject.kOverwrite)
                n_synced += 1

            elif obj.InheritsFrom("PlotUtils::MnvH2D"):
                cv = ROOT.TH2D(obj)
                cv.SetDirectory(0)

                for bandname in obj.GetVertErrorBandNames():
                    band = obj.GetVertErrorBand(str(bandname))
                    if band:
                        cv.Copy(band)

                obj.Write(name, ROOT.TObject.kOverwrite)
                n_synced += 1

            else:
                n_skipped += 1

        except Exception as e:
            print("[WARN] Could not sync {}: {}".format(name, e))
            n_skipped += 1

    f.Close()
    print("[DONE] Synced {} MnvH* objects in {}".format(n_synced, path))
    print("[DONE] Skipped {} objects".format(n_skipped))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sync_band_cvs_file.py file.root")
        sys.exit(1)

    sync_mnv_band_cvs_in_file(sys.argv[1])