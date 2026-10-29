import os
import sys
import re
import multiprocessing as mpi
import subprocess
import ROOT
import PlotUtils
from itertools import chain
from config.AnalysisConfig import AnalysisConfig
from config.SignalDef import SIGNAL_DEFINITION
from tools import Utilities
import shutil
import math

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SelectionFilesRegex = "(kin|truth)_dist_(data|mc)(.+)_"+ AnalysisConfig.selection_tag+"_"+AnalysisConfig.ntuple_tag+"(_[0-9]+)?\.root"

# # Match only files for the current playlist in the shared output directory.
# # Example matches:
# #   kin_dist_mcle1_p6_TestCuts_UntilEtheta_MAD_10000.root
# #   kin_dist_datale13C_p6_TestCuts_UntilEtheta_MAD_2000.root
# playlist_chunk = AnalysisConfig.playlist
# if not playlist_chunk.endswith("_p6"):
#     playlist_chunk = playlist_chunk + "_p6"
# SelectionFilesRegex = (
#     rf"(kin|truth)_dist_(data|mc)({re.escape(playlist_chunk)})_"
#     rf"{re.escape(AnalysisConfig.selection_tag)}_"
#     rf"{re.escape(AnalysisConfig.ntuple_tag)}"
#     rf"(_[0-9]+)?\.root"
# )


def is_valid_root_file(filepath):
    try:
        print(f"🔍 Checking file: {filepath}")
        f = ROOT.TFile.Open(filepath, "READ")
        if not f:
            print(f"⚠️ Could not open file: {filepath}")
            return False
        if f.IsZombie():
            print(f"⚠️ Zombie file: {filepath}")
            f.Close()
            return False
        if f.TestBit(ROOT.TFile.kRecovered):
            print(f"⚠️ Recovered file (incomplete write): {filepath}")
            f.Close()
            return False
        keys = f.GetListOfKeys()
        if not keys or keys.GetSize() == 0:
            print(f"⚠️ File has no histograms: {filepath}")
            f.Close()
            return False
        f.Close()
        return True
    except Exception as e:
        print(f"⚠️ Exception while checking {filepath}: {e}")
        return False


def MakeNo2p2hName(path):
    base, ext = os.path.splitext(path)
    return base + "_no2p2h" + ext
def MakeOnly2p2hName(path):
    base, ext = os.path.splitext(path)
    return base + "_only2p2h" + ext


def AddOneFile(input_string,output_string, pot_scale,bigGenie=False):
    input_file = ROOT.TFile.Open(input_string)
    output_file = ROOT.TFile.Open(output_string,"UPDATE")
    keylist = input_file.GetListOfKeys()
    for key in keylist:
        hist = input_file.Get(key.GetName())
        if isinstance(hist,ROOT.TTree):
            continue
        hist.Scale(pot_scale)
        output_hist = output_file.Get(hist.GetName())
        if output_hist and not bigGenie:
            output_hist.Add(hist)
        elif output_hist and bigGenie and any(signal in key.GetName() for signal in SIGNAL_DEFINITION):
            print("cloning {} from BigGenie selection".format(key.GetName()))
            output_hist = hist.Clone()
        elif output_hist:
            output_hist.Add(hist)
        else:
            output_hist = hist.Clone()
        output_hist.Write("",ROOT.TObject.kOverwrite)
        del hist

    keylist.Clear()
    del keylist
    input_file.Clear()
    input_file.Close()
    del input_file
    print("done a file")



def MergeTH3Histograms(input_files, output_path):
    """
    Merge CV-only TH3D histograms from ROOT files.

    MnvH1D/MnvH2D objects are handled separately by madd.
    This function overwrites only TH3D objects in the merged file.
    """
    merged = {}
    expected_names = None

    input_files = list(input_files)

    if not input_files:
        raise RuntimeError("[TH3D] No input files provided")

    for path in input_files:
        f = ROOT.TFile.Open(path, "READ")

        if not f or f.IsZombie():
            raise RuntimeError(
                "[TH3D] Cannot open input: {}".format(path)
            )

        try:
            names = {
                key.GetName()
                for key in f.GetListOfKeys()
                if key.GetClassName() == "TH3D"
            }

            if expected_names is None:
                expected_names = names

                if not names:
                    raise RuntimeError(
                        "[TH3D] No TH3D histograms in {}".format(path)
                    )

            elif names != expected_names:
                raise RuntimeError(
                    "[TH3D] Histogram mismatch in {}\n"
                    "Missing: {}\n"
                    "Extra: {}".format(
                        path,
                        sorted(expected_names - names),
                        sorted(names - expected_names),
                    )
                )

            for name in sorted(names):
                h = f.Get(name)

                if not h or not h.InheritsFrom("TH3D"):
                    raise RuntimeError(
                        "[TH3D] Cannot read {} from {}".format(
                            name, path
                        )
                    )

                if name not in merged:
                    merged[name] = h.Clone(name)
                    merged[name].SetDirectory(0)
                    merged[name].Sumw2()
                    continue

                target = merged[name]

                for axis_target, axis_input in (
                    (target.GetXaxis(), h.GetXaxis()),
                    (target.GetYaxis(), h.GetYaxis()),
                    (target.GetZaxis(), h.GetZaxis()),
                ):
                    if axis_target.GetNbins() != axis_input.GetNbins():
                        raise RuntimeError(
                            "[TH3D] Bin count mismatch: {} in {}".format(
                                name, path
                            )
                        )

                    for ibin in range(1, axis_target.GetNbins() + 2):
                        edge_a = axis_target.GetBinLowEdge(ibin)
                        edge_b = axis_input.GetBinLowEdge(ibin)

                        if not math.isclose(
                            edge_a,
                            edge_b,
                            rel_tol=1e-10,
                            abs_tol=1e-12,
                        ):
                            raise RuntimeError(
                                "[TH3D] Bin edge mismatch: "
                                "{} in {}".format(name, path)
                            )

                if not target.Add(h):
                    raise RuntimeError(
                        "[TH3D] ROOT Add failed: {} in {}".format(
                            name, path
                        )
                    )

        finally:
            f.Close()

    fout = ROOT.TFile.Open(output_path, "UPDATE")

    if not fout or fout.IsZombie():
        raise RuntimeError(
            "[TH3D] Cannot open output: {}".format(output_path)
        )

    try:
        fout.cd()

        for name, h in merged.items():
            if h.Write(name, ROOT.TObject.kOverwrite) <= 0:
                raise RuntimeError(
                    "[TH3D] Failed to write {}".format(name)
                )

            print(
                "[TH3D] {}: entries={}, integral={:.12g}".format(
                    name,
                    h.GetEntries(),
                    h.Integral(),
                ),
                flush=True,
            )

    finally:
        fout.Close()

    print(
        "[TH3D] Successfully merged {} histograms from {} files".format(
            len(merged),
            len(input_files),
        ),
        flush=True,
    )



# def MaddWrapper(output_playlist, input_files, is_data):
#     valid_files = [f for f in input_files if is_valid_root_file(f)]
#     print(f"✅ Found {len(valid_files)} valid input files for playlist {output_playlist}")

#     if not valid_files:
#         print(f"❌ No valid files to combine for {output_playlist} — skipping.")
#         return

#     # Optional: sort to make sure the first file is reliable
#     valid_files.sort()

#     args = ["madd", AnalysisConfig.SelectionHistoPath(output_playlist, is_data)]
#     args.extend(valid_files)

#     print(f"📦 Running madd on output: {args[1]}")
#     subprocess.run(args, stdout=subprocess.DEVNULL)
#     print("✅ Finished madd.")


def MaddWrapper(output_playlist, input_files, is_data, output_file=None):
    input_files = list(input_files)

    valid_files = [
        f for f in input_files
        if is_valid_root_file(f)
    ]

    if len(valid_files) != len(input_files):
        raise RuntimeError(
            "[MADD] {} invalid input files".format(
                len(input_files) - len(valid_files)
            )
        )

    if not valid_files:
        raise RuntimeError(
            "[MADD] No files for playlist {}".format(output_playlist)
        )

    valid_files.sort()

    out = (
        output_file
        or AnalysisConfig.SelectionHistoPath(output_playlist, is_data)
    )

    if os.path.exists(out):
        os.remove(out)

    args = ["madd", out] + valid_files

    print(
        "[MADD] Combining {} files into {}".format(
            len(valid_files), out
        ),
        flush=True,
    )

    subprocess.run(args, check=True)

    if not os.path.exists(out):
        raise RuntimeError(
            "[MADD] Output was not created: {}".format(out)
        )

    MergeTH3Histograms(valid_files, out)

    print(
        "[MADD] Completed playlist {}".format(output_playlist),
        flush=True,
    )





# def MergeHistograms():
#     for sample_type in dict_of_files:
#         #smaple_type is data or mc
#         if sample_type not in AnalysisConfig.data_types:
#            continue
#         if sample_type == "data":
#             MaddWrapper(AnalysisConfig.playlist,chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),True)
#         elif sample_type == "mc":
#             MaddWrapper(AnalysisConfig.playlist,chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),False)
#     for special_sample in dict_of_special_mc_samples:
#         if "BigGenie" in list(dict_of_special_mc_samples[special_sample].keys())[0]:
#             MaddWrapper(AnalysisConfig.playlist+str(special_sample),chain.from_iterable(iter(dict_of_special_mc_samples[special_sample].values())),False)
#             MergeTuples(AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,False),AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist+str(special_sample),False),None,None,True)
#         else:
#             MaddWrapper(AnalysisConfig.playlist+str(special_sample),chain.from_iterable(iter(dict_of_special_mc_samples[special_sample].values())),False)
#             MergeTuples(AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,False),AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist+str(special_sample),False))



def MergeHistograms():
    nominal_mc = None
    has_special_sample = len(dict_of_special_mc_samples) > 0

    # First produce the standard merged outputs
    for sample_type in dict_of_files:
        if sample_type not in AnalysisConfig.data_types:
            continue

        if sample_type == "data":
            MaddWrapper(
                AnalysisConfig.playlist,
                chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),
                True
            )

        elif sample_type == "mc":
            MaddWrapper(
                AnalysisConfig.playlist,
                chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),
                False
            )

            nominal_mc = AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist, False)

            if has_special_sample:
                mc_no2p2h = MakeNo2p2hName(nominal_mc)

                if os.path.exists(nominal_mc):
                    shutil.copy2(nominal_mc, mc_no2p2h)
                    print(f"Saved nominal-only MC as: {mc_no2p2h}")
                else:
                    raise FileNotFoundError(f"Nominal MC file not found: {nominal_mc}")

    if nominal_mc is None:
        return

    if not has_special_sample:
        return

    # Then build each special-sample-only file and add it into the usual MC file
    for i, special_sample in enumerate(dict_of_special_mc_samples):
        special_files = chain.from_iterable(
            iter(dict_of_special_mc_samples[special_sample].values())
        )

        if len(dict_of_special_mc_samples) == 1:
            special_output_file = MakeOnly2p2hName(nominal_mc)
        else:
            base, ext = os.path.splitext(MakeOnly2p2hName(nominal_mc))
            special_output_file = f"{base}_{i}{ext}"

        MaddWrapper(
            AnalysisConfig.playlist,
            special_files,
            False,
            output_file=special_output_file
        )

        first_key = list(dict_of_special_mc_samples[special_sample].keys())[0]
        if "BigGenie" in first_key:
            MergeTuples(nominal_mc, special_output_file, None, None, True)
        else:
            MergeTuples(nominal_mc, special_output_file)

# def MergeHistograms():
#     nominal_mc = None

#     # First produce the standard merged outputs
#     for sample_type in dict_of_files:
#         if sample_type not in AnalysisConfig.data_types:
#             continue

#         if sample_type == "data":
#             MaddWrapper(
#                 AnalysisConfig.playlist,
#                 chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),
#                 True
#             )

#         elif sample_type == "mc":
#             MaddWrapper(
#                 AnalysisConfig.playlist,
#                 chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),
#                 False
#             )

#             nominal_mc = AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist, False)
#             mc_no2p2h = MakeNo2p2hName(nominal_mc)

#             if os.path.exists(nominal_mc):
#                 shutil.copy2(nominal_mc, mc_no2p2h)
#                 print(f"Saved nominal-only MC as: {mc_no2p2h}")
#             else:
#                 raise FileNotFoundError(f"Nominal MC file not found: {nominal_mc}")

#     if nominal_mc is None:
#         return

#     # Then build each special-sample-only file and add it into the usual MC file
#     for i, special_sample in enumerate(dict_of_special_mc_samples):
#         special_files = chain.from_iterable(
#             iter(dict_of_special_mc_samples[special_sample].values())
#         )

#         if len(dict_of_special_mc_samples) == 1:
#             special_output_file = MakeOnly2p2hName(nominal_mc)
#         else:
#             base, ext = os.path.splitext(MakeOnly2p2hName(nominal_mc))
#             special_output_file = f"{base}_{i}{ext}"

#         MaddWrapper(
#             AnalysisConfig.playlist,
#             special_files,
#             False,
#             output_file=special_output_file
#         )

#         first_key = list(dict_of_special_mc_samples[special_sample].keys())[0]
#         if "BigGenie" in first_key:
#             MergeTuples(nominal_mc, special_output_file, None, None, True)
#         else:
#             MergeTuples(nominal_mc, special_output_file)


def MergeTuples(tuple1,tuple2,pot1=None,pot2=None,bigGenie=False):
    """
    Wanna merge tuple two tuples, by direct merging or scale by reading Meta tree.
    """
    #tuple2 is special sample
    #tuple1 is OG mc sample
    pot1 = pot1 or Utilities.getPOTFromFile(tuple1)
    pot2 = pot2 or Utilities.getPOTFromFile(tuple2,bigGenie)
    print("tuple1: ",tuple1," with POT ",pot1)
    print("tuple2: ",tuple2," with POT ",pot2)
    print("scaling the special sample with potscale={:.2f}".format(pot1/pot2))
    AddOneFile(tuple2,tuple1,pot1/pot2,bigGenie)

def AddRegexMatchedFiles(dir_path,f = None):
    if f is None:
        files = os.listdir(dir_path)
    else:
        files = [f]
    filesmap = {}
    count = 0
    for string in files:
        match = re.match(SelectionFilesRegex,string)
        #print match.group(1,2,3)
        if match is not None:
            filesmap.setdefault(match.group(2),{}).setdefault(match.group(1),{}).setdefault(match.group(3),[]).append(dir_path+"/"+match.string)
            count+=1
    print(("added {} files".format(count)))
    return filesmap

if __name__ == '__main__':

    dict_of_files = AddRegexMatchedFiles(AnalysisConfig.input_dir)
    for sample_type in dict_of_files:
        if sample_type == "data":

            print("found data")
    dict_of_special_mc_samples={}
    i=0
    while (True):
        print("Do you want to add more special MC sample?")
        print("Type in the path or directory if yes.")
        print("Press Enter if no.")
        path = input()
        if len(path)==0:
            break
        elif os.path.isdir(path):
            tmp = AddRegexMatchedFiles(path)
        elif os.path.isfile(path):
            tmp = AddRegexMatchedFiles("/".join(path.split("/")[:-1]),path.split("/")[-1])
        else:
            tmp = {}
        if tmp:
            dict_of_special_mc_samples[i] = tmp["mc"]["kin"]
            i+=1
        else:
            print("I can't find the directory or file you typed.")
            print(("Your input is {}").format(path))

    MergeHistograms()