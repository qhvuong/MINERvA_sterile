import os
import sys
import re
import multiprocessing as mpi
import subprocess
import ROOT
import PlotUtils
from itertools import chain
from config.AnalysisConfig import AnalysisConfig
from config.SignalDef import SIGNAL_DEFINATION
from tools import Utilities

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SelectionFilesRegex = "(kin|truth)_dist_(data|mc)(.+)_"+ AnalysisConfig.selection_tag+"_"+AnalysisConfig.ntuple_tag+"(_[0-9]+)?\.root"


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
        elif output_hist and bigGenie and any(signal in key.GetName() for signal in SIGNAL_DEFINATION):
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


def MaddWrapper(output_playlist, input_files, is_data):
    valid_files = [f for f in input_files if is_valid_root_file(f)]
    print(f"✅ Found {len(valid_files)} valid input files for playlist {output_playlist}")

    if not valid_files:
        print(f"❌ No valid files to combine for {output_playlist} — skipping.")
        return

    # Optional: sort to make sure the first file is reliable
    valid_files.sort()

    args = ["madd", AnalysisConfig.SelectionHistoPath(output_playlist, is_data)]
    args.extend(valid_files)

    print(f"📦 Running madd on output: {args[1]}")
    subprocess.run(args, stdout=subprocess.DEVNULL)
    print("✅ Finished madd.")


def MergeHistograms():
    for sample_type in dict_of_files:
        #smaple_type is data or mc
        if sample_type not in AnalysisConfig.data_types:
           continue
        if sample_type == "data":
            MaddWrapper(AnalysisConfig.playlist,chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),True)
        elif sample_type == "mc":
            MaddWrapper(AnalysisConfig.playlist,chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),False)
    for special_sample in dict_of_special_mc_samples:
        if "BigGenie" in list(dict_of_special_mc_samples[special_sample].keys())[0]:
            MaddWrapper(AnalysisConfig.playlist+str(special_sample),chain.from_iterable(iter(dict_of_special_mc_samples[special_sample].values())),False)
            MergeTuples(AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,False),AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist+str(special_sample),False),None,None,True)
        else:
            MaddWrapper(AnalysisConfig.playlist+str(special_sample),chain.from_iterable(iter(dict_of_special_mc_samples[special_sample].values())),False)
            MergeTuples(AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist,False),AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist+str(special_sample),False))



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



import argparse

parser = argparse.ArgumentParser(description="Combine ROOT files for a playlist")
parser.add_argument("--i", dest="input_dirs", action="append", required=True,help="Input directory (repeat for multiple)")
parser.add_argument("--playlist", type=str, required=True)
parser.add_argument("--ntuple_tag", type=str, required=True)
parser.add_argument("--selection_tag", type=str, required=True)
parser.add_argument("--mc_only", action="store_true")
parser.add_argument("--data_only", action="store_true")
parser.add_argument("--cal_POT", action="store_true")
args = parser.parse_args()





if __name__ == '__main__':

    
    AnalysisConfig.input_dirs = args.input_dirs
    AnalysisConfig.playlist = args.playlist
    AnalysisConfig.ntuple_tag = args.ntuple_tag
    AnalysisConfig.selection_tag = args.selection_tag

    # define regex here after ntuple_tag and selection_tag are known
    SelectionFilesRegex = (
        r"(kin|truth)_dist_(data|mc)(.+)_" + args.selection_tag + "_" + args.ntuple_tag + r"(_[0-9]+)?\.root"
    )

    # dict_of_files = AddRegexMatchedFiles(AnalysisConfig.input_dir)
    dict_of_files = {}
    for dir_path in AnalysisConfig.input_dirs:
        new_files = AddRegexMatchedFiles(dir_path)
        for sample_type in new_files:
            for reco_type in new_files[sample_type]:
                for sample in new_files[sample_type][reco_type]:
                    dict_of_files.setdefault(sample_type, {}).setdefault(reco_type, {}).setdefault(sample, []).extend(
                        new_files[sample_type][reco_type][sample]
                    )
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