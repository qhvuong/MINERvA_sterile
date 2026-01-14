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

def MaddWrapper(output_playlist,input_files,is_data):
    args = ["madd", AnalysisConfig.SelectionHistoPath(output_playlist,is_data)]
    args.extend(input_files)
    #cmd = "madd {} {}".format(AnalysisConfig.SelectionHistoPath(output_playlist,is_data)," ".join(input_files))
    #print (cmd)
    #os.system(cmd)
    print(args)
    subprocess.run(args,stdout=subprocess.DEVNULL)
    print ("done")

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
