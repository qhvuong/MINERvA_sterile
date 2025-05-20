import ROOT
import os
import json
import PlotUtils
import config
import inspect
MACRO_ROOT = os.path.dirname(os.path.abspath(__file__))+"/../"
POT_FILE = MACRO_ROOT+"configs/POT.json"
#SCALE_FILE = MACRO_ROOT+"background_fit/scale.json"


def loadJSON(path):
    with open(path,"r") as pl:
        data = json.load(pl)
    return data

def writeJSON(path,data):
    with open(path,"w") as pl:
        json.dump(data, pl, indent=4)

def guessPlaylistName(playlist,st,nickname):
    return "file_option/playlist_ccqenue_{}{}_{}.txt".format(st,playlist,nickname)

def composeTChain(input_txt,tree, start=None, count=None):
    chain = PlotUtils.ChainWrapper(tree)
    if start is None and count is None:
        #process the playlist at once.
        chain.AddPlaylist(input_txt)
    elif start is not None and count is not None:
        #process the playlist at slices for grid parallelization
        with open(input_txt) as playlist:
            counter = 0
            for i, line in enumerate(playlist):
                if  (counter > 0) or i == start:
                    chain.Add(line.rstrip("\n"))
                    print(i,line.rstrip("\n"))
                    counter += 1
                if counter >= count:
                    break
    else:
        #Only recieve start or count varible, cant proceed
        raise ValueError("I don't know how to proceed b/c only one of start/count is given")
    return chain

def findPlaylistFile(playlist,st,nickname):
    data = loadJSON(POT_FILE)
    tree = "NuECCQE"
    try:
        input_txt= data[playlist][st][nickname]["playlist_location"]
        if "tree" in data[playlist][st][nickname]:
            tree = data[playlist][st][nickname]["tree"]
    except KeyError:
        print(("Requested ntuple not found: ",playlist, st, nickname))
        guessing = guessPlaylistName(playlist,st,nickname)
        data.setdefault(playlist,{}).setdefault(st,{})[nickname]={
            "playlist_location" : guessing
        }
        if os.path.isfile(MACRO_ROOT+guessing):
            print(("I guess the playlist file located at {}".format(guessing)))
            writeJSON(POT_FILE,data)
            return MACRO_ROOT+guessing,tree
        else:
            raise KeyError
    return MACRO_ROOT+input_txt,tree

def fileChain(playlist, st, nickname ,tree=None, start = None,count = None):
    """
    Create TChain from name of playlist and name of TTree.
    The Location of root files are contained in a playlist file, which is inferred from playlist name.
    """
    path,main_tree = findPlaylistFile(playlist,st,nickname)
    return composeTChain(path,tree if (tree is not None) else main_tree ,start,count)


def countLines(playlist,st,nickname):
    with open(findPlaylistFile(playlist,st,nickname)[0]) as f:
        nlines = sum(1 for _ in f)
    return nlines

def calcPOT(mytree,getTotal=False):
    sumPOT = 0.0
    print('POT calculation ------')
    for i in range(mytree.GetEntries()):
        sumPOT += mytree.GetValue("POT_Total" if getTotal else "POT_Used",i)
    return sumPOT

def getPOT(playlists,st,nickname, start = None, count=None, cal_POT=False):
    if not isinstance(playlists, list): playlists = [playlists]
    data = loadJSON(POT_FILE)
    playlists_used_POT = 0.0
    local_cal_POT = cal_POT or count is not None 
    for playlist in playlists:
        try:
            used_POT=data[playlist][st][nickname]["used_POT"]
        except KeyError:
            print("Failed to load POT, will calculate")
            local_cal_POT= True
        if local_cal_POT:
            used_POT = calcPOT(fileChain(playlist,st,nickname,"Meta",start,count))
            if count is None:
                data[playlist][st][nickname]["used_POT"]=used_POT
                writeJSON(POT_FILE,data)
        playlists_used_POT+=float(used_POT)
    return playlists_used_POT

def getPOTFromFile(filename,bigGenie=False):
    metatree = ROOT.TChain("Meta")
    if metatree.Add(filename,-1):
        if not bigGenie:
            return ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
        else:
            return 2e17 * ROOT.PlotUtils.POTCounter().getPOTfromTChain(metatree)
    else:
        return None


def ArachneLink(event):
    head = "https://minerva05.fnal.gov/Arachne/arachne.html"
    #head = "https://minerva05.fnal.gov/rodriges/Arachne/arachne.html"
    is_data = False
    if hasattr(event, "mc"):
        is_data = not event.mc
        run = event.ev_run if is_data else event.mc_run
    else:
        try:
            run = event.mc_run
        except AttributeError as RuntimeError:
            #no mc_run value, is data
            run = event.ev_run
            is_data = True

    if is_data:
        subrun = event.ev_subrun 
        gate = event.ev_gate
        slices = event.slice_numbers if isinstance(event.slice_numbers,float) else event.slice_numbers[0]
        head += "?det=MV&recoVer=v21r1p1"
    else:
        subrun = event.mc_subrun
        gate = event.mc_nthEvtInFile+1
        slices = 1
        head +=  "?det=SIM_minerva&recoVer=v21r1p1"

    print((head+"&run=%d&subrun=%d&gate=%d&slice=%d" % (run,subrun,gate,slices)))


def PlaylistLookup(run):
    if (run>=110000 and run<=110149) or  (run>=118000 and run<=118999):
        return "minervame1a"
    elif run>=111000 and run<=111029:
        return "minervame1b"

    elif run>=111030 and run<=111099:
        return "minervame1c"

    elif run>=111100 and run<=111324:
        return "minervame1d"

    elif run>=111325 and run<=111489:
        return "minervame1e"

    elif run>=111490 and run<=111754:
        return "minervame1f"

    elif run>=110150 and run<=110374:
        return "minervame1g"

    elif run>=113000 and run<=113019:
        return "minervame1l"

    elif run>=113020 and run<=113269:
        return "minervame1m"

    elif run>=113270 and run<=113374:
        return "minervame1n"

    elif run>=113375 and run<=113424:
        return "minervame1o"

    elif run>=112000 and run<=112074:
        return "minervame1p"

    #RHC playlists
    elif run>=123000 and run<=123099:
        return "minervame5a"

    elif run>=122000 and run<=122239:
        return "minervame6a"

    elif run>=122240 and run<=122389:
        return "minervame6b"

    elif run>=122390 and run<=122564:
        return "minervame6c"

    elif run>=122565 and run<=122729:
        return "minervame6d"

    elif run>=122730 and run<=122879:
        return "minervame6e"

    elif run>=122880 and run<=122999:
        return "minervame6f"

    elif run>=123100 and run<=123141:
        return "minervame6g"

    elif run>=123250 and run<=123300:
        return "minervame6h"

    elif run>=123425 and run<=123465:
        return "minervame6i"

    elif run>=123467 and run<=123606:
        return "minervame6J"

    elif run>=131000 and run<=131019:
        return "minervame3a"

    #Low Energy
    # neutrino mode (FHC)
    elif run >= 10200 and run < 10250:
        return "minerva1"
    elif run >= 10250 and run < 10255:
        return "minerva7"
    elif run >= 10255 and run < 10260:
        return "minerva9"
    elif run >= 13200 and run < 13300:
        return "minerva13c"
    # anti-neutrino mode (RHC)
    elif run >= 50200 and run < 50250:
        return "minerva5"
    else:
        return None


def DelNestedIterable(inobject):
    if isinstance(inobject,list):
        for _ in inobject:
            DelNestedIterable(_)
    elif isinstance(inobject,dict):
        for _,__ in inobject.items():
            DelNestedIterable(__)
            del _
    else:
        del inobject

def GetHistogram(ifile,plot):
    hist = ifile.Get(plot)
    ROOT.SetOwnership(hist,True)
    if not hist:
        print ("histogram not found: "+plot)
        return None
    return hist

def PopulateMnvH2DByMnvH1D(hist_out_tempalte,hist_in,populateXAxis = None):
    """
    Duplicate value in hist in to a higher dimention histogram like hist_out_template
    PopulatedXAxis = None : There is only one bin in input hist, just duplicate it everybin
    PopulatedXAxis = True : The Y projection of output histogram is NXbins times input histogram
    """
    hist_out = hist_out_template.Clone()
    hist_out.ClearAllErrorBands()
    hist_out.Reset();
    if populateXAxis is not None:
        if hist_in.GetNbins() != (hist_out.GetNbinsY() if populateXAxis else hist_out.GetNbinsX()):
            print("Error, populate histogram with different number of bins")
            exit(1)

    #First CV:
    PopulateTH2DByTH1D(hist_out,hist_in,populateXAxis)
    #Second Fill ErrorBand With CV:
    hist_out.AddMissingErrorBandsWithFillWithCV(hist_in)
    # then universes:
    for bandname in hist_in.GetErrorBandNames():
        h2d_errorband = hist_out.GetVertErrorBand(bandname)
        h1d_errorband = hist_in.GetVertErrorBand(bandname)
        for i in range(h2d_errorband.GetNHists()):
            PopulateTH2DByTH1D(h2d_errorband.GetHist(i),h1d_errorband.GetHist(i),populateXAxis)

    return hist_out

def PopulateTH2DByTH1D(hist_out,hist_in,populateXAxis = None):
    if poplulateXAxis is None:
        for i in range(1,hist_out.GetNbinsX()+1):
            for j in range(1,hist_out.GetNbinsY()+1):
                hist_out.SetBinContent(i,j,hist_in.GetBinContent(1))
                hist_out.SetBinError(i,j,hist_in.GetBinError(1))
    elif populateXAixs:
        for i in range(1,hist_out.GetNbinsX()+1):
            for j in range(1,hist_out.GetNbinsY()+1):
                hist_out.SetBinContent(i,j,hist_in.GetBinContent(j))
                hist_out.SetBinError(i,j,hist_in.GetBinError(j))
    else:
        for i in range(1,hist_out.GetNbinsX()+1):
            for j in range(1,hist_out.GetNbinsY()+1):
                hist_out.SetBinContent(i,j,hist_in.GetBinContent(i))
                hist_out.SetBinError(i,j,hist_in.GetBinError(i))

def decorator_PrintNReturn(func):
    def f(*args,**kwargs):
        r = func(*args,**kwargs)
        #print("returning: {}".format(r))
        return r
    return f

def decorator_ReLU(func,minimal=0):
    def Rectifier(*args,**kwargs):
        return max(minimal,func(*args,**kwargs))
    return Rectifier

def decorator_Cap(func,cap):
    def Capped(*args,**kwargs):
        return min(cap,func(*args,**kwargs))
    return Capped

def decorator_RestrictByBinning(func,binning):
    return decorator_Cap(decorator_ReLU(func,binning[0]),binning[-2])

def getFilesAndPOTScale(playlist, type_path_map, ntuple_tag,raw_pot = False):
    pots = [None,None]
    files = [None,None]
    for i,t in enumerate(["data","mc"]):
        try:
            path = type_path_map[t]
        except KeyError:
            continue
        pots[i]= getPOTFromFile(path) or getPOT(playlist,t,ntuple_tag)
        files[i]=ROOT.TFile(path) or None
    print(pots[0], pots[1], 'Data MC ') 
    pot_scale = pots[0]/pots[1] if pots.count(None) == 0 else 1.0
    
    if raw_pot:
        return files[0],files[1],pot_scale,pots[0],pots[1]
    else:
        return files[0],files[1],pot_scale
