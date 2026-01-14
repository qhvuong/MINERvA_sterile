import os
import math
import ROOT
import PlotUtils
from PlotUtils.HistWrapper import HistWrapper
from array import array
from functools import partial
import copy
from tools.PlotLibrary import TranslateSettings,VariantPlotsNamingScheme
from config.SignalDef import TRUTH_CATEGORIES,EXTRA_OTHER
from config.CutConfig import SAMPLE_CUTS,KINEMATICS_CUTS
from tools.CutLibrary import CUTS

class HistWrapper1D(HistWrapper):
    def __init__(self,title,bins):
        # Note: bins is an array of lower edges of each bin
        super(HistWrapper1D,self).__init__(title,len(bins)-1,array('d',bins),{})

    @property
    def name(self):
        return self.hist.GetName()

    @name.setter
    def name(self,name):
        self.hist.SetName(name)

    # def AddErrorBands(self,cv_universes):
    #     self.AddUniverses(cv_universes)

    # def Fill(self, value, wgt, universe):
    #     self.FillUniverse(universe,value,wgt)
    def Clone(self,name):
        clone = copy.copy(self)
        clone.hist = self.hist.Clone(name)
        clone.eventToUnivMap = {}
        return clone

    def Write(self):
        self.SyncCVHistos()
        self.hist.Write()
        del self.hist

class HistWrapper2D(HistWrapper):
    def __init__(self,title,Xbins,Ybins):
        # Note: bins is an array of lower edges of each bin
        super(HistWrapper2D,self).__init__(title,len(Xbins)-1,array('d',Xbins),len(Ybins)-1,array('d',Ybins),{})

    @property
    def name(self):
        return self.hist.GetName()

    @name.setter
    def name(self,name):
        self.hist.SetName(name)

    def Clone(self,name):
        clone = copy.copy(self)
        clone.hist = self.hist.Clone(name)
        clone.eventToUnivMap = {}
        return clone

    def Write(self):
        self.SyncCVHistos()
        self.hist.Write()
        del self.hist

class MnvResponseWrapper(object):
    def __init__(self,title, bins):
        self.reco_hist = HistWrapper2D(title,bins[0],bins[1])
        self.truth_hist = HistWrapper2D(title,bins[2],bins[3])

        # +2 for overflow/underflow bins of reco/truth histogram
        migrationNbinsX = (self.reco_hist.hist.GetNbinsX()+2)*(self.reco_hist.hist.GetNbinsY()+2)
        migrationNbinsY = (self.truth_hist.hist.GetNbinsX()+2)*(self.truth_hist.hist.GetNbinsY()+2)
        self.migration_hist = HistWrapper2D(title, list(range(migrationNbinsX+1)),list(range(migrationNbinsY+1)))

    @property
    def name(self):
        return self.migration_hist.name

    @name.setter
    def name(self,name):
        self.reco_hist.name = name+"_reco"
        self.truth_hist.name = name+"_truth"
        self.migration_hist.name = name

    def AddUniverses(self,cv_universes):
        self.reco_hist.AddUniverses(cv_universes)
        self.truth_hist.AddUniverses(cv_universes)
        self.migration_hist.AddUniverses(cv_universes)

    def FillUniverse(self,universe, reco_x,reco_y,truth_x,truth_y, wgt):
        mig_x = self.migration_hist.hist.GetXaxis().GetBinCenter(super(PlotUtils.MnvH2D, self.reco_hist.hist).FindBin(reco_x,reco_y))
        mig_y = self.migration_hist.hist.GetYaxis().GetBinCenter(super(PlotUtils.MnvH2D, self.truth_hist.hist).FindBin(truth_x,truth_y))

        #print(mig_x,mig_y)

        self.reco_hist.FillUniverse(universe,reco_x,reco_y,wgt)
        self.truth_hist.FillUniverse(universe,truth_x,truth_y,wgt)
        self.migration_hist.FillUniverse(universe,mig_x,mig_y,wgt)

    def Write(self):
        self.reco_hist.Write()
        self.truth_hist.Write()
        self.migration_hist.Write()


class PlotProcessor():
    def __init__(self,name,title,binning,value_getter = None, cuts = None, weight_function = None,**kwargs):
        
        if value_getter is None:
            print("You didn't tell me what to fill for plot: " + name)
            exit(1)
        
        if len(binning)!=len(value_getter):
            print("inconsistent number of binning and value getter.")
            exit(1)

        for l in binning:
            if len(l)<3 or l[-1]<=l[0]:
                print("wrong binning of plot {}: {}".format(name,l))
                exit(1)

        if len(value_getter) == 1:
            self.histwrapper = HistWrapper1D(title,binning[0])
        elif len (value_getter) == 2:
            self.histwrapper = HistWrapper2D(title,binning[0],binning[1])

        elif len(value_getter) == 4:
            self.histwrapper = MnvResponseWrapper(title, binning)
        else :
            print("Do not support more than 2D histograms yet.")
            exit(1)
        self.value_getter = value_getter
        self.cuts = []
        if cuts is not None:
            for i in cuts:
                self.AddCut(i)

        self.weight_function =  weight_function
        self.histwrapper.name = name

    def Clone(self, name):
        #shalow copy of attrs
        clone = copy.copy(self)
        #deep copy of histwrapper and cuts.
        clone.histwrapper = self.histwrapper.Clone(name)
        clone.cuts = list(self.cuts)
        return clone

    def AddErrorBands(self,universes):
        self.histwrapper.AddUniverses(universes)

    def Process(self,universe):

        if all(cut(universe) for cut in self.cuts[::-1]):
            try:
                value = [_(universe) for _ in self.value_getter]
            except Exception as e:
                print(self.histwrapper.name)
                print(("Error", e))
                return None
        else :
            return None

        wgt = self.weight_function(universe) if self.weight_function else universe.GetWeight()
        if isinstance(value[0],list):
            if not isinstance(wgt,list):
                wgt = len(value[0])*[wgt]
            #filling the multiple entries per event.
            for v in map(lambda wgt,*args: (list(args),wgt), wgt, *value):
                self.FillHist(universe,*v)
        else:
            self.FillHist(universe,value,wgt)

    def FillHist(self,universe,value,wgt):
        if value.count(None) == 0:
            args = value+[wgt]
            self.histwrapper.FillUniverse(universe, *args)

    def Finalize(self):
        self.histwrapper.Write()

    def AddCut(self, func):
        self.cuts.append(func)

    

def PlotProcessorProliferater(plot, name_func):
    plots = []
    for k,v  in name_func.items():
        plots.append(plot.Clone(VariantPlotsNamingScheme(plot.histwrapper.name,k)))
        plots[-1].AddCut(v)
    return plots

# complecated function for interpret tags of a plot. Got to revisit after sometime and find out how to organize multiple tags.
def MakePlotProcessors(**kwargs):
    settings = TranslateSettings(kwargs["key"])
    #print settings
    plots = []
    tags = settings["tags"]
    settings.pop("tags")
    # tags = settings.get("tags", [])
    # settings.pop("tags", None)

    
    if "suffix" in kwargs["key"]:
        settings["name"]+=kwargs["key"]["suffix"]

    if "mc_only" in tags and not kwargs["mc"]:
        return plots
    else:
        plots.append( PlotProcessor(**settings))

    if "ignore_selection" not in tags:
        name_func = {}
        if "sideband" in tags:
            for region in kwargs["region"]:
                name_func[VariantPlotsNamingScheme(region)] = MakePlotProcessors.func.setdefault(region,partial(lambda universe,_: universe.classifier.side_band == _, _=region))

        tmp = []
        for plot in plots:
            tmp.extend(PlotProcessorProliferater(plot,name_func))
            plot.AddCut(MakePlotProcessors.func.setdefault("Signal",lambda universe: universe.classifier.side_band =="Signal"))
        plots.extend(tmp)

    if "signal_only" in tags:
        if kwargs["mc"]:
            for plot in plots:
                plot.AddCut(lambda universe: universe.classifier.is_true_signal)
                plot.weight_function = lambda event:event.GetWeight(False)
        else:
            print("cant't make signal only plots for data")
            return plots

    if "truth_class" in tags and kwargs["mc"]:
        name_func = {}
        for cate in list(TRUTH_CATEGORIES.keys())+["Other"]:
            if cate in EXTRA_OTHER:
                continue
            name_func[cate] = MakePlotProcessors.func.setdefault(cate,partial(lambda universe,_: universe.classifier.truth_class == _, _=cate))


        tmp = []
        for plot in plots:
            tmp.extend(PlotProcessorProliferater(plot,name_func))
        plots.extend(tmp)

    return plots

MakePlotProcessors.func = {}

def MakeSelectionPlotsByCut(**kwargs):
    settings = TranslateSettings(kwargs["key"])
    plots = []
    cuts = [ partial(lambda event, cut : event.classifier.reco_cuts_passed[cut] , cut = _ ) for _ in SAMPLE_CUTS["Signal"]]
    cuts.extend( [ partial(lambda event, cut : event.classifier.reco_cuts_passed[cut] , cut = "Reco{}".format(_) ) for _ in KINEMATICS_CUTS] )

    for i in range(len(cuts)):
        local_settings = settings.copy()
        local_settings["name"] = "{}_{}{}".format(settings["name"],i,"cut")
        local_settings["cuts"] = cuts[:i+1]
        plots.append(PlotProcessor(**local_settings))
        local_settings["name"] = "{}_{}{}".format(settings["name"],i,"cut_signal")
        local_settings["cuts"] = [lambda event:event.classifier.is_true_signal]+cuts[:i+1]
        plots.append(PlotProcessor(**local_settings))

    return plots


def MakeCutVariablePlot(cutname,value_getter=None):
    plots = []
    settings={
        "name": cutname,
        "title":";{};NEvents".format(cutname),
    }

    print (binning,value_getter,name)

    if value_getter is None:
        settings["value_getter"]=[CUTS[cutname].Values]
    else :
        settings["value_getter"]=[value_getter]

    settings["binning"]=[CUTS[cutname]._variable_range]

    settings["cuts"] = [lambda event: all(event.classifier.reco_cuts_passed[_] for _ in SAMPLE_CUTS["Signal"] if _ != cutname)]
    settings["cuts"].extend([lambda event: all(event.classifier.reco_cuts_passed["Reco{}".format(_)] for _ in KINEMATICS_CUTS)])

    plots.append(PlotProcessor(**settings))
    settings["cuts"] = settings["cuts"]+[lambda event: event.classifier.is_true_signal]
    settings["name"] = "{}_signal".format(cutname)
    plots.append(PlotProcessor(**settings))

    return plots

