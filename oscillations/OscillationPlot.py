import os
import sys
import ROOT
import PlotUtils

import numpy as np
import math
from array import array
import argparse

# Parse plotting-only args first, before AnalysisConfig sees sys.argv.
_plot_parser = argparse.ArgumentParser(add_help=False)
_plot_parser.add_argument("--plot-input-dir", default="chi2s")
_plot_parser.add_argument("--plot-mode", default="noFluxProfile")
_plot_parser.add_argument("--plot-no-fc", dest="plot_use_fc", action="store_false", default=True)
_plot_parser.add_argument("--critical-only", action="store_true", default=False)

_plot_args, _remaining_argv = _plot_parser.parse_known_args()
sys.argv = [sys.argv[0]] + _remaining_argv

from tools.PlotLibrary import HistHolder
from tools.StitchedHistogram import *
from tools.Fitters import *
from config.AnalysisConfig import AnalysisConfig
ccnueroot = os.environ.get('CCNUEROOT')

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FixedLocator
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patheffects as patheffects

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

str_to_latex = {
        "dm2"   : r"$\Delta m^{2}$ [ $eV^{2}$ ]",
        "ue4"   : r"$|U_{e4}|^{2}$",
        "umu4"  : r"$|U_{\mu4}|^{2}$",
        "utau4" : r"$|U_{{\tau 4}}|^2"
        }

str_to_axis = {
        "dm2"   : np.logspace(-1,2,100),
        "ue4"   : 0.15*np.logspace(-4,0,100),
        "umu4"  : 0.41*np.logspace(-5,0,100)
        }

str_to_index = {
        "dm2"   : 0,
        "ue4"   : 2,
        "umu4"  : 1
        }

str_to_indices = {
        "dm2"   : [20,27,32,39,44,51,61,71,92],
        "ue4"   : [0,53,70,78,79,80,82,84,92],
        "umu4"  : [0,50,58,62,65,67,69,72,87]
        }

str_to_zoom = {
        "dm2"   : 1e0,
        "ue4"   : 5e-3,
        "umu4"  : 5e-4
        }

fill_exclusions = False

hatch_styles = ["//",r"\\","+","oo","**"]

class ExperimentContour:
    def __init__(self,title,fname,fill,color="black",h_index=0):
        self.title = title

        self.color = color
        self.patch = Line2D([0], [0], color=self.color)
        self.style = "solid"
        self.fill = fill
        self.fname = fname
        self.hatch = hatch_styles[h_index]

        self.data = {}
        self.LoadData()

    def interpolate(self,inArr,size=100):
        ret = []
        for i in range(len(inArr)-1):
            x = [inArr[i],inArr[i+1]]
            new = np.linspace(x[0],x[1],size)
            ret.append(new)
        ret = np.concatenate(ret,axis=None)
        return(ret)

    def LoadData(self,header=1):
        if isinstance(self.fname,str):
            head = open(self.fname).readline().strip('\n')
            xaxis = head.split(',')[0]
            yaxis = head.split(',')[1]
            result = np.loadtxt(self.fname,delimiter=',',skiprows=header)
            for i,axis in enumerate([xaxis,yaxis]):
                if "theta" in axis:
                    if "mue" in axis:
                        data = result[:,i]/4
                    else:
                        data = np.array([np.roots([-4,4,-coeff])[1] for coeff in result[:,i]])
                else:
                    data = result[:,i]

                data = self.interpolate(data)
                self.data[axis] = data

        elif isinstance(self.fname,list):
            head = open(self.fname[0]).readline().strip('\n')
            xaxis = head.split(',')[0]
            yaxis = head.split(',')[1]
            self.data[xaxis] = []
            self.data[yaxis] = []
            for name in self.fname:
                result = np.loadtxt(name,delimiter=',',skiprows=header)
                for i,axis in enumerate([xaxis,yaxis]):
                    if "theta" in axis:
                        if "mue" in axis:
                            data = result[:,i]
                        else:
                            data = np.array([np.roots([-4,4,-coeff])[1] for coeff in result[:,i]])
                    else:
                        data = result[:,i]

                    data = self.interpolate(data)
                    self.data[axis].append(data)

        for axis in [xaxis,yaxis]:
            if "e4" in axis:
                new_axis = "ue4"
            elif "mu4" in axis:
                new_axis = "umu4"
            elif "mue" in axis:
                new_axis = "ue4_umu4"
            else:
                continue

            self.data[new_axis] = self.data[axis]
            del self.data[axis]

    def SetColor(self,color):
        self.color = color

    def SetLineStyle(self,style):
        self.style = style

    def SetPatch(self,graphic):
        graphic = graphic.lower()
        if graphic == "line":
            self.patch = Line2D([0], [0], color=self.color, linestyle=self.style, alpha=0.4)
        elif graphic == "patch":
            self.patch = Patch(color=self.color) 
        elif graphic == "hatch":
            self.patch = Patch(fill=False,hatch=self.hatch, edgecolor='black')
        else:
            raise ValueError("Graphic type not supported for legend")

    def GetPatch(self):
        return(self.patch)

    def GetTitle(self):
        return(self.title)

    def GetIntersects(self,data_ref,data_comp,axis,line,intName):
        if not isinstance(self.fname,list):
            intersections = []
            if line > data_comp[0]:
                intersections.append(data_ref[0])

            for i in range(1,len(data_ref)-1):
                if line > data_comp[i] and line < data_comp[i-1]:
                    intersections.append(data_ref[i])
                elif line < data_comp[i] and line > data_comp[i-1]:
                    intersections.append(data_ref[i])

            if line > data_comp[-1]:
                if "MINOS" in self.fname:
                    intersections.append(str_to_axis[intName][-1])
                else:
                    intersections.append(data_ref[-1])

            box1 = []
            box2 = []
            if len(intersections) > 0:
                for i in range(0,len(intersections),2):
                    point1 = intersections[i]
                    point2 = intersections[i+1]

                    box1.append([axis[0],axis[0],axis[-1],axis[-1]])
                    box2.append([point1,point2,point2,point1])

            return(box1,box2)
        else:
            intersections = []

            for i in range(1,len(data_ref)-1):
                if line > data_comp[i] and line < data_comp[i-1]:
                    intersections.append(data_ref[i])
                elif line < data_comp[i] and line > data_comp[i-1]:
                    intersections.append(data_ref[i])

            box1 = []
            box2 = []
            if len(intersections) > 1:
                for i in range(0,len(intersections),2):
                    point1 = intersections[i]
                    point2 = intersections[i+1]

                    box1.append([axis[0],axis[0],axis[-1],axis[-1]])
                    box2.append([point1,point2,point2,point1])

            return(box1,box2)

    def TranslateMixing(self,panel,line):
        if panel == 'dm2':
            intersect = 0
            for i in range(1,len(self.data['dm2'])):
                if self.data['dm2'][i] > line and self.data['dm2'][i-1] < line:
                    intersect = self.data['ue4_umu4'][i] # want largest intersect for limits, so don't break
            ue4_0 = intersect / (4*str_to_axis["umu4"][0])
            ue4_1 = intersect / (4*str_to_axis["umu4"][-1])
            umu4_0 = intersect / (4*str_to_axis["ue4"][0])
            umu4_1 = intersect / (4*str_to_axis["ue4"][-1])
            return(np.array([ue4_0,ue4_1]),np.array([str_to_axis["umu4"][0],str_to_axis["umu4"][-1]]))
        else:
            x = []
            y = []
            for i,dm2 in enumerate(self.data['dm2']):
                mixing = self.data['ue4_umu4'][i]
                x.append(mixing / (4*line))
                y.append(dm2)
            return(np.array(x),np.array(y))

    def PlotBox(self,axis,boxx,boxy):
        if self.fill:
            self.SetPatch('Patch')
            for i in range(len(boxx)):
                axis.fill(boxx[i],boxy[i],alpha=0.4,color=self.color)
        else:
            self.SetPatch('Hatch')
            for i in range(len(boxx)):
                axis.fill(boxx[i],boxy[i],fill=False,hatch=self.hatch,alpha=0.2,color='black')

    def Draw(self,axis,x,y):
        if self.fill:
            self.SetPatch("patch")
            axis.fill(x,y,color=self.color,alpha=0.4)
        else:
            self.SetPatch("line")
            axis.plot(x,y,color=self.color,linestyle=self.style,alpha=0.4)

    def Plot(self,axis,xaxis,yaxis,line,panel):
        axes = self.data.keys()
        t = self.data[list(self.data.keys())[0]] 
        
        if isinstance(t,list):
            for i in range(len(t)):
                if xaxis in axes and yaxis in axes:
                    x = self.data[xaxis][i]
                    y = self.data[yaxis][i]
                    self.Draw(axis,x,y)
                elif "ue4_umu4" in axes:
                    x,y = self.TranslateMixing(panel,line)
                    self.Draw(axis,x,y)
                else: # need to recreate an axis
                    if xaxis not in axes:
                        ydata = self.data[yaxis][i]
                        pdata = self.data[panel][i]
                        boxx,boxy = self.GetIntersects(ydata,pdata,str_to_axis[xaxis],line,yaxis)
                        self.PlotBox(axis,boxx,boxy)
                    elif yaxis not in axes:
                        xdata = self.data[xaxis][i]
                        pdata = self.data[panel][i]
                        boxx,boxy = self.GetIntersects(xdata,pdata,str_to_axis[yaxis],line,xaxis)
                        self.PlotBox(axis,boxy,boxx)
        else:
            if xaxis in axes and yaxis in axes:
                x = self.data[xaxis]
                y = self.data[yaxis]
                self.Draw(axis,x,y)
            elif "ue4_umu4" in axes:
                x,y = self.TranslateMixing(panel,line)
                self.Draw(axis,x,y)
            else: # need to recreate an axis
                if xaxis not in axes:
                    ydata = self.data[yaxis]
                    pdata = self.data[panel]
                    boxx,boxy = self.GetIntersects(ydata,pdata,str_to_axis[xaxis],line,yaxis)
                    self.PlotBox(axis,boxx,boxy)
                elif yaxis not in axes:
                    xdata = self.data[xaxis]
                    pdata = self.data[panel]
                    boxx,boxy = self.GetIntersects(xdata,pdata,str_to_axis[yaxis],line,xaxis)
                    self.PlotBox(axis,boxy,boxx)

class PanelPlot:
    def __init__(self,title,xaxis,yaxis,panel,name=""):
        self.title = title
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.panel = panel
        self.name = name

        self.x = str_to_axis[self.xaxis]
        self.y = str_to_axis[self.yaxis]
        self.p = str_to_axis[self.panel]

        self.zoomx = self.x[0]
        self.zoomy = self.y[0]

        self.indices = str_to_indices[panel]

        self.exclusion_results = []
        self.fit_results = []
        self.artists = []
        self.texts = []
        self.sens = None
        self.excl = None

    def SetName(self,name):
        self.name = name

    def Zoom(self):
        self.zoomx = str_to_zoom[self.xaxis]
        self.zoomy = str_to_zoom[self.yaxis]
        self.SetName(self.name[:-4] + "_zoomed.png")
        for txt in self.texts:
            txt.remove()

        textprops = dict(facecolor='white',edgecolor='white', alpha=0.8)
        for i,ax in enumerate(self.axes.flatten()):
            ax.set_xlim((self.zoomx,self.x[-1]))
            ax.set_ylim((self.zoomy,self.y[-1]))

            rx = self.x[7] / self.x[0]
            ry = self.y[7] / self.y[0]

            textx = rx * self.zoomx
            texty = ry * self.zoomy

            text = str_to_latex[self.panel]
            plabel = self.p[self.indices[i]]
            if self.panel == "dm2":
                text += "= {:.1f}".format(plabel)
                if i == 6:
                    texty=0.1
            elif self.panel == "ue4":
                text += "= {:.3f}".format(plabel)
            elif self.panel == "umu4":
                text += "= {:.4f}".format(plabel)
            self.texts.append(ax.text(textx,texty,s=text,bbox=textprops,c="black",fontweight="bold"))

        self.Save()
        self.texts = []

    def SetTitle(self,title):
        self.title = title

    def SetXaxis(self,axis):
        self.xaxis = axis
        self.x = str_to_axis[axis]

    def SetYaxis(self,axis):
        self.yaxis = axis
        self.y = str_to_axis[axis]

    def SetPanel(self,axis):
        self.panel = axis
        self.p = str_to_axis[axis]
        self.indices = str_to_indices[self.panel]

    def AddExclusions(self,exp):
        for e in exp:
            self.exclusion_results.append(e)

    def AddAlloweds(self,exp):
        for e in exp:
            self.fit_results.append(e)

    def CreateAxis(self):
        nrows = math.ceil(math.sqrt(len(self.indices)))
        ncols = nrows
        self.fig, self.axes = plt.subplots(nrows,ncols,sharex=True,sharey=True,figsize=(12,8),gridspec_kw={'wspace':0, 'hspace':0})
        art1 = self.fig.suptitle(self.title,y=0.97,size=20)
        art2 = self.fig.supxlabel(str_to_latex[self.xaxis],y=0.02,size=18)
        art3 = self.fig.supylabel(str_to_latex[self.yaxis],x=.07,size=18)

        self.artists.extend([art1,art2,art3])
        textprops = dict(facecolor='white',edgecolor='white', alpha=0.8)

        for i,ax in enumerate(self.axes.flatten()):
            ax.label_outer()
            ax.set_aspect("auto")
            ax.set_xscale("log")
            ax.set_yscale("log")
            if self.yaxis == "dm2":
                ax.yaxis.set_major_locator(FixedLocator([1e2,1e1,1e0]))
            elif self.yaxis == "umu4":
                ax.yaxis.set_major_locator(FixedLocator([1e-1,1e-2,1e-3,1e-4,1e-5]))
            elif self.yaxis == "ue4":
                ax.yaxis.set_major_locator(FixedLocator([1e-1,1e-2,1e-3,1e-4]))
            ax.set_xlim((self.x[0],self.x[-1]))
            ax.set_ylim((self.y[0],self.y[-1]))
                
            text = str_to_latex[self.panel]
            plabel = self.p[self.indices[i]]
            if self.panel == "dm2":
                text += "= {:.1f}".format(plabel)
            elif self.panel == "ue4":
                text += "= {:.3f}".format(plabel)
            elif self.panel == "umu4":
                text += "= {:.4f}".format(plabel)

            self.texts.append(ax.text(self.x[7],self.y[7],s=text,bbox=textprops,c="black",fontweight="bold"))

    def PlotFeldmanCousins(self, FC_excl, FC_sens, limits):
        contour_labels = [str(i)+'%' for i in limits]
        contour_colors = ['red','blue']

        X,Y = np.meshgrid(self.x,self.y)

        for i,ax in enumerate(self.axes.flatten()):
            self.sens = ax.contour(X,Y,FC_sens[i],levels=limits,colors=contour_colors,origin="lower",linestyles='dashed')
            if fill_exclusions:
                self.sens.set(path_effects=[patheffects.Stroke(linewidth=3, foreground='black'), patheffects.Normal()])
                self.excl = ax.contourf(X,Y,FC_excl[i],levels=[95,100],alpha=0.6,colors=contour_colors,origin="lower")
            else:
                self.excl = ax.contour(X,Y,FC_excl[i],levels=limits,colors=contour_colors,origin="lower")

    def PlotFeldmanCousinsLists(self, FC_excl, FC_sens, colors, limits):
        contour_labels = [str(i)+'%' for i in limits]

        X,Y = np.meshgrid(self.x,self.y)

        for i,ax in enumerate(self.axes.flatten()):
            for j in range(len(colors)):
                self.sens = ax.contour(X,Y,FC_sens[j][i],levels=limits,colors=colors[j],origin="lower",linestyles='dashed')
                self.excl = ax.contour(X,Y,FC_excl[j][i],levels=limits,colors=colors[j],origin="lower")
                #self.excl.set(path_effects=[patheffects.withTickedStroke()])

    def PlotExclusions(self):
        for i,ax in enumerate(self.axes.flatten()):
            panel = self.p[self.indices[i]]
            for exp in self.exclusion_results:
                exp.Plot(ax,self.xaxis,self.yaxis,panel,self.panel)

    def PlotAlloweds(self):
        for i,ax in enumerate(self.axes.flatten()):
            panel = self.p[self.indices[i]]
            for exp in self.fit_results:
                exp.Plot(ax,self.xaxis,self.yaxis,panel,self.panel)

    def PlotLegend(self,limits):
        handles = []
        contour_labels = [str(i)+"%" for i in limits]
        contour_colors = ['red','blue','green']

        for c in range(len(contour_labels)):
            if fill_exclusions:
                line = plt.Line2D([0,0], [0,0], color=contour_colors[c],linestyle='dashed',path_effects=[patheffects.Stroke(linewidth=3, foreground='black'), patheffects.Normal()])
            else:
                line = plt.Line2D([0,0], [0,0], color=contour_colors[c],linestyle='dashed')

            handles.append(line)
        for c in range(len(contour_labels)):
            if fill_exclusions:
                line = Patch(color=contour_colors[c])
            else:
                line = plt.Line2D([0,0], [0,0], color=contour_colors[c])
            handles.append(line)
            
        ph = [plt.plot([],marker="", ls="")[0]]*2
        contour_labels+=contour_labels

        half = len(limits)
        handles = ph[:1] + handles[:half] + ph[1:] + handles[half:]
        contour_labels = ["Sensitivity:"] + contour_labels[:half] + ["Exclusion:"] + contour_labels[half:]

        leg = self.fig.legend(handles,contour_labels,ncol=2,title='MINERvA Feldman Cousins',bbox_to_anchor=(1.2, 0.4),loc='upper right',fontsize=14)
        leg.get_title().set_fontsize(14)

        for vpack in leg._legend_handle_box.get_children():
            for hpack in vpack.get_children()[:1]:
                hpack.get_children()[0].set_width(0)
                
        excl_handles = [exp.GetPatch() for exp in self.exclusion_results]
        excl_handles.extend([exp.GetPatch() for exp in self.fit_results])

        excl_labels = [exp.GetTitle() for exp in self.exclusion_results]
        excl_labels.extend([exp.GetTitle() for exp in self.fit_results])
        leg2 = self.fig.legend(excl_handles,excl_labels,title="External Experiments",bbox_to_anchor=(1.2, 0.75),loc='upper right',fontsize=14)
        leg2.get_title().set_fontsize(14)

        self.artists.extend([leg,leg2])

    def PlotListLegend(self,limits,titles,colors):
        handles = []
        contour_labels = [title for title in titles]
        contour_colors = colors

        for c in range(len(contour_labels)):
            line = plt.Line2D([0,0], [0,0], color=contour_colors[c], linestyle='dashed')
            handles.append(line)
        for c in range(len(contour_labels)):
            line = plt.Line2D([0,0], [0,0], color=contour_colors[c])
            handles.append(line)
            
        ph = [plt.plot([],marker="", ls="")[0]]*2
        contour_labels+=contour_labels

        half = len(titles)
        handles = ph[:1] + handles[:half] + ph[1:] + handles[half:]
        contour_labels = ["95% Sensitivity:"] + contour_labels[:half] + ["95% Exclusion:"] + contour_labels[half:]

        leg = self.fig.legend(handles,contour_labels,ncol=2,title='MINERvA Feldman Cousins',bbox_to_anchor=(1.2, 0.4),loc='upper right',fontsize=12)
        leg.get_title().set_fontsize(14)

        for vpack in leg._legend_handle_box.get_children():
            for hpack in vpack.get_children()[:1]:
                hpack.get_children()[0].set_width(0)
                
        excl_handles = [exp.GetPatch() for exp in self.exclusion_results]
        excl_handles.extend([exp.GetPatch() for exp in self.fit_results])

        excl_labels = [exp.GetTitle() for exp in self.exclusion_results]
        excl_labels.extend([exp.GetTitle() for exp in self.fit_results])
        leg2 = self.fig.legend(excl_handles,excl_labels,title="External Experiments",bbox_to_anchor=(1.2, 0.75),loc='upper right',fontsize=14)
        leg2.get_title().set_fontsize(14)

        self.artists.extend([leg,leg2])

    def Save(self):
        print("saving figure {}".format(self.name))
        plt.savefig(self.name,bbox_extra_artists=self.artists,bbox_inches='tight')

    def MakePlot(self,FC_excl,FC_sens,limits):
        self.CreateAxis()
        self.PlotFeldmanCousins(FC_excl,FC_sens,limits)
        self.PlotExclusions()
        self.PlotAlloweds()
        self.PlotLegend(limits)
        self.Save()

    def MakeCompPlot(self,excls,sens,titles,limits):
        colors = ["red","blue"]
        self.CreateAxis()
        self.PlotFeldmanCousinsLists(excls,sens,colors,limits)
        self.PlotExclusions()
        self.PlotAlloweds()
        self.PlotListLegend(limits,titles,colors)
        self.Save()

    def Animate(self,FC_excl,FC_sens,limits,index):
        self.fig, self.axes = plt.subplots(figsize=(12,8),gridspec_kw={'wspace':0, 'hspace':0})
        self.axes.set_title(self.title,size=20)
        self.axes.set_xlabel(str_to_latex[self.xaxis],size=18)
        self.axes.set_ylabel(str_to_latex[self.yaxis],size=18)

        textprops = dict(facecolor='white',edgecolor='white', alpha=0.8)

        self.axes.label_outer()
        self.axes.set_aspect("auto")
        self.axes.set_xscale("log")
        self.axes.set_yscale("log")
            
        self.axes.set_xlim((self.x[0],self.x[-1]))
        self.axes.set_ylim((self.y[0],self.y[-1]))

        self.axes.text(self.x[7],self.y[7],s=str_to_latex[self.panel]+ "= {:.5f}".format(self.p[index]),bbox=textprops,c="black",fontweight="bold")

        contour_labels = [str(i)+'%' for i in limits]
        contour_colors = ['red','blue']

        X,Y = np.meshgrid(self.x,self.y)

        self.axes.contour(X,Y,FC_sens,levels=limits,colors=contour_colors,origin="lower",linestyles='dashed')
        self.axes.contour(X,Y,FC_excl,levels=limits,colors=contour_colors,origin="lower")

        panel = self.p[index]
        for exp in self.exclusion_results:
            exp.Plot(self.axes,self.xaxis,self.yaxis,panel,self.panel)
        panel = self.p[index]
        for exp in self.fit_results:
            exp.Plot(self.axes,self.xaxis,self.yaxis,panel,self.panel)
        self.PlotLegend(limits)
        self.Save()
        plt.close()

    def PlotRawSurfacePanel(
        self,
        surface3d,
        zlabel=r"$\chi^2$",
        contour_levels=None,
        contour_color="white",
        contour_label_fmt="%.2f",
    ):
        """
        Plot raw 2D surface values on the same 3x3 panel layout used by FC plots.

        surface3d is indexed as:
          [dm2, umu4, ue4]

        self.panel determines which axis is fixed:
          panel="dm2"  -> x=ue4,  y=umu4, fixed dm2
          panel="umu4" -> x=ue4,  y=dm2,  fixed umu4
          panel="ue4"  -> x=umu4, y=dm2,  fixed ue4
        """

        self.CreateAxis()

        X, Y = np.meshgrid(self.x, self.y)

        z_slices = []

        for idx in self.indices:
            if self.panel == "dm2":
                z = surface3d[idx, :, :]
            elif self.panel == "umu4":
                z = surface3d[:, idx, :]
            elif self.panel == "ue4":
                z = surface3d[:, :, idx]
            else:
                raise ValueError("Unknown panel: {}".format(self.panel))

            z_slices.append(z)

        zmin = min(np.nanmin(z) for z in z_slices)
        zmax = max(np.nanmax(z) for z in z_slices)

        pcm = None

        for i, ax in enumerate(self.axes.flatten()):
            if i >= len(z_slices):
                ax.axis("off")
                continue

            z = z_slices[i]

            pcm = ax.pcolormesh(
                X,
                Y,
                z,
                shading="auto",
                vmin=zmin,
                vmax=zmax,
            )

            if contour_levels is not None:
                good_levels = [
                    lev for lev in contour_levels
                    if lev > np.nanmin(z) and lev < np.nanmax(z)
                ]

                if len(good_levels) > 0:
                    cs = ax.contour(
                        X,
                        Y,
                        z,
                        levels=good_levels,
                        colors=contour_color,
                        linewidths=1.1,
                    )
                    ax.clabel(cs, fmt=contour_label_fmt, fontsize=7)

        cbar = self.fig.colorbar(
            pcm,
            ax=self.axes.ravel().tolist(),
            shrink=0.90,
            pad=0.02,
        )
        cbar.set_label(zlabel, fontsize=14)

        self.artists.append(cbar.ax)

        self.Save()

def PlotAsimovRawChi2Panels(asimov_chi2s, mode, title):
    """
    Make original-style 3x3 panel plots for raw Asimov chi2.
    Uses the same str_to_indices slices as the FC contour plots.
    """

    outbase = "/exp/minerva/data/users/qvuong/surfaces/plots"

    pplot = PanelPlot(title + "\nAsimov raw $\\chi^2$", "ue4", "umu4", "dm2")

    # 1. fixed dm2: ue4 vs umu4
    pplot.SetXaxis("ue4")
    pplot.SetYaxis("umu4")
    pplot.SetPanel("dm2")
    pplot.SetName("{}/Asimov_raw_chi2_ue4_vs_umu4_{}.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(asimov_chi2s, zlabel=r"Asimov raw $\chi^2$")
    pplot.Zoom()

    # 2. fixed umu4: ue4 vs dm2
    pplot = PanelPlot(title + "\nAsimov raw $\\chi^2$", "ue4", "dm2", "umu4")
    pplot.SetName("{}/Asimov_raw_chi2_ue4_vs_dm2_{}.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(asimov_chi2s, zlabel=r"Asimov raw $\chi^2$")
    pplot.Zoom()

    # 3. fixed ue4: umu4 vs dm2
    pplot = PanelPlot(title + "\nAsimov raw $\\chi^2$", "umu4", "dm2", "ue4")
    pplot.SetName("{}/Asimov_raw_chi2_umu4_vs_dm2_{}.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(asimov_chi2s, zlabel=r"Asimov raw $\chi^2$")
    pplot.Zoom()

def PlotAsimovDChi2Panels(asimov_dchi2s, mode, title, fc_levels=None):
    """
    Make original-style 3x3 panel plots for Asimov delta-chi2.
    Uses the same str_to_indices slices as the FC contour plots.
    """

    outbase = "/exp/minerva/data/users/qvuong/surfaces/plots"

    # 1. fixed dm2: ue4 vs umu4
    pplot = PanelPlot(title + "\nAsimov $\\Delta\\chi^2$", "ue4", "umu4", "dm2")
    pplot.SetName("{}/Asimov_dchi2_ue4_vs_umu4_{}_95CL.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(
        asimov_dchi2s,
        zlabel=r"Asimov $\Delta\chi^2$",
        contour_levels=fc_levels,
    )
    pplot.Zoom()

    # 2. fixed umu4: ue4 vs dm2
    pplot = PanelPlot(title + "\nAsimov $\\Delta\\chi^2$", "ue4", "dm2", "umu4")
    pplot.SetName("{}/Asimov_dchi2_ue4_vs_dm2_{}_95CL.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(
        asimov_dchi2s,
        zlabel=r"Asimov $\Delta\chi^2$",
        contour_levels=fc_levels,
    )
    pplot.Zoom()

    # 3. fixed ue4: umu4 vs dm2
    pplot = PanelPlot(title + "\nAsimov $\\Delta\\chi^2$", "umu4", "dm2", "ue4")
    pplot.SetName("{}/Asimov_dchi2_umu4_vs_dm2_{}_95CL.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(
        asimov_dchi2s,
        zlabel=r"Asimov $\Delta\chi^2$",
        contour_levels=fc_levels,
    )
    pplot.Zoom()

def PlotDataDChi2Panels(data_dchi2s, mode, title, fc_levels=None):
    """
    Make original-style 3x3 panel plots for data delta-chi2.
    Uses the same str_to_indices slices as the FC contour plots.
    Draws FC critical contour levels if fc_levels is provided.
    """

    outbase = "/exp/minerva/data/users/qvuong/surfaces/plots"

    # 1. fixed dm2: ue4 vs umu4
    pplot = PanelPlot(title + "\nData $\\Delta\\chi^2$", "ue4", "umu4", "dm2")
    pplot.SetName("{}/Data_dchi2_ue4_vs_umu4_{}_95CL.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(
        data_dchi2s,
        zlabel=r"Data $\Delta\chi^2$",
        contour_levels=fc_levels,
    )
    pplot.Zoom()

    # 2. fixed umu4: ue4 vs dm2
    pplot = PanelPlot(title + "\nData $\\Delta\\chi^2$", "ue4", "dm2", "umu4")
    pplot.SetName("{}/Data_dchi2_ue4_vs_dm2_{}_95CL.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(
        data_dchi2s,
        zlabel=r"Data $\Delta\chi^2$",
        contour_levels=fc_levels,
    )
    pplot.Zoom()

    # 3. fixed ue4: umu4 vs dm2
    pplot = PanelPlot(title + "\nData $\\Delta\\chi^2$", "umu4", "dm2", "ue4")
    pplot.SetName("{}/Data_dchi2_umu4_vs_dm2_{}_95CL.png".format(outbase, mode))
    pplot.PlotRawSurfacePanel(
        data_dchi2s,
        zlabel=r"Data $\Delta\chi^2$",
        contour_levels=fc_levels,
    )
    pplot.Zoom()

def GetFCSlices(dchi2s,achi2s,results,panel):
    indices = str_to_indices[panel]
    slc = str_to_index[panel]
    sens_list = []
    excl_list = []
    for i in indices:
        FC_sens = achi2s.take(indices=i,axis=slc)
        FC_excl = dchi2s.take(indices=i,axis=slc)
        for iy, ix in np.ndindex(FC_sens.shape):
            FC_sens[iy,ix] = 100*(results[FC_sens[iy,ix] > results].shape[0]/results.shape[0])
            FC_excl[iy,ix] = 100*(results[FC_excl[iy,ix] > results].shape[0]/results.shape[0])

        sens_list.append(FC_sens)
        excl_list.append(FC_excl)
    return(sens_list,excl_list)

def GetFCSlice(dchi2s,achi2s,results,index,slc):
    FC_sens = achi2s.take(indices=index,axis=slc)
    FC_excl = dchi2s.take(indices=index,axis=slc)
    for iy, ix in np.ndindex(FC_sens.shape):
        FC_sens[iy,ix] = 100*(results[FC_sens[iy,ix] > results].shape[0]/results.shape[0])
        FC_excl[iy,ix] = 100*(results[FC_excl[iy,ix] > results].shape[0]/results.shape[0])

    return(FC_sens,FC_excl)


def GetContinuousBestFitChi2(mode):
    """
    Continuous best-fit chi2 values from the nominal no-throw fits.

    These are the values analogous to Ryan's hard-coded best_fit values,
    but for our updated configurations.
    """

    best_fit_chi2 = {
        # prodNueel_noRatio, --profile-flux
        "prodNueel_noRatio_profiledFlux": 35.139265,

        # prodNueel, --profile-flux --exclude ratio
        # The merged surface mode name does not include excludeRatio unless
        # you added it to the gridSurfaces/makeSurface output naming.
        "prodNueel_profiledFlux": 38.769660,

        # Optional alias if you later encode exclude ratio in the mode name
        "prodNueel_profiledFlux_excluderatio": 38.769660,
        "prodNueel_profiledFlux_excludeRatio": 38.769660,
    }

    return best_fit_chi2.get(mode, None)


def LoadMergedSurfaces(input_dir, mode, use_fc=True):
    """
    Load merged outputs from merge_2D_quinn.py.

    Expected files:
      data_chi2s_<mode>.npy
      asimov_chi2s_<mode>.npy
      data_penalties_<mode>.npy
      asimov_penalties_<mode>.npy
      delta_m_values_<mode>.npy

    Optional FC file:
      asimov_deltachi2s.npy
    """

    data_file = os.path.join(input_dir, "data_chi2s_{}.npy".format(mode))
    asimov_file = os.path.join(input_dir, "asimov_chi2s_{}.npy".format(mode))
    data_penalty_file = os.path.join(input_dir, "data_penalties_{}.npy".format(mode))
    asimov_penalty_file = os.path.join(input_dir, "asimov_penalties_{}.npy".format(mode))
    dm_file = os.path.join(input_dir, "delta_m_values_{}.npy".format(mode))

    for f in [data_file, asimov_file, data_penalty_file, asimov_penalty_file, dm_file]:
        if not os.path.isfile(f):
            raise IOError("Missing required merged file: {}".format(f))

    data_chi2s = np.load(data_file)
    asimov_chi2s = np.load(asimov_file)
    data_penalties = np.load(data_penalty_file)
    asimov_penalties = np.load(asimov_penalty_file)
    delta_m_values = np.load(dm_file)

    print("Loaded mode:", mode)
    print("  data_chi2s shape       =", data_chi2s.shape)
    print("  asimov_chi2s shape     =", asimov_chi2s.shape)
    print("  data_penalties max     =", np.nanmax(data_penalties))
    print("  asimov_penalties max   =", np.nanmax(asimov_penalties))
    print("  delta_m_values shape   =", delta_m_values.shape)

    # Data delta-chi2 relative to the scanned grid minimum.
    data_grid_min = np.nanmin(data_chi2s)
    data_grid_idx = np.unravel_index(np.nanargmin(data_chi2s), data_chi2s.shape)
    data_dchi2s = data_chi2s - data_grid_min

    # Asimov delta-chi2 relative to the scanned grid minimum.
    asimov_grid_min = np.nanmin(asimov_chi2s)
    asimov_grid_idx = np.unravel_index(np.nanargmin(asimov_chi2s), asimov_chi2s.shape)
    asimov_dchi2s = asimov_chi2s - asimov_grid_min

    print("  data grid min chi2    =", data_grid_min)
    print("  data grid min index   =", data_grid_idx)
    print("  data grid min dm2     =", delta_m_values[data_grid_idx[0]])
    print("  data grid min umu4    =", str_to_axis["umu4"][data_grid_idx[1]])
    print("  data grid min ue4     =", str_to_axis["ue4"][data_grid_idx[2]])
    print("  data dchi2 min/max    =", np.nanmin(data_dchi2s), np.nanmax(data_dchi2s))

    print("  asimov grid min chi2  =", asimov_grid_min)
    print("  asimov grid min index =", asimov_grid_idx)
    print("  asimov grid min dm2   =", delta_m_values[asimov_grid_idx[0]])
    print("  asimov grid min umu4  =", str_to_axis["umu4"][asimov_grid_idx[1]])
    print("  asimov grid min ue4   =", str_to_axis["ue4"][asimov_grid_idx[2]])
    print("  asimov dchi2 min/max  =", np.nanmin(asimov_dchi2s), np.nanmax(asimov_dchi2s))

    results = None
    if use_fc:
        fc_file = os.path.join(input_dir, "asimov_deltachi2s_{}.npy".format(mode))
        if not os.path.isfile(fc_file):
            raise IOError(
                "Missing FC toy file: {}. "
                "Run MergeAsimovs first, or use --plot-no-fc for fixed chi2/surface plots.".format(fc_file)
            )
        results = np.load(fc_file)
        print("  FC results shape       =", results.shape)
        print("  FC results min/max     =", np.nanmin(results), np.nanmax(results))

    return data_dchi2s, asimov_dchi2s, results, delta_m_values


def PlotOneSurface(surface3d, panel="dm2", index=0, name="surface.png",
                   title="Chi2 surface", levels=None):
    """
    Plot one 2D slice from a merged 3D chi2 surface.

    surface3d is indexed as:
      [dm2, umu4, ue4]
    """

    # if levels is None:
    #     levels = [1, 2.71, 3.84, 6.63, 9.0, 11.83]

    if panel == "dm2":
        # fixed dm2 slice: y=umu4, x=ue4
        z = surface3d[index, :, :]
        xaxis = "ue4"
        yaxis = "umu4"
        fixed_value = str_to_axis["dm2"][index]

    elif panel == "umu4":
        # fixed umu4 slice: y=dm2, x=ue4
        z = surface3d[:, index, :]
        xaxis = "ue4"
        yaxis = "dm2"
        fixed_value = str_to_axis["umu4"][index]

    elif panel == "ue4":
        # fixed ue4 slice: y=dm2, x=umu4
        z = surface3d[:, :, index]
        xaxis = "umu4"
        yaxis = "dm2"
        fixed_value = str_to_axis["ue4"][index]

    else:
        raise ValueError("panel must be one of: dm2, umu4, ue4")

    x = str_to_axis[xaxis]
    y = str_to_axis[yaxis]
    X, Y = np.meshgrid(x, y)

    print("Plotting:", name)
    print("  panel =", panel)
    print("  index =", index)
    print("  fixed value =", fixed_value)
    print("  z shape =", z.shape)
    print("  z min/max =", np.nanmin(z), np.nanmax(z))

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(str_to_latex[xaxis], fontsize=14)
    ax.set_ylabel(str_to_latex[yaxis], fontsize=14)

    label = "{} = {:.4g}".format(panel, fixed_value)
    ax.set_title("{}\n{}".format(title, label), fontsize=14)

    # Filled color map for the surface
    pcm = ax.pcolormesh(X, Y, z, shading="auto")
    cbar = fig.colorbar(pcm, ax=ax)
    if levels is None:
        cbar.set_label(r"$\chi^2$", fontsize=14)
    else:
        cbar.set_label(r"$\Delta\chi^2$", fontsize=14)

    # Contour lines
    if levels is not None:
        good_levels = [lev for lev in levels if lev > np.nanmin(z) and lev < np.nanmax(z)]
        if len(good_levels) > 0:
            print("  Drawing contour levels:", good_levels)
            cs = ax.contour(X, Y, z, levels=good_levels, colors="white", linewidths=1.2)
            ax.clabel(cs, fmt="%.2f", fontsize=8)
        else:
            print("  No contour levels inside z range.")

    os.makedirs(os.path.dirname(name), exist_ok=True)
    fig.tight_layout()
    fig.savefig(name)
    plt.close(fig)

def PlotFCContours(data_chi2s, asimov_chi2s, results, mode, title):
    """
    Make the original-style FC contour panel plots.

    data_chi2s should already be data chi2 - best_fit.
    asimov_chi2s is the null/Asimov chi2 surface.
    results is the FC toy dchi2 distribution.
    """

    limits = [95]

    pplot = PanelPlot(title, "ue4", "dm2", "umu4")

    stereo = ExperimentContour("STEREO 95% Excl.", "exp_results/stereo_2Dexcl.csv", False, "black", h_index=0)
    minos = ExperimentContour("MINOS 90% Excl.", "exp_results/MINOS.csv", False, "black", h_index=1)
    prospect = ExperimentContour("PROSPECT 95% Excl.", "exp_results/prospect.csv", False, "blue", h_index=2)
    katrin = ExperimentContour("KATRIN 95% Excl.", "exp_results/katrin_noMass.csv", False, "purple", h_index=3)
    microboone = ExperimentContour("MicroBooNE 95% Excl.", "exp_results/microboone.csv", False, "teal", h_index=4)
    neutrino4 = ExperimentContour(
        "Neutrino-4 $2\sigma$ Conf.",
        ["exp_results/n4_c1.csv", "exp_results/n4_c2.csv", "exp_results/n4_c3.csv", "exp_results/n4_c4.csv"],
        True,
        "pink",
        h_index=2
    )
    raa = ExperimentContour("RAA 90% Allowed", "exp_results/RAA.csv", True, "gray", h_index=3)

    neutrino4.SetPatch("Patch")
    raa.SetPatch("Patch")

    pplot.AddExclusions([stereo, minos, prospect, katrin, microboone])

    if not fill_exclusions:
        pplot.AddAlloweds([neutrino4, raa])

    # os.makedirs("plots", exist_ok=True)

    pplot.SetXaxis("ue4")
    pplot.SetYaxis("umu4")
    pplot.SetPanel("dm2")
    sens_list, excl_list = GetFCSlices(data_chi2s, asimov_chi2s, results, "dm2")
    pplot.SetName("/exp/minerva/data/users/qvuong/surfaces/plots/FC_ue4_vs_umu4_{}.png".format(mode))
    pplot.MakePlot(excl_list, sens_list, limits)
    pplot.Zoom()

    pplot.SetXaxis("ue4")
    pplot.SetYaxis("dm2")
    pplot.SetPanel("umu4")
    sens_list, excl_list = GetFCSlices(data_chi2s, asimov_chi2s, results, "umu4")
    pplot.SetName("/exp/minerva/data/users/qvuong/surfaces/plots/FC_ue4_vs_dm2_{}.png".format(mode))
    pplot.MakePlot(excl_list, sens_list, limits)
    pplot.Zoom()

    pplot.SetXaxis("umu4")
    pplot.SetYaxis("dm2")
    pplot.SetPanel("ue4")
    sens_list, excl_list = GetFCSlices(data_chi2s, asimov_chi2s, results, "ue4")
    pplot.SetName("/exp/minerva/data/users/qvuong/surfaces/plots/FC_umu4_vs_dm2_{}.png".format(mode))
    pplot.MakePlot(excl_list, sens_list, limits)
    pplot.Zoom()

def PlotFCCriticalDistribution(input_dir, mode):
    fc_file = os.path.join(input_dir, "asimov_deltachi2s_{}.npy".format(mode))

    if not os.path.isfile(fc_file):
        raise IOError("Missing FC toy file: {}".format(fc_file))

    results = np.load(fc_file)
    results = np.asarray(results).ravel()
    results = results[np.isfinite(results)]

    print("\nLoaded FC toy dchi2 file:")
    print("  file   =", fc_file)
    print("  ntoys  =", results.size)
    print("  min    =", np.min(results))
    print("  max    =", np.max(results))
    print("  mean   =", np.mean(results))
    print("  median =", np.median(results))

    cls = [68, 90, 95, 99]
    crit = {cl: np.percentile(results, cl) for cl in cls}

    print("\nFC critical dchi2 levels:")
    for cl in cls:
        print("  {}% = {}".format(cl, crit[cl]))

    outdir = os.path.join(input_dir, "plots")
    os.makedirs(outdir, exist_ok=True)

    outname = os.path.join(
        outdir,
        "FC_dchi2_distribution_{}.png".format(mode)
    )

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.hist(
        results,
        bins=60,
        histtype="stepfilled",
        alpha=0.45,
        density=True,
    )

    ax.set_xlim(0, 20)

    ymax = ax.get_ylim()[1]

    line_styles = {
        68: "--",
        90: "-.",
        95: "-",
        99: ":",
    }

    for cl in cls:
        x = crit[cl]
        ax.axvline(x, linestyle=line_styles[cl], linewidth=2)
        ax.text(
            x,
            0.92 * ymax,
            "{}% = {:.3g}".format(cl, x),
            rotation=90,
            va="top",
            ha="right",
            fontsize=10,
        )

    # plot_max = 300.0
    # plot_results = results[results <= plot_max]
    # n_overflow = results.size - plot_results.size

    # fig, ax = plt.subplots(figsize=(8, 6))

    # ax.hist(
    #     plot_results,
    #     bins=60,
    #     range=(0, plot_max),
    #     histtype="stepfilled",
    #     alpha=0.45,
    #     density=False,
    # )

    # ax.set_xlim(-20, plot_max)

    # ymax = ax.get_ylim()[1]

    # line_styles = {
    #     68: "--",
    #     90: "-.",
    #     95: "-",
    #     99: ":",
    # }

    # for cl in cls:
    #     x = crit[cl]

    #     if x <= plot_max:
    #         ax.axvline(x, linestyle=line_styles[cl], linewidth=2)
    #         ax.text(
    #             x,
    #             0.92 * ymax,
    #             "{}% = {:.3g}".format(cl, x),
    #             rotation=90,
    #             va="top",
    #             ha="right",
    #             fontsize=10,
    #         )
    #     else:
    #         ax.text(
    #             0.98 * plot_max,
    #             0.92 * ymax,
    #             "{}% = {:.3g} > {}".format(cl, x, plot_max),
    #             rotation=90,
    #             va="top",
    #             ha="right",
    #             fontsize=10,
    #         )

    # if n_overflow > 0:
    #     ax.text(
    #         0.98,
    #         0.85,
    #         "Overflow > {}: {} toys".format(plot_max, n_overflow),
    #         transform=ax.transAxes,
    #         ha="right",
    #         va="top",
    #         fontsize=10,
    #     )

    ax.set_xlabel(r"Toy $\Delta\chi^2$")
    # ax.set_ylabel("Number of toys")
    ax.set_ylabel("Probability density")
    ax.set_title("FC toy distribution\n{}".format(mode))
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(outname)
    plt.close(fig)

    print("\nSaved FC dchi2 distribution plot:")
    print("  {}".format(outname))

    return results, crit




if __name__ == "__main__":

    if _plot_args.critical_only:
        PlotFCCriticalDistribution(
            _plot_args.plot_input_dir,
            _plot_args.plot_mode,
        )
        sys.exit(0)

    title = "MINERvA Sterile Neutrino Search"
    if _plot_args.plot_mode == "noFluxProfile":
        title += "\nNo flux profiling"
    elif _plot_args.plot_mode == "profiledFlux":
        title += "\nFlux profiled"

    data_dchi2s, asimov_dchi2s, results, delta_m_values = LoadMergedSurfaces(
        _plot_args.plot_input_dir,
        _plot_args.plot_mode,
        use_fc=_plot_args.plot_use_fc
    )

    fc_levels = None
    fc1sigma = np.percentile(results, 68.27)
    fc2sigma = np.percentile(results, 95.45)
    fc3sigma = np.percentile(results, 99.73)

    # fc_levels = [
    #     fc1sigma,
    #     fc2sigma,
    #     fc3sigma,
    # ]

    fc_levels = [
        # fc1sigma,
        fc2sigma,
        # fc3sigma,
    ]

    print("\nFC sigma-equivalent critical dchi2 levels:")
    print("  1 sigma (68.27%) =", fc1sigma)
    print("  2 sigma (95.45%) =", fc2sigma)
    print("  3 sigma (99.73%) =", fc3sigma)

    # --------------------------------------------------
    # 1. Raw surface plots
    # --------------------------------------------------
    outdir = "/exp/minerva/data/users/qvuong/surfaces/plots/p8/surfaces_{}".format(_plot_args.plot_mode)
    os.makedirs(outdir, exist_ok=True)

    # # --------------------------------------------------
    # # 1a. Fixed dm2 slices: ue4 vs umu4
    # # surface index: [dm2, umu4, ue4]
    # # --------------------------------------------------
    # for i in range(data_dchi2s.shape[0]):
    #     dm2_val = delta_m_values[i]

    #     PlotOneSurface(
    #         data_dchi2s,
    #         panel="dm2",
    #         index=i,
    #         name="{}/data_surface_ue4_vs_umu4_dm2idx{:03d}.png".format(outdir, i),
    #         title="Data surface, {}: dm2 = {:.4g}".format(_plot_args.plot_mode, dm2_val),
    #         levels=fc_levels
    #     )

    #     PlotOneSurface(
    #         asimov_dchi2s,
    #         panel="dm2",
    #         index=i,
    #         name="{}/asimov_surface_ue4_vs_umu4_dm2idx{:03d}.png".format(outdir, i),
    #         title="Asimov raw χ² surface, {}: dm2 = {:.4g}".format(_plot_args.plot_mode, dm2_val),
    #         levels=None
    #     )

    # # --------------------------------------------------
    # # 1b. Fixed umu4 slices: ue4 vs dm2
    # # surface index: [dm2, umu4, ue4]
    # # --------------------------------------------------
    # for i in range(data_dchi2s.shape[1]):
    #     umu4_val = str_to_axis["umu4"][i]

    #     PlotOneSurface(
    #         data_dchi2s,
    #         panel="umu4",
    #         index=i,
    #         name="{}/data_surface_ue4_vs_dm2_umu4idx{:03d}.png".format(outdir, i),
    #         title="Data surface, {}: umu4 = {:.4g}".format(_plot_args.plot_mode, umu4_val),
    #         levels=fc_levels
    #     )

    #     PlotOneSurface(
    #         asimov_dchi2s,
    #         panel="umu4",
    #         index=i,
    #         name="{}/asimov_surface_ue4_vs_dm2_umu4idx{:03d}.png".format(outdir, i),
    #         title="Asimov raw χ² surface, {}: umu4 = {:.4g}".format(_plot_args.plot_mode, umu4_val),
    #         levels=None
    #     )

    # # --------------------------------------------------
    # # 1c. Fixed ue4 slices: umu4 vs dm2
    # # surface index: [dm2, umu4, ue4]
    # # --------------------------------------------------
    # for i in range(data_dchi2s.shape[2]):
    #     ue4_val = str_to_axis["ue4"][i]

    #     PlotOneSurface(
    #         data_dchi2s,
    #         panel="ue4",
    #         index=i,
    #         name="{}/data_surface_umu4_vs_dm2_ue4idx{:03d}.png".format(outdir, i),
    #         title="Data surface, {}: ue4 = {:.4g}".format(_plot_args.plot_mode, ue4_val),
    #         levels=fc_levels
    #     )

    #     PlotOneSurface(
    #         asimov_dchi2s,
    #         panel="ue4",
    #         index=i,
    #         name="{}/asimov_surface_umu4_vs_dm2_ue4idx{:03d}.png".format(outdir, i),
    #         title="Asimov raw χ² surface, {}: ue4 = {:.4g}".format(_plot_args.plot_mode, ue4_val),
    #         levels=None
    #     )

    # --------------------------------------------------
    # 2. FC contour plots
    # --------------------------------------------------
    if _plot_args.plot_use_fc:
        PlotFCContours(
            data_dchi2s,
            asimov_dchi2s,
            results,
            _plot_args.plot_mode,
            title
        )

    PlotAsimovDChi2Panels(
        asimov_dchi2s,
        _plot_args.plot_mode,
        title,
        fc_levels=fc_levels,
    )

    PlotDataDChi2Panels(
        data_dchi2s,
        _plot_args.plot_mode,
        title,
        fc_levels=fc_levels,
    )
