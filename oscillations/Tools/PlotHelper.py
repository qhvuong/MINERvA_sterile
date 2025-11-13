import ROOT
import ctypes

def RatioAxis(hist,MNVPLOTTER):
    #MNVPLOTTER.ApplyAxisStyle(hist)

    hist.SetTitle("")
    hist.GetYaxis().SetNdivisions(304) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    hist.GetYaxis().SetTitleOffset(2.5)
    hist.GetYaxis().SetTitleFont(43)
    hist.GetYaxis().SetTitleSize(26)
    hist.GetYaxis().SetLabelFont(43)
    hist.GetYaxis().SetLabelSize(34)
    hist.GetYaxis().CenterTitle()
    hist.GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.
    hist.GetXaxis().SetTitleOffset(2.5)
    hist.GetXaxis().SetTitleFont(43)
    hist.GetXaxis().SetTitleSize(26)
    hist.GetXaxis().SetLabelFont(43)
    hist.GetXaxis().SetLabelSize(34)
    hist.GetXaxis().CenterTitle()

    #hist.GetYaxis().SetLabelSize(.13)
    #hist.GetYaxis().SetTitleSize(0.1)
    #hist.GetYaxis().SetTitleOffset(0.6)
    #hist.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    #hist.GetXaxis().SetTitleSize(0.16)
    #hist.GetXaxis().SetTitleOffset(0.9)
    #hist.GetXaxis().SetLabelSize(.15)

def CanvasPartition(C, Nx, Ny, lMargin, rMargin, bMargin, tMargin):
    vSpacing = 0.0
    vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny

    hSpacing = 0.0
    hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx

    vposd, vposu, vmard, vmaru, vfactor = 0, 0, 0, 0, 0
    hposl, hposr, hmarl, hmarr, hfactor = 0, 0, 0, 0, 0

    for i in range(Nx):
        if (i == 0):
            hposl = 0.0
            hposr = lMargin + hStep
            hfactor = hposr - hposl
            hmarl = lMargin / hfactor
            hmarr = 0.0
        elif (i == Nx - 1):
            hposl = hposr + hSpacing
            hposr = hposl + hStep + rMargin
            hfactor = hposr - hposl
            hmarl = 0.0
            hmarr = rMargin / (hposr - hposl)
        else:
            hposl = hposr + hSpacing
            hposr = hposl + hStep
            hfactor = hposr - hposl
            hmarl = 0.0
            hmarr = 0.0

        for j in range(Ny):
            if (j == 0):
                vposd = 0.0
                vposu = bMargin + vStep
                vfactor = vposu - vposd
                vmard = bMargin / vfactor
                vmaru = 0.0
            elif (j == Ny - 1):
                vposd = vposu + vSpacing
                vposu = vposd + vStep + tMargin
                vfactor = vposu - vposd
                vmard = 0.0
                vmaru = tMargin / (vposu - vposd)
            else:
                vposd = vposu + vSpacing
                vposu = vposd + vStep
                vfactor = vposu - vposd
                vmard = 0.0
                vmaru = 0.0

            C.cd(0)

            name = "pad_{}_{}".format(i, j)
            pad = ROOT.TPad(name, "", hposl, vposd, hposr, vposu)
            pad.SetLeftMargin(hmarl)
            pad.SetRightMargin(hmarr)
            pad.SetBottomMargin(vmard)
            pad.SetTopMargin(vmaru)

            pad.SetFrameBorderMode(0)
            pad.SetBorderMode(0)
            pad.SetBorderSize(0)

            pad.Draw()
 
def XtoPad(x):
   xl = ctypes.c_double(0)
   yl = ctypes.c_double(0)
   xu = ctypes.c_double(0)
   yu = ctypes.c_double(0)
   ROOT.gPad.GetPadPar(xl, yl, xu, yu)
   pw = ctypes.c_double(xu.value - xl.value)
   lm = ROOT.gPad.GetLeftMargin()
   rm = ROOT.gPad.GetRightMargin()
   fw = pw.value - pw.value * lm - pw.value * rm
   return (x * fw + pw.value * lm) / pw.value
 
def YtoPad(y):
   xl = ctypes.c_double(0)
   yl = ctypes.c_double(0)
   xu = ctypes.c_double(0)
   yu = ctypes.c_double(0)
   ROOT.gPad.GetPadPar(xl, yl, xu, yu)
   ph = ctypes.c_double(yu.value - yl.value)
   tm = ROOT.gPad.GetTopMargin()
   bm = ROOT.gPad.GetBottomMargin()
   fh = ph.value - ph.value * bm - ph.value * tm
   return (y * fh + bm * ph.value) / ph.value
