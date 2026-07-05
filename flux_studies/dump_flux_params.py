#!/usr/bin/env python
import ROOT
import sys

path = sys.argv[1]
f = ROOT.TFile.Open(path)

print("[FILE]", path)
for key in f.GetListOfKeys():
    name = key.GetName()
    obj = key.ReadObj()

    cls = obj.ClassName()
    if "TParameter" in cls:
        try:
            print(name, cls, obj.GetVal())
        except Exception:
            print(name, cls, "<no GetVal>")
    elif cls == "TObjString":
        print(name, cls, obj.GetString().Data())