
import ROOT

path = (
    "root://fndcadoor.fnal.gov:1094/"
    "pnfs/fnal.gov/usr/minerva/persistent/DataPreservation/"
    "Low_Energy_Era/p8/StandardMC/Playlist_minerva1/"
    "MasterAnaDev_mc_AnaTuple_run00010200_Playlist.root"
)

f = ROOT.TFile.Open(path)
if not f or f.IsZombie():
    raise RuntimeError("Cannot open input ROOT file")

print("Available TTrees:")
for key in f.GetListOfKeys():
    obj = key.ReadObj()
    if obj.InheritsFrom("TTree"):
        print(" ", key.GetName(), obj.GetEntries())

# Locate the MasterAnaDev event tree automatically.
trees = []
for key in f.GetListOfKeys():
    obj = key.ReadObj()
    if obj.InheritsFrom("TTree"):
        if obj.GetBranch("mc_charm") and obj.GetBranch("mc_intType"):
            trees.append(obj)

if not trees:
    raise RuntimeError("No tree containing mc_charm and mc_intType")

tree = max(trees, key=lambda t: t.GetEntries())
print(f"\nUsing tree: {tree.GetName()}")

base = "mc_current==1 && mc_incoming==14"
qe = base + " && mc_intType==1"

n_cc = tree.GetEntries(base)
n_qe = tree.GetEntries(qe)
n_charm = tree.GetEntries(qe + " && mc_charm==1")
n_noncharm = tree.GetEntries(qe + " && mc_charm!=1")

print("\n=== RAW EVENT COUNTS ===")
print("CC numu:", n_cc)
print("QE (all):", n_qe)
print("QE charm==1:", n_charm)
print("QE charm!=1:", n_noncharm)

if n_qe:
    print(f"Charm fraction of QE: {100*n_charm/n_qe:.4f}%")

print("\n=== mc_charm values within QE ===")
tree.Draw("mc_charm", qe, "goff")

values = {}
for i in range(tree.GetSelectedRows()):
    x = int(tree.GetV1()[i])
    values[x] = values.get(x, 0) + 1

for value, count in sorted(values.items()):
    print(f"mc_charm={value}: {count}")

f.Close()
