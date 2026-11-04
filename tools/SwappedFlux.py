"""Generated-flux conversion to accompany the mapped-PDG CV flux weight.

For an ordinary generated nue event, the existing mapped-PDG CV call gives
Phi_mu_CV / Phi_mu_generated. This helper supplies Phi_mu_generated /
Phi_e_generated. Do not also apply a CV mu/e flux ratio.
"""

import math
from functools import lru_cache
from pathlib import Path

import ROOT


_SWAP_PDG = {12: 14, -12: -14, 14: 12, -14: -12}
_LE_FLUX_TAG = {
    "minerva1": "minerva1",
    "minerva7": "minerva1",
    "minerva9": "minerva1",
    "minerva13": "minerva13",
    "minerva13c": "minerva13",
    "minerva13e": "minerva13",
    "minerva5": "minerva5",
}
_DEFAULT_FLUX_DIR = Path(__file__).resolve().parents[1] / "custom_plotutils/data/flux"


@lru_cache(maxsize=None)
def _load_generated_flux(flux_dir, playlist_tag, pdg):
    path = Path(flux_dir) / "flux-g4numiv5-pdg{}-{}.root".format(pdg, playlist_tag)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Cannot open generated flux: {}".format(path))
    try:
        h = f.Get("flux_E_unweighted")
        if not h or h.GetDimension() != 1:
            raise RuntimeError("Missing 1D flux_E_unweighted in {}".format(path))
        clone = h.Clone("swap_generated_{}_{}".format(playlist_tag, pdg))
        clone.SetDirectory(0)
        return clone
    finally:
        f.Close()


def _flux_value(hist, energy):
    axis = hist.GetXaxis()
    if not axis.GetXmin() <= energy < axis.GetXmax():
        raise ValueError("True energy {} GeV outside generated-flux range [{}, {})".format(
            energy, axis.GetXmin(), axis.GetXmax()))
    # Match GetFluxCVWeight: interpolate, except above 75 GeV.
    if energy > 75.0:
        return float(hist.GetBinContent(axis.FindFixBin(energy)))
    return float(hist.Interpolate(energy))


def get_swapped_flux_ratio(universe, flux_dir=None):
    """Return Phi_target_generated / Phi_source_generated at true energy.

    Source PDG is the ORIGINAL mc_incoming, before the Python PDG swap.
    Playlist comes from universe.GetPlaylist(); event sign determines nu/nubar.
    Thus wrong-sign events use their own sign within the same beam playlist.
    Files are cached; no absolute-flux threshold, ratio cap, or silent fallback.
    This helper supports the configured LE playlists only. The histogram files
    must use the same absolute flux units and correspond to the MC generation.
    """
    source_pdg = int(universe.mc_incoming)
    if source_pdg not in _SWAP_PDG:
        raise ValueError("Cannot swap incoming PDG {}".format(source_pdg))
    playlist = str(universe.GetPlaylist()).strip().lower()
    if playlist not in _LE_FLUX_TAG:
        raise ValueError("Unsupported LE flux playlist: {!r}".format(playlist))
    energy = float(universe.mc_incomingE) * 1e-3
    if not math.isfinite(energy):
        raise ValueError("Nonfinite true neutrino energy")
    directory = str(Path(flux_dir).expanduser().resolve() if flux_dir is not None
                    else _DEFAULT_FLUX_DIR)
    tag = _LE_FLUX_TAG[playlist]
    source = _load_generated_flux(directory, tag, source_pdg)
    target = _load_generated_flux(directory, tag, _SWAP_PDG[source_pdg])
    denominator = _flux_value(source, energy)
    numerator = _flux_value(target, energy)
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("Invalid generated source flux: playlist={}, PDG={}, E={}, flux={}".format(
            tag, source_pdg, energy, denominator))
    if not math.isfinite(numerator) or numerator < 0:
        raise ValueError("Invalid generated target flux: playlist={}, PDG={}, E={}, flux={}".format(
            tag, _SWAP_PDG[source_pdg], energy, numerator))
    ratio = numerator / denominator
    if not math.isfinite(ratio):
        raise ValueError("Nonfinite generated-flux ratio")
    return ratio
