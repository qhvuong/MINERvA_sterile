#pragma once

#include <string>
#include <vector>

class TChain;
class TEventList;

namespace PlotUtils
{
    class MnvH1D;
}

class FluxCalculatorLoop
{
public:
    FluxCalculatorLoop() {};
    virtual ~FluxCalculatorLoop() {};

    void EventLoop(
        TChain* chain,
        const TEventList* evtList,
        PlotUtils::MnvH1D* histogram,
        std::string branchName,
        double additionalWeight = 1.0,
        bool cvWeighted = false,
        unsigned int beamUniverseOffset = 0,
        int beamFocusMode = 0,
        const std::vector<double>& flatBeamFocusWeights =
            std::vector<double>(),
        const std::vector<double>& beamFocusEnergyEdges =
            std::vector<double>(),
        const std::vector<double>& beamFocusResponseWeights =
            std::vector<double>()
    );
};