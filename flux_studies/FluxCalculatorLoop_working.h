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
        bool useFlatBeamFocus = false,
        const std::vector<double>& flatBeamFocusWeights =
            std::vector<double>()
    );
};