#pragma once

#include <string>
class TChain;
class TEventList;
namespace PlotUtils
{
        class MnvH1D;
}
class FluxCalculatorLoop
{
        public:
                FluxCalculatorLoop()  {};
                virtual ~FluxCalculatorLoop() {};
                
                void EventLoop(TChain * chain, const TEventList * evtList, PlotUtils::MnvH1D * histogram, std::string branchName, double additionalWeight=1, bool cvWeighted=false);
        
        private:
};