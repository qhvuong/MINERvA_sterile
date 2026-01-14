#pragma once

#include <string>
// #include <map>
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
                
                void EventLoop(TChain * chain, const TEventList * evtList, std::vector<PlotUtils::MnvH1D*>& parentHistos, std::string branchName, double additionalWeight=1, bool cvWeighted=false);
        
        private:
};