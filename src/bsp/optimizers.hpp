#pragma once

#include "common.hpp"
#include "utils/logging.hpp"
#include "sco/optimizers.hpp"

namespace BSP {

class BSPTrustRegionSQP : public BasicTrustRegionSQP {
public:
	BSPTrustRegionSQP();
	BSPTrustRegionSQP(OptProbPtr prob);
	void addMeritDoneCallback(const Callback& cb);
	void addAlphaDoneCallback(const Callback& cb);
	OptStatus optimize();
	OptStatus optimize2();
	int max_alpha_increases_;
protected:
	vector<Callback> merit_done_callbacks_;
	vector<Callback> alpha_done_callbacks_;
	void callMeritDoneCallbacks(DblVec& x);
	void callAlphaDoneCallbacks(DblVec& x);;
};

}
