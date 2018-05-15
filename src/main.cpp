#include "common.hpp"
#include "WrightFisher.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[])
{
    int N = atoi(argv[1]);

    WF::Matrix wf_ref = WF::Single(N, WF::NEITHER, 0, 0.5, 1e-9, 1e-9, 0);
    wf_ref.Q.save_market("ref.mtx");
    // cout << wf_ref.Q.dense() << endl << endl;

    WF::Matrix wf_ext = WF::Single(N, WF::EXTINCTION, 0, 0.5, 1e-9, 1e-9, 0);
    wf_ext.Q.save_market("ext.mtx");
    // cout << wf_ext.Q.dense() << endl << endl;
    // cout << wf_ext.R << endl << endl;

    WF::Matrix wf_fix = WF::Single(N, WF::FIXATION, 0, 0.5, 1e-9, 1e-9, 0);
    wf_fix.Q.save_market("fix.mtx");
    // cout << wf_fix.Q.dense() << endl << endl;
    // cout << wf_fix.R << endl << endl;

    WF::Matrix wf_abs = WF::Single(N, WF::BOTH, 0, 0.5, 1e-9, 1e-9, 0);
    wf_abs.Q.save_market("abs.mtx");
    // cout << wf_abs.Q.dense() << endl << endl;
    // cout << wf_abs.R << endl << endl;

    return 0;
}