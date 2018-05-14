#include "common.hpp"
#include "WrightFisher.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[])
{
	WF::Row(0, 10, 10);

	cout << r.Q << endl;

	// dmat x(2, 2);

	// dmat z(x);

	// cout << z << endl;
	WF::Matrix wf_ref = WF::Single(2, WF::NEITHER);
	cout << "about to print" << endl;
	wf_ref.Q.debug_print();
	cout << wf_ref.Q.dense() << endl;

	// WF::Matrix wf_ext = WF::Single(2, WF::EXTINCTION);
	// wf_ext.Q.debug_print();
	// cout << wf_ext.Q.dense() << endl;


	return 0;
}