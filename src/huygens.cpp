// huygens.cpp // main
// an odd kind of sympathy

// user headers
#include "tests.h" // gets soundmath

#include "audio.h"
#include "midi.h"

using namespace std;

static bool running = true;
void interrupt(int ignore)
{
	cout << " [keyboard interrupt, exiting ...]" << endl;
	running = false;
}


int process(const float* in, float* out)
{
	for (int i = 0; i < bsize; i++)
	{
		// fmtest(out, i);	
		// mirrortest(in, out, i);
		// gravitytest(in, out, i);
		// resontest(in, out, i);
		// particletest(in, out, i);
		// bouncertest(in, out, i);
		gravtest(in, out, i);
	}

	return 0;
}

static Audio A = Audio(process);
static MidiIn MI = MidiIn();
static MidiOut MO = MidiOut();

int main()
{
	// particletest2();
	// staticbounctest();
	// return 0;

	// bind keyboard interrupt to program exit
	signal(SIGINT, interrupt);

	gravinit();
	// bunc.reset(100);


	A.startup(); // startup audio engine

	// Pa_Sleep(500);

	// bunc.reset(200);

	// Pa_Sleep(500);

	// bunc.reset(300);

	while (true)
	{
		Pa_Sleep(1000);
		poly.print();
	}

	MI.startup(); // startup midi engine
	MO.startup();
	MI.getports();
	MO.getports();

	MI.open("Launchpad Pro MK3 LPProMK3 MIDI");
	MO.open("Launchpad Pro MK3 LPProMK3 MIDI");

	vector<unsigned char> message;
	int nBytes, i;
	double stamp;

	while (running)
	{
		stamp = MI.get(&message);
		nBytes = message.size();
		for (i = 0; i < nBytes; i++)
			cout << "Byte " << i << " = " << (int)message[i] << ", ";
		if (nBytes > 0)
			cout << "stamp = " << stamp << endl;

		Pa_Sleep(1);
	}

	MI.shutdown(); // shutdown midi engine
	MO.shutdown();

	A.shutdown(); // shutdown audio engine
	return 0;


}