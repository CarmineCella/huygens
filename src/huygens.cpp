// huygens.cpp // main
// an odd kind of sympathy

// standard libraries
#include <stdexcept>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <signal.h>

// user headers
#include "soundmath.h"
#include "audio.h"
#include "midi.h"

using namespace std;

static bool running = true;
void interrupt(int ignore)
{
	cout << " [keyboard interrupt] ...\n" << endl;
	running = false;
}

static Soundmath<double> sm = Soundmath<double>();
static Synth<double> osc1 = Synth<double>(&sm.cycle, 330, 0, 0.0001);
static Synth<double> osc2 = Synth<double>(&sm.cycle, 165, 0, 0.0001);
static Synth<double> osc3 = Synth<double>(&sm.square, 20, 0, 0.0001);
static Synth<double> osc4 = Synth<double>(&sm.cycle, 0.05, 0, 0.0001);


int process(float* in, float* out)
{
	for (int i = 0; i < bsize; i++)
	{
		// out[2 * i] = 0.5 * (float)sm.cycle(osc1());
		// out[2 * i + 1] = 0.5 * (float)sm.cycle(osc1());
		out[2 * i] = 0.5 * (float)osc1();
		out[2 * i + 1] = 0.5 * (float)osc1();
		
		osc1.freqmod(330 + (165 + 165 * osc3()) * osc2());
		osc3.freqmod(10 + 10 * osc4());
		
		osc1.tick();
		osc2.tick();
		osc3.tick();
		osc4.tick();
	}

	return 0;
}

static Audio A = Audio(process);
static MidiIn MI = MidiIn();
static MidiOut MO = MidiOut();

int main()
{
	// bind keyboard interrupt to program exit
	signal(SIGINT, interrupt);

	A.startup(); // startup audio engine

	while (running)
		Pa_Sleep(1000);

	A.shutdown(); // shutdown audio engine
	return 0;

	MI.startup(); // startup midi engine
	MO.startup();
	MI.getports();
	MO.getports();

	MI.open("Launchpad Pro MK3 LPProMK3 MIDI");
	MO.open("Launchpad Pro MK3 LPProMK3 MIDI");

	std::vector<unsigned char> message;
	int nBytes, i;
	double stamp;

	while (running)
	{
		stamp = MI.get(&message);
		nBytes = message.size();
		for (i = 0; i < nBytes; i++)
		{
			std::cout << "Byte " << i << " = " << (int)message[i] << ", ";
		}
		if (nBytes > 0)
			std::cout << "stamp = " << stamp << std::endl;

		Pa_Sleep(1);
	}

	MI.shutdown(); // shutdown midi engine
	MO.shutdown();

}