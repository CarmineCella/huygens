// huygens.cpp // main
// an odd kind of sympathy

// user headers
#include "tests.h" // gets soundmath

#include "audio.h"
#include "midi.h"
#include "launchpad.h"

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
		// polytest(in, out, i);
		// phasetest(in, out, i);
		// fmtest(out, i);
		// noisetest(in, out, i);
		mirrortest(in, out, i);

	}

	return 0;
}

static Audio A = Audio(process);
static MidiIn MI = MidiIn();
static MidiOut MO = MidiOut();
static Launchpad L = Launchpad();

void welcome()
{
	int notes[] = {48, 58, 65, 64, 69, 74, 67};
	int n = 7;

	for (int i = 0; i < n; i++)
	{
		poly.makenote(notes[i], 127);
		Pa_Sleep(300);
	}

	Pa_Sleep(2000);

	for (int i = 0; i < n; i++)
	{
		poly.endnote(notes[i]);
		Pa_Sleep(300);
	}
}

int main()
{
	// bind keyboard interrupt to program exit
	signal(SIGINT, interrupt);

	// double max = -1;
	// double min = 1;
	// for (int i = 0; i < TABSIZE; i++)
	// {
	// 	double sample = sm.noise.table[i];
	// 	if (sample > max)
	// 		max = sample;
	// 	if (sample < min)
	// 		min = sample;
	// }

	// cout << max << ", " << min << endl;
	// return 0;
	// Filter<complex<double>> F(1, {-1, 1i, -1i}, {0.7 + 0.7i, 0.7 - 0.7i});
	// double sample = 1;
	// F(sample);
	// for (int i = 0; i < 1000; i++)
	// {
	// 	cout << F(0) << endl;
	// 	F.tick();
	// }

	// vector<complex<double>> coeffs = F.coefficients({-1, 1i, -1i});
	// for (int i = 0; i < 4; i++)
	// 	cout << coeffs[i] << endl;


	// return 0;

	A.startup(); // startup audio engine
	MI.startup(); // startup midi engine
	MI.ignore();
	MO.startup();
	MI.getports();
	MO.getports();

	MI.open("IAC Driver Bus 1");

	vector<unsigned char> message;
	int nBytes;
	double stamp;

	// welcome();

	while (running)
	{
		stamp = MI.get(&message);
		nBytes = message.size();

		if (nBytes && (Status)message[0] == note_on)
		{
			int pitch = (int)message[1];
			int velocity = (int)message[2];
			if (velocity)
				poly.makenote(pitch, velocity);
			else
				poly.endnote(pitch);
		}

		Pa_Sleep(1);
	}

	MI.shutdown(); // shutdown midi engine
	MO.shutdown();

	A.shutdown(); // shutdown audio engine

	return 0;
}