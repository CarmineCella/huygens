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
		polytest(in, out, i);
		// phasetest(in, out, i);
		// fmtest(out, i);
		// noisetest(in, out, i);
		// mirrortest(in, out, i);
		// snaptest(in, out, i);
		// bandpasstest(in, out, i);
		// restest(in, out, i);
		// subtest(in, out, i);
		// addtest(in, out, i);
		// fbanktest(in, out, i);
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

	// for (int j = 0; j < 2; j++)
	// {
	// 	double frequency = 400 * (j + 1);
	// 	pair<vector<double>, vector<double>> params = FilterType<double>::resonant(0.01, frequency, 0.9999);
	// 	fbank.coefficients(params.first, params.second, j);
	// 	for (int i = 0; i < 3; i++)
	// 		cout << params.first[i] << ", ";
	// 	cout << endl;
	// 	for (int i = 0; i < 3; i++)
	// 		cout << params.second[i] << ", ";
	// 	cout << endl << endl;
	// }

	A.startup(); // startup audio engine

	MO.startup(); // startup midi out
	MO.getports();

	MI.startup(); // startup midi in
	MI.ignore();
	MI.getports();
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
				poly.makenote(pitch, dbtoa(-8 * (1 - (double)velocity / 127)));
			else
				poly.endnote(pitch);
		}

		Pa_Sleep(1);
	}

	MO.shutdown(); // shutdown midi engine
	MI.shutdown();

	A.shutdown(); // shutdown audio engine

	return 0;
}