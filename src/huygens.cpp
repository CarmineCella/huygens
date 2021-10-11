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
	}

	return 0;
}

static Audio A = Audio(process);
static MidiIn MI = MidiIn();
static MidiOut MO = MidiOut();
static Launchpad L = Launchpad();

int main()
{
	// bind keyboard interrupt to program exit
	signal(SIGINT, interrupt);

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

	int pitches[poly.voices];
	for (int i = 0; i < poly.voices; i++)
		pitches[i] = 0;

	while (running)
	{
		stamp = MI.get(&message);
		nBytes = message.size();


		if (nBytes && (int)message[0] == 144)
		{
			// int pitch = 36 + L.note(message[1]);
			int pitch = (int)message[1];
			if ((int)message[2])
				pitches[poly.request(mtof(pitch), dbtoa(-8 * (1 - (double)message[2] / 127)))] = pitch;
			else
			{
				for (int j = 0; j < poly.voices; j++)
					if (pitches[j] == pitch)
						poly.release(j);
			}
		}

		Pa_Sleep(10);
	}

	MI.shutdown(); // shutdown midi engine
	MO.shutdown();

	A.shutdown(); // shutdown audio engine

	return 0;
}