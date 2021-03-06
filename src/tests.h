// tests.h
#include "soundmath.h"

Synth<double> osc1 = Synth<double>(&cycle, 330);
Synth<double> osc2 = Synth<double>(&cycle, 165);
Synth<double> osc3 = Synth<double>(&cycle, 1);
Synth<double> osc4 = Synth<double>(&square, 0.05);
int fmtest(float* out, int i)
{
	out[2 * i] = 0.5 * (float)osc1();
	out[2 * i + 1] = 0.5 * (float)osc1();
	// osc1.freqmod(330 + (165 + 165 * osc3()) * osc2());
	osc3.freqmod(10 + 10 * osc4());
	osc1.freqmod(330 + (165 + 165 * osc3()));
	osc1.tick();
	// osc2.tick();
	osc3.tick();
	osc4.tick();

	return 0;
}

//////////////////////////////////////////////////////////

// Delay<double> delay(10, 2 * SR); // 10 delay lines, buffer size 2 seconds

// void mirrorinit()
// {
// 	delay.coefficients({{0,1}}, {{20000,0.5}, {30000,0.5}});
// }


// int mirrortest(const float* in, float* out, int i)
// {
// 	float sample = (float)(delay((in[i])));
// 	out[2 * i] = sample;
// 	out[2 * i + 1] = sample;

// 	delay.tick();
// 	// lowpass.tick();

// 	return 0;
// }

//////////////////////////////////////////////////////////

Reson<double> big = Reson<double>(&cycle, 100, 8, 0.9);
Synth<double> bigosc = Synth<double>(&cycle, 50);
Synth<double> biglfo = Synth<double>(&cycle, 0.25);
int resontest(const float* in, float* out, int i)
{
	float sample = (float)big();
	out[2 * i] = sample;
	out[2 * i + 1] = sample;

	// big.fundmod(100 + 50 * bigosc());
	big.decaymod(0.8 + 0.1 * biglfo());

	big.tick();
	bigosc.tick();
	biglfo.tick();

	return 0;
}

//////////////////////////////////////////////////////////

// Polyphon<double> poly = Polyphon<double>(&cycle, 16, 13, 0.5);
// Polyphon<double> poly = Polyphon<double>(&saw, 16, 9, 0.4);
// Polyphon<double> poly = Polyphon<double>(&cycle, 16, 11, 0.7);
// Polyphon<double> poly(&saw, 10, 11, 0.7);
Polyphon<double> poly(&cycle, 7, 11, 0.7);
// Polyres<double> poly(7, 7, 0.7);

double amplification = 4;

int polytest(const float* in, float* out, int i)
{
	float sample = 2 * poly();
	// float sample = poly(amplification * 2 * (double)rand() / RAND_MAX - 1);
	// float sample = poly(amplification * (((1 + lfo()) / 2) * in[i] + (1 - (1 + lfo()) / 2) * 2 * (double)rand() / RAND_MAX - 1) / 2);
	// float sample = poly(in[i]);
	out[2 * i] = sample;
	out[2 * i + 1] = sample;

	if (i % 8 == 0)
	{
		poly.physics();
	}

	poly.tick();

	return 0;
}

//////////////////////////////////////////////////////////

Synth<double> oscA(&cycle, 100, 0, 0);
Synth<double> oscB(&cycle, 200, 0, 0);
Synth<double> oscC(&cycle, 300, 0, 0);
Synth<double> oscD(&cycle, 400, 0, 0);
Synth<double> oscE(&cycle, 500, 0, 0);
Synth<double> oscF(&cycle, 600, 0, 0);
Synth<double> oscG(&cycle, 700, 0, 0);

Synth<double> oscUp(&cycle, 300, 0, 0);
Synth<double> oscDown(&cycle, 225, 0, 0);
Synth<double> toggle(&square, 0.05);
Filter<double> smoother({0.0001}, {0, -0.9999});

Synth<double> lfo(&cycle, 0.05);
double sensitivity = 0;

int n = 5;
Synth<double>* synths[] = { &oscA, &oscB, &oscC, &oscD, &oscE, &oscF, &oscG };

int snaptest(const float* in, float* out, int i)
{
	out[2 * i] = oscUp();
	out[2 * i + 1] = oscDown();

	double distance = oscUp.phase - oscDown.phase;
	oscUp.phasemod(-distance * 0.01 * smoother(1 + toggle()) / 2);
	oscDown.phasemod(distance * 0.01 * smoother(1 + toggle()) / 2);

	oscUp.tick();
	oscDown.tick();
	toggle.tick();
	smoother.tick();
	return 0;
}

int phasetest(const float* in, float* out, int i)
{
	double sample = 0;
	for (int j = 0; j < n; j++)
		sample += (*synths[j])();
	sample /= n;

	out[2 * i] = (float)sample;
	out[2 * i + 1] = (float)sample;

	sensitivity = (1 + lfo()) / 120;
	for (int j = 0; j < n; j++)
		for (int k = j + 1; k < n; k++)
		{
			double distance = synths[j]->phase - synths[k]->phase;
			synths[k]->phasemod(distance * sensitivity);
			synths[j]->phasemod(-distance * sensitivity);
		}
	
	for (int j = 0; j < n; j++)
	{
		synths[j]->tick();
	}

	lfo.tick();

	return 0;
}

//////////////////////////////////////////////////////////

Synth<double> noise1(&noise, 0.05);

int noisetest(const float* in, float* out, int i)
{
	float sample = noise1();
	noise1.phasemod(sample);

	if (sample > 1 || sample < -1)
		cout << "failure " << noise1.phase << endl;

	// sample = lowpass(sample);
	out[2 * i] = sample;
	out[2 * i + 1] = sample;

	noise1.tick();
	// lowpass.tick();
	return 0;
}


// Filter<double> bp({0.001,0,-0.001}, {1, -2 * 0.99999 * cos(2 * PI * 300 / SR), 0.99999 * 0.99999});
// Filter<double> bp2({0.001,0,-0.001}, {1, -2 * 0.99999 * cos(2 * PI * 450 / SR), 0.99999 * 0.99999});
// Filter<double> bp3({0.001,0,-0.001}, {1, -2 * 0.99999 * cos(2 * PI * 750 / SR), 0.99999 * 0.99999});
Filter<double> bp({0, 0, 0}, {0, 0, 0});
Filter<double> bp2({0, 0, 0}, {0, 0, 0});
Filter<double> bp3({0, 0, 0}, {0, 0, 0});
Synth<double> bplfo(&cycle, 1);

void bandpasstest(const float* in, float* out, int i)
{
	double random = (double)rand() / RAND_MAX * 2 - 1;
	float sample = (bp(random) + bp2(random) + bp3(random)) / 3;
	if (sample > 1 || sample < -1)
		cout << "failure" << endl;

	out[2 * i] = sample;
	out[2 * i + 1] = sample;

	bp.resonant(300 + 3 * bplfo(), 0.99999);
	bp2.resonant(450 + 4.5 * bplfo(), 0.99999);
	bp3.resonant(750 + 7.5 * bplfo(), 0.99999);

	bp.tick();
	bp2.tick();
	bp3.tick();
	bplfo.tick();
}

//////////////////////////////////////////////////////////

Additive<double> addi(&cycle, 10, 7, 0.7);

void addtest(const float* in, float* out, int i)
{
	float sample = addi();

	out[2 * i] = sample;
	out[2 * i + 1] = sample;

	if (i % 8 == 0)
		addi.physics();

	addi.tick();
}