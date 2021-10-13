// tests.h

#include "soundmath.h"

Soundmath<double> sm = Soundmath<double>();

Synth<double> osc1 = Synth<double>(&sm.cycle, 330);
Synth<double> osc2 = Synth<double>(&sm.cycle, 165);
Synth<double> osc3 = Synth<double>(&sm.square, 1);
Synth<double> osc4 = Synth<double>(&sm.square, 0.05);
int fmtest(float* out, int i)
{
	out[2 * i] = 0.5 * (float)osc1();
	out[2 * i + 1] = 0.5 * (float)osc1();
	osc1.freqmod(330 + (165 + 165 * osc3()) * osc2());
	osc3.freqmod(10 + 10 * osc4());
	// osc1.freqmod(330 + (165 + 165 * osc3()));
	osc1.tick();
	osc2.tick();
	osc3.tick();
	// osc4.tick();

	return 0;
}

//////////////////////////////////////////////////////////

// distortion
int mirrortest(const float* in, float* out, int i)
{
	out[2 * i] = sm.cycle((1 + in[i]) * 2) / 4;
	out[2 * i + 1] = sm.cycle((1 + in[i]) * 2) / 4;

	return 0;
}

//////////////////////////////////////////////////////////

double x1 = 330;
double v1 = 0;
double x2 = 495;
double v2 = 0;
double m = 10;

Synth<double> cycle1 = Synth<double>(&sm.triangle, x1);
Synth<double> cycle2 = Synth<double>(&sm.cycle, x2);

double force(double f1, double f2)
{
	double distance = abs(f1 - f2);
	double cubed = distance * distance * distance;
	double epsilon = 0.0001;
	return 1 / (epsilon + cubed);
}

int gravitytest(const float* in, float* out, int i)
{
	out[2 * i] = 0.5 * (float)(cycle1() + cycle2());
	out[2 * i + 1] = 0.5 * (float)(cycle1() + cycle2());

	v1 += (x2 - x1) * force(x1, x2) / m;
	v2 += (x1 - x2) * force(x1, x2) / m;

	x1 += v1;
	x2 += v2;

	cycle1.freqmod(x1);
	cycle2.freqmod(x2);

	cycle1.tick();
	cycle2.tick();

	return 0;
}

//////////////////////////////////////////////////////////

Reson<double> big = Reson<double>(&sm.cycle, 100, 8, 0.9);
Synth<double> bigosc = Synth<double>(&sm.cycle, 50);
Synth<double> biglfo = Synth<double>(&sm.cycle, 0.25);
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

// Polyphon<double> poly = Polyphon<double>(&sm.cycle, 16, 13, 0.5);
// Polyphon<double> poly = Polyphon<double>(&sm.saw, 16, 9, 0.4);
// Polyphon<double> poly = Polyphon<double>(&sm.cycle, 16, 11, 0.7);
Polyphon<double> poly(&sm.cycle, 10, 11, 0.7);
// Polyphon<double> poly = Polyphon<double>(&sm.cycle, 7, 11, 0.7);

// Filter<double> lowpass({1,0}, {0,0.8});
// Filter<complex<double>> F(1, {-1, 1i, -1i}, {0.7 + 0.7i, 0.7 - 0.7i});


int polytest(const float* in, float* out, int i)
{
	float sample = (float)(poly());
	out[2 * i] = sample;
	out[2 * i + 1] = sample;

	if (i %  8 == 0)
		poly.physics();

	poly.tick();
	return 0;
}

//////////////////////////////////////////////////////////

Synth<double> oscA(&sm.cycle, 100, 0, 0);
Synth<double> oscB(&sm.cycle, 200, 0, 0);
Synth<double> oscC(&sm.cycle, 300, 0, 0);
Synth<double> oscD(&sm.cycle, 400, 0, 0);
Synth<double> oscE(&sm.cycle, 500, 0, 0);
Synth<double> oscF(&sm.cycle, 600, 0, 0);
Synth<double> oscG(&sm.cycle, 700, 0, 0);
Synth<double> lfo(&sm.triangle, 0.05);
double sensitivity = 0;

int n = 7;
Synth<double>* synths[] = { &oscA, &oscB, &oscC, &oscD, &oscE, &oscF, &oscG };

int phasetest(const float* in, float* out, int i)
{
	float sample = 0;
	for (int j = 0; j < n; j++)
		sample += (float)((*synths[j])());
	sample /= n;

	out[2 * i] = sample;
	out[2 * i + 1] = sample;

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
