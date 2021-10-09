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

// Particle m1 = Particle(1, 500, 0, 0, 10);
// Particle m2 = Particle(1, 300, 0, 0, 10);
// Spring s1 = Spring(&m1, &m2, 10, 400);
// Synth<double> mass1 = Synth<double>(&sm.cycle, 500);
// Synth<double> mass2 = Synth<double>(&sm.cycle, 300);
// void particletest(const float* in, float* out, int i)
// {
// 	float sample = (float)(mass1() + mass2());
// 	out[2 * i] = sample;
// 	out[2 * i + 1] = sample;

// 	mass1.freqmod(m1.position);
// 	mass2.freqmod(m2.position);
	
// 	mass1.tick();
// 	mass2.tick();

// 	m1.tick();
// 	m2.tick();

// 	m1.prepare();
// 	m2.prepare();
// 	s1.tick();
// }

//////////////////////////////////////////////////////////

// Bouncer<double> bunc = Bouncer<double>(&sm.cycle, 100, 12, 0.7);
// Synth<double> bunclfo = Synth<double>(&sm.square, 0.1);
// int bouncertest(const float* in, float* out, int i)
// {
// 	float sample = (float)bunc();
// 	out[2 * i] = sample;
// 	out[2 * i + 1] = sample;

// 	// bunc.fundmod(100 + 50 * bunclfo());

// 	bunc.tick();
// 	// bunclfo.tick();

// 	return 0;
// }

// void staticounctest()
// {
// 	bunc.fundmod(200);

// 	for (int i = 0; i < 10; i++)
// 	{
// 		bunc.tick();
// 		bunclfo.tick();
// 	}

// }

// Particle m3 = Particle(1, 200, 0, 0, 10);
// Spring s2 = Spring(&m2, &m3, 10, 100);
// void particletest2()
// {
// 	for (int i = 0; i < 10000; i++)
// 	{
// 		if (i % 1000 == 0)
// 			cout << m1.position << ", " << m2.position << ", " << m3.position << endl;

// 		m1.tick();
// 		m2.tick();
// 		m3.tick();

// 		m1.prepare();
// 		m2.prepare();
// 		m3.prepare();

// 		s1.tick();
// 		s2.tick();
// 	}
// }

//////////////////////////////////////////////////////////

Polyphon<double> poly = Polyphon<double>(&sm.cycle, 2, 1, 0.7);

int gravinit()
{
	poly.request(mtof(68));
	poly.request(mtof(68.5));

	poly.print();
	// poly.bouncers[0].reset(100);
	// poly.bouncers[1].reset(200);
	return 0;
}

int gravtest(const float* in, float* out, int i)
{
	float sample = (float)poly();
	out[2 * i] = sample;
	out[2 * i + 1] = sample;

	poly.prepare();
	poly.tick();
	return 0;
}
