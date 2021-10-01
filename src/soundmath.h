// soundmath.h

#include <vector>
#include <cmath>

typedef unsigned long ulong;

const int TABSIZE = 65536;
const double PI = 3.14159265359;
const int SR = 44100;

template <typename T> class Wave
{
public:
	T table[TABSIZE];

	Wave() { }

	Wave(std::function<T(double)> shape)
	{
		for (int i = 0; i < TABSIZE; i++)
		{
			table[i] = shape((double)i / TABSIZE);
		}
	}

	T lookup(double phase)
	{
		int center = (int)(phase * TABSIZE) % TABSIZE;
		int before = (center - 1) % TABSIZE;
		int after = (center + 1) % TABSIZE;
		int later = (center + 2) % TABSIZE;

		double disp = (phase * TABSIZE - center);
		disp -= int(disp);

		// cubic interpolation
		T sample = table[before] * ((disp - 0) * (disp - 1) * (disp - 2)) / ((-1 - 0) * (-1 - 1) * (-1 - 2)) + 
				   table[center] * ((disp + 1) * (disp - 1) * (disp - 2)) / (( 0 + 1) * ( 0 - 1) * ( 0 - 2)) + 
				   table[after]  * ((disp + 1) * (disp + 0) * (disp - 2)) / (( 1 + 1) * ( 1 - 0) * ( 1 - 2)) + 
				   table[later]  * ((disp + 1) * (disp + 0) * (disp - 1)) / (( 2 + 1) * ( 2 - 0) * ( 2 - 1));

		return sample;
	}

	T operator()(double phase)
	{
		return lookup(phase);
	}

	Wave<T> operator+(const Wave<T>& other)
	{
		Wave<T> sum;
		for (int i = 0; i < TABSIZE; i++)
		{
			sum.table[i] = this->table[i] + other.table[i];
		}
		return sum;
	}
};

// container for holding Waves
template <typename T> class Soundmath
{
public:
	Wave<T> cycle = Wave<T>([] (double phase) -> T { return sin(2 * PI * phase); });
	Wave<T> saw = Wave<T>([] (double phase) -> T { return 2 * phase - 1; });
	Wave<T> triangle = Wave<T>([] (double phase) -> T { return abs(std::fmod(4 * phase + 3, 4.0) - 2) - 1; });
	Wave<T> square = Wave<T>([] (double phase) -> T { return phase > 0.5 ? 1 : (phase < 0.5 ? -1 : 0); });
	Wave<T> phasor = Wave<T>([] (double phase) -> T { return phase; });
};


// An Oscillator has an associated frequency and phase. The frequency is used to update the phase each sample.
// This update requires a call to tick.
// Phase and frequency can be modulated; such modulation is smoothed according to a "drag coefficient"
template <typename T> class Oscillator
{
public:
	// k is relaxation time in seconds
	Oscillator(double f, double phi, double k)
	{
		frequency = abs(f);
		target_freq = frequency;
		phase = fmax(0, phi);
		target_phase = phase;
		if (k == 0)
			stiffness = 0;
		else {
			static double epsilon = std::numeric_limits<double>::epsilon();
			static double order = log2(epsilon);
			stiffness = pow(2.0, order / (fmax(0, k) * SR));
		}
	}

	// needs to be called once per sample
	void tick()
	{
		phase += frequency / SR;
		target_phase += frequency / SR;

		frequency = target_freq * (1 - stiffness) + frequency * stiffness;
		phase = target_phase * (1 - stiffness) + phase * stiffness;

		phase -= int(phase);
		target_phase -= int(target_phase);
	}

	double lookup()
	{ return phase; }

	T operator()()
	{ return phase; }

	void freqmod(double target)
	{ target_freq = target; }

	void phasemod(double offset)
	{
		target_phase = phase + offset;
		target_phase -= int(target_phase);
	}

private:
	double frequency;
	double target_freq;
	double phase;
	double target_phase;
	double stiffness;
};


template <typename T> class Synth : public Oscillator<T>
{
public:
	Synth(Wave<T>* form, double f, double phi = 0, double k = 0.05) : Oscillator<T>(f, phi, k)
	{ waveform = form; }

	T operator()()
	{ return (*waveform)(this->lookup()); }

private:
	Wave<T>* waveform;
};