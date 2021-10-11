// soundmath.h

#include "includes.h"
#include "particle.h"

const int TABSIZE = 65536;

template <typename T> class Wave
{
public:
	T table[TABSIZE];

	Wave() { }

	Wave(function<T(double)> shape)
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
	Wave<T> triangle = Wave<T>([] (double phase) -> T { return abs(fmod(4 * phase + 3, 4.0) - 2) - 1; });
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
	Oscillator(double f = 0, double phi = 0, double k = 2.0 / SR)
	{
		frequency = abs(f);
		target_freq = frequency;
		phase = fmax(0, phi);
		target_phase = phase;

		stiffness = relaxation(k);
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

	void reset(double f)
	{
		frequency = target_freq = f;
		phase = target_phase = 0;
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
	Synth(Wave<T>* form, double f, double phi = 0, double k = 2.0 / SR) : Oscillator<T>(f, phi, k)
	{ waveform = form; }

	T operator()()
	{ return (*waveform)(this->lookup()); }

private:
	Wave<T>* waveform;
};


// a bank of synthesizers, with optional harmonicity.
template <typename T> class Reson
{
public:
	Reson(Wave<T>* form, double fundamental, uint overtones, double decay, double harmonicity = 1, double k = 2.0 / SR)
	{
		this->fundamental = this->target_fundamental = fundamental;
		this->harmonicity = this->target_harmonicity = harmonicity;
		this->overtones = overtones;
		this->decay = this->target_decay = decay;

		for (int i = 0; i < overtones; i++)
		{
			Synth<T> synth(form, fundamental * pow(i + 1, harmonicity));
			resonators.push_back(synth);
		}

		normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;

		stiffness = relaxation(k);
	}

	T operator()()
	{
		T sample = 0;
		for (int i = 0; i < overtones; i++)
		{
			sample += pow(decay, i) * resonators[i]() / normalization;
		}
		return (T)sample;
	}

	void tick()
	{
		fundamental = target_fundamental * (1 - stiffness) + fundamental * stiffness;
		decay = target_decay * (1 - stiffness) + decay * stiffness;
		harmonicity = target_harmonicity * (1 - stiffness) + harmonicity * stiffness;

		for (int i = 0; i < overtones; i++)
		{
			resonators[i].freqmod(fundamental * pow(i + 1, harmonicity));
			resonators[i].tick();
		}

		normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;
	}

	void decaymod(double target)
	{ target_decay = target; }

	void harmmod(double target)
	{ target_harmonicity = target; }

	void fundmod(double target)
	{ target_fundamental = target; }

private:
	uint overtones;
	double fundamental;
	double target_fundamental;
	double decay;
	double target_decay;
	double normalization;
	double harmonicity;
	double target_harmonicity;
	double stiffness;
	vector<Synth<T>> resonators;
};


template <typename T> class Polyphon
{
public:
	uint voices; // number of harmonic series
	uint overtones; // number of overtones per series
	double decay; // rate of decay of overtones
	double normalization; // to accommodate different decay rates
	double harmonicity; // harmonicity

	Wave<T>* waveform;
	Oscillator<T>* oscillators;
	Particle* particles;
	Particle* guides;
	Spring* springs;
	Gravity* gravities;

	double* active;
	double* amplitudes;

	~Polyphon()
	{
		delete [] oscillators;
		delete [] particles;
		delete [] guides;
		delete [] springs;
		delete [] gravities;
		delete [] active;
		delete [] amplitudes;
	}

	Polyphon(Wave<T>* waveform, uint voices, uint overtones, double decay, double harmonicity = 1)
	{
		this->waveform = waveform;
		this->voices = voices;
		this->overtones = overtones;
		this->decay = decay;
		this->normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;
		this->harmonicity = harmonicity;

		oscillators = new Oscillator<T>[voices * overtones];
		particles = new Particle[voices * overtones];
		guides = new Particle[voices];
		springs = new Spring[voices * overtones];
		gravities = new Gravity[voices * (voices - 1) * overtones * overtones / 2];
		active = new double[voices];
		amplitudes = new double[voices];

		for (int i = 0; i < voices; i++)
		{
			active[i] = 0;
			amplitudes[i] = 0;
		}

		for (int i = 0; i < voices; i++)
		{
			springs[i * overtones].bind(&guides[i], &particles[i * overtones]);
			springs[i * overtones].strength(1000000);
			for (int j = 1; j < overtones; j++)
			{
				springs[i * overtones + j].bind(&particles[i * overtones + j - 1], &particles[i * overtones + j]);
				springs[i * overtones + j].strength(500000);
			}

		}

		int count = 0;
		for (int i = 0; i < voices; i++)
			for (int j = i + 1; j < voices; j++)
				for (int k = 0; k < overtones; k++)
					for (int m = 0; m < overtones; m++)
					{
						gravities[count].bind(&particles[i * overtones + k], &particles[j * overtones + m]);
						gravities[count].strength(100000);
						count++;
					}
	}

	T operator()()
	{
		T sample = 0;
		for (int i = 0; i < voices; i++)
			if (amplitudes[i])
				for (int j = 0; j < overtones; j++)
					sample += amplitudes[i] * pow(decay, j) * (*waveform)(oscillators[i * overtones + j]()) / (2 * voices * normalization);

		return (T)sample;
	}

	void physics()
	{
		// zero the forces on particles
		for (int i = 0; i < voices; i++)
			if (active[i])
			{
				guides[i].prepare();
				guides[i].mass = active[i] * pow(decay, i);
				for (int j = 0; j < overtones; j++)
				{
					particles[i * overtones + j].prepare();
					particles[i].mass = active[i];
				}
			}

		// apply spring forces
		for (int i = 0; i < voices; i++)
			if (active[i])
				for (int j = 0; j < overtones; j++)
					springs[i * overtones + j].tick();

		// apply forces from gravity
		int count = 0;
		for (int i = 0; i < voices; i++)
		{
			if (active[i]) // if one voice is active
				for (int j = i + 1; j < voices; j++)
				{
					if (active[j]) // as is another
						for (int k = 0; k < overtones; k++)
							for (int m = 0; m < overtones; m++)
							{
								gravities[count].tick(); // allow their overtones to interact
								count++;
							}
					else
						count += overtones * overtones;
				}
			else
				count += overtones * overtones * (voices - i - 1);
		}

		// tick the particles
		for (int i = 0; i < voices; i++)
			if (active[i])
				for (int j = 0; j < overtones; j++)
					particles[i * overtones + j].tick();
	}

	void tick()
	{
		for (int i = 0; i < voices; i++)
			amplitudes[i] = 0.005 * active[i] + 0.995 * amplitudes[i];

		// update oscillator frequencies; tick oscillators
		for (int i = 0; i < voices; i++)
			if (active[i] || amplitudes[i])
				for (int j = 0; j < overtones; j++)
				{
					oscillators[i * overtones + j].freqmod(mtof(particles[i * overtones + j]()));
					oscillators[i * overtones + j].tick();
				}
	}

	// request new voice at a given frequency; return voice number
	int request(double fundamental, double amplitude = 0)
	{
		int voice = -1;
		for (int i = 0; i < voices; i++)
		{
			if (!active[i])
			{
				voice = i;
				break;
			}
		}

		if (voice >= 0)
		{
			double frequency = fundamental;
			double previous;
			active[voice] = amplitude;
			guides[voice].initialize(1, ftom(fundamental));
			for (int j = 0; j < overtones; j++)
			{
				previous = frequency;
				frequency = fundamental * pow(j + 1, harmonicity);

				particles[voice * overtones + j].initialize(1, ftom(frequency));
				oscillators[voice * overtones + j].reset(frequency);
				springs[voice * overtones + j].target(ftom(frequency) - ftom(previous));
			}
			// cout << "allocated voice " << voice << endl; 
		}

		return voice;
	}

	void release(int voice)
	{
		if (voice >= 0)
		{
			active[voice] = 0;
			return;
		}

		for (int i = 0; i < voices; i++)
		{
			active[i] = 0;
		}
	}

	void print()
	{
		for (int i = 0; i < voices; i++)
			if (active[i])
			{
				cout << "voice " << i << " at ";
				for (int j = 0; j < overtones; j++)
					cout << mtof(particles[i * overtones + j]()) << " ";
				cout << endl;
			}
		cout << endl;
	}
};
