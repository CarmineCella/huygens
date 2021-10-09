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
	Oscillator(double f, double phi, double k)
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


template <typename T> class Bouncer
{
public:
	Bouncer(Wave<T>* form, double fundamental, uint overtones, double decay, double harmonicity = 1, double k = 2.0 / SR)
		: guide(1, ftom(fundamental))
	{
		this->fundamental = this->target_fundamental = fundamental;
		this->harmonicity = this->target_harmonicity = harmonicity;
		this->overtones = overtones;
		this->decay = this->target_decay = decay;

		resonators.reserve(overtones);
		particles.reserve(overtones);
		springs.reserve(overtones);

		for (int i = 0; i < overtones; i++)
		{
			Particle particle(1, ftom(fundamental * pow(i + 1, harmonicity)), 0, 0, 0.01);
			particles.push_back(particle);

			Synth<T> synth(form, fundamental * pow(i + 1, harmonicity));
			resonators.push_back(synth);
			if (i > 0)
			{
				Spring spring(&particles[i-1], &particles[i], 0, particles[i].position - particles[i-1].position);
				springs.push_back(spring);
				// cout << "paired particle " << i-1 << " at position " << particles[i-1].position << " with particle " << i << " at position " << particles[i].position << endl;
				// cout << "  equilibrium " << abs(particles[i-1].position - particles[i].position) << endl;
			}
			else
			{
				Spring spring(&guide, &particles[0], 0, 0);
				springs.push_back(spring);
				// cout << "paired guide particle at position " << guide.position << " with particle 0 at position " << particles[0].position << endl;
				// cout << springs[0].first << " " << springs[0].second << " " << springs[0].k;
			}
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

	void prepare()
	{
		for (int i = 0; i < overtones; i++)
		{
			particles[i].prepare();
			guide.prepare();
		}
	}

	void tick()
	{
		fundamental = target_fundamental; // * (1 - stiffness) + fundamental * stiffness;
		decay = target_decay * (1 - stiffness) + decay * stiffness;
		harmonicity = target_harmonicity * (1 - stiffness) + harmonicity * stiffness;
		normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;

		guide.position = ftom(fundamental);

		for (int i = 0; i < overtones; i++)
			particles[i].tick();

		for (int i = 0; i < overtones; i++)
		{
			resonators[i].freqmod(mtof(particles[i].position));
			resonators[i].tick();
		}

		for (int i = 0; i < overtones; i++)
			springs[i].tick();
	}

	// assign a new fundamental, reset all phases (no interpolation)
	void reset(double fundamental)
	{
		this->fundamental = this->target_fundamental = fundamental;

		guide.position = ftom(fundamental);
		for (int i = 0; i < overtones; i++)
		{
			double frequency = fundamental * pow(i + 1, harmonicity);
			particles[i].reset(ftom(frequency));
			resonators[i].reset(frequency);
		}
	}

	// modulate the decay parameter
	void decaymod(double target)
	{ target_decay = target; }

	// modulate the harmonicity parameter
	void harmmod(double target)
	{ target_harmonicity = target; }

	// modulate the fundamental
	void fundmod(double target)
	{ target_fundamental = target; }


private:
	uint overtones;
	double target_fundamental;
	double decay;
	double target_decay;
	double normalization;
	double harmonicity;
	double target_harmonicity;
	double stiffness;
	vector<Synth<T>> resonators;
	Particle guide;
	vector<Spring> springs;

public:
	vector<Particle> particles;
	double fundamental;
};


// template <typename T> class Polyphon
// {
// public:
// 	Polyphon(Wave<T>* form, uint voices, uint overtones, double decay, double harmonicity = 1, double k = 2.0 / SR)
// 	{
// 		bouncers.reserve(voices);
// 		active.reserve(voices);
// 		for (int i = 0; i < voices; i++)
// 		{
// 			Bouncer bounce(form, 1, overtones, decay, harmonicity, k);
// 			bouncers.push_back(bounce);
// 			active.push_back(false);
// 		}

// 		this->decay = decay;
// 		this->harmonicity = harmonicity;
// 		this->voices = voices;
// 		this->overtones = overtones;

// 		gravities.reserve(voices * (voices - 1) * overtones * (overtones + 1) / 4);
// 		for (int i = 0; i < voices; i++)
// 		{
// 			for (int j = i + 1; j < voices; j++)
// 			{
// 				for (int k = 0; k < overtones; k++)
// 				{
// 					for (int m = k; m < overtones; m++)
// 					{
// 						Gravity gravity(&bouncers[i].particles[k], &bouncers[j].particles[m], 1000);
// 						gravities.push_back(gravity);
// 					}
// 				}
// 			}
// 		}
// 	}

// 	T operator()()
// 	{
// 		T sample = 0;
// 		for (int i = 0; i < voices; i++)
// 		{
// 			if (active[i])
// 			{
// 				sample += bouncers[i]();
// 			}
// 		}
// 		return (T)sample;
// 	}

// 	void prepare()
// 	{
// 		for (int i = 0; i < voices; i++)
// 		{
// 			bouncers[i].prepare();
// 		}
// 	}

// 	void tick()
// 	{
// 		int count = 0;
// 		for (int k = 0; k < overtones; k++)
// 		{
// 			for (int m = k; m < overtones; m++)
// 			{
// 				for (int i = 0; i < voices; i++)
// 				{
// 					for (int j = i + 1; j < voices; j++)
// 					{
// 						if (active[i] && active[j])
// 						{
// 							gravities[count].tick();
// 						}
// 						count++;
// 					}
// 				}
// 			}
// 		}

// 		for (int i = 0; i < voices; i++)
// 		{
// 			if (active[i])
// 			{
// 				bouncers[i].tick();
// 			}
// 		}

// 	}

// 	// request new voice at a given frequency; return voice number
// 	int request(double fundamental)
// 	{
// 		int voice = -1;
// 		for (int i = 0; i < voices; i++)
// 		{
// 			if (!active[i])
// 			{
// 				voice = i;
// 				break;
// 			}
// 		}

// 		if (voice >= 0)
// 		{
// 			active[voice] = true;
// 			bouncers[voice].reset(fundamental);
// 			cout << "allocated voice " << voice << endl; 
// 		}

// 		return voice;
// 	}

// 	void release(int voice)
// 	{
// 		if (voice >= 0)
// 		{
// 			active[voice] = false;
// 			return;
// 		}

// 		for (int i = 0; i < voices; i++)
// 		{
// 			voices[i] = false;
// 		}
// 	}

// 	// modulate the decay parameter
// 	// void decaymod(double target)
// 	// { target_decay = target; }

// 	// modulate the harmonicity parameter
// 	// void harmmod(double target)
// 	// { target_harmonicity = target; }

// 	// modulate the fundamental
// 	// void fundmod(double target)
// 	// { target_fundamental = target; }

// 	void print()
// 	{
// 		for (int i = 0; i < voices; i++)
// 			if (active[i])
// 				cout << "voice " << i << " active at frequency " << bouncers[i].fundamental << "\n";
// 	}


// private:
// 	uint overtones;
// 	uint voices;
// 	double decay;
// 	double harmonicity;

// public:
// 	vector<Bouncer<T>> bouncers;
// 	vector<bool> active;
// 	vector<Gravity> gravities;
// };

template <typename T> class Polyphon
{
public:
	uint voices; // number of harmonic series
	uint overtones; // number of overtones per series
	double decay; // rate of decay of overtones
	double harmonicity; // harmonicity

	Wave<T>* waveform; // base waveform
	// vector<Oscillator<T>> oscillators; // phasors freq-modulated by
	// vector<Particle> particles; // particles, storing frequencies
	// vector<Particle> guides;
	// vector<Spring> springs; // bind overtones together
	// vector<Gravity> gravities; // pull modes toward one another across voices
	// Spring** springs;
	// Gravity** gravities;

	Wave<T>* waveform;
	Oscillator<T>* oscillators;
	Particle* particles;
	Particle* guides;
	Spring* springs;
	Gravity* gravities;

	~Polyphon()
	{
		delete [] oscillators;
		delete [] particles;
		delete [] guides;
		delete [] springs;
		delete [] gravities;

		// int count = 0;
		// for (int i = 0; i < voices; i++)
		// {
		// 	for (int j = 0; j < overtones; j++)
		// 	{
		// 		delete oscillators[count];
		// 		count++;
		// 	}
		// }
		// delete oscillators;

		// count = 0;
		// for (int i = 0; i < voices; i++)
		// {
		// 	for (int j = 0; j < overtones; j++)
		// 	{
		// 		delete particles[count];
		// 		count++
		// 	}
		// }
		// delete particles;

		// count = 0;
		// for (int i = 0; i < voices; i++)
		// {
		// 	for (int j = i + 1; j < voices; j++)
		// 	{
		// 		for (int k = 0; k < overtones; k++)
		// 		{
		// 			for (int m = 0; m < overtones; m++)
		// 			{
		// 				delete gravities[count];
		// 				count++;
		// 			}
		// 		}
		// 	}
		// }

		// delete gravities;


	}

	Polyphon(Wave<T>* waveform, uint voices, uint overtones, double decay, double harmonicity = 1, double k = 2.0 / SR)
	{
		this->waveform = form;
		this->voices = voices;
		this->overtones = overtones;
		this->decay = decay;
		this->harmonicity = harmonicity;

		oscillators = new Oscillator[voices * overtones];
		particles = new Particle[voices * overtones];
		guides = new Particle[voices];
		springs = new Spring[voices * overtones];
		gravities = new Gravity[voices * (voices - 1) * overtones * overtones / 2];

		int count = 0;
		for (int i = 0; i < voices; i++)
		{
			for (int j = i + 1; j < voices; j++)
			{
				for (int k = 0; k < overtones; k++)
				{
					for (int m = 0; m < overtones; m++)
					{
						gravities[count].initialize(&particles[i * overtones + k], &particles[j * overtones + m], 1000);
						count++;
					}
				}
			}
		}
	}

	T operator()()
	{
		T sample = 0;
		for (int i = 0; i < voices; i++)
		{
			if (active[i])
			{
				for (int j = 0; j < overtones; j++)
				{
					sample += (*waveform)(oscillators[i * overtones + j]());
				}
			}
		}
		return (T)sample;
	}

	// zero out forces on particles
	void prepare()
	{
		for (int i = 0; i < voices * overtones; i++)
		{
			particles[i].prepare();
		}
	}

	void tick()
	{
		int count = 0;
		for (int k = 0; k < overtones; k++)
		{
			for (int m = k; m < overtones; m++)
			{
				for (int i = 0; i < voices; i++)
				{
					for (int j = i + 1; j < voices; j++)
					{
						if (active[i] && active[j])
						{
							gravities[count].tick();
						}
						count++;
					}
				}
			}
		}

		for (int i = 0; i < voices; i++)
		{
			if (active[i])
			{
				bouncers[i].tick();
			}
		}

	}

	// request new voice at a given frequency; return voice number
	int request(double fundamental)
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
			active[voice] = true;
			guides[voice].initialize(1, fundamental);
			for (int j = 0; j < overtones; j++)
			{
				double frequency = fundamental * pow(j + 1, harmonicity);
				particles[voice * overtones + j].initialize(1, ftom(frequency));
				oscillators[voice * overtones + j].initialize(1, frequency);
			}
			bouncers[voice].reset(fundamental);
			cout << "allocated voice " << voice << endl; 
		}

		return voice;
	}

	void release(int voice)
	{
		if (voice >= 0)
		{
			active[voice] = false;
			return;
		}

		for (int i = 0; i < voices; i++)
		{
			voices[i] = false;
		}
	}

	void print()
	{
		for (int i = 0; i < voices; i++)
			if (active[i])
				cout << "voice " << i << " active at frequency " << bouncers[i].fundamental << "\n";
	}
};