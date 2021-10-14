// soundmath.h

#include "includes.h"
#include "particle.h"

#define FORCE 50000

const int TABSIZE = 65536;

enum Interp
{
	none = 0, linear = 1, quadratic = 2, cubic = 3
};

template <typename T> class Wave
{
public:
	T table[TABSIZE];

	Wave() { }

	Wave(function<T(double)> shape, Interp interp = Interp::cubic)
	{
		this->interp = interp;
		for (int i = 0; i < TABSIZE; i++)
		{
			table[i] = shape((double)i / TABSIZE);
		}
	}

	T none(int center)
	{
		return table[center];
	}

	T linear(int center, int after, double disp)
	{
		return table[center] * (1 - disp) + table[after] * disp;
	}

	T quadratic(int before, int center, int after, double disp)
	{
		return table[before] * ((disp - 0) * (disp - 1)) / ((-1 - 0) * (-1 - 1)) + 
			   table[center] * ((disp + 1) * (disp - 1)) / (( 0 + 1) * ( 0 - 1)) + 
			   table[after]  * ((disp + 1) * (disp + 0)) / (( 1 + 1) * ( 1 - 0));
	}

	T cubic(int before, int center, int after, int later, double disp)
	{
		return table[before] * ((disp - 0) * (disp - 1) * (disp - 2)) / ((-1 - 0) * (-1 - 1) * (-1 - 2)) + 
			   table[center] * ((disp + 1) * (disp - 1) * (disp - 2)) / (( 0 + 1) * ( 0 - 1) * ( 0 - 2)) + 
			   table[after]  * ((disp + 1) * (disp + 0) * (disp - 2)) / (( 1 + 1) * ( 1 - 0) * ( 1 - 2)) + 
			   table[later]  * ((disp + 1) * (disp + 0) * (disp - 1)) / (( 2 + 1) * ( 2 - 0) * ( 2 - 1));
	}

	T lookup(double phase)
	{
		phase += 1;
		phase -= int(phase);
		int center = (int)(phase * TABSIZE) % TABSIZE;
		int before = (center - 1 + TABSIZE) % TABSIZE;
		int after = (center + 1) % TABSIZE;
		int later = (center + 2) % TABSIZE;

		double disp = (phase * TABSIZE - center);
		disp -= int(disp);

		// interpolation
		switch (interp)
		{
			case Interp::none : return none(center);
			case Interp::linear : return linear(center, after, disp);
			case Interp::quadratic : return quadratic(before, center, after, disp);
			case Interp::cubic : return cubic(before, center, after, later, disp);
		}
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

private:
	Interp interp;
};

// container for holding Waves
// template <typename T> class Soundmath
// {
// public:
// 	Wave<T> cycle = Wave<T>([] (double phase) -> T { return sin(2 * PI * phase); });
	
// };

Wave<double> saw([] (double phase) -> double { return 2 * phase - 1; });
Wave<double> triangle([] (double phase) -> double { return abs(fmod(4 * phase + 3, 4.0) - 2) - 1; });
Wave<double> square([] (double phase) -> double { return phase > 0.5 ? 1 : (phase < 0.5 ? -1 : 0); });
Wave<double> phasor([] (double phase) -> double { return phase; });
Wave<double> noise([] (double phase) -> double { return  2 * ((double)rand() / RAND_MAX) - 1; }, Interp::linear);
Wave<double> cycle([] (double phase) -> double { return sin(2 * PI * phase); });

// An Oscillator has an associated frequency and phase. The frequency is used to update the phase each sample.
// This update requires a call to tick.
// Phase and frequency can be modulated; such modulation is smoothed according to a "drag coefficient"
template <typename T> class Oscillator
{
public:
	double phase;

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
		double weight = (1 - stiffness) * cycle(2 * abs(target_phase - phase) + 0.25); // deals with ambiguity of phasemod(0.5)
		phase = weight * target_phase + (1 - weight) * phase;

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
		// ensure target_phase is in [0,1)
		target_phase += offset;
		target_phase -= int(target_phase);
		target_phase += 1;
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
	double* pitches;

	~Polyphon()
	{
		delete [] oscillators;
		delete [] particles;
		delete [] guides;
		delete [] springs;
		delete [] gravities;
		delete [] active;
		delete [] amplitudes;
		delete [] pitches;
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
		pitches = new double[voices];

		for (int i = 0; i < voices; i++)
		{
			active[i] = 0;
			amplitudes[i] = 0;
		}

		for (int i = 0; i < voices; i++)
		{
			springs[i * overtones].bind(&guides[i], &particles[i * overtones]);
			springs[i * overtones].strength(2 * FORCE);
			for (int j = 1; j < overtones; j++)
			{
				springs[i * overtones + j].bind(&particles[i * overtones + j - 1], &particles[i * overtones + j]);
				springs[i * overtones + j].strength(1 * FORCE);
			}

		}

		int count = 0;
		for (int i = 0; i < voices; i++)
			for (int j = i + 1; j < voices; j++)
				for (int k = 0; k < overtones; k++)
					for (int m = 0; m < overtones; m++)
					{
						gravities[count].bind(&particles[i * overtones + k], &particles[j * overtones + m]);
						gravities[count].strength(1 * FORCE);
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
				guides[i].mass = active[i];
				for (int j = 0; j < overtones; j++)
				{
					particles[i * overtones + j].prepare();
					// particles[i * overtones + j].mass = active[i] * pow(decay, j);
					particles[i * overtones + j].mass = active[i];
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

	void makenote(double pitch, double velocity)
	{
		pitches[request(mtof(pitch), dbtoa(-8 * (1 - (double)velocity / 127)))] = pitch;
	}

	void endnote(double pitch)
	{
		for (int j = 0; j < voices; j++)
			if (pitches[j] == pitch)
				release(j);
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


template <typename T> class Phasepusher
{
public:
	Phasepusher() { }
	~Phasepusher()
	{
		delete [] oscillators;
		delete [] active;
		delete [] amplitudes;
	}

	uint voices; // number of harmonic series
	uint overtones; // number of overtones per series
	double decay; // rate of decay of overtones
	double normalization; // to accommodate different decay rates
	double harmonicity; // harmonicity

	Wave<T>* waveform;
	Oscillator<T>* oscillators;

	double* active;
	double* amplitudes;

	Phasepusher(Wave<T>* waveform, uint voices, uint overtones, double decay, double harmonicity = 1)
	{
		this->waveform = waveform;
		this->voices = voices;
		this->overtones = overtones;
		this->decay = decay;
		this->normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;
		this->harmonicity = harmonicity;

		oscillators = new Oscillator<T>[voices * overtones];
		active = new double[voices];
		amplitudes = new double[voices];

		for (int i = 0; i < voices; i++)
		{
			active[i] = 0;
			amplitudes[i] = 0;
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
								// gravities[count].tick(); // allow their overtones to interact
								count++;
							}
					else
						count += overtones * overtones;
				}
			else
				count += overtones * overtones * (voices - i - 1);
		}
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
					// oscillators[i * overtones + j].freqmod(mtof(particles[i * overtones + j]()));
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
			for (int j = 0; j < overtones; j++)
			{
				previous = frequency;
				frequency = fundamental * pow(j + 1, harmonicity);

				oscillators[voice * overtones + j].reset(frequency);
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

};


template <typename T> class Filter
{
public:
	Filter() { }
	~Filter() { }


	// initialize a filter with its feedforward and feedback coefficients
	Filter(const vector<T>& feedforward, const vector<T>& feedback)
	{
		forward = feedforward;
		back = feedback;
		back[0] = 0;
		order = max(forward.size(), back.size());

		forward.resize(order, 0);
		back.resize(order, 0);

		forget();
	}

	// initialize a filter with its gain, and the poles and zeros of its transfer function
	Filter(T gain, const vector<T>& zeros, const vector<T>& poles)
	{
		order = max(zeros.size(), poles.size());

		forward = coefficients(zeros);
		reverse(forward.begin(), forward.end());
		for (int i = 0; i < order; i ++)
			forward[i] *= gain;

		back = coefficients(poles);
		reverse(back.begin(), back.end());
		back[0] = 0;

		forward.resize(order, 0);
		back.resize(order, 0);

		forget();
	}

	void forget()
	{
		output = vector<T>(order, 0);
		history = vector<T>(order, 0);
		origin = 0;
		computed = false;
	}

	// get coefficients of a polynomial multiplication (z - a1)(z - a2) ... given roots
	static vector<T> coefficients(const vector<T>& zeros, int start = 0)
	{
		int degree = zeros.size() - start;
		if (degree == 0)
			return {0};
		if (degree == 1)
			return {zeros[start], 1};

		vector<T> total(degree + 1, 0);
		T first = zeros[start];

		vector<T> shifted(coefficients(zeros, start + 1));

		for (int i = 0; i < degree; i++)
		{
			total[i] += first * shifted[i];
			total[i + 1] += shifted[i];
		}

		return total;
	}

	void tick()
	{
		origin++;
		origin %= order;
		computed = false;
	}

	T operator()(T sample)
	{
		if (!computed)
			history[origin] = sample;

			T out = 0;
			int index;
			for (int i = 0; i < order; i++)
			{
				index = (origin - i + order) % order;
				out += forward[i] * history[index] - back[i] * output[index];
			}

			output[origin] = out;
			computed = true;

		return output[origin];
	}

private:
	int order;
	int origin;
	bool computed;

	vector<T> forward;
	vector<T> back;
	vector<T> output;
	vector<T> history;
};

template <typename T> class Delay
{
public:
	Delay() { }
	~Delay() { }

	// a "sparse" filter
	Delay(const vector<pair<uint, T>>& feedforward, const vector<pair<uint, T>>& feedback)
	{
		forward = feedforward;
		back = feedback;

		order = 0;
		for (int i = 0; i < forward.size(); i++)
			if (forward[i].first > order)
				order = forward[i].first;

		for (int i = 0; i < back.size(); i++)
		{
			if (back[i].first > order)
				order = back[i].first;

			if (back[i].first == 0)
				back[i].second = 0;
		}

		forget();
	}

	void forget()
	{
		output = vector<T>(order, 0);
		history = vector<T>(order, 0);
		origin = 0;
		computed = false;
	}

	void tick()
	{
		origin++;
		origin %= order;
		computed = false;
	}

	T operator()(T sample)
	{
		if (!computed)
			history[origin] = sample;

			T out = 0;
			uint index, delay;
			T attenuation;

			for (int i = 0; i < forward.size(); i++)
			{
				delay = forward[i].first;
				attenuation = forward[i].second;
				index = (origin - delay + order) % order;
				out += attenuation * history[index];
			}

			for (int i = 0; i < back.size(); i++)
			{
				delay = back[i].first;
				attenuation = back[i].second;
				index = (origin - delay + order) % order;
				out -= attenuation * output[index];
			}

			output[origin] = out;
			computed = true;

		return output[origin];
	}




private:
	int order;
	int origin;
	bool computed;

	vector<pair<uint, T>> forward;
	vector<pair<uint, T>> back;
	vector<T> output;
	vector<T> history;
};