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
	Wave() { }
	~Wave() { }

	Wave(function<T(double)> shape, Interp interp = Interp::cubic, double left = 0, double right = 1, bool periodic = true)
	{
		this->interp = interp;
		this->left = left;
		this->right = right;
		this->periodic = periodic;

		for (int i = 0; i < TABSIZE; i++)
		{
			double phase = (double) i / TABSIZE;
			table[i] = shape((1 - phase) * left + phase * right);
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

	T lookup(double input)
	{
		double phase = (input - left) / (right - left);
		
		// get value at endpoint if input is out of bounds
		if (!periodic && (phase < 0 || phase >= 1))
		{
			int center = (phase < 0) ? 0 : TABSIZE - 1;
			return none(center);
		}
		else
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

protected:
	T table[TABSIZE];

private:
	Interp interp;
	double left; // endpoints of range
	double right;
	bool periodic;
};


Wave<double> saw([] (double phase) -> double { return 2 * phase - 1; });
Wave<double> triangle([] (double phase) -> double { return abs(fmod(4 * phase + 3, 4.0) - 2) - 1; });
Wave<double> square([] (double phase) -> double { return phase > 0.5 ? 1 : (phase < 0.5 ? -1 : 0); });
Wave<double> phasor([] (double phase) -> double { return phase; });
Wave<double> noise([] (double phase) -> double { return  2 * ((double)rand() / RAND_MAX) - 1; }, Interp::linear);
Wave<double> cycle([] (double phase) -> double { return sin(2 * PI * phase); });
Wave<double> click([] (double phase) -> double { return (int)(phase == 0); }, Interp::linear);


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

protected:
	T lookup()
	{ return phase; }

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

	void makenote(double pitch, double amplitude)
	{
		pitches[request(mtof(pitch), amplitude)] = pitch;
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
		order = max(forward.size(), back.size()) - 1;

		forward.resize(order + 1, 0);
		back.resize(order + 1, 0);

		reset();
	}

	// initialize a filter with its gain, and the poles and zeros of its transfer function
	Filter(T gain, const vector<T>& zeros, const vector<T>& poles)
	{
		order = max(zeros.size(), poles.size());

		forward = coefficients(zeros);
		reverse(forward.begin(), forward.end());
		for (int i = 0; i < order + 1; i ++)
			forward[i] *= gain;

		back = coefficients(poles);
		reverse(back.begin(), back.end());
		back[0] = 0;

		reset();
	}

	void initialize(uint order)
	{
		forward = vector<T>(order + 1, 0);
		back = vector<T>(order + 1, 0);

		this->order = order;

		reset();
	}

	// adds some padding to the ring buffer
	void reset()
	{
		buffsize = order + 2;
		output = vector<T>(buffsize, 0);
		history = vector<T>(buffsize, 0);
		origin = 0;
		computed = false;
	}

	void forget()
	{
		output = vector<T>(buffsize, 0);
		computed = false;
	}

	// get coefficients of a monic polynomial given its roots
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
		origin %= buffsize;
		computed = false;
	}

	T operator()(T sample)
	{
		if (!computed)
		{
			history[origin] = sample;

			T out = 0;
			int index;
			for (int i = 0; i < order + 1; i++)
			{
				index = (origin - i + buffsize) % buffsize;
				out += forward[i] * history[index] - back[i] * output[index];
			}

			output[origin] = out;
			computed = true;
		}

		return output[origin];
	}

	void resonant(T frequency, T Q)
	{
		T cosine = cos(2 * PI * frequency / SR);

		complex<T> cosine2(cos(4 * PI * frequency / SR), 0);
		complex<T> sine2(sin(4 * PI * frequency / SR), 0);
		
		complex<T> maximum = 1.0 / (Q - 1) - 1.0 / (Q - cosine2 - 1.0i * sine2);
		double amplitude = 1 / sqrt(abs(maximum));

		forward = vector<T>({amplitude, 0, -amplitude});
		back = vector<T>({0, -2 * Q * cosine, Q * Q});
		order = 2;
	}

private:
	int order;
	int origin;
	int buffsize;
	bool computed;

	vector<T> forward;
	vector<T> back;
	vector<T> output;
	vector<T> history;
};


template <typename T> class Polyres
{
public:
	uint voices; // number of harmonic series
	uint overtones; // number of overtones per series
	double decay; // rate of decay of overtones
	double normalization; // to accommodate different decay rates
	double harmonicity; // harmonicity

	Filter<T>* filters;
	Particle* particles;
	Particle* guides;
	Spring* springs;
	Gravity* gravities;

	double* active;
	double* amplitudes;
	double* pitches;
	
	double Q;

	~Polyres()
	{
		delete [] filters;
		delete [] particles;
		delete [] guides;
		delete [] springs;
		delete [] gravities;

		delete [] active;
		delete [] amplitudes;
		delete [] pitches;
	}

	Polyres(uint voices, uint overtones, double decay, double Q = 0.99999, double harmonicity = 1)
	{
		this->voices = voices;
		this->overtones = overtones;
		this->decay = decay;
		this->normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;
		this->harmonicity = harmonicity;
		this->Q = Q;

		filters = new Filter<T>[voices * overtones];
		particles = new Particle[voices * overtones];
		guides = new Particle[voices];
		springs = new Spring[voices * overtones];
		gravities = new Gravity[voices * (voices - 1) * overtones * overtones / 2];

		active = new double[voices];
		amplitudes = new double[voices];
		pitches = new double[voices];

		for (int i = 0; i < voices; i++)
			for (int j = 0; j < overtones; j++)
			{
				filters[i * overtones + j].initialize(3);
			}

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

	T operator()(T input)
	{
		T sample = 0;
		for (int i = 0; i < voices; i++)
			if (amplitudes[i])
				for (int j = 0; j < overtones; j++)
					sample += amplitudes[i] * pow(decay, j) * filters[i * overtones + j](input) / (2 * voices * normalization);

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
					filters[i * overtones + j].resonant(mtof(particles[i * overtones + j]()), Q);
					filters[i * overtones + j].tick();
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
				// filters[voice * overtones + j].forget();
				// filters[voice * overtones + j].resonant(frequency, Q);
				springs[voice * overtones + j].target(ftom(frequency) - ftom(previous));
			}
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

	void makenote(double pitch, double amplitude)
	{
		pitches[request(mtof(pitch), amplitude)] = pitch;
	}

	void endnote(double pitch)
	{
		for (int j = 0; j < voices; j++)
			if (pitches[j] == pitch)
				release(j);
	}
};

// circular buffer object
template <typename T> class Buffer
{
public:
	~Buffer()
	{
		delete [] data;
	}

	Buffer(uint size = 0)
	{
		size += (size == 0) ? 1 : 0; // disallow size zero
		this->size = size;
		origin = 0;
		// data = vector<T>(size, 0);
		data = new T[size];
		memset(data, 0, size);
	}

	void tick()
	{
		origin++;
		origin %= size;
	}

	// lookup at some past position
	T operator()(uint position = 0)
	{
		return data[(origin - position + size) % size];
	}

	void write(T value)
	{
		data[origin] = value;
	}

	void accum(T value)
	{
		data[origin] += value;
	}


private:
	// vector<T> data;
	T* data;
	uint size;
	uint origin;
};


// many filters in parallel
template <typename T> class Filterbank
{
public:
	Filterbank() { }
	~Filterbank()
	{
		delete [] forwards;
		delete [] backs;
	}

	// initialize N filters of a given order, reading and working with circular buffers
	Filterbank(uint order, uint N = 1) : input(order + 2)
	{
		outputs.reserve(N);
		for (int i = 0; i < N; i++)
			outputs.push_back(Buffer<T>(order + 2));

		forwards = new T[(order + 1) * N];
		backs = new T[(order + 1) * N];

		for (int i = 0; i < order + 1; i++)
			for (int j = 0; j < N; j++)
			{
				forwards[i * N + j] = 0; // jth filter, ith coefficient
				backs[i * N + j] = 0;
			}

		this->order = order;
		this->N = N;
		computed = false;
	}

	// initalize the nth filter's coefficients
	void coefficients(const vector<T>& forward, const vector<T>& back, uint n = 0)
	{
		// assert(order * N + n);

		for (int i = 0; i < order + 1; i++)
		{
			forwards[i * N + n] = forward[i];
			backs[i * N + n] = back[i];
		}
		backs[n] = 0;
	}

	// get the result of filter applied to a sample
	T operator()(T sample)
	{
		if (!computed)
		{
			T out = 0;
			input.write(sample);

			for (int j = 0; j < N; j++)     		// j is filter number
			{
				outputs[j].write(0);
				for (int i = 0; i < order + 1; i++) // i is delay time
					outputs[j].accum(forwards[i * N + j] * input(i) - backs[i * N + j] * outputs[j](i));
				out += outputs[j]();
			}

			computed = true;
			return out;
		}

		T out = 0;
		for (int j = 0; j < N; j++)
			out += outputs[j]();

		return out;
	}

	// timestep
	void tick()
	{
		input.tick();
		for (int j = 0; j < N; j++)
			outputs[j].tick();
		computed = false;
	}

private:
	Buffer<T> input; // circular buffers of inputs and outputs
	vector<Buffer<T>> outputs; 
	uint N; // number of filters
	uint order;

	T* forwards; // feedforward coefficient lists
	T* backs; // feedback coefficient lists

	bool computed; // flag in case of repeated calls to operator()
};


// many delays
template <typename T> class Delay
{
public:
	Delay() { }
	~Delay()
	{
		delete [] forwards;
		delete [] backs;
	}

	// initialize N delays of given sparsity and maximum time
	Delay(uint sparsity, uint time) : input(time), output(time)
	{
		forwards = new pair<uint, T>[sparsity]; 
		backs = new pair<uint, T>[sparsity];

		for (int i = 0; i < sparsity; i++)
		{
			forwards[i] = {0, 0};
			backs[i] = {0, 0};
		}

		this->sparsity = sparsity;
		computed = false;
	}

	// initalize coefficients
	void coefficients(const vector<pair<uint, T>>& forward, const vector<pair<uint, T>>& back)
	{
		uint order = min(sparsity, (uint)forward.size());
		for (int i = 0; i < order; i++)
			forwards[i] = forward[i];
		for (int i = order; i < sparsity; i++)
			forwards[i] = {0, 0}; // zero out trailing coefficients

		order = min(sparsity, (uint)back.size());
		for (int i = 0; i < min(sparsity, (uint)back.size()); i++)
		{
			if (back[i].first != 0) // don't allow zero-time feedback
				backs[i] = back[i];
			else
				backs[i] = {0, 0};
		}
		for (int i = order; i < sparsity; i++)
			backs[i] = {0, 0}; // zero out trailing coefficients
	}

	// modulate the feedforward coeffs of nth delay
	void modulate_forward(uint n, const pair<uint, T>& forward)
	{ forwards[n] = forward; }

	// modulate the feedback coeffs of nth delay
	void modulate_back(uint n, const pair<uint, T>& back)
	{
		if (back.first == 0)
			backs[n] = {0, 0};
		else
			backs[n] = back;
	}

	// get the result of filter applied to a sample
	T operator()(T sample)
	{
		if (!computed)
		{
			input.write(sample);
			output.write(0);

			for (int i = 0; i < sparsity; i++) // i is delay time
			{
				pair<uint, T> forward = forwards[i];
				pair<uint, T> back = backs[i];
				output.accum(forward.second * input(forward.first) - back.second * output(back.first));
			}

			computed = true;
		}

		return output(0);
	}

	// timestep
	void tick()
	{
		input.tick();
		output.tick();
		computed = false;
	}

private:
	Buffer<T> input; // circular buffers of inputs and outputs
	Buffer<T> output; 
	uint sparsity;

	pair<uint, T>* forwards; // feedforward times and coefficients
	pair<uint, T>* backs; // feedback times and coefficients

	bool computed; // flag in case of repeated calls to operator()
};


template <typename T> class FilterType
{
public:
	FilterType() { }
	~FilterType() { }

	static pair<vector<T>, vector<T>> resonant(T gain, T frequency, T modulus)
	{
		T cosine = cos(2 * PI * frequency / SR);

		complex<T> cosine2(cos(4 * PI * frequency / SR), 0);
		complex<T> sine2(sin(4 * PI * frequency / SR), 0);
		
		complex<T> maximum = 1.0 / (modulus - 1) - 1.0 / (modulus - cosine2 - 1.0i * sine2);
		T amplitude = gain / sqrt(abs(maximum));

		vector<T> forward({amplitude, 0, -amplitude});
		vector<T> back({0, -2 * modulus * cosine, modulus * modulus});

		return {forward, back};
	}
};


template <typename T> class Compressor
{
public:
	Compressor() { }
	~Compressor() { }

	Compressor(T threshold, T ratio, T attack, T release, T knee = 0, bool makeup = false)
	{
		this->threshold = threshold;
		this->correction = atodb(ratio);

		up_stiffness = relaxation(attack);
		down_stiffness = relaxation(release);
		amplitude = 0;
	}

	void tick()
	{
		if (amplitude < threshold)
			gain = 0;
		else if (amplitude > threshold)
			gain = correction;

		if (gain > target_gain)
			gain = target_gain * (1 - down_stiffness) + gain * down_stiffness;
		else if (gain < target_gain)
			gain = target_gain * (1 - up_stiffness) + gain * up_stiffness;
	}

	T operator()(T sample)
	{
		amplitude = abs(sample);
		return dbtoa(gain) * sample;
	}

private:
	T gain;
	T target_gain;
	T up_stiffness;
	T down_stiffness;
	T threshold;
	T correction;
	T amplitude;
};


// an n-shot Synth<T>
template <typename T> class Envelope
{
public:
	Envelope() { }
	~Envelope() { }

	Envelope(const vector<function<T(double)>>& functions, Interp interp = Interp::cubic)
	{
		stages = functions.size();
		shapes.reserve(stages);
		for (int i = 0; i < stages; i++)
			shapes.push_back(Wave<T>(functions[i], interp, 0, 1, false));

		stage = 0;
		active = false;
	}

	void trigger(T time, int stage = -1)
	{
		active = true;
		if (stage >= 0)
			this->stage = (uint)stage;
		else
			this->stage++;

		this->stage %= stages;
		phase = 0;
		this->rate = time / SR;
	}

	void tick()
	{
		if (active)
		{
			phase += rate;
			if (phase > 1)
			{
				phase = 1;
				active = false;
			}
		}
	}

	T operator()()
	{
		return shapes[stage](phase);
	}

private:
	vector<Wave<T>> shapes;
	uint stages;
	uint stage;
	T rate;
	T phase;

	bool active;
};

Envelope<double> hann({[] (double phase) -> double { return 0.5 * (1 - cos(2 * PI * phase)); }});


template <typename T> class Granulator
{
public:
	Granulator() { }
	~Granulator() { }

	Granulator(Envelope<T>* window)
	{

	}
};



// receives requests for particle systems; runs physics on them
template <typename T> class Minimizer
{
public:
	Minimizer() { }
	~Minimizer()
	{
		delete [] springs;
		delete [] gravities;
		delete [] guides;
		delete [] particles;		
		delete [] active;
		delete [] pitches;
	}

	Minimizer(uint voices, uint overtones, T decay, T harmonicity) : 
		voices(voices), overtones(overtones), decay(decay), harmonicity(harmonicity)
	{
		particles = new Particle[voices * overtones];
		guides = new Particle[voices];
		springs = new Spring[voices * overtones];
		gravities = new Gravity[voices * (voices - 1) * overtones * overtones / 2];

		active = new T[voices];
		pitches = new T[voices];

		for (int i = 0; i < voices; i++)
			active[i] = 0;

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

	// request a new note, return allocated voice number
	int request(T fundamental, T amplitude = 0)
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
			T frequency = fundamental;
			T previous;
			active[voice] = amplitude;
			guides[voice].initialize(1, ftom(fundamental)); // mass depends on decay / amplitude?
			for (int j = 0; j < overtones; j++)
			{
				previous = frequency;
				frequency = fundamental * pow(j + 1, harmonicity);

				particles[voice * overtones + j].initialize(1, ftom(frequency));
				springs[voice * overtones + j].target(ftom(frequency) - ftom(previous));
			}
		}

		return voice;
	}

	// release a given voice, or all voices
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

	void makenote(T pitch, T amplitude)
	{
		int voice = request(mtof(pitch), amplitude);
		if (voice >= 0)
			pitches[voice] = pitch;
	}

	void endnote(T pitch)
	{
		for (int j = 0; j < voices; j++)
			if (pitches[j] == pitch)
				release(j);
	}

protected:
	uint voices;
	uint overtones;
	T decay;
	T harmonicity;

	Particle* particles; // allows subclasses to access particle positions
	T* active; // and turn voices on / off
	T* pitches;

private:
	Particle* guides;
	Spring* springs;
	Gravity* gravities;
};


template <typename T> class Additive : public Minimizer<T>
{
public:
	using Minimizer<T>::voices, Minimizer<T>::overtones, Minimizer<T>::decay, Minimizer<T>::harmonicity, 
		  Minimizer<T>::particles, Minimizer<T>::active, Minimizer<T>::pitches;
	
	Additive() { }
	~Additive()
	{
		delete [] oscillators;
		delete [] amplitudes;
	}

	Additive(Wave<T>* waveform, uint voices, uint overtones, T decay, T harmonicity = 1.0, T k = 0.1) : 
		Minimizer<T>(voices, overtones, decay, harmonicity), waveform(waveform), attack(relaxation(k))
	{
		normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;
		
		oscillators = new Oscillator<T>[voices * overtones];
		amplitudes = new T[voices];

		for (int i = 0; i < voices; i++)
			amplitudes[i] = 0;

		for (int i = 0; i < voices * overtones; i++)
			oscillators[i] = Oscillator<T>();
	}

	void tick()
	{
		for (int i = 0; i < voices; i++)
			amplitudes[i] = (1 - attack) * active[i] + attack * amplitudes[i];

		// update oscillator frequencies; tick oscillators
		for (int i = 0; i < voices; i++)
			if (active[i] || amplitudes[i])
				for (int j = 0; j < overtones; j++)
				{
					oscillators[i * overtones + j].freqmod(mtof(particles[i * overtones + j]()));
					oscillators[i * overtones + j].tick();
				}
	}

	T operator()()
	{
		T sample = 0;
		for (int i = 0; i < voices; i++)
			if (amplitudes[i])
				for (int j = 0; j < overtones; j++)
					sample += amplitudes[i] * pow(decay, j) * (*waveform)(oscillators[i * overtones + j]()) / (voices * normalization);

		return (T)sample;		
	}

private:
	Wave<T>* waveform;
	Oscillator<T>* oscillators;
	T* amplitudes;
	T attack;
	T normalization;
};

template <typename T> class Subtractive : public Minimizer<T>
{
public:
	using Minimizer<T>::voices, Minimizer<T>::overtones, Minimizer<T>::decay, Minimizer<T>::harmonicity,
		  Minimizer<T>::particles, Minimizer<T>::active, Minimizer<T>::pitches;

	Subtractive() { }
	~Subtractive()
	{
		delete [] amplitudes;
	}

	Subtractive(uint voices, uint overtones, T decay, T harmonicity = 1.0, T k = 0.1) : 
		Minimizer<T>(voices, overtones, decay, harmonicity), attack(relaxation(k))
	{
		normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;
		amplitudes = new T[voices];

		for (int i = 0; i < voices; i++)
			amplitudes[i] = 0;

		filters = Filterbank<T>(2, voices * overtones);
	}

	void tick()
	{
		for (int i = 0; i < voices; i++)
			amplitudes[i] = (1 - attack) * active[i] + attack * amplitudes[i];

		// update and tick filters
		for (int i = 0; i < voices; i++)
			for (int j = 0; j < overtones; j++)
			{
				T frequency = mtof(particles[i * overtones + j]());
				pair<vector<T>, vector<T>> params = FilterType<T>::resonant(pow(decay, j) * amplitudes[i] / (voices * normalization), frequency, 0.99999);
				filters.coefficients(params.first, params.second, i * overtones + j);
			}

		filters.tick();
	}

	T operator()(T sample)
	{
		return filters(sample);
	}

private:
	Filterbank<T> filters;
	T* amplitudes;
	T normalization;
	T attack;
};
