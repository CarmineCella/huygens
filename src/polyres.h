// polyres.h
#include "includes.h"

#ifndef POLYRES

namespace soundmath
{
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
	};
}

#define POLYRES
#endif