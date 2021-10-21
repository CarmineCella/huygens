// sinusoids.h
#include "includes.h"

#ifndef SINUSOIDS

namespace soundmath
{
	// a bank of synthesizers, with optional harmonicity.
	template <typename T> class Sinusoids
	{
	public:
		Sinusoids() { }
		~Sinusoids() { }

		Sinusoids(Wave<T>* form, double fundamental, uint overtones, double decay, double harmonicity = 1, double k = 2.0 / SR)
		{
			this->fundamental = this->target_fundamental = fundamental;
			this->harmonicity = this->target_harmonicity = harmonicity;
			this->overtones = overtones;
			this->decay = this->target_decay = decay;

			for (int i = 0; i < overtones; i++)
			{
				Synth<T> synth(form, fundamental * pow(i + 1, harmonicity));
				synths.push_back(synth);
			}

			normalization = decay != 1 ? (1 - pow(decay, overtones)) / (1 - decay) : overtones;

			stiffness = relaxation(k);
		}

		T operator()()
		{
			T sample = 0;
			for (int i = 0; i < overtones; i++)
			{
				sample += pow(decay, i) * synths[i]() / normalization;
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
				synths[i].freqmod(fundamental * pow(i + 1, harmonicity));
				synths[i].tick();
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
		vector<Synth<T>> synths;
	};
}

#define SINUSOIDS
#endif