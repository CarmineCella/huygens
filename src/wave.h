// wave.h
#include "includes.h"

#ifndef WAVE

namespace soundmath
{
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
}

#define WAVE
#endif