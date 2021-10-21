// buffer.h
#include "includes.h"

#ifndef BUFFER

namespace soundmath
{
	// circular buffer object
	template <typename T> class Buffer
	{
	public:
		Buffer() { }
		~Buffer() { }

		Buffer(uint size)
		{
			this->size = size;
			origin = 0;
			data = vector<T>(size, 0);
		}

		void tick()
		{
			origin++;
			origin %= size;
		}

		// lookup at some past position
		T operator()(uint position)
		{
			return data[(origin - position + size) % size];
		}

		void write(T value)
		{
			data[origin] = value;
		}


	private:
		vector<T> data;
		uint size;
		uint origin;
	};
}

#define BUFFER
#endif