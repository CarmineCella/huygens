// delay.h
#include "includes.h"

#ifndef DELAY

namespace soundmath
{
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

			order++; // make room in the ring buffer

			reset();
		}

		void reset()
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
			{
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
			}

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
}

#define DELAY
#endif