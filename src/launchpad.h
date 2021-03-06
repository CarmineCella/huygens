#include "midi.h"

static int notes[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
					   0,  1,  2,  3,  4,  5,  6,  7, -1, -1,
					   5,  6,  7,  8,  9, 10, 11, 12, -1, -1,
					  10, 11, 12, 13, 14, 15, 16, 17, -1, -1,
					  15, 16, 17, 18, 19, 20, 21, 22, -1, -1,
					  20, 21, 22, 23, 24, 25, 26, 27, -1, -1,
					  25, 26, 27, 28, 29, 30, 31, 32, -1, -1,
					  30, 31, 32, 33, 34, 35, 36, 37, -1, -1,
					  35, 36, 37, 38, 39, 40, 41, 42, -1, -1,
					  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
					  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
					  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
					  -1, -1, -1, -1, -1, -1, -1 };

class Launchpad
{
public:
	Launchpad()
	{
	}

	int note(unsigned char message)
	{
		return notes[message];
	}

};