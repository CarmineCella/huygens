// includes.h

#ifndef INCLUDED

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <signal.h>

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

#include <portaudio.h>
#include "RtMidi.h"

using namespace std;

typedef unsigned long ulong;

const double PI = 3.14159265359;
const int SR = 44100;

const static double epsilon = numeric_limits<double>::epsilon();
const static double order = log2(epsilon);

// return the stiffness coefficient so that relaxation
// occurs in k seconds (for interpolation)
double relaxation(double k)
{
	if (k == 0)
		return 0; 
	return pow(2.0, order / (fmax(0, k) * SR));
}

// midi to frequency
double mtof(double midi)
{
	return 440 * pow(2, (midi - 69) / 12);
}

double ftom(double frequency)
{
	return 69 + log2(frequency / 440) * 12;
}

double atodb(double amplitude)
{
	return 20 * log(amplitude);
}

double dbtoa(double db)
{
	return pow(10, db / 20);
}

#endif
#define INCLUDED