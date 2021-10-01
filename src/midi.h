// midi.h

#include "RtMidi.h"
#include <iostream>
#include <vector>

using namespace std;


class MidiIn
{
private:
	RtMidiIn *midi = 0;
	vector<string> names;
	int port = -1;

public:
	MidiIn()
	{

	}

	void startup()
	{
		try { midi = new RtMidiIn(); }
		catch (RtMidiError &error)
		{
			error.printMessage();
			exit(EXIT_FAILURE);
		}
	}

	void shutdown()
	{
		delete midi;
	}

	void getports(bool verbose = false)
	{
		names.clear();
		int n = midi->getPortCount();
		string name;
		
		for (int i = 0; i < n; i++)
		{
			try { name = midi->getPortName(i); }
			catch (RtMidiError &error)
			{
				error.printMessage();
				continue;
			}
			if (verbose)
				cout << "  input port #" << i+1 << ": " << name << '\n';
			names.push_back(name);
		}
	}

	int open(int number)
	{
		try
		{
			midi->openPort(number);
		}
		catch (RtMidiError &error)
		{
			error.printMessage();
			return 1;
		}
		port = number;
		return 0;
	}


	int open(string name)
	{
		for (int i = 0; i < names.size(); i++)
			if (names[i] == name)
				return open(i);
		return 1;
	}

	void ignore(bool sysex = true, bool time = true, bool sense = true)
	{
		midi->ignoreTypes(sysex, time, sense);
	}

	double get(vector<unsigned char>* message)
	{
		return midi->getMessage(message);
	}

};


class MidiOut
{
private:
	RtMidiOut *midi = 0;
	vector<string> names;
	int port = -1;

public:
	MidiOut()
	{

	}

	void startup()
	{
		try { midi = new RtMidiOut(); }
		catch (RtMidiError &error)
		{
			error.printMessage();
			exit(EXIT_FAILURE);
		}
	}

	void shutdown()
	{
		delete midi;
	}

	void getports(bool verbose = false)
	{
		names.clear();
		int n = midi->getPortCount();
		string name;
		
		for (int i = 0; i < n; i++)
		{
			try { name = midi->getPortName(i); }
			catch (RtMidiError &error)
			{
				error.printMessage();
				continue;
			}
			if (verbose)
				cout << "  output port #" << i+1 << ": " << name << '\n';
			names.push_back(name);
		}
	}

	int open(int number)
	{
		try
		{
			midi->openPort(number);
		}
		catch (RtMidiError &error)
		{
			error.printMessage();
			return 1;
		}
		port = number;
		return 0;
	}


	int open(string name)
	{
		for (int i = 0; i < names.size(); i++)
			if (names[i] == name)
				return open(i);
		return 1;
	}
};
