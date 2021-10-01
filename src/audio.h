// audio.h

#include <portaudio.h>
#include <iostream>

const int bsize = 64;

int callback(const void*, void*, unsigned long, const PaStreamCallbackTimeInfo*, PaStreamCallbackFlags, void*);

class Audio
{
public:
	PaStream* stream;
	int (*process)(float*, float*);

	Audio(int (*processor)(float* in, float* out))
	{
		process = processor;
	}

	static PaError clean(PaError err)
	{
		if (err != paNoError)
		{
			Pa_Terminate();
			std::cout << "PortAudio error: " << Pa_GetErrorText(err);
		}
		return err;
	}

	void startup()
	{
		clean(Pa_Initialize());
		clean(Pa_OpenDefaultStream(&stream, 0, 2, paFloat32, SR, bsize, callback, this));
		clean(Pa_StartStream(stream));
	}

	void shutdown()
	{
		clean(Pa_StopStream(stream));
		clean(Pa_CloseStream(stream));
		clean(Pa_Terminate());
	}
};

int callback(const void* inputBuffer, void* outputBuffer, 
			  unsigned long framesPerBuffer, 
			  const PaStreamCallbackTimeInfo* timeInfo, 
			  PaStreamCallbackFlags statusFlags,
			  void* userData)
{
	Audio* A = (Audio*)userData;
	A->process((float*) inputBuffer, (float*) outputBuffer);
	return 0;
}
