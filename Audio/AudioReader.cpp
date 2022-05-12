#include <iostream>
#include "AudioFile.h"

#pragma once
void audioReader(double *vSignal, double *rSignal, int szSignal, int szTest, char *signal)
{
    // 1. Set a file path to an audio file on your machine
    std::string inputFilePath = std::string(signal);
    
    // 2. Create an AudioFile object and load the audio file
    AudioFile<float> a;
    bool loadedOK = a.load (inputFilePath);
    
    /** If you hit this assert then the file path above
     probably doesn't refer to a valid audio file */
    assert (loadedOK);
	
    for (int i = 0; i < szSignal; i++){
        vSignal[i] = a.samples[0][i*10];
    }
	
    for (int i = 0; i < szTest; i++){
        rSignal[i] = a.samples[0][i];
    }
    
    return;
}
