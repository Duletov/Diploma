#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>

#include "Audio/AudioReader.cpp"

#include "Dictionaries/Dictionary.cpp"
#include "Dictionaries/GaborDictionary.cpp"
#include "Dictionaries/DctDictionary.cpp"
#include "Dictionaries/SplineDictionary.cpp"
#include "Dictionaries/HypSplineDictionary.cpp"
#include "Dictionaries/TrigSplineDictionary.cpp"
#include "Dictionaries/MaxTrigSplineDictionary.cpp"
#include "Dictionaries/MaxSplineDictionary.cpp"
#include "Dictionaries/MaxHypSplineDictionary.cpp"

#include "Algorithms/Algorithm.cpp"
#include "Algorithms/OMP.cpp"
#include "Algorithms/OMP_str.cpp"
#include "Algorithms/Lasso.cpp"
#include "Algorithms/StOMP.cpp"
#include "Algorithms/SOMP.cpp"

int main(int argc, char** argv) {
	
	if (!(argc == 9))	{
		return -1;
	}
	
	int nAtoms =  std::stoi(argv[1]);
	int szSignal = std::stoi(argv[2]);
	int szTest = std::stoi(argv[3]);
	int rightBorder = std::stoi(argv[4]);
	char dictType = argv[5][0];
	char algoType = argv[6][0];
	char signalType = argv[7][0];
	char *signal = argv[8];
	
	int temp=0;
	if(dictType == 's' or dictType == 't' or dictType == 'h'){
		for(int i=4; i<nAtoms;i=i*2)
			temp+=i;
	}
	
	
	double *vSignal = new double[szSignal];
	double *rSignal = new double[szTest];
	double *mDictionary = new double[(nAtoms + temp) * szSignal];
	double *fullDictionary = new double[(nAtoms + temp) * szTest];
	
	switch(dictType){
		case 'd':{
			DctDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary, fullDictionary);
			break;
		}
		case 'g':{
			GaborDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary, fullDictionary);
			break;
		}
		case 's':{
			SplineDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary, fullDictionary);
			break;
		}
		case 'h':{
			HypSplineDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary, fullDictionary);
			break;
		}
		case 't':{
			TrigSplineDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary, fullDictionary);
			break;
		}
		case 'a':{
			MaxSplineDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary,fullDictionary);
			break;
		}
		case 'b':{
			MaxTrigSplineDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary,fullDictionary);
			break;
		}
		case 'c':{
			MaxHypSplineDictionary generator(nAtoms, szSignal, szTest, rightBorder);
			generator.CreateDictionary(mDictionary,fullDictionary);
			break;
		}
		default:{
			std::cout << "dictType is not defined";
			return 1;
		}
	}
	
	
	nAtoms = nAtoms + temp;
	
	if(signalType == 'a'){
		audioReader(vSignal, rSignal, szSignal, szTest, signal);
	}
	else{
		if(signalType == 's'){
			switch(signal[0]){
				case 'x':{
					for(int i=0;i<szSignal;i++){
						vSignal[i] = (i*1.0/szSignal)*(i*1.0/szSignal);
					}
					for(int i=0;i<szTest;i++){
						rSignal[i] = (i*1.0/szTest)*(i*1.0/szTest);
					}
					break;
				}
				case 'c':{
					for(int i=0;i<szSignal;i++){
						vSignal[i] = 0.2*cosh(i*5.0/szSignal);
					}
					for(int i=0;i<szTest;i++){
						rSignal[i] = 0.2*cosh(i*5.0/szTest);
					}
					break;
				}
				case 'a':{
					for(int i=0;i<szSignal;i++){
						vSignal[i] = atan(i * 10.0 / szSignal);
					}
					for(int i=0;i<szTest;i++){
						rSignal[i] = atan(i*10.0/szTest);
					}
					break;
				}
				default:{
					std::cout << "signal is not defined";
					return 1;
				}
			}
		}
	}
	
	auto start = clock();
	
	switch(algoType){
		case 'o':{
			OMPAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 'l':{
			LassoAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 's':{
			StOMPAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 't':{
			OMPsAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 'p':{
			SplineOMPAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		default:{
			std::cout << "dictType is not defined";
			return 1;
		}
	}
	auto end = clock();
	
	std::cout << "\n" << (end - start) / CLOCKS_PER_SEC << " seconds";
	
	delete [] mDictionary;
	delete [] fullDictionary;
	return 0;
}
