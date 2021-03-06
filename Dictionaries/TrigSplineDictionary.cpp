#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#define sqr(x) ((x)*(x))

#include "Dictionary.cpp"

class TrigSplineDictionary : public Dictionary {
	
	public:
		TrigSplineDictionary(int atomsCount, int signalSize, int testSize, int rightBorder) : Dictionary(atomsCount, signalSize, testSize, rightBorder) {}

		void CreateDictionary(double* atoms, double* tests) override {
			int counter = 0;
			for(int iter=4;iter<=atomsCount;iter=iter*2){
				double delta_x = double(rightBorder)/(iter - 2.0);
				double xi[iter+4];
				xi[0] = -2.0*delta_x;
				for(int i=0;i<iter+3;i++){
					//std::cout << xi[i] << ' ';
					xi[i+1]=xi[i]+delta_x;
				}
				//std::cout << xi[atomsCount+3] << std::endl << std::endl;
				for(int i=0;i<iter;counter++,i++){
					for(int j=0;j<signalSize;j++){
						double cur_position = (double(rightBorder)/(signalSize-1)*j);
						//std::cout << cur_position << ' ';
						if((cur_position >= xi[i]) and (cur_position < xi[i+1])){
							atoms[counter*signalSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((cur_position - xi[i])/2)))/((sin((xi[i+2]-xi[i])/2))*(sin((xi[i+1]-xi[i])/2)));
							//std::cout << 0 << ' ' << 1 << ' ' << cur_position << ' ' << xi[i] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+1] << std::endl;
						}
						else if((cur_position >= xi[i+1]) and (cur_position < xi[i+2])){
							atoms[counter*signalSize+j] = ((cos((xi[i+2]-xi[i+1])/2))/(sin((xi[i+1]-xi[i])/2)))*(((sqr(sin((cur_position-xi[i])/2)))/(sin((xi[i+2]-xi[i])/2)))-((sin((xi[i+3]-xi[i])/2))*(sqr(sin((cur_position-xi[i+1])/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+2]-xi[i+1])/2)))));
							//std::cout << 1 << ' ' << 2 << ' ' << cur_position << ' ' << xi[i+1] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+2] << std::endl;
						}
						else if((cur_position >= xi[i+2]) and (cur_position < xi[i+3])){
							atoms[counter*signalSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((xi[i+3] - cur_position)/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+3]-xi[i+2])/2)));
							//std::cout << 2 << ' ' << 3 << ' ' << cur_position << ' ' << xi[i+2] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+3] << std::endl;
						}
						else{
							atoms[counter*signalSize+j]=0.0;
						}
						//std::cout << atoms[i*signalSize+j] << ' ';
					}
					//std::cout << std::endl;
					for(int j=0;j<testSize;j++){
						double cur_position = (double(rightBorder)/(testSize-1)*j);
						//std::cout << cur_position << ' ';
						if((cur_position >= xi[i]) and (cur_position < xi[i+1])){
							tests[counter*testSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((cur_position - xi[i])/2)))/((sin((xi[i+2]-xi[i])/2))*(sin((xi[i+1]-xi[i])/2)));
							//std::cout << 0 << ' ' << 1 << ' ' << cur_position << ' ' << xi[i] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+1] << std::endl;
						}
						else if((cur_position >= xi[i+1]) and (cur_position < xi[i+2])){
							tests[counter*testSize+j] = ((cos((xi[i+2]-xi[i+1])/2))/(sin((xi[i+1]-xi[i])/2)))*(((sqr(sin((cur_position-xi[i])/2)))/(sin((xi[i+2]-xi[i])/2)))-((sin((xi[i+3]-xi[i])/2))*(sqr(sin((cur_position-xi[i+1])/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+2]-xi[i+1])/2)))));
							//std::cout << 1 << ' ' << 2 << ' ' << cur_position << ' ' << xi[i+1] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+2] << std::endl;
						}
						else if((cur_position >= xi[i+2]) and (cur_position < xi[i+3])){
							tests[counter*testSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((xi[i+3] - cur_position)/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+3]-xi[i+2])/2)));
							//std::cout << 2 << ' ' << 3 << ' ' << cur_position << ' ' << xi[i+2] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+3] << std::endl;
						}
						else{
							tests[counter*testSize+j]=0.0;
						}
						//std::cout << atoms[i*signalSize+j] << ' ';
					}
					//std::cout << std::endl;
				}
			}
			int iter = atomsCount;
			double delta_x = double(rightBorder)/(iter - 2.0);
			double xi[iter+4];
			xi[0] = -1.0*delta_x;
			for(int i=0;i<iter+3;i++){
				//std::cout << xi[i] << ' ';
				xi[i+1]=xi[i]+delta_x;
			}
			for(int i=0;i<iter;counter++,i++){
				for(int j=0;j<signalSize;j++){
					double cur_position = (double(rightBorder)/(signalSize-1)*j);
					//std::cout << cur_position << ' ';
					if((cur_position >= xi[i]) and (cur_position < xi[i+1])){
						atoms[counter*signalSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((cur_position - xi[i])/2)))/((sin((xi[i+2]-xi[i])/2))*(sin((xi[i+1]-xi[i])/2)));
						//std::cout << 0 << ' ' << 1 << ' ' << cur_position << ' ' << xi[i] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+1] << std::endl;
					}
					else if((cur_position >= xi[i+1]) and (cur_position < xi[i+2])){
						atoms[counter*signalSize+j] = ((cos((xi[i+2]-xi[i+1])/2))/(sin((xi[i+1]-xi[i])/2)))*(((sqr(sin((cur_position-xi[i])/2)))/(sin((xi[i+2]-xi[i])/2)))-((sin((xi[i+3]-xi[i])/2))*(sqr(sin((cur_position-xi[i+1])/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+2]-xi[i+1])/2)))));
						//std::cout << 1 << ' ' << 2 << ' ' << cur_position << ' ' << xi[i+1] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+2] << std::endl;
					}
					else if((cur_position >= xi[i+2]) and (cur_position < xi[i+3])){
						atoms[counter*signalSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((xi[i+3] - cur_position)/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+3]-xi[i+2])/2)));
						//std::cout << 2 << ' ' << 3 << ' ' << cur_position << ' ' << xi[i+2] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+3] << std::endl;
					}
					else{
						atoms[counter*signalSize+j]=0.0;
					}
					//std::cout << atoms[i*signalSize+j] << ' ';
				}
				//std::cout << std::endl;
				for(int j=0;j<testSize;j++){
					double cur_position = (double(rightBorder)/(testSize-1)*j);
					//std::cout << cur_position << ' ';
					if((cur_position >= xi[i]) and (cur_position < xi[i+1])){
						tests[counter*testSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((cur_position - xi[i])/2)))/((sin((xi[i+2]-xi[i])/2))*(sin((xi[i+1]-xi[i])/2)));
						//std::cout << 0 << ' ' << 1 << ' ' << cur_position << ' ' << xi[i] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+1] << std::endl;
					}
					else if((cur_position >= xi[i+1]) and (cur_position < xi[i+2])){
						tests[counter*testSize+j] = ((cos((xi[i+2]-xi[i+1])/2))/(sin((xi[i+1]-xi[i])/2)))*(((sqr(sin((cur_position-xi[i])/2)))/(sin((xi[i+2]-xi[i])/2)))-((sin((xi[i+3]-xi[i])/2))*(sqr(sin((cur_position-xi[i+1])/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+2]-xi[i+1])/2)))));
						//std::cout << 1 << ' ' << 2 << ' ' << cur_position << ' ' << xi[i+1] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+2] << std::endl;
					}
					else if((cur_position >= xi[i+2]) and (cur_position < xi[i+3])){
						tests[counter*testSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((xi[i+3] - cur_position)/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+3]-xi[i+2])/2)));
						//std::cout << 2 << ' ' << 3 << ' ' << cur_position << ' ' << xi[i+2] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+3] << std::endl;
					}
					else{
						tests[counter*testSize+j]=0.0;
					}
					//std::cout << atoms[i*signalSize+j] << ' ';
				}
				//std::cout << std::endl;
			}
			return;
		}
};
