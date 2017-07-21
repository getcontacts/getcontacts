// Author: Anthony Kai Kwang Ma
// Email: akma327@stanford.edu
// Date: February 18, 2016

// Process Noncovalent Interaction by Time Frame for MOR Trajectories 
// C++ implementation that will use a binary search tree data structure 
// In the case of water mediated hydrogen bonds, where the dictionary size
// Is very large - Anthony Ma 02/22/16

// Usage:
// ./processTimeFrameData <Stitch Interaction Input Path> <Input File Name> <Binary Dictionary Output Path> <Output File Name> <interaction type><binarize flag>
// Example:
// StitchedInteractionPath="/scratch/PI/rondror/akma327/noncovalent_Interaction_Scripts/mOR-InteractionOutput/MOR_active_waters/mor_active_refine150_agonist_noNanobody_noNTerm/rep_2"
// INPUT_FILE_NAME="salt_bridge_result_stitch.txt"
// OUTPUT_PATH="/scratch/PI/rondror/akma327/noncovalent_Interaction_Scripts/mOR-InteractionOutputDictionary/MOR_active_waters/mor_active_refine150_agonist_noNanobody_noNTerm/rep_2"
// OUTPUT_FILE_NAME="salt_bridge_result_dict.txt"
// ./processTimeFrameData $StitchedInteractionPath $INPUT_FILE_NAME $OUTPUT_PATH $OUTPUT_FILE_NAME -sb -b

#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <time.h>
using namespace std;

std::map<string, string> FLAG_INTERACTION;

int stringToInt(string s){
	int x;
	stringstream convert(s);
	convert >> x;
	return x;
}

// C++ Implementation of python's split implementation that converts a string
// to a vector containing the tokens split upon a desired delimeter
vector<string> split(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }
  return internal;
}


int getChainId(string line){
	string atom1 = split(line, '-')[1];
	if(atom1.find("_") != std::string::npos){
		return stringToInt(split(atom1, '_')[1]);
	}
	return 0;
}

string retrieveDataFileReader(char* argv[]){
	cout << "retrieveDataFileReader()" <<endl;
	string input_path = argv[1];
	string input_filename = argv[2];
	string filename = input_path + "/" + input_filename;
	cout << "Creating Binary Dictionary For: " << filename <<endl;
	cout << "Final filename: " << filename << endl;
	return filename;
}

map<int, std::map<string, std::vector<int> > > createTimeFrameDict(char* argv[], int &stride, std::vector<string> &TrajectoryPaths, string &TopologyPath, int &totalFrames){
	cout << "createTimeFrameDict()" <<endl;
	int numLinesRead = 0;
	std::map<int, std::map<string, std::vector<int> > > timelapse_interaction_by_chain_dict;
	string filename = retrieveDataFileReader(argv);
	ifstream f; 
	f.open(filename.c_str());
	if(!f){
		cout << "Interaction Results Not Found" << endl;
	}else{
		int currFrame = 0;
		string currProtein = "";
		std::string line;
		while(std::getline(f, line)){
			numLinesRead ++;
			// if(numLinesRead >1000000){
			// 	break;
			// } 
			if(numLinesRead %10000 == 0){
				cout << "Line: " << numLinesRead << " ** " << line << endl;
			}
			if(line.find("Stride:") != std::string::npos){
				stride = stringToInt(split(line, ':')[1]);
			}else if(line.find("TrajectoryPath:") != std::string::npos){
				TrajectoryPaths.push_back(split(line, ':')[1]);
			}else if(line.find("TopologyPath:") != std::string::npos){
				TopologyPath = split(line, ':')[1];
			}else if(line.find("nFrames:") != std::string::npos){
				continue;
			}else if(line.find("Computing Time:") != std::string::npos){
				continue;
			}else if(line.find("Frame:") != std::string::npos){
				currFrame = stringToInt(split(line, ':')[1]);
			}else if(line.find("Salt Bridges:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Pi-Cation:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Pi-Stacking:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("T-Stacking:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Van Der Waals:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Hydrogen Bond-Water Mediated:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Backbone Backbone Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Sidechain Backbond Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Sidechain Sidechain Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Residue Water Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Water Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Extended Water Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}
			// Ligand Based Interactions
			else if(line.find("Ligand Backbone Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Ligand Sidechain Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Ligand Water Hydrogen Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Ligand Water Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}else if(line.find("Ligand Extended Water Bonds:") != std::string::npos){
				currProtein = split(line, ':')[1];
			}
			else if(!line.empty()){
				string key = line;
				int chain_id = getChainId(line);
				if(timelapse_interaction_by_chain_dict.find(chain_id) == timelapse_interaction_by_chain_dict.end()){
					std::map<string, std::vector<int> > internalDict; //interaction keys to value = [1,0,0,0,1...]
					std::vector<int> presentTimeFrames;
					presentTimeFrames.push_back(currFrame);
					internalDict.insert(std::pair<string, std::vector<int> >(key, presentTimeFrames));
					timelapse_interaction_by_chain_dict.insert(std::pair<int, std::map<string, std::vector<int> > >(chain_id, internalDict));
				}else{
					if(timelapse_interaction_by_chain_dict[chain_id].find(key) == timelapse_interaction_by_chain_dict[chain_id].end()){
						std::vector<int> presentTimeFrames;
						presentTimeFrames.push_back(currFrame);
						timelapse_interaction_by_chain_dict[chain_id].insert(std::pair<string, std::vector<int > >(key, presentTimeFrames));
					}else{
						timelapse_interaction_by_chain_dict[chain_id][key].push_back(currFrame);
					}
				}
			}
		}
		totalFrames = currFrame + 1;
	}
	return timelapse_interaction_by_chain_dict;
}

string writeRootPath(string stitchedInteractionPath){
	size_t len = strlen("mOR-InteractionOutput");
	size_t found = stitchedInteractionPath.find("mOR-InteractionOutput");
	if(found != std::string::npos){
		stitchedInteractionPath.replace(found, len, "mOR-InteractionOutputDictionary");
	}
	return stitchedInteractionPath;
}


void writeFrequencyMap(char*argv[], int stride, std::vector<string> &TrajectoryPaths, string &TopologyPath, 
						std::map<int,std::map<string, std::vector<int> > > & timelapse_interaction_by_chain_dict, int &totalFrames){
	cout << "writeFrequencyMap()" << endl;
	int numKeysWritten = 0;
	string output_path = argv[3];
	string output_filename = argv[4];
	string interaction_selection = argv[5];
	string filename = output_path + "/" + output_filename;
	string mkcmd = "mkdir " + output_path;
	system(mkcmd.c_str());
	string buffer;
	ofstream f;
	f.open(filename.c_str());
	f << "Stride:" << stride << endl;
	for(int i = 0; i<TrajectoryPaths.size(); i++){
		string TrajectoryPath = TrajectoryPaths[i];
		f <<"TrajectoryPath:" << TrajectoryPath <<endl;
	}
	f << "TopologyPath:" << TopologyPath << endl;
	f << FLAG_INTERACTION[interaction_selection] << " Dictionary Heat Map " <<endl;
	f << "TotalFrames: " << totalFrames << endl;
	for(std::map<int, std::map<string, std::vector<int> > >::iterator it1 = timelapse_interaction_by_chain_dict.begin(); it1 != timelapse_interaction_by_chain_dict.end(); ++it1){
		int chain_id = it1->first;
		f << "Dictionary for Chain: " << chain_id << endl;
		std::map<string, std::vector<int> > timeFrameDict = it1->second;
		for(std::map<string, std::vector<int> >::iterator it2 = timeFrameDict.begin(); it2 != timeFrameDict.end(); ++it2){
			if(numKeysWritten % 1000 == 0){
				cout <<"Key #: " << numKeysWritten << endl;
			}
			string interaction_key = it2->first;
			f << interaction_key << "~[";
			std::vector<int> bitVector = it2->second;
			int bitVecSize = bitVector.size();
			for (int i = 0; i<bitVecSize; i++){
				if(i < bitVecSize-1){
					f << bitVector[i] << ",";
				}
			}
			f<<bitVector[bitVecSize-1] << "]" <<endl;
			numKeysWritten++;
		}
	}
	cout << "Wrote results to: " << filename << endl;
}

//Should probably write directly to file rather than saving in intermediary map. Memory storage too high 
std::vector<int> binarizeTime(std::vector<int> timeFrames, int totalFrames){
	vector<int> binaryVec;
	for (int i = 0; i<totalFrames; i++){
		binaryVec.push_back(0);
	}
	for (int j = 0; j<timeFrames.size(); j++){
		int frameIndex = timeFrames[j];
		binaryVec[frameIndex] = 1;
	}
	return binaryVec;
}


void binarizeAndWriteTimePoints(char*argv[], int stride, std::vector<string> &TrajectoryPaths, string &TopologyPath, 
																				std::map<int,std::map<string, std::vector<int> > > timelapse_interaction_by_chain_dict, 
																				int totalFrames){
	cout << "binarizeTimePoints()" << endl;
	int numKeysWritten = 0;
	string output_path = argv[3];
	string output_filename = argv[4];
	string interaction_selection = argv[5];
	string filename = output_path + "/" + output_filename;
	string buffer;
	ofstream f;
	f.open(filename.c_str());
	f << "Stride:" << stride << endl;
	for(int i = 0; i<TrajectoryPaths.size(); i++){
		string TrajectoryPath = TrajectoryPaths[i];
		f <<"TrajectoryPath:" << TrajectoryPath <<endl;
	}
	f << "TopologyPath:" << TopologyPath << endl;
	f << FLAG_INTERACTION[interaction_selection] << " Dictionary Heat Map " <<endl;
	f << "TotalFrames: " << totalFrames << endl;
	for(std::map<int, std::map<string, std::vector<int> > >::iterator it1 = timelapse_interaction_by_chain_dict.begin(); it1 != timelapse_interaction_by_chain_dict.end(); ++it1){
		int chain_id = it1->first;
		f << "Dictionary for Chain: " << chain_id << endl;
		std::map<string, std::vector<int> > timeFrameDict = it1->second;
		for(std::map<string, std::vector<int> >::iterator it2 = timeFrameDict.begin(); it2 != timeFrameDict.end(); ++it2){
			if(numKeysWritten % 1000 == 0){
				cout <<"Key #: " << numKeysWritten << endl;
			}
			string interaction_key = it2->first;
			std::vector<int> timeFrames = it2->second;
			std::vector<int> binaryVec = binarizeTime(timeFrames, totalFrames);
			f << interaction_key << "~[";
			int bitVecSize = binaryVec.size();
			for(int i = 0; i<bitVecSize;i++){
				if(i<bitVecSize-1){
					f << binaryVec[i] << ",";
				}
			}
			f<<binaryVec[bitVecSize-1] <<"]" <<endl;
			numKeysWritten ++;
		}
	}
	cout << "Wrote results to: " << filename << endl;
}


void initialize_flags(){
	FLAG_INTERACTION.insert(std::pair<string, string> ("-sb", "Salt Bridges"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-pc", "Pi-Cation"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-ps", "Pi-Stacking"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-ts", "T_Stacking"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-vdw", "Van Der Waals"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-hb", "Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-hbw", "Hydrogen Bond-Water Mediated"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-hbbb", "Backbone Backbone Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-hbsb", "Sidechain Backbond Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-hbss", "Sidechain Sidechain Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-rw", "Residue Water Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-wb", "Water Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-wb2", "Extended Water Bonds"));
	// Ligand Based Interactions
	FLAG_INTERACTION.insert(std::pair<string, string> ("-hlb", "Ligand Backbone Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-hls", "Ligand Sidechain Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-lw", "Ligand Water Hydrogen Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-lwb", "Ligand Water Bonds"));
	FLAG_INTERACTION.insert(std::pair<string, string> ("-lwb2", "Ligand Extended Water Bonds"));

}

int main(int argc, char* argv[]){
	initialize_flags();
	int stride = 1;
	std::vector<string> TrajectoryPaths;
	string TopologyPath = "";
	int totalFrames = 0;
	std::map<int,std::map<string, std::vector<int> > > timelapse_interaction_by_chain_dict =  createTimeFrameDict(argv, stride, TrajectoryPaths, TopologyPath, totalFrames);
	cout<< "Number of Distinct Keys: " << timelapse_interaction_by_chain_dict[0].size() << endl;
	if(argc > 6){
		string binarizeFlag = argv[6];
		if(binarizeFlag.find("-b") != std::string::npos){
			binarizeAndWriteTimePoints(argv, stride, TrajectoryPaths, TopologyPath, timelapse_interaction_by_chain_dict, totalFrames);
		}
	}else{
		writeFrequencyMap(argv, stride, TrajectoryPaths, TopologyPath, timelapse_interaction_by_chain_dict, totalFrames);
	}
	return 0;
}

