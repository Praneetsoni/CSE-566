#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <unordered_set>

#include "rlcsa.h"

using namespace std;
using namespace CSA;
//usint max_string =0;
string test;
typedef unsigned char uchar;
typedef unsigned long usint;
typedef struct{
	string header;
	string seq;
} sequence_t;

vector<sequence_t> seq_list;

int get_data(const string& contigFileName, 
	uchar*& data, 
	vector<sequence_t>& seq_list,
	uint64_t& char_count
	)
{
	try 
	{
		ifstream contigFile;
		contigFile.open(contigFileName);
	   	char_count = 0;
	   	string line;
		
		sequence_t temp;
	   	string genome;
	   	// counting total number of reads and their length
	   	while (getline(contigFile , line))
	   	{
	   		line = line.substr(2,line.length() - 2); // removing the leading ">m"	   			   		
	   		temp.header = line;
			
	   		getline(contigFile, line); // this is the actual sequence
			char_count += line.length() + 1; // +1 for the '\0' separator
			temp.seq = line;
			test = line;
	   		seq_list.push_back(temp);
	   	}

	   	contigFile.close();

	   	cout << "Input file " << contigFileName << " contains " << seq_list.size() << " sequences." << endl;
	   	cout << "The temporary data array will have size " << char_count/(double)1000000 << "MB." << endl;

	   	uint64_t i = 0;
		data = new uchar[char_count];
		for (auto temp : seq_list)
		{
	    	for (uint64_t j = 0; j < temp.seq.length(); j++)
	    	{
	    		data[i] = (uchar)temp.seq[j];
	    		i++;
	    		
	    	}
	    	data[i] = '\0';
	    	i++;		
		}
	   	cout << "Created the data array" << endl;

	   	return EXIT_SUCCESS;
	} catch (exception& error)
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

}
string extractAccession(usint idx){
	string s = seq_list[idx].header;
	return s.substr(2, s.find(" ")-2);
	
}
void create_output(RLCSA& rlcsa){
	fstream fs;
	map<usint, vector<pair<usint, usint> > > cover;
	
	usint index;
	string pattern;
	usint *pos;
	char track[rlcsa.getNumberOfSequences()];
	fs.open("redundant.csv",fstream::out);
	
	for(usint j=0; j<rlcsa.getNumberOfSequences() ;j++){
		//string pattern((char *)rlcsa.display(j));
		
		pattern = seq_list[j].seq;
		pair_type result = rlcsa.count(pattern);
		if(length(result)>1){
			
			usint check =1,first_occ = 0;
			pos = rlcsa.locate(result);
			first_occ = rlcsa.getSequenceForPosition(pos[0]);
			for(usint p=0; p<length(result); p++){
				index = rlcsa.getSequenceForPosition(pos[p]);
				if(index != j){
					fs <<extractAccession(j)<< ",";
					fs << extractAccession(index)<<endl;
					
					cover[index].push_back(make_pair(pos[p], pattern.size()));
				}	
				
				//find exact matches
				//cout<<"sequecne range "<<length(rlcsa.getSequenceRange(index)) << "string length " << pattern.size()<<endl; 
				if(length(rlcsa.getSequenceRange(index)) != pattern.size() ){
					//cout<<"size not match" << rlcsa.display(index) << "    :   "<<pattern << endl;
					check = 0;
				}
			}
			
			track[j] = 'n';
			// if check ==1 then all matches are perfect and we need to include the first of them in the list
			if(check == 1){
				track[first_occ] = 'y';
			}
			
			/*if(strcmp((seq_list[rlcsa.getSequenceForPosition(pos[1])].seq).c_str(), (char *)rlcsa.display(rlcsa.getSequenceForPosition(pos[1])) )){
				cout<<"Not a match!!!"<< endl;
				cout<< "position number "<< rlcsa.getSequenceForPosition(pos[1]) <<endl;
				cout<< rlcsa.display(rlcsa.getSequenceForPosition(pos[1])) <<endl;
			}*/
		}
		else{
			track[j] = 'y';
			
		}
		pattern.clear();
	}
	fs.close();
	

	//Create non redundant csv
	fs.open("Non_redundant.fasta",fstream::out);
	int cnt = 0;
	for(usint i =0; i< rlcsa.getNumberOfSequences(); i++){
		
		if(track[i] == 'y'){
			cnt++;
			fs << ">m" << seq_list[i].header<<endl;
			fs << seq_list[i].seq<<endl;
		}
	}
	cout<<"The count of non-redundant string si "<<cnt<<endl;
	
	fs.close();
	
	//Create Coverage file
	fs.open("Coverage.csv", fstream::out);
	cnt =0;
	//for(usint i =0; i< rlcsa.getNumberOfSequences(); i++){
	for (map<usint, vector<pair<usint, usint>>>::iterator it=cover.begin(); it!=cover.end(); ++it){	
		usint i = it->first;
		
		if(track[i]=='y'){
			cnt++;
			uchar* superstr = rlcsa.display(i); 
			usint size = (it->second).size();
			
			pair_type info = rlcsa.getSequenceRange(i);
			usint *count = new usint[length(info)];
			bzero(count,length(info));
			for(usint j=0; j< size ;j++){
				pair<usint,usint> curr = (it->second)[j];
				
				usint offset = curr.first - info.first;
				for(usint m=  offset ; m< offset+curr.second ;m++ ){
					count[m]+= 1;
				}
			}
			fs<<superstr<<",";
			for(usint j =0 ; j< length(info) ;j++){
				fs << count[j] <<",";
			}
			fs<<endl;
			delete[] count;
		}
	}
	cout<<"The coverage count is "<<cnt<<endl;
	fs.close();
	
}
int main(int argc, char** argv){
	uint64_t char_count;
    uchar *data = NULL;
    string FileName;
	

 	double startTime = readTimer();

	if (argc <= 1) 
	{
		cerr << "Usage: ./compress <contig file> " << endl;
	 	return 0;
	}

	FileName = argv[1];

	cout << FileName <<	 endl;

	if (EXIT_FAILURE == get_data(FileName, data, seq_list, char_count))
	{
		return EXIT_FAILURE;
	}
	
	//RLCSABuilder builder(32,0,char_count/(double)1000000,1);
	// Build RLCSA and report some information.
	RLCSA rlcsa(data, char_count, 32, 128, 1, false); // parameter with value '1' is number of threads; available is compiles with muti-thread support
	//data = 0; // The constructor deleted the data.

	if ( !(rlcsa.isOk()) ) 
	{
		return EXIT_FAILURE;
	}
	rlcsa.printInfo();
	rlcsa.reportSize(true);
	rlcsa.writeTo("out.txt");
	
	
	// we sample internal nodes and check their abundances
 	startTime = readTimer();
	create_output(rlcsa);
	
 	cout << "Time for printing maximal contigs: " << readTimer() - startTime << "sec" << endl;
 	

	return EXIT_SUCCESS;
}
	
