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

string test;
typedef unsigned char uchar;
typedef unsigned long usint;

typedef struct{
	string header;
	string seq;
} sequence_t;

vector<sequence_t> seq_list;		//list of sequences and headers read from FASTA file
vector<sequence_t> seq_list_in;
vector<sequence_t> seq_list_test;

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
//Extracts the sequence description and removes the leading">m"
string extractAccession(usint idx){
	string s = seq_list[idx].header;
	return s.substr(2, s.find(" ")-2);
	
}

// Generate Output Files for the case of the Exact match
void exactMatch(RLCSA& rlcsa){
	fstream fs;
	map<usint, vector<pair<usint, usint> > > cover;				// For each sequence, we have a list containing <start,size> of its sub-strings
	
	usint index;
	string pattern;
	usint *pos;
	char track[rlcsa.getNumberOfSequences()];				// to identify the position of unique sequences from the entire set
	fs.open("data/redundant.csv",fstream::out);
	
	for(usint j=0; j<rlcsa.getNumberOfSequences() ;j++){
		//string pattern((char *)rlcsa.display(j));
		
		pattern = seq_list[j].seq;
		pair_type result = rlcsa.count(pattern);			//count returns the suffix array range of matches, which is essentially the number of matches
		
		if(length(result)>1){
			
			usint check =1,first_occ = 0;
			pos = rlcsa.locate(result);						// locate gives the string position value stored at a position of suffix array
			first_occ = rlcsa.getSequenceForPosition(pos[0]);
			
			for(usint p=0; p<length(result); p++){
				index = rlcsa.getSequenceForPosition(pos[p]);	// used to easily extract the sequence at position p of the complete string
				
				// Generating the file containing FASTA headers of each sequence and it's corresponding matches
				if(index != j){
					fs <<extractAccession(j)<< ",";
					fs << extractAccession(index)<<endl;
					
					// creating a data structure which stores the positions of its sub-strings, found in the dataset
					cover[index].push_back(make_pair(pos[p], pattern.size()));
				}	
				
				//find exact matches of same length				
				if(length(rlcsa.getSequenceRange(index)) != pattern.size() ){					
					check = 0;
				}
			}
			
			//since this string has matches in the dataset, we are marking it to be redundant
			track[j] = 'n';
			
			// if check ==1 then all matches are perfect and we need to include only one of them in the list
			if(check == 1){
				track[first_occ] = 'y';
			}
		}
		else{			
			// count didn't return any other matches, so marking this string to be unique
			track[j] = 'y';
			
		}
		pattern.clear();
	}
	fs.close();
	

	//Create file containing list of non-redundant sequences in FASTA format
	fs.open("data/Non_redundant.fasta",fstream::out);
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
	fs.open("data/Coverage.csv", fstream::out);
	cnt =0;
	//for(usint i =0; i< rlcsa.getNumberOfSequences(); i++){
	for (map<usint, vector<pair<usint, usint>>>::iterator it=cover.begin(); it!=cover.end(); ++it){	
		usint i = it->first;
		
		if(track[i]=='y'){
			cnt++;
			uchar* superstr = rlcsa.display(i); 			// fetches the sequence number i from rlcsa data structure
			usint size = (it->second).size();				// Number of sub-string matches, given by the length of the vector
			
			pair_type info = rlcsa.getSequenceRange(i);
			usint *count = new usint[length(info)];
			bzero(count,length(info));
			for(usint j=0; j< size ;j++){
				pair<usint,usint> curr = (it->second)[j];	// start and length info about the sub-string match
				
				usint offset = curr.first - info.first;		// Finding the match location of sub-string within the super string and adding coverage data
				for(usint m=  offset ; m< offset+curr.second ;m++ ){
					count[m]+= 1;
				}
			}
			
			//Populating the coverage file
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


// Generate Output files for the case of Inexact Matches
void approxMatch(RLCSA& rlcsa){
	
	char track[rlcsa.getNumberOfSequences()];
	usint gino =0;	
	map<usint, vector<pair<usint, usint>>> cover;			// sama data structure as before
	fstream fs;
	fs.open("data/redundant_inexact.csv",fstream::out);
	
	for(usint j=0; j<rlcsa.getNumberOfSequences() ;j++){
		//string pattern((char *)rlcsa.display(j));
		
		track[j] = 'y';
		string pattern = seq_list_in[j].seq;
		usint half = (usint)pattern.size()/2;
		
		//Left half search
		pair_type candidate = rlcsa.count(pattern.substr(0,half));			// Search the left half of the sequence to find possible matches and then search in their neighbourhood
		usint size = length(candidate);
		if(size > 1){
			usint* px = rlcsa.locate(candidate);
			gino++;
			for(usint i=0; i< size ;i++){
				
				usint off = (usint) pattern.size();
				pair_type range = (rlcsa.getSequenceRangeForPosition(px[i]));
				usint starting = px[i] - range.first;
				usint remainder_size = range.second -(px[i]);    //check  add +1
				usint index = rlcsa.getSequenceForPosition(px[i]);
				//uchar* test = new uchar[length(rlcsa.getSequenceRangeForPosition(px[i]))];
				
				
				uchar* test = (rlcsa.display(index));				// fetch the candidate string to search for match
				// conditions in "if" try to filter out infeasible matches to reduce operations
				if(remainder_size >= pattern.size() && index != j && length(range)>pattern.size()){
					int err=0;
					for(usint m=half; m< (usint)pattern.size(); m++){
											
						if(test[starting +m] != pattern[m]){
							err++;
							
							// Allow for one mismatch, if more than one break the search
							if(err>1)
								break;
						}
					}
					// Add to the redudant file if only one error
					if(err<=1 ){
						fs << extractAccession(j)<< ",";
						fs << extractAccession(index)<<endl;
						track[j]='n';
						cover[index].push_back(make_pair(px[i], off));			// if match, add to the list of super strings
					}
				}
			}			
		}
		
		//Do the same search for Right Half of the pattern
		candidate = rlcsa.count(pattern.substr(half));
		size = length(candidate);
		if(size > 1){
			usint* px = rlcsa.locate(candidate);
			gino++;
			
			for(usint i=0; i< size ;i++){
				
				usint off = (usint) pattern.size();
				pair_type range = (rlcsa.getSequenceRangeForPosition(px[i]));
				usint starting = px[i] - half;
				usint remainder_size = (px[i]) - range.first;    //check
				usint index = rlcsa.getSequenceForPosition(px[i]);
				
				//uchar* test = new uchar[length(rlcsa.getSequenceRangeForPosition(px[i]))];
				uchar* test = (rlcsa.display(index));
				if(remainder_size >= half && index != j && length(range)>pattern.size()){
					int err=0;
					for(usint m=half-1; m>=0 ; m--){
											
						if(test[starting +m] != pattern[m]){
							err++;
						
							if(err>1)
								break;
						}
					}
					
					//Add to the file if only one mismatch
					if(err<=1 ){
						fs << extractAccession(j)<< ",";
						fs << extractAccession(index)<<endl;
						track[j]='n';
						cover[index].push_back(make_pair(px[i], off));
					}
				}
			}			
		}
	}
	
	fs.close();
	
	
	//Create non redundant csv
	fs.open("data/Non_redundant_inexact.fasta",fstream::out);
	int cnt = 0;
	for(usint i =0; i< rlcsa.getNumberOfSequences(); i++){
		
		if(track[i] == 'y'){
			cnt++;
			fs << ">m" << seq_list_in[i].header<<endl;
			fs << seq_list_in[i].seq<<endl;
		}
	}
	cout<<"The count of non-redundant string si "<<cnt<<endl;
	cout<< "Total iterations "<<gino<<endl;
	fs.close();
	
	
	
	//Create Coverage file
	fs.open("data/Coverage_inexact.csv", fstream::out);
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
				
				usint offset = curr.first - info.first;					// Find the match position in sequence and add coverage info
				for(usint m=  offset ; m< offset+curr.second ;m++ ){
					count[m]+= 1;
				}
			}
			
			// Add to the coverage file
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
	
	
	// Create the Meta Data for Exact matches
 	startTime = readTimer();
	exactMatch(rlcsa);
		
 	cout << "Time for printing exact matches: " << readTimer() - startTime << "sec" << endl;
	
	// For the case of inexact Match, use the List of Non-redundant sequence generated by previous method
	if (EXIT_FAILURE == get_data("data/Non_redundant.fasta", data, seq_list_in, char_count))
	{
		return EXIT_FAILURE;
	}
	// RLCSA data structure built on the list of unique sequences
	RLCSA rlcsa_in(data, char_count, 32, 128, 1, false); 
	if ( !(rlcsa_in.isOk()) ) 
	{
		return EXIT_FAILURE;
	}
	rlcsa.printInfo();
	rlcsa.reportSize(true);
	rlcsa.writeTo("out.txt");
	
	//Generate Output files for Inexact Matches
	startTime = readTimer();
	approxMatch(rlcsa_in);
	cout << "Time for Approx matching: " << readTimer() - startTime << "sec" << endl;
	
	
	//Testing the unique sequences of generated Data
	
	if (EXIT_FAILURE == get_data("data/Non_redundant_inexact.fasta", data, seq_list_test, char_count))
	{
		return EXIT_FAILURE;
	}
	RLCSA rlcsa_test(data, char_count, 32, 128, 1, false); 
	int gino=0;
	
	// If gino is > 0 at the end of loop, then there are repeats in the non-redundant list generated
	for(usint j=0; j<rlcsa_test.getNumberOfSequences() ;j++){
		string pattern = seq_list_test[j].seq;
		pair_type result = rlcsa_test.count(pattern);
		if(length(result)>1){
			gino++;
			cout<<"There's a repeat!!"<<endl;
		}
	}
	//cout<<"Total count"<<gino<<endl;
	gino=0;
	for(usint j=0; j<rlcsa_test.getNumberOfSequences() ;j++){
		string pattern = seq_list_test[j].seq;
		pair_type result = rlcsa.count(pattern);
		if(length(result)==0){
			gino++;
			cout<<"There's an error!!"<<endl;
		}
	}
	
	return EXIT_SUCCESS;
}
	
