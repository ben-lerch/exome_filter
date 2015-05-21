// exomefilter v1
//
// Simulate whole-exome sequencing study to identify mendelian disease loci using 1000 genomes reference data
//
// Author: Ben Lerch
// Last Update: Aug 2014

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <time.h>

using namespace std;



time_t T = time(0);
boost::random::mt11213b gen(T);

int numReps;
int numRefIndv;
int numVariants;

int bernoulliTrial(){
boost::random::uniform_int_distribution<> bernoulli(0,1);
return bernoulli(gen);
}

int pickCrossoverPoint(){
boost::random::exponential_distribution<> switchPoint(.00000001);
return switchPoint(gen);
}


int myfindString(vector <string> v, string keyword){
int result = find(v.begin(), v.end(), keyword) - v.begin(); 
return result;
}


// Read the pedigree file into vectors representing each column

vector <string> FAMID, ID, FATID, MOTID;
//male=0,female=1; if indv was sequenced set to 1, if the individual is a founder, set to 1, if the individual has the disease, set to 1, otherwise set to 0
vector <int> DISEASE, SEX, SEQUENCE, FOUNDER;

vector <string> geno;
vector <string> chrom;
vector <int> pos;
vector <string> info;
string line;
string out;

string genofile;
string chromfile;
string posfile;

vector <string> filters;
vector <double> filterInts;
vector <int> passTotal;
vector <string> results;

//af is 1000G allele freq
double maxAF;
vector <double> af;

//gp is GERP score col
//the larger the score, the more conserved the site
double minGERP;
vector <double> GERP;

//pc is phastCons score col
//range of 0 to 1, with 1 being the most conserved
double minphastCons;
vector <double> phastCons;

//py is phyloP score col
//the larger the score, the more conserved the site
double minphyloP;
vector <double> phyloP;

//pp is PolyPhen2 score col
//The score ranges from 0 to 1 with 1 being the most damaging
double minPolyPhen2;
vector <double> PolyPhen2;

//sf is SIFT score col
//Ranges from 0 to 1. The amino acid substitution is predicted damaging if the score is <= 0.05
double maxSIFT;
vector <double> SIFT;

//cd is raw CADD score col
// the annotation profile for a given variant suggests that that variant is likely to be "observed" (negative values) vs "simulated" (positive values)
double minCADD;
vector <double> CADD;

//mt is MutationTaster score col
//higher scores are more damaging
double minMutationTasterScore;
vector <double> MutationTasterScore;

//MutationTaster pred col
//categories are "A" ("disease_causing_automatic"), "D" ("disease_causing"), "N" ("polymorphism") or "P" ("polymorphism_automatic")
vector <string> MutationTasterPred;


//MutationAssessor score col
//higher scores are more damaging
double minMutationAssessorScore;
vector <double> MutationAssessorScore;

//MutationAssessor pred col
//categories are high ("H") or medium ("M"), or predicted non-functional, i.e. low ("L") or neutral ("N")
vector <string> MutationAssessorPred;



//lr is LRT score col
//The score ranges from 0 to 1 and a larger score signifies that the codon is more constrained or a NS is more likely to be deleterious.
double minLRT;
vector <double> LRT;


string AFfile = "ref/AF";
string Annofile = "ref/Anno";
string GERPfile = "ref/dbNSFP_GERP++_RS";
string phastConsfile = "ref/dbNSFP_phastCons100way_vertebrate";
string phyloPfile = "ref/dbNSFP_phyloP100way_vertebrate";
string PolyPhen2file = "ref/dbNSFP_Polyphen2_HDIV_score";
string SIFTfile = "ref/dbNSFP_SIFT_score";
string CADDfile = "ref/dbNSFP_CADD_raw";
string MutationTasterScorefile = "ref/dbNSFP_MutationTaster_score";
string MutationTasterPredfile = "ref/dbNSFP_MutationTaster_pred";
string MutationAssessorScorefile = "ref/dbNSFP_MutationAssessor_score";
string MutationAssessorPredfile = "ref/dbNSFP_MutationAssessor_pred";
string LRTfile = "ref/dbNSFP_LRT_score";





//anno is column saying if variant is synonymous, nonsynonymous, etc
vector <string> anno;

void readPedigree(string pedigree){


string line;
ifstream pedigreeFile (pedigree.c_str());

if(pedigreeFile.is_open()){
while( getline(pedigreeFile,line) ){

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
boost::char_separator<char> sep("\t");
tokenizer tokens(line, sep);
tokenizer::iterator tok_iter = tokens.begin();

FAMID.push_back(*tok_iter); ++tok_iter; ID.push_back(*tok_iter); ++tok_iter; FATID.push_back(*tok_iter); ++tok_iter; MOTID.push_back(*tok_iter); ++tok_iter;
SEX.push_back(boost::lexical_cast<int>(*tok_iter)); ++tok_iter; DISEASE.push_back(boost::lexical_cast<int>(*tok_iter)); ++tok_iter; SEQUENCE.push_back(boost::lexical_cast<int>(*tok_iter));
if(FATID.back()=="." && MOTID.back()=="."){ FOUNDER.push_back(1);}
else{FOUNDER.push_back(0);}

}
}
else{cout << "ERROR: Unable to open pedigree file." << endl;}
}

//read command line arguments

void readArguments(int argc, char** argv){

if(argc == 1){

cout << "\n\nstudyNoise v1\nBen Lerch\n\n" << endl;

cout << "description: Compute distribution of number of variants remaining after filtering a vcf containing genotypes assigned to provided list of individuals." << endl;
cout << "Genotype assignments for gen 1 of pedigree are randomly drawn from 1000G phase 1, and children are assigned genotypes by drawing from parent genotypes, switching haplotype based on exp(100 Million).\n\n" << endl;

cout << "-rep [Int]: number of times to sample from ref vcf" << endl;
cout << "-numRefIndv [Int]: number of individuals in reference vcf" << endl;
cout << "-numVariants [Int]: number of variants in ref vcf" << endl;
cout << "-genoFile [file]: genotype fields from ref VCF" << endl;
cout << "-chromFile [file]: chrom column from ref VCF" << endl;
cout << "-posFile [file]: pos column from ref VCF" << endl;


cout << "PEDIGREE ARGUMENTS" << endl;
cout << "-unrelatedStudy [Int]: construct vcf of Int unrelated individuals" << endl;
cout << "-pedigree [file]: construct study based on pedigree file\n" << endl;

cout << "FILTERING ARGUMENTS" << endl;

cout << "Frequency\n" << endl;
cout << "-maxAF [Double] [file]: filter out variants that are not rare in 1000G\n" << endl;
cout << "-novelindbSNP [file]: filter out variants that are present in dbSNP\n" << endl;
cout << "-novelin1000G: filter out variants that are present in 1000 Genomes with 1000 Genomes Allele Count >= Allele Count of founders in VCF\n" << endl;


cout << "-nonSyn [file]: only keep nonsynonymous variants\n\n" << endl; 


cout << "Conservation\n" << endl;
cout << "-minGERP [Double] [file]: filter out variants that are not highly conserved based on GERP score\n" << endl;
cout << "-minphastCons [Double] [file]: filter out variants that are not highly conserved based on phastCons score\n" << endl;
cout << "-minphyloP [Double] [file]: filter out variants that are not highly conserved based on phyloP score\n" << endl;

cout << "Impact\n" << endl;
cout << "-minPolyPhen2 [Double] [file]: filter out variants that are not predicted to be damaging based on PolyPhen2 score\n" << endl;
cout << "-maxSIFT [Double] [file]: filter out variants that are not predicted to be damaging based on SIFT score\n" << endl;
cout << "-minCADD [Double] [file]: filter out variants that are not predicted to be damaging based on CADD score\n" << endl;
cout << "-minMutationTaster [Double] [file]: filter out variants that are not predicted to be damaging based on MutationTaster score\n" << endl;
cout << "-predMutationTaster [file]: filter out variants that are not predicted to be deleterious by MutationTaster (N or P)\n" << endl;
cout << "-minMutationAssessor [Double] [file]: filter out variants that are not predicted to be damaging based on MutationAssessor score\n" << endl;
cout << "-predMutationAssessor [String] [file]: filter out variants that are not predicted to be deleterious by MutationAssessor\n" << endl;
cout << "-minLRT [Double] [file]: filter out variants that are not predicted to be damaging based on LRT score\n" << endl;


}

for(int i = 1; i < argc; i++){

if(string(argv[i])=="-rep"){numReps = boost::lexical_cast<int>(argv[i+1]);}
if(string(argv[i])=="numRefIndv"){numRefIndv = boost::lexical_cast<int>(argv[i+1]);}
if(string(argv[i])=="numVariants"){numVariants = boost::lexical_cast<int>(argv[i+1]);}
if(string(argv[i])=="-genoFile"){genofile = argv[i+1];}
if(string(argv[i])=="-chromFile"){chromfile = argv[i+1];}
if(string(argv[i])=="-posFile"){posfile = argv[i+1];}


if(string(argv[i])=="-unrelatedStudy"){ int numunrelatedStudy = boost::lexical_cast<int>(argv[i+1]); 
for(int i=0;i<numunrelatedStudy;i++){ FAMID.push_back(string("FAM")+boost::lexical_cast<string>(i)); ID.push_back(string("ID")+boost::lexical_cast<string>(i)); 
FATID.push_back("."); MOTID.push_back("."); SEX.push_back(bernoulliTrial()); DISEASE.push_back(bernoulliTrial()); SEQUENCE.push_back(1); FOUNDER.push_back(1);} 
} 

if(string(argv[i]) == "-pedigree"){ string pedfile = argv[i+1]; readPedigree(pedfile); }

if(string(argv[i]) == "-out" ){ out = argv[i+1];}

if(string(argv[i]) == "-Nonsynonymous"){filters.push_back(argv[i]); filterInts.push_back(0); passTotal.push_back(0); results.push_back(""); ifstream n (Annofile.c_str()); if(n.is_open()){ while( getline(n,line) ){ anno.push_back(line);}}}

if(string(argv[i]) == "-AlleleFrequency"){filters.push_back(argv[i]); maxAF = boost::lexical_cast<double>(argv[i+1]); filterInts.push_back(boost::lexical_cast<double>(argv[i+1]));  passTotal.push_back(0); results.push_back("");
ifstream a (AFfile.c_str()); if(a.is_open()){ while( getline(a,line) ){ af.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-novelin1000G"){filters.push_back(argv[i]);}


if(string(argv[i]) == "-GERP"){filters.push_back(argv[i]); minGERP = boost::lexical_cast<double>(argv[i+1]); filterInts.push_back(boost::lexical_cast<double>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (GERPfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ GERP.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-phastCons"){filters.push_back(argv[i]); minphastCons = boost::lexical_cast<double>(argv[i+1]); filterInts.push_back(boost::lexical_cast<double>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (phastConsfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ phastCons.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-phyloP"){filters.push_back(argv[i]); minphyloP = boost::lexical_cast<double>(argv[i+1]); filterInts.push_back(boost::lexical_cast<double>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (phyloPfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ phyloP.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-PolyPhen2"){filters.push_back(argv[i]); minPolyPhen2 = boost::lexical_cast<double>(argv[i+1]); filterInts.push_back(boost::lexical_cast<double>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (PolyPhen2file.c_str()); if(g.is_open()){ while( getline(g,line) ){ PolyPhen2.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-SIFT"){filters.push_back(argv[i]); maxSIFT = boost::lexical_cast<int>(argv[i+1]); filterInts.push_back(boost::lexical_cast<int>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (SIFTfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ SIFT.push_back(boost::lexical_cast<int>(line));}}}

if(string(argv[i]) == "-CADD"){filters.push_back(argv[i]); minCADD = boost::lexical_cast<int>(argv[i+1]); filterInts.push_back(boost::lexical_cast<int>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (CADDfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ CADD.push_back(boost::lexical_cast<int>(line));}}}

if(string(argv[i]) == "-MutationTasterScore"){filters.push_back(argv[i]); minMutationTasterScore = boost::lexical_cast<double>(argv[i+1]); filterInts.push_back(boost::lexical_cast<double>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (MutationTasterScorefile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationTasterScore.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-MutationTasterPred"){filters.push_back(argv[i]); filterInts.push_back(0); passTotal.push_back(0); results.push_back("");
ifstream g (MutationTasterPredfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationTasterPred.push_back(line);}}}

if(string(argv[i]) == "-MutationAssessorScore"){filters.push_back(argv[i]); minMutationAssessorScore = boost::lexical_cast<double>(argv[i+1]); filterInts.push_back(boost::lexical_cast<double>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (MutationAssessorScorefile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationAssessorScore.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-MutationAssessorPred"){filters.push_back(argv[i]); filterInts.push_back(0); passTotal.push_back(0); results.push_back("");
ifstream g (MutationAssessorPredfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationAssessorPred.push_back(line);}}}

if(string(argv[i]) == "-LRT"){filters.push_back(argv[i]); minLRT = boost::lexical_cast<int>(argv[i+1]); filterInts.push_back(boost::lexical_cast<int>(argv[i+1])); passTotal.push_back(0); results.push_back("");
ifstream g (LRTfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ LRT.push_back(boost::lexical_cast<int>(line));}}}




if(string(argv[i]) == "-segregatesDom"){filters.push_back(argv[i]); filterInts.push_back(0); passTotal.push_back(0); results.push_back("");}
if(string(argv[i]) == "-segregatesRec"){filters.push_back(argv[i]); filterInts.push_back(0); passTotal.push_back(0); results.push_back("");}


}

}

int pickRefSample(){
boost::random::uniform_int_distribution<> refSample(1,numRefIndv);
return refSample(gen);
}

//apply Filters

void applyFilters(string myVCFgeno, vector <string> inVCF,  vector <int> DISEASEinVCF, vector <int> SEQUENCEinVCF, int k){

int filterPass = 0;

for(int i=0; i<filters.size(); i++){


//consider gene-level versions of these

if(filters[i]=="-segregatesRec"){ 
int pass = 1;
for(int i=0; i<inVCF.size(); i++){ if( (SEQUENCEinVCF[i]==1 && DISEASEinVCF[i]==0 && myVCFgeno.substr(i*2,2)=="11") || (SEQUENCEinVCF[i]==1 && DISEASEinVCF[i]==1 && myVCFgeno.substr(i*2,2)!="11") ){pass=0; break;}}
if(pass==1){filterPass++;} 
}

if(filters[i]=="-segregatesDom"){ 
int pass = 1; 
for(int i=0; i<inVCF.size(); i++){ if( (SEQUENCEinVCF[i]==1 && DISEASEinVCF[i]==0 && myVCFgeno.substr(i*2,2)!="00" ) || (SEQUENCEinVCF[i]==1 && DISEASEinVCF[i]==1 && myVCFgeno.substr(i*2,2)=="00") ){pass=0; break;}} 
if(pass==1){filterPass++;} 
}

if(filters[i]=="-Nonsynonymous"){if(anno[k] == "Nonsynonymous"){filterPass++;}else{break;}}

if(filters[i]=="-AlleleFrequency"){ if( af[k] < maxAF){filterPass++;}else{break;}}

if(filters[i]=="-GERP"){ if( GERP[k] > minGERP ){filterPass++;}else{break;}}

if(filters[i]=="-phastCons"){ if( phastCons[k] > minphastCons ){filterPass++;}else{break;}}

if(filters[i]=="-phyloP"){ if( phyloP[k] > minphyloP){filterPass++;}else{break;}}

if(filters[i]=="-PolyPhen2"){ if( PolyPhen2[k] > minPolyPhen2 ){filterPass++;}else{break;}}

if(filters[i]=="-SIFT"){ if( SIFT[k] < maxSIFT ){filterPass++;}else{break;}}

if(filters[i]=="-CADD"){ if( CADD[k] > minCADD){filterPass++;}else{break;}}

if(filters[i]=="-MutationTasterScore"){ if( MutationTasterScore[k] > minMutationTasterScore ){filterPass++;}else{break;}}

if(filters[i]=="-MutationTasterPred"){if(MutationTasterPred[k]=="D" || MutationTasterPred[k]=="A"){filterPass++;}else{break;}}

if(filters[i]=="-LRT"){ if( LRT[k] > minLRT){filterPass++;}else{break;}}



}

for(int i=0; i<filterPass; i++){passTotal[i]++;}

}







//convert input values to a random vcf


void pedToVCF(vector <string> chrom, vector <int> pos, vector <string> info, vector <string> geno){

vector <string> inVCF, FAMIDinVCF;
vector <int> colMOT, colFAT, oriMOT, oriFAT, upToMOT, upToFAT, FOUNDERinVCF, SEXinVCF, SEQUENCEinVCF, DISEASEinVCF;


string chr = "1";

//initialize values for sampling from VCF
while(inVCF.size() < FAMID.size()){ 
for(int i = 0; i < FAMID.size(); i++){
if(find(inVCF.begin(), inVCF.end(), ID[i]) == inVCF.end()){
if(FOUNDER[i]==1){ inVCF.push_back(ID[i]); FOUNDERinVCF.push_back(FOUNDER[i]); SEXinVCF.push_back(SEX[i]); SEQUENCEinVCF.push_back(SEQUENCE[i]); DISEASEinVCF.push_back(DISEASE[i]); FAMIDinVCF.push_back(FAMID[i]);
colMOT.push_back(pickRefSample()); colFAT.push_back(pickRefSample()); oriMOT.push_back(bernoulliTrial()); oriFAT.push_back(bernoulliTrial()); upToMOT.push_back(pickCrossoverPoint()); upToFAT.push_back(pickCrossoverPoint());}

if(find(inVCF.begin(), inVCF.end(), FATID[i]) != inVCF.end() && find(inVCF.begin(), inVCF.end(), MOTID[i]) != inVCF.end() ){ 
inVCF.push_back(ID[i]); FOUNDERinVCF.push_back(FOUNDER[i]); SEXinVCF.push_back(SEX[i]); SEQUENCEinVCF.push_back(SEQUENCE[i]); DISEASEinVCF.push_back(DISEASE[i]); FAMIDinVCF.push_back(FAMID[i]);
colMOT.push_back(myfindString(inVCF,MOTID[i])); colFAT.push_back(myfindString(inVCF,FATID[i]));
oriMOT.push_back(bernoulliTrial()); oriFAT.push_back(bernoulliTrial()); upToMOT.push_back(pickCrossoverPoint()); upToFAT.push_back(pickCrossoverPoint()); }
}}
}


//construct fake VCF line-by-line

int iter;
//start by reading in ref VCF and creating first 9 fields of my VCF
for(int k=0; k<numVariants;k++){
string myVCFgeno;
iter = k;

//genotype each individual and add to my VCF line
for(int i=0;i<inVCF.size();i++){


if(FOUNDERinVCF[i]==1){
if(chr!= chrom[k]){ chr=chrom[k]; oriMOT[i]=bernoulliTrial(); oriFAT[i]=bernoulliTrial(); upToMOT[i]=pickCrossoverPoint(); upToFAT[i]=pickCrossoverPoint();}
if(upToMOT[i] < pos[k]){oriMOT[i] = (oriMOT[i]+1) % 2; upToMOT[i]=upToMOT[i]+pickCrossoverPoint(); }
if(upToFAT[i] < pos[k]){oriFAT[i] = (oriFAT[i]+1) % 2; upToFAT[i]=upToFAT[i]+pickCrossoverPoint(); }
if(chr!="X"){myVCFgeno += geno[k].substr(colMOT[i]*2+oriMOT[i],1); myVCFgeno += geno[k].substr(colFAT[i]*2+oriFAT[i],1);}
else if(SEXinVCF[i]==1){upToFAT[i]=500000000; myVCFgeno += geno[k].substr(colMOT[i]*2+oriMOT[i],1); myVCFgeno += geno[k].substr(colFAT[i]*2,1);}
else if(SEXinVCF[i]==0){myVCFgeno += geno[k].substr(colMOT[i]*2+oriMOT[i],1)+string(".");}
}



if(FOUNDERinVCF[i]==0){
string m = myVCFgeno.substr(colMOT[i]*2+oriMOT[i],1);
string f = myVCFgeno.substr(colFAT[i]*2+oriFAT[i],1);

if(chr!= chrom[k]){ chr=chrom[k]; oriMOT[i]=bernoulliTrial(); oriFAT[i]=bernoulliTrial(); upToMOT[i]=pickCrossoverPoint(); upToFAT[i]=pickCrossoverPoint();}
if(upToMOT[i] < pos[k]){oriMOT[i] = (oriMOT[i]+1) % 2; upToMOT[i]=upToMOT[i]+pickCrossoverPoint(); }
if(upToFAT[i] < pos[k]){oriFAT[i] = (oriFAT[i]+1) % 2; upToFAT[i]=upToFAT[i]+pickCrossoverPoint(); }
if(chr!="X"){myVCFgeno += m; myVCFgeno += f;}
else if(SEXinVCF[i]==1){upToFAT[i]=500000000; myVCFgeno += m; f = myVCFgeno.substr(colFAT[i]*2,1); myVCFgeno += f;}
else if(SEXinVCF[i]==0){myVCFgeno += m; myVCFgeno += ".";}
}

}

//define AC for this variant
int myAC = 0;
for(int i=0; i<inVCF.size(); i++){ 
if(SEQUENCEinVCF[i]==1 && (myVCFgeno.substr(i,2) == ".1" || myVCFgeno.substr(i,2) == "1." || myVCFgeno.substr(i,2) == "01" || myVCFgeno.substr(i,2) == "10")){myAC++;}
if(SEQUENCEinVCF[i]==1 && myVCFgeno.substr(i,2) == "11"){myAC++; myAC++;}
}


//remove if monomorphic
if(myAC==0){continue;}

applyFilters(myVCFgeno, inVCF, DISEASEinVCF, SEQUENCEinVCF, iter); 

} 

}




int main(int argc, char** argv){

readArguments(argc, argv);

ifstream g (genofile.c_str());
if(g.is_open()){
while( getline(g,line) ){
geno.push_back(line);
}}

ifstream c (chromfile.c_str());
if(c.is_open()){
while( getline(c,line) ){
chrom.push_back(line);
}}

ifstream p (posfile.c_str());
if(p.is_open()){
while( getline(p,line) ){
pos.push_back(boost::lexical_cast<int>(line));
}}

int counter = 1;
for(int i=0;i<numReps;i++){ 
pedToVCF(chrom, pos, info, geno);
for(int i=0; i<filters.size(); i++){ results[i]=results[i]+string(",")+boost::lexical_cast<string>(passTotal[i]);}
for(int i=0; i<passTotal.size(); i++){passTotal[i]=0;}
cout << "Iteration " << counter << " complete." << endl;
counter++;
}



string rfile = string("output/") + out + string("/plots.R");

ofstream myfile;
myfile.open (rfile.c_str());

for(int i=0; i<filters.size();i++){ myfile << filters[i].substr(1,string::npos) << "<-c(" << results[i].substr(1,string::npos) << ");"; }

for(int i=0; i<filters.size();i++){ myfile << "png(filename = 'output/" << out << "/" << filters[i].substr(1,string::npos) << ".png'); max<-max("<< filters[i].substr(1,string::npos)<<"); min<-min("
    << filters[i].substr(1,string::npos) << ");mybreaks<-max/(max-min)*100; hist(" << filters[i].substr(1,string::npos) << 
    ",breaks=mybreaks,col=rgb(2/255, 112/255, 200/255, 0.90),xlab='Number of Variants',ylab='Iterations');dev.off();";}


myfile.close();


}

