#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sstream>
#include <climits>
#include <string>
#include <vector>


namespace _Cide{
	extern std::vector< std::vector<int> > graphT;
    extern std::vector<int> attentionQuota;
}

typedef std::pair<int,double> infPair; // node - monetary marginal gain
typedef long long int64;
typedef unsigned long long uint64;



#ifdef _WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

using namespace std;

float getCurrentMemoryUsage();
float getRunningTime(time_t startTime);

void stringTokenizer(string& str, float *tokens, int size, string& delimiters);
void stringTokenizer(string& str, double *tokens, int size, string& delimiters);

void rtrim(char *str);
void ltrim(char *str);
void trim(char *str);

inline unsigned int strToInt(string s) {
  unsigned int i;
  istringstream myStream(s);

  if(myStream >> i) 
    return i;
  
  else {
    std::cout << "String " << s << " is not a number." << endl;
    return atoi(s.c_str());
  }
  
  return i;
}

inline unsigned int strToInt(const char* s) {
  unsigned int i;
  istringstream myStream(s);

  if (myStream >> i)
    return i;
  
  else {
    cout << "String " << s << " is not a number." << endl;
    return atoi(s);
    }
  
  return i;
}

inline float strToFloat(const char* s) {
  return atof(s);  
}

inline float strToFloat(string s) {
  return atof(s.c_str());
}

inline string floatToStr(float f) {
  stringstream ss;
  ss << f;
  return ss.str();
}

inline int64_t strToInt64(string s) {
  int64_t i;
  istringstream myStream(s);

  if (myStream >> i)
    return i;
  
  else {
    cout << "String " << s << " is not a number." << endl;
    exit(1);    
  }
  
  return i;
}

inline string intToStr(int i) {
  stringstream ss;
  ss << i;
  return ss.str();  
}

inline double strToDouble(string s) {	
// 	return std::stod(s);
	double a = 0; 
	stringstream ss;
	ss << s;
	ss >> a;
	return a;
	
	
}

float getCurrentMemoryUsage();

#endif
