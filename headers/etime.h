#include <ctime>
#include <string>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <cmath>

#ifndef ETIME_H
#define ETIME_H

/*
***** BEGIN LICENSE BLOCK *****
*	
*	The content of this code is released under
*	the (http://goo.gl/zNe9hw) MIT License.
*
*	Developer: Vagner Bessa (vagner.fisica@gmail.com)
*	
***** END LICENSE BLOCK *******
*/

/*
*******************************
*	etime C++ class's
*	- Purpose: The main purpose of this library is
*			   to compute the elapsed time between 
*			   computations
*
*	- Info: This class can be found at
*		    https://github.com/vagner-fisica/ioput.git
*			Please check README to see more details.
*******************************
*/

#define HOUR 3600
#define MINUTE 60

using std::stringstream;
using std::string;
using std::setw;
using std::setfill;
using std::cout;

class etime{
	
	public:
		
//		default constructor
		etime();

		void start();

		void end(string);
		double end();
		
	private:

		clock_t clock_;
		double etime_;
		
		void start_clock();
		void end_clock();

		string get_formated_time(double);		

};
etime::etime (){

}
/*		
	register early CPU clock
*/
void etime::start_clock (){
	clock_ = clock();
}
/*		
	register late CPU clock,
	compute difference between previous
	and set etime_ member as the elapsed
	time in seconds between previous call
	of start() and forward call of end();
*/
void etime::end_clock (){
	clock_ = clock() - clock_;
	etime_ = ((double)clock_)/CLOCKS_PER_SEC;
}
/*		
	convert etime_ (in seconds) into
	a formated string of the type hh:mm:ss
	and return it.
*/
string etime::get_formated_time (double time_target){
	int hour, minute, second;
	stringstream aux;
	string str;
	
	if (time_target >= 1.0) {
		hour = (int)round(time_target/HOUR);
		second = (int)round(fmod(time_target,HOUR));
		minute = (int)round(second/MINUTE);
		second = (int)round(fmod(second,MINUTE));
	
		aux << setw(2) << setfill('0') << hour << ":" << setw(2) << setfill('0') << minute << ":" << setw(2) << setfill('0') << second;
		aux >> str;
		return str;
	} else {
		aux << time_target << "s";
		aux >> str;
		return str;
	}

}
/*
	calls start_clock()
*/
void etime::start(){
	start_clock();
}
/*
	calls end_clock() and prints
	the string, and get_formated_time(etime_)
	with an ending line.
*/
void etime::end(string s){
	end_clock();
	cout << s;
	cout << get_formated_time(etime_);
	cout << std::endl;
}
/*
	return etime_ in seconds.
*/
double etime::end(){
	end_clock();
	return etime_;
}
#endif
