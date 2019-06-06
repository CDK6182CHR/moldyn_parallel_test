#ifndef ETIME_H
#define ETIME_H

#include <ctime>
#include <string>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <cmath>



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
        void register_time(string);
		double end();
		
	private:

		clock_t clock_;
		double etime_;
		
		void start_clock();
		void end_clock();

		string get_formated_time(double);		

};
#endif
