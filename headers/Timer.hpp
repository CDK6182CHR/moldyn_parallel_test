#ifndef Timer_H
#define Timer_H

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
*	Developer: Vagner Bessa (bessavagner@gmail.com)
*	
***** END LICENSE BLOCK *******
*/

/*
*******************************
*	Timer C++ class's
*	- Purpose: The main purpose of this class is
*			   to compute the elapsed time between 
*			   computations
*
*	- Info: This class can be found at
*		    https://github.com/vagner-fisica/ioput.git
*			Please check README to see more details.
*******************************
*/

class Timer{

	private:
        constexpr static int m_hour = 3600;
        constexpr static int m_minute = 60;
        
		clock_t m_clock;
		double m_timer;
		
		void start_clock();
		void end_clock();

		std::string get_formated_time(double);

	public:
		
//		Default and only constructor.
		Timer();

		void start_timer();
		
		void end_timer(std::string);
		double end_timer();
		
        void register_time(std::string);
				
};
#endif
