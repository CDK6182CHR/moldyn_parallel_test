#include "Timer.hpp"

Timer::Timer (){

}
/*		
	Private function that start timer
*/
void Timer::start_clock (){
	m_clock = clock();
}
//  Private function that
//	register current CPU clock,
//	compute difference between previous
//	and set m_timer_ member as the elapsed
//	time in seconds between previous call
//	of start() and forward call of end();
void Timer::end_clock (){
	m_clock = clock() - m_clock;
	m_timer = ((double)m_clock)/CLOCKS_PER_SEC;
}
//  Private function that
//	converts time_target = m_timer (in seconds) into
//	a formated string of the type hh:mm:ss
//	if it is higher than 1s	and return it as a std::string
std::string Timer::get_formated_time (double time_target){
	int hour_digit, minute_digit, second_digit;
	std::stringstream aux;
	std::string str;
	
	if (time_target >= 1.0) {
		hour_digit = (int)round(time_target/m_hour);
		second_digit = (int)round(fmod(time_target, m_hour));
		minute_digit = (int)round(second_digit/m_minute);
		second_digit = (int)round(fmod(second_digit, m_minute));
	
		aux << std::setw(2) << std::setfill('0') << hour_digit << ":";
		aux << std::setw(2) << std::setfill('0') << minute_digit << ":";
		aux << std::setw(2) << std::setfill('0') << second_digit;
		
		aux >> str;
		
		return str;
		
	} else {
		aux << time_target << "s";
		
		aux >> str;
		
		return str;
	}
}
//  Public function that
//	register current CPU clock,
//	compute difference between previous one
//	and set m_timer member as the elapsed
//	time in seconds between previous call
//	of start_timer() and forward call of end_timer(),
//  and prints te massagem followed by
//  the return of get_formated_time(m_timer)
void Timer::register_time(std::string message){
	m_timer = ((double)(clock() - m_clock))/CLOCKS_PER_SEC;
	std::cout << message;
	std::cout << get_formated_time(m_timer);
	std::cout << std::endl;	
}
//	Public function that
//	calls start_clock()
void Timer::start_timer(){
	start_clock();
}
//	Public function that
//	calls end_clock(), prints
//	message and calls private function
//  get_formated_time(m_timer) with an ending line.
void Timer::end_timer(std::string message){
	end_clock();
	std::cout << message;
	std::cout << get_formated_time(m_timer);
	std::cout << std::endl;
}
//  Overloaded end_timer() public
//  function that returns m_timer
//  in seconds.
double Timer::end_timer(){
	end_clock();
	return m_timer;
}
