/* Copyright <2015> LLNL
 * Author	:	poudel1
 * Date		:	06/19/2015
 */
#ifndef _NANO_TIME_HPP_INCLUDED
#define _NANO_TIME_HPP_INCLUDED

#include <iostream>
#include <sstream>
#include <string>

class nano_time {
 public:
	uint32_t seconds;
	uint32_t nanoseconds;

	nano_time() : seconds(0), nanoseconds(0) { }
	nano_time(uint32_t sec, uint32_t nano_sec) :
				seconds(sec), nanoseconds(nano_sec) { }

	explicit nano_time(std::string str) {
		int pos = str.find(".");
		seconds = std::stoi(str.substr(0, pos));

		std::stringstream nanoseconds_str;
		nanoseconds_str << str.substr(pos + 1);
		if(str.length() - pos - 1 < 9) {
			for(unsigned int i = 0; i < (9 - str.length() + pos + 1); i++)
				nanoseconds_str << "0";
		}
		nanoseconds = std::stoi(nanoseconds_str.str());
	}

	std::string to_string() {
		std::stringstream str;
		str << seconds << "." << nanoseconds;
		return str.str();
	}

	friend nano_time operator+(nano_time lhs, nano_time rhs) {
		uint32_t _seconds = lhs.seconds + rhs.seconds;
		uint32_t _nanoseconds = lhs.nanoseconds + rhs.nanoseconds;
		if(_nanoseconds >= 1000000000) {
			_nanoseconds = _nanoseconds % 1000000000;
			_seconds++;
		}
		return nano_time(_seconds, _nanoseconds);
	}

	friend nano_time operator-(nano_time lhs, nano_time rhs) {
		if (lhs.nanoseconds < rhs.nanoseconds)
			return nano_time(
				lhs.seconds - rhs.seconds - 1, lhs.nanoseconds * 10 - rhs.nanoseconds);
		else
			return nano_time(
				lhs.seconds - rhs.seconds, lhs.nanoseconds - rhs.nanoseconds);
	}

	friend bool operator>(const nano_time& lhs, const nano_time& rhs) {
		if (lhs.seconds > rhs.seconds)
			return true;
		else if (lhs.seconds == rhs.seconds)
			return lhs.nanoseconds > rhs.nanoseconds;
		else
			return false;
	}
} __attribute__ ((packed));

#endif  //  _NANO_TIME_HPP_INCLUDED
