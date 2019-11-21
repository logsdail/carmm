#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include <sstream>
#include <string>
#include <sstream>
#include <iostream>

template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
{
	std::istringstream iss(s);
	return !(iss >> f >> t).fail();
}

template <class T>
bool from_string(T& t, const std::string& s)
{
	std::istringstream iss(s);
	return !(iss >> std::dec >> t).fail();
}

template <class T>
std::string to_string(T& t)
{
	std::string retString;
	std::stringstream os;
	os << t;
	os >> retString;
	return retString;
}

#endif
