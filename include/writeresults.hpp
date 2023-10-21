#ifndef _WRITERESULTS_HPP_DEFINED_
#define _WRITERESULTS_HPP_DEFINED_

#include <iomanip>

#pragma once
class writeresults {
	public:
	void output( double , double , double , double , double , double , double, double, double, double, double);
	
	template < typename A >
        std::string to_string_with_precision(A , int );
};
#endif
