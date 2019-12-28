#ifndef OUTPUT_H
#define OUTPUT_H
#include <cstdlib>
#include <string>
#include <typedefs.h>

class Output{
	
	public:
		
		Output();
														
		~Output();
		
		void write_view(const std::string &fn, ScalarField::HostMirror data);

		FILE* file;
		
		size_t frame = 0;

};

#endif // OUTPUT_H
