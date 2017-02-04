#include "VML.h"

#include <sstream>

namespace VML
{
	std::string Vector::toString()const
	{
		std::ostringstream ss;
		ss << "(" << getX() << ", " << getY() << ", " << getZ() << ", " << getW() << ")";
		return ss.str();
	}
}