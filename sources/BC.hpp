#ifndef BC_HPP
#define BC_HPP

#include <map>
#include "typedefs.hpp"
#include "inlet.hpp"
#include "outlet.hpp"
#include "wall.hpp"

using namespace std;

map<string, bCondition> bcList = {condition("inlet", inlet),
				  condition("outlet", outlet),
				  condition("wall", wall)};


#endif
