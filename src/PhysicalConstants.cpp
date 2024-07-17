#include "PhysicalConstants.hpp"

PhysicalConstants& PhysicalConstants::getInstance() {
        static PhysicalConstants instance;
        return instance;
    }