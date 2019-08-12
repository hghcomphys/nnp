

#include "logger.h"
#include <iostream>

/* ----------------------------------------------------------------------
   setup for class Logger  
------------------------------------------------------------------------- */

LOG::LOG() {}

LOG::LOG(LOG_t type) 
{
    msglevel = type;
    if(false) {
        operator << ("["+getLabel(type)+"]");
    }
}

LOG::~LOG() {
    if(opened) {
        std::cout << std::endl;
    }
    opened = false;
}

LOG& LOG::operator<< (const std::string& msg) 
{
    if( msglevel >= verbosity ) {
        std::cout << msg;
        opened = true;
    }
    return *this;
}

inline std::string getLabel(LOG_t type) 
{
    std::string label;
    switch(type) {
        case DEBUG: label = "DEBUG"; break;
        case INFO:  label = "INFO "; break;
        case WARN:  label = "WARN "; break;
        case ERROR: label = "ERROR"; break;
    }
    return label;
}
