

#include "logger.h"
#include <iostream>
#include <stdexcept>

/* ----------------------------------------------------------------------
   setup for class Logger  
------------------------------------------------------------------------- */

LOG::LOG(): LOG(INFO) {}

LOG::LOG(LOG_t level):  isOpened(false), messageLevel(DEBUG)
{
    messageLevel = level;
    if( isHeader ) {
        operator << ("["+getLabel(level)+"] ");
    }
}

LOG::~LOG() 
{
    if(isOpened) {

        buffer << "\n";
        switch ( messageLevel )
        {
        case ERROR:
            std::cerr << buffer.str();
            break;

        default:
            std::cout << buffer.str();
        }      
    }
    isOpened = false;
}

LOG& LOG::operator<< (const std::string& message) 
{
    if( messageLevel >= verbosity ) {
        buffer << message;
        isOpened = true;
    }
    return *this;
}

inline std::string LOG::getLabel(LOG_t level) 
{
    std::string label;
    switch(level) {
        case DEBUG: label = "DEBUG"; break;
        case INFO:  label = "INFO "; break;
        case WARN:  label = "WARN "; break;
        case ERROR: label = "ERROR"; break;
    }
    return label;
}

const std::string LOG::toString() const {
    return buffer.str();
}