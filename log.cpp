

#include "log.h"
#include <iostream>
#include <stdexcept>

/* ----------------------------------------------------------------------
   setup for class Logger  
------------------------------------------------------------------------- */

Log::Log(): Log(INFO) {}

Log::Log(LOG_t level):  isOpened(false), messageLevel(DEBUG)
{
    messageLevel = level;
    if( isHeader ) {
        operator << ("["+getLabel(level)+"] ");
    }
}

Log::~Log() 
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

Log& Log::operator<< (const std::string& message) 
{
    if( messageLevel >= verbosity ) {
        buffer << message;
        isOpened = true;
    }
    return *this;
}

inline std::string Log::getLabel(LOG_t level) 
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

const std::string Log::toString() const {
    return buffer.str();
}