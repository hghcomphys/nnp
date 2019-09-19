

#include "logger.h"
#include <stdexcept>
#include <string>
#include <iostream>

/* ----------------------------------------------------------------------
   setup for class Logger  
------------------------------------------------------------------------- */
Log::Log(): Log(DEBUG) {}

Log::Log(LOG_t level):  isOpened(false), messageLevel(level) { setLable(); }

Log::~Log() 
{
    if ( isOpened ) {
        endl(); // go next line
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

Log& Log::operator<< (const char * message) 
{
    if( messageLevel >= verbosity ) {
        buffer << message;
        isOpened = true;
    }
    return *this;
}

Log& Log::operator<< (int message) { 
    operator<< (std::to_string(message)); 
    return *this;
} 

Log& Log::operator<< (double message) { 
    operator<< (std::to_string(message)); 
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

void Log::setLable() {
    if( isHeader ) {
        operator << ("["+getLabel(messageLevel)+"] ");
    }
}

const std::string Log::toString() const { return buffer.str(); }

void Log::clear() { 
    buffer.clear(); 
    setLable();
}

void Log::endl() { buffer << "\n"; }