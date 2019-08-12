#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>

enum LOG_t { 
    DEBUG, 
    INFO, 
    WARN, 
    ERROR 
};

class LOG {
public:
    LOG();
    LOG(LOG_t type);
    ~LOG();
    LOG& operator<< (const std::string& msg);

private:
    bool opened = false;
    LOG_t msglevel = LOG_t::DEBUG;
    std::string getLabel(LOG_t type);
    // TODO: set outside the class
    LOG_t verbosity = LOG_t::INFO;
};

#endif 