#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <sstream>  

enum LOG_t { 
    DEBUG, 
    INFO, 
    WARN, 
    ERROR 
};

class LOG {
public:
    LOG();
    LOG(LOG_t level);
    ~LOG();
    LOG& operator<< (const std::string& message);
    const std::string toString() const;

private:
    bool isOpened;
    LOG_t messageLevel;
    std::string getLabel(LOG_t level);
    std::stringstream buffer;
    // TODO: set outside the class
    LOG_t verbosity = INFO;
    bool isHeader = true;
};

#endif 