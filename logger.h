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

class Log {
public:
    Log();
    Log(LOG_t level);
    ~Log();
    Log& operator<< (const std::string& message);
    Log& operator<< (const char * message);
    Log& operator<< (int message);
    Log& operator<< (double message);
    const std::string toString() const;
    void clear();
    void endl();

private:
    bool isOpened;
    LOG_t messageLevel;
    std::stringstream buffer;
    std::string getLabel(LOG_t level);
    void setLable();
    bool isBlank();
    // TODO: set outside the class
    LOG_t verbosity = INFO;
    bool isHeader = true;
};

#endif 