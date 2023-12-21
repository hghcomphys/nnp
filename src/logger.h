/*
  logger.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Hossein Ghorbanfekr

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <sstream>

enum LOG_t
{
    DEBUG,
    INFO,
    WARN,
    ERROR
};

class Log
{
public:
    Log();
    Log(LOG_t level);
    ~Log();
    Log &operator<<(const std::string &message);
    Log &operator<<(const char *message);
    Log &operator<<(int message);
    Log &operator<<(double message);
    const std::string toString() const;
    void clear();
    void endl();

private:
    bool isOpened;
    LOG_t messageLevel;
    std::stringstream buffer;
    std::string getLabel(LOG_t level);
    void setLabel();
    bool isBlank();
    // TODO: set outside the class
    LOG_t verbosity = INFO;
    bool isHeader = true;
};

#endif