/*
 logger.cpp: This file is part of Free Molecular Dynamics

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

#include "logger.h"
#include <stdexcept>
#include <string>
#include <iostream>

/* ----------------------------------------------------------------------
   setup for class Logger
------------------------------------------------------------------------- */
Log::Log() : Log(DEBUG) {}

Log::Log(LOG_t level) : isOpened(false), messageLevel(level) { setLabel(); }

Log::~Log()
{
    if (isOpened)
    {
        endl(); // go next line
        switch (messageLevel)
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

Log &Log::operator<<(const std::string &message)
{
    if (messageLevel >= verbosity)
    {
        buffer << message;
        isOpened = true;
    }
    return *this;
}

Log &Log::operator<<(const char *message)
{
    if (messageLevel >= verbosity)
    {
        buffer << message;
        isOpened = true;
    }
    return *this;
}

Log &Log::operator<<(int message)
{
    operator<<(std::to_string(message));
    return *this;
}

Log &Log::operator<<(double message)
{
    operator<<(std::to_string(message));
    return *this;
}

inline std::string Log::getLabel(LOG_t level)
{
    std::string label;
    switch (level)
    {
    case DEBUG:
        label = "DEBUG";
        break;
    case INFO:
        label = "INFO ";
        break;
    case WARN:
        label = "WARN ";
        break;
    case ERROR:
        label = "ERROR";
        break;
    }
    return label;
}

void Log::setLabel()
{
    if (isHeader)
    {
        operator<<("[" + getLabel(messageLevel) + "] ");
    }
}

const std::string Log::toString() const { return buffer.str(); }

void Log::clear()
{
    buffer.clear();
    setLabel();
}

void Log::endl() { buffer << "\n"; }