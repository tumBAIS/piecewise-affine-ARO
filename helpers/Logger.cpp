//
// Created by simon on 12.08.21.
//

#include "Logger.h"
#include "helpers.h"

namespace helpers{

std::string time_stamp() {
    std::time_t time_now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm* now = std::localtime(&time_now);

    std::stringstream s;
    s << now->tm_year + 1900 << "_" << now->tm_mon + 1 << "_" << now->tm_mday << "_"
    << now->tm_hour << "_" << now->tm_min << "_" << now->tm_sec;
    return s.str();
}

Logger::Logger(std::string const& filename) : _filename(filename) {}

void Logger::warning(std::string const& msg) {
    if(_warnings){
        *this << "WARNING: " + msg;
    }
}

void Logger::error(std::string const& msg) {
    *this << "ERROR: " + msg;
    throw MyException(msg);
}

void Logger::set_logfile(std::string const& filename) {
    _filename = filename;
}

void Logger::set_warnings(bool warnings) {
    _warnings = warnings;
}

Logger global_logger;

}
