//
// Created by simon on 08.05.21.
//

#include <sstream>
#include <numeric>
#include "helpers.h"
#include "ctime"

namespace helpers {

MyException::MyException(std::string message) : _message(std::move(message)) {}

const char* MyException::what() const throw() {
    return _message.c_str();
}

void exception_throw(std::string const& msg, Logger& logger) {
    exception_check(false, msg, logger);
}

void exception_check(const bool test, const std::string& msg, Logger & logger) {
    if (not test) {
        logger.error(msg);
    }
}

void warning_throw(std::string const& msg, Logger& logger) {
    warning_check(false, msg, logger);
}

void warning_check(bool const test, std::string const& msg, Logger& logger) {
    if (not test) {
        logger.warning(msg);
    }
}

double mean(std::vector<double> const& data) {
    return std::reduce(data.begin(), data.end())/double(data.size());
}

double variance(std::vector<double> const& data) {
    double const m = mean(data);
    return std::reduce(data.begin(), data.end(), 0,
                       [m](auto a, auto x){return a + (x-m)*(x-m);})/double(data.size());
}

double standard_deviation(std::vector<double> const& data) {
    return std::sqrt(variance(data));
}

}
