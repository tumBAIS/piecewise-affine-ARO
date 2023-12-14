//
// Created by simon on 12.08.21.
//

#ifndef ROBUSTOPTIMIZATION_LOGGER_H
#define ROBUSTOPTIMIZATION_LOGGER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <mutex>

namespace helpers {

std::string time_stamp();

class Logger {
public:
    explicit Logger(std::string const& filename = "");

    template<class T>
    void operator<<(T const& t) {
        _mutex.lock();
        if (not _filename.empty()) {
            std::ofstream out(_filename, std::ios_base::app);
            out << time_stamp() << ": " << t << std::endl;
        }
        std::cout << time_stamp() << ": " << t << std::endl;
        _mutex.unlock();
    }

    void warning(const std::string& msg);

    void error(const std::string& msg);

    void set_logfile(std::string const& filename);

    void set_warnings(bool warnings);

private:
    std::string _filename;
    std::mutex _mutex;
    bool _warnings = true;
};

extern Logger global_logger;

}

#endif //ROBUSTOPTIMIZATION_LOGGER_H
