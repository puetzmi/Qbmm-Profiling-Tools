/**
 * @file io_utils.hpp
 * @author M. Puetz
 * @brief This file contains useful functions for IO-operations, in particular utilities to parse command line arguments and input files.
 * @date 2022-10-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef IO_UTILS_HPP
#define IO_UTILS_HPP

#include <string>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <iomanip>
#include <chrono>


/**
 * @brief Generic function to parse command line arguments of the form 'key=value' and return value of argument with the given name.
 * 
 * @tparam T Type of the argument.
 * @param argc Number of command line arguments.
 * @param argv Command line arguments
 * @param argName Name of the command line argument whose value is requested.
 * @return T Value of the command line argument.
 */
template <class T>
T parseArgument(int argc, char *argv[], const std::string &argName);


/**
 * @brief Parse string command line arguments of the form 'key=value' and return value of argument with the given name.
 * 
 * @param argc Number of command line arguments.
 * @param argv Command line arguments
 * @param argName Name of the command line argument whose value is requested.
 * @return std::string Value of the command line argument.
 */
template <>
std::string parseArgument<std::string>(int argc, char *argv[], const std::string &argName)
{

    if (argc < 2) {
        std::string errorMessage = "No arguments provided. Exiting.";
        throw std::runtime_error(errorMessage);
    }

    for (int i=1; i<argc; i++) {

        std::string arg(argv[i]);
        size_t pos = arg.find("=");

        if (pos == std::string::npos) {
            std::string errorMessage = "Got invalid argument '" + arg + "'."
                "Parameters must have the form `KEY=VALUE`. Exiting.";
            throw std::runtime_error(errorMessage);
        }

        std::string key = arg.substr(0, pos);
        for (unsigned int i=0; i<key.size(); i++) {
            key[i] = std::tolower(key[i]);
        }

        if (key == argName) {
            return arg.substr(pos + 1);
        }
    }

    std::string errorMessage = "'" + argName + "' not found in command line arguments. Exiting.";
    throw std::runtime_error(errorMessage);

}


/**
 * @brief Parse integer command line arguments of the form 'key=value' and return value of argument with the given name.
 * 
 * @param argc Number of command line arguments.
 * @param argv Command line arguments
 * @param argName Name of the command line argument whose value is requested.
 * @return int of the command line argument.
 */
template <>
int parseArgument<int>(int argc, char *argv[], const std::string &argName)
{
    return std::stoi(parseArgument<std::string>(argc, argv, argName));
}


/**
 * @brief Parse double command line arguments of the form 'key=value' and return value of argument with the given name.
 * 
 * @param argc Number of command line arguments.
 * @param argv Command line arguments
 * @param argName Name of the command line argument whose value is requested.
 * @return double of the command line argument.
 */
template <>
double parseArgument<double>(int argc, char *argv[], const std::string &argName)
{
    return std::stod(parseArgument<std::string>(argc, argv, argName));
}


/**
 * @brief Parse boolean command line arguments of the form 'key=value' and return value of argument with the given name.
 * 
 * @param argc Number of command line arguments.
 * @param argv Command line arguments
 * @param argName Name of the command line argument whose value is requested.
 * @return bool of the command line argument.
 */
template <>
bool parseArgument<bool>(int argc, char *argv[], const std::string &argName)
{
    std::string value = parseArgument<std::string>(argc, argv, argName);
    for (unsigned int i=0; i<value.size(); i++) {
        value[i] = std::tolower(value[i]);
    }

    if (value == "true" || value == "1") {
        return true;
    }
    else if (value == "false" || value == "0") {
        return false;
    }

    std::string errMsg = "Boolean argument '" + argName + "' has in valid value '" + value + "'.";
    throw std::runtime_error(errMsg);
}


/**
 * @brief Parse input file with entries of the form 'KEY=VALUE;', where 'VALUE' may be a list.
 * 
 * @param filename Full path of the file to be read.
 * @return std::map<std::string, std::vector<std::string> > Map containing all found entries.
 */
std::map<std::string, std::vector<std::string> > parseSetupFile(const std::string &filename)
{

    std::ifstream inputFile(filename);
    std::vector<std::string> entries;

    // Check if input file exists
    if (!inputFile) {
        std::string errorMessage = "Got invalid input filename: The file '" 
            + filename + "' does not exist. Exiting.";
        throw std::runtime_error(errorMessage);
    }

    // Read file
    std::stringstream ss;
    ss << inputFile.rdbuf();
    std::string str;
    std::string tmp;
    
    // Remove comment lines
    while (std::getline(ss, tmp)) {
        size_t pos = tmp.find("//");
        size_t end = (pos == std::string::npos)  ? tmp.size() : pos;
        str += tmp.substr(0, end);
    }

    ss.clear();
    ss << str;

    // Parse entries removing whitespace
    while (std::getline(ss, tmp, ';')) {
        tmp.erase(std::remove_if(tmp.begin(), tmp.end(), ::isspace), tmp.end());
        if (tmp.size() > 0) {
            entries.push_back(tmp);
        }
    }
    inputFile.close();

    // Create map to be returned eventually
    std::map<std::string, std::vector<std::string> > setup;

    for (auto& entry : entries) {
        size_t pos = entry.find('=');

        // Check if it is a valid entry and if so what type of entry it is
        if (pos == std::string::npos) {
            std::string errorMessage = "Got invalid setup entry '" + entry + "'."
                "Entries must have the form `KEY=VALUE`.";
            throw std::runtime_error(errorMessage);
        }

        // The key is the string before '='
        std::string key = entry.substr(0, pos);
        setup[key] = std::vector<std::string>();

        // Remove brackets
        std::string value = entry.substr(pos+1, entry.size());
        size_t begin = value.find("{");
        size_t end = value.find("}");
        if (begin != std::string::npos && end != std::string::npos) {   // if opening and closing brackets were found
            value = value.substr(begin+1, end-1);
        }
        else if(begin != std::string::npos && end == std::string::npos) { // if only opening bracket was found
            std::string errorMessage = "Unmatched '{' in entry " + entry + " in input file.";
            throw std::runtime_error(errorMessage);
        }
        else if(begin == std::string::npos && end != std::string::npos) { // if only closing bracket was found
            std::string errorMessage = "Unmatched '{' in entry " + entry + " in input file.";
            throw std::runtime_error(errorMessage);
        }

        std::stringstream ss(value);
        std::string element;

        // Read all elements of entry removing whitespace and double quotes
        while (std::getline(ss, element, ',')) {
            element.erase(
                            std::remove_if(element.begin(), element.end(), 
                                [](unsigned char c){return std::isspace(c) || c == '"';}
                            ), 
                            element.end()
                        );
            if (element.size() > 0) {
                setup[key].push_back(element);
            }
        }
    }

    return setup;
}


/**
 * @brief Read numeric array from file.
 * 
 * @param[in] filename Full path of the file to be read.
 * @param[in] nRows Number of rows to be read.
 * @param[in] nCols Number of columns to be read.
 * @param[out] arrayOut Output array with allocated memory for nRows*nCols doubles.
 */
void readArrayFromFile(const std::string &filename, int nRows, int nCols, double* arrayOut) 
{

    std::ifstream infile(filename);

    if (!infile) {
        std::string errMsg = "The file with the given name '" + filename + "' does not exist.";
        throw std::runtime_error(errMsg);
    }

    std::string line;
    int i = 0;
    while (i < nRows) {
        if (!(std::getline(infile, line))) {
            std::string errMsg("Parameter `nRows` exceeds number of rows in input file.");
            throw std::runtime_error(errMsg);
        }

        std::size_t begin = line.find_first_not_of(' ');
        bool readLine = (begin != std::string::npos);       // check if line is empty
        readLine = readLine && line[begin] != '#';          // check if line is Python-style comment
        readLine = readLine && line[begin] != '/'           // check if line is C-style comment
                    && line.size() > begin + 1 && line[begin+1] != '/';

        if (readLine) {
            std::stringstream iss(line);
            for (int j=0; j<nCols; j++) {
                if (!(iss >> arrayOut[i*nCols + j])) {
                    std::string errMsg("Parameter `nCols` exceeds number of columns in input file.");
                    throw std::runtime_error(errMsg);
                }
            }
            i++;
        }
    }

    infile.close();
}


/**
 * @brief Get timestamp (date, time, timezone) as string using the system clock.
 * 
 * @return std::string Timestamp as string.
 */
std::string getTimestamp() 
{

    std::stringstream ss;
    time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ss << std::put_time(std::localtime(&now), "%F %T %Z");

    return ss.str();
}


/**
 * @brief Open output file and write header.
 * 
 * @param appName Name of the application calling this function.
 * @param outputFilename Name of the output file.
 * @param buildId The application's build ID.
 * @return std::FILE* Pointer to the output file.
 */
std::FILE * openOutputFile(const std::string &outputFilename, const std::string &appName,
                            const char *buildId)
{

    FILE * outputFile = std::fopen(outputFilename.c_str(), "w+");
    std::fprintf(outputFile, "# This file was generated by the application '%s'"
        "(build %s). [%s]\n#\n", appName.c_str(), buildId, getTimestamp().c_str());

    return outputFile;
}


/**
 * @brief Close given output file after writing footer.
 * 
 * @param outputFile Pointer to the output file
 * @return int See `std::fclose`.
 */
int closeOutputFile(std::FILE *outputFile)
{
    std::fprintf(outputFile, "#\n# End of file. [%s]", getTimestamp().c_str());
    return std::fclose(outputFile);
}

#endif // IO_UTILS_HPP