// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <assert.h>
#include <stdexcept>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <math.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>


#if __i386 || _M_IX86
#define BITCT 32
#define ONEUL 1UL
#elif __x86_64__ || _M_X64
#define BITCT 64
#define ONEUL 1ULL
#else
#error "Error: Cannot detect if it is 32bit or 64 bit system"
#endif


namespace misc
{
// my own codes

template <class T>
class vec2d
{
public:
    vec2d() {}
    vec2d(size_t row, size_t col, T def)
    {
        if (row == 0 || col == 0)
        { throw std::invalid_argument("Dimension of 2d vector must be >0"); }
        m_storage.resize(row * col, def);
        m_row = row;
        m_col = col;
    }
    vec2d(size_t row, size_t col)
    {
        if (row == 0 || col == 0)
        { throw std::invalid_argument("Dimension of 2d vector must be >0"); }
        m_storage.resize(row * col);
        m_row = row;
        m_col = col;
    }
    T operator()(size_t row, size_t col) const
    {
        if (row > m_row || col > m_col)
            throw std::out_of_range("2d vector out of range!");
        return m_storage[row * m_col + col];
    }
    T& operator()(size_t row, size_t col)
    {
        if (row > m_row || col > m_col)
            throw std::out_of_range("2d vector out of range!");
        return m_storage[row * m_col + col];
    }
    void clear() { m_storage.clear(); }
    size_t rows() const { return m_row; }
    size_t cols() const { return m_col; }

private:
    size_t m_row = 0;
    size_t m_col = 0;
    std::vector<T> m_storage;
};


inline bool to_bool(const std::string& input)
{
    std::string str = input;
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    if (str.compare("T") == 0 || str.compare("TRUE") == 0)
        return true;
    else if (str.compare("F") == 0 || str.compare("FALSE") == 0)
        return false;
    else
    {
        std::string error_message = "Input is not True/False values: " + str;
        throw std::runtime_error(error_message);
    }
}

// function from John D.Cook
// https://www.johndcook.com/blog/standard_deviation/
class RunningStat
{
public:
    RunningStat() {}
    void clear()
    {
        n = 0;
        M1 = M2 = M3 = M4 = 0.0;
    }
    void push(double x)
    {
        double delta, delta_n, delta_n2, term1;

        size_t n1 = n;
        n++;
        delta = x - M1;
        assert(n > 0);
        delta_n = delta / n;
        delta_n2 = delta_n * delta_n;
        term1 = delta * delta_n * n1;
        M1 += delta_n;
        M4 += term1 * delta_n2 * (n * n - 3 * n + 3) + 6 * delta_n2 * M2
              - 4 * delta_n * M3;
        M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
        M2 += term1;
    }
    size_t get_n() const { return n; }

    double mean() const { return M1; }

    double var() const { return M2 / ((double) n - 1.0); }

    double sd() const { return sqrt(var()); }

private:
    size_t n = 0;
    double M1 = 0, M2 = 0, M3 = 0, M4 = 0;
};


// Functions from R
double dnorm(double x, double mu = 0.0, double sigma = 1.0, bool log = false);
double qnorm(double p, double mu = 0.0, double sigma = 1.0,
             bool lower_tail = true, bool log_p = false);

// codes from stackoverflow
inline std::vector<std::string> split(const std::string& seq,
                                      const std::string& separators = "\t ")
{
    std::size_t prev = 0, pos;
    std::vector<std::string> result;
    while ((pos = seq.find_first_of(separators, prev)) != std::string::npos)
    {
        if (pos > prev) result.emplace_back(seq.substr(prev, pos - prev));
        prev = pos + 1;
    }
    if (prev < seq.length())
        result.emplace_back(seq.substr(prev, std::string::npos));
    return result;
}

inline void split(std::vector<std::string>& result, const std::string& seq,
                  const std::string& separators = "\t ")
{
    std::size_t prev = 0, pos;
    result.clear();
    while ((pos = seq.find_first_of(separators, prev)) != std::string::npos)
    {
        if (pos > prev) result.emplace_back(seq.substr(prev, pos - prev));
        prev = pos + 1;
    }
    if (prev < seq.length())
        result.emplace_back(seq.substr(prev, std::string::npos));
}

template <typename T>
inline T convert(const std::string& str)
{
    std::istringstream iss(str);
    T obj;
    iss >> obj;

    if (!iss.eof() || iss.fail())
    { throw std::runtime_error("Unable to convert the input"); }
    return obj;
}

template <typename T>
inline std::string to_string(T value)
{
    std::stringstream out;
    out << value;
    return out.str();
}

// trim from start (in place)
inline void ltrim(std::string& s)
{
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))));
}
// trim from end (in place)
inline void rtrim(std::string& s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace)))
                .base(),
            s.end());
}
// trim from both ends (in place)
inline void trim(std::string& s)
{
    ltrim(s);
    rtrim(s);
}
// trim from start (copying)
inline std::string ltrimmed(std::string s)
{
    ltrim(s);
    return s;
}
// trim from end (copying)
inline std::string rtrimmed(std::string s)
{
    rtrim(s);
    return s;
}
// trim from both ends (copying)
inline std::string trimmed(std::string s)
{
    trim(s);
    return s;
}
// From http://stackoverflow.com/a/24386991/1441789
template <class T>
inline T base_name(T const& path, T const& delims = "/\\")
{
    return path.substr(path.find_last_of(delims) + 1);
}

template <class T>
inline T remove_extension(T const& filename)
{
    typename T::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != T::npos ? filename.substr(0, p) : filename;
}

inline void replace_substring(std::string& s, const std::string& search,
                              const std::string& replace)
{
    for (size_t pos = 0;; pos += replace.length())
    {
        // Locate the substring to replace
        pos = s.find(search, pos);
        if (pos == std::string::npos) break;
        // Replace by erasing and inserting
        s.erase(pos, search.length());
        s.insert(pos, replace);
    }
}

/*!
 * \brief Function to check if two double are equal from
 *        https://stackoverflow.com/a/4010279/1441789
 * \param a the first double
 * \param b the second double
 * \param error_factor level of error, should be of no concern to us at the
 *        moment
 * \return True if two double are equal
 */
inline bool logically_equal(double a, double b, double error_factor = 1.0)
{
    return ((a == b)
            || (std::abs(a - b) < std::abs(std::min(a, b))
                                      * std::numeric_limits<double>::epsilon()
                                      * error_factor));
}

}
