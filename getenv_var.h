#ifndef GET_ENV_VAR_H_
#define GET_ENV_VAR_H_

#include <cerrno>
#include <cstdlib>
#include <cctype>
#include <cstring>

#include <stdexcept>
#include <vector>
#include <utility>

#include <iostream>

template <typename T>
struct str_parser {
    static T parse(const char *str, char **end_ptr);
};

template <typename T>
T getenv_var(const char *var_name, const T default_value = T())
{
    const char *var = std::getenv(var_name);

    if (var == nullptr)
        return default_value;
    else {
        char *end_ptr;
        T ret{str_parser<T>::parse(var, &end_ptr)};
        if (end_ptr == var)
            throw std::runtime_error("getenv_var go invalid var content");
        
        return ret;
    }
}

template <>
struct str_parser<std::string> {
    static std::string parse(const char *str, char **end_ptr)
    {
        const std::size_t len = std::strlen(str);
        *end_ptr = (char*)str + len;
        return std::string(str);
    }
};

template <>
struct str_parser<int> {
    static int parse(const char *str, char **end_ptr)
    {
        return std::strtol(str, end_ptr, 10);
    }
};


template <>
struct str_parser<float> {
    static float parse(const char *str, char **end_ptr)
    {
        return std::strtof(str, end_ptr);
    }
};

template <>
struct str_parser<double> {
    static double parse(const char *str, char **end_ptr)
    {
        return std::strtod(str, end_ptr);
    }
};

//  str is supposed to be like v1, v2, v3,
template <typename T>
struct str_parser<std::vector<T> > {
    static std::vector<T> parse(const char *str, char **end_ptr)
    {
        std::vector<T> ret;

        const char *cur = str;
        const char *end;

        for(;;) {
            char *end;

            try {
                T value = str_parser<T>::parse(cur, &end);

                if (end == cur) {
                    *end_ptr = end;
                    return ret;
                }
                else {
                    ret.push_back(value);
                    cur = end;
                    while (std::isspace(*cur) || *cur == ',')
                        cur++;
                }
            }
            catch(...) {
                *end_ptr = end;
                return ret;
            }
        }

    }
};

//  str is supposed to be like (v1, v2)
template <typename T1, typename T2>
struct str_parser<std::pair<T1, T2> > {
    static std::pair<T1, T2> parse(const char *str, char **end_ptr)
    {
        const char *cur = str;
        *end_ptr = (char*)str; // useful if fail


        while (std::isspace(*cur))
            cur++;   //  skip space

        //  skip '('
        if (*cur != '(')
            return {};
        cur++;

        char *end;
        
        const T1 first = str_parser<T1>::parse(cur, &end);
        if (cur == end)
            return {};
        cur = end;

        while (std::isspace(*cur) || *cur == ',')
            cur++;

        const T2 second = str_parser<T2>::parse(cur, &end);
        if (cur == end)
            return {};
        cur = end;

        while (std::isspace(*cur))
            cur++;   //  skip space

        if (*cur == ')') {
            *end_ptr = (char*)cur + 1;
            return std::make_pair(first, second);
        }   
        else {
            return {};
        }
    }
};


#endif