
#include <cstring>

#include "Laplacian2DPara.h"
#include "getenv_var.h"



#define CASE_ENUM(match, enum_case)                 \
    {                                               \
        const auto len = std::strlen(match);        \
        if (std::strncmp(str, match, len) == 0) {   \
            *end = (char*)str + len;                \
            return enum_case;                       \
        }                                           \
    }

template <>
struct str_parser<Laplacian2D::CL> {
    static Laplacian2D::CL parse(const char *str, char **end)
    {
        CASE_ENUM("dirichlet", Laplacian2D::CL::DIRICHLET)
        CASE_ENUM("Dirichlet", Laplacian2D::CL::DIRICHLET)

        CASE_ENUM("neumann_non_constant", Laplacian2D::CL::NEUMANN_NON_CONSTANT)
        CASE_ENUM("Neumann_non_constant", Laplacian2D::CL::NEUMANN_NON_CONSTANT)

        CASE_ENUM("neumann", Laplacian2D::CL::NEUMANN)
        CASE_ENUM("Neumann", Laplacian2D::CL::NEUMANN)

        throw std::runtime_error("Got Invalid CL");
    }
};

template <>
struct str_parser<Laplacian2D::Source> {
    static Laplacian2D::Source parse(const char *str, char **end)
    {
        CASE_ENUM("non", Laplacian2D::Source::NON)
        CASE_ENUM("Non", Laplacian2D::Source::NON)

        CASE_ENUM("polynomial", Laplacian2D::Source::POLYNOMIAL)
        CASE_ENUM("Polynomial", Laplacian2D::Source::POLYNOMIAL)

        CASE_ENUM("trigonometrique", Laplacian2D::Source::TRIGONOMETRIQUE)
        CASE_ENUM("Trigonometrique", Laplacian2D::Source::TRIGONOMETRIQUE)
    
        CASE_ENUM("instationnaire", Laplacian2D::Source::INSTATIONNAIRE)
        CASE_ENUM("Instationnaire", Laplacian2D::Source::INSTATIONNAIRE)

        throw std::runtime_error("Got Invalid Source");
    }
};