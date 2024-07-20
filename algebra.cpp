#include <stdexcept>
#include "algebra.h"

// Default constructor
Vector2D::Vector2D() : _x{ 0.0 }, _y{ 0.0 }
{
    // std::cout << "Vector2D()" << std::endl;
}

// Parameterized constructor
Vector2D::Vector2D(double x, double y) : _x{ x }, _y{ y }
{
    // std::cout << "Vector2D(double x, double y)" << std::endl;
}

// Copy constructor
Vector2D::Vector2D(const Vector2D& other) : _x{ other._x }, _y{ other._y }
{
    // std::cout << "Vector2D(const Vector2D& other)" << std::endl;
}

// Move constructor
Vector2D::Vector2D(Vector2D&& other) noexcept : _x{ other._x }, _y{ other._y }
{
    other._x = 0.0;
    other._y = 0.0;
    // std::cout << "Vector2D(Vector2D&& other) noexcept" << std::endl;
}

// Destructor
Vector2D::~Vector2D()
{
    // std::cout << "~Vector2D()" << std::endl;
}

// Operator overloads
Vector2D Vector2D::operator+(const Vector2D& other) const
{
    return Vector2D(_x + other._x, _y + other._y);
}

Vector2D Vector2D::operator-(const Vector2D& other) const
{
    return Vector2D(_x - other._x, _y - other._y);
}

Vector2D Vector2D::operator*(double scalar) const
{
    return Vector2D(_x * scalar, _y * scalar);
}

Vector2D& Vector2D::operator=(const Vector2D& other)
{
    if (this != &other) {
        _x = other._x;
        _y = other._y;
    }
    return *this;
}

Vector2D& Vector2D::operator=(Vector2D&& other) noexcept
{
    if (this != &other) {
        _x = other._x;
        _y = other._y;
        other._x = 0.0;
        other._y = 0.0;
    }
    return *this;
}

Vector2D& Vector2D::operator+=(const Vector2D& other)
{
    _x += other._x;
    _y += other._y;
    return *this;
}

Vector2D& Vector2D::operator-=(const Vector2D& other)
{
    _x -= other._x;
    _y -= other._y;
    return *this;
}

Vector2D& Vector2D::operator*=(double scalar)
{
    _x *= scalar;
    _y *= scalar;
    return *this;
}

double& Vector2D::operator[](int index)
{
    if (index == 0)
        return _x;
    else if (index == 1)
        return _y;
    else 
        throw std::out_of_range("Index out of range");
}

const double& Vector2D::operator[](int index) const
{
    if (index == 0)
        return _x;
    else if (index == 1)
        return _y;
    else 
        throw std::out_of_range("Index out of range");
}

// Overload operator* to multiply a scalar with a Vector2D object
Vector2D operator*(double scalar, const Vector2D& vector)
{
    return Vector2D(scalar * vector[0], scalar * vector[1]);
}

void Vector2D::setZero()
{
    _x = 0.0;
    _y = 0.0;
}

// Overload operator<< to print a Vector2D object
std::ostream& operator<<(std::ostream& os, const Vector2D& vector)
{
    os << "(" << vector[0] << ", " << vector[1] << ")";
    return os;
}