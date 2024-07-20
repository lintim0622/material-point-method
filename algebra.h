#pragma once
#include <iostream>

class Vector2D {
public:
    // Default constructor
    Vector2D();

    // Parameterized constructor
    Vector2D(double x, double y);

    // Copy constructor
    Vector2D(const Vector2D& other);

    // Move constructor
    Vector2D(Vector2D&& other) noexcept;

    // Destructor
    ~Vector2D();

    // Operator overloads
    Vector2D operator+(const Vector2D& other) const;
    Vector2D operator-(const Vector2D& other) const;
    Vector2D operator*(double scalar) const;
    Vector2D& operator=(const Vector2D& other);
    Vector2D& operator=(Vector2D&& other) noexcept;

    Vector2D& operator+=(const Vector2D& other);
    Vector2D& operator-=(const Vector2D& other);
    Vector2D& operator*=(double scalar);
    double& operator[](int index);
    const double& operator[](int index) const;

    // reset vector
    void setZero();

private:
    double _x;
    double _y;
};

Vector2D operator*(double scalar, const Vector2D& vector);
std::ostream& operator<<(std::ostream& os, const Vector2D& vector);