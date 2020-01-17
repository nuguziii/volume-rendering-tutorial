#pragma once
class Vec {

public:

    union {
        float data[3];
        struct {
            float x;
            float y;
            float z;
        };
    };

    // Constructors

    // Vectors default to 0, 0, 0.
    Vec() {
        x = 0;
        y = 0;
        z = 0;
    }

    // Construct with values, 3D
    Vec(float ax, float ay, float az) {
        x = ax;
        y = ay;
        z = az;
    }

    // Construct with values, 2D
    Vec(float ax, float ay) {
        x = ax;
        y = ay;
        z = 0;
    }

    // Copy constructor
    Vec(const Vec& o) {
        x = o.x;
        y = o.y;
        z = o.z;
    }

    // Addition

    Vec operator+(const Vec& o) {
        return Vec(x + o.x, y + o.y, z + o.z);
    }

    Vec& operator+=(const Vec& o) {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }

    Vec operator+(const float s) {
        return Vec(x + s, y + s, z + s);
    }

    // Subtraction

    Vec operator-() {
        return Vec(-x, -y, -z);
    }

    Vec operator-(const Vec o) {
        return Vec(x - o.x, y - o.y, z - o.z);
    }

    Vec& operator-=(const Vec o) {
        x -= o.x;
        y -= o.y;
        z -= o.z;
        return *this;
    }

    // Multiplication
    Vec operator*(const Vec o) {
        return Vec(x * o.x, y * o.y, z * o.z);
    }

    Vec& operator*=(const Vec o) {
        x *= o.x;
        y *= o.y;
        z *= o.z;
        return *this;
    }

    // Multiplication by scalars

    Vec operator*(const float s) {
        return Vec(x * s, y * s, z * s);
    }

    Vec& operator*=(const float s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    // Division by scalars

    Vec operator/(const float s) {
        return Vec(x / s, y / s, z / s);
    }

    Vec& operator/=(const float s) {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }

    // Dot product

    float operator^(const Vec o) {
        return (x * o.x) + (y * o.y) + (z * o.z);
    }


    // An in-place dot product does not exist because
    // the result is not a vector.

    // Cross product

    //Vec operator^(const Vec o) {
    //    float nx = y * o.z - o.y * z;
    //    float ny = z * o.x - o.z * x;
    //    float nz = x * o.y - o.x * y;
    //    return Vec(nx, ny, nz);
    //}

    //Vec& operator^=(const Vec o) {
    //    float nx = y * o.z - o.y * z;
    //    float ny = z * o.x - o.z * x;
    //    float nz = x * o.y - o.x * y;
    //    x = nx;
    //    y = ny;
    //    z = nz;
    //    return *this;
    //}

    // Other functions

    // Length of vector
    float magnitude() {
        return sqrt(magnitude_sqr());
    }

    // Length of vector squared
    float magnitude_sqr() {
        return (x * x) + (y * y) + (z * z);
    }

    // Returns a normalised copy of the vector
    // Will break if it's length is 0
    Vec normalised() {
        return Vec(*this) / magnitude();
    }

    // Modified the vector so it becomes normalised
    Vec& normalise() {
        (*this) /= magnitude();
        return *this;
    }

};