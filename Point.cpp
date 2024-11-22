#include <bits/stdc++.h>


template<typename Type>
class TPoint {
 public:
    Type x, y;
    int id;

    // the permissible error of calculations
    static constexpr Type eps = static_cast<Type>(1e-9);

    // initialization
    constexpr TPoint() : x(0), y(0), id(-1) {}
    template<typename T1, typename T2>
    TPoint(const T1 &x, const T2 &y) :
    x(static_cast<Type>(x)), y(static_cast<Type>(y)), id(-1) {}
    template<typename T1, typename T2, typename T3>
    TPoint(const T1 &x, const T2 &y, const T3 &id) :
    x(static_cast<Type>(x)), y(static_cast<Type>(y)), id(static_cast<int>(id)) {}

    // 90º rotation
    void rotate_clockwise() { swap(x, y); y *= -1; }
    void rotate_counterclockwise() { swap(x, y); x *= -1; }

    // the square of the length (can be useful with integer types)
    [[nodiscard]] auto len2() const {
        if constexpr (std::is_integral_v<Type>) {
            return static_cast<long long>(x) * x + static_cast<long long>(y) * y;
        } else {
            return static_cast<long double>(x) * x + static_cast<long double>(y) * y;
        }
    }
    [[nodiscard]] long double len() const { return sqrt(len2()); }

    [[nodiscard]] long double dist_to_point(const TPoint &other) const { return (*this - other).len(); }
    [[nodiscard]] long double dist_to_line(const TPoint &p1, const TPoint &p2) const {
#ifdef _GLIBCXX_DEBUG
        if (p1 == p2) throw std::runtime_error("A line cannot be defined by a single point\n");
#endif
        return abs((p1 - *this) ^ (p2 - *this)) / (p1 - p2).len();
    }
    [[nodiscard]] long double dist_to_ray(const TPoint &start, const TPoint &other) const {
#ifdef _GLIBCXX_DEBUG
        if (start == other) throw std::runtime_error("A ray cannot be defined by a single point\n");
#endif
        if ((*this - start) * (other - start) >= 0)
            return dist_to_line(start, other);
        return dist_to_point(start);
    }
    [[nodiscard]] long double dist_to_section(const TPoint &p1, const TPoint &p2) const {
        if ((*this - p1) * (p2 - p1) >= 0 && (*this - p2) * (p1 - p2) >= 0) {
            if (p1 == p2) return dist_to_point(p1);
            return dist_to_line(p1, p2);
        }
        return min(dist_to_point(p1), dist_to_point(p2));
    }

    // check for being in the upper half-plane
    [[nodiscard]] bool is_upper() const { return y > eps || (abs(y) <= eps && x > eps); }

    // polar angle comparison
    // -1 - this point angle is smaller
    //  0 - they are equal
    //  1 - this point angle is bigger
    [[nodiscard]] int cmp_polar(const TPoint &other) const {
#ifdef _GLIBCXX_DEBUG
        if (!*this || !other) throw std::runtime_error("The angle is undefined for the point (0; 0)\n");
#endif
        if (const bool a = is_upper(), b = other.is_upper(); a != b) { return a ? -1 : 1; }
        Type cross_prod = *this ^ other;
        return cross_prod > eps ? -1 : cross_prod < -eps ? 1 : 0;
    }

    std::tuple<Type, Type, int> operator()() const { return make_tuple(x, y, id); }

    template<typename T> TPoint &operator +=(const TPoint<T> &other) {
        x += static_cast<Type>(other.x); y += static_cast<Type>(other.y); return *this;
    }
    template<typename T> TPoint &operator -=(const TPoint<T> &other) {
        x -= static_cast<Type>(other.x); y -= static_cast<Type>(other.y); return *this;
    }
    template<typename T> TPoint &operator *=(const T &value) {
        x *= static_cast<Type>(value); y *= static_cast<Type>(value); return *this;
    }
    template<typename T> TPoint &operator /=(const T &value) {
        x /= static_cast<Type>(value); y /= static_cast<Type>(value); return *this;
    }
    TPoint operator -() const { return TPoint(-x, -y, id); }
    bool operator !() const { return abs(x) <= eps && abs(y) <= eps; }

    template<typename VT> TPoint operator *(const VT &value) const { return TPoint(x * value, y * value); }
    template<typename VT> TPoint operator /(const VT &value) const { return TPoint(x / value, y / value); }

    auto operator *(const TPoint &other) const {
        if constexpr (std::is_integral_v<Type>) {
            return static_cast<long long>(x) * other.x + static_cast<long long>(y) * other.y;
        } else {
            return static_cast<long double>(x) * other.x + static_cast<long double>(y) * other.y;
        }
    }
    auto operator %(const TPoint &other) const {
        if constexpr (std::is_integral_v<Type>) {
            return static_cast<long long>(x) * other.y - static_cast<long long>(y) * other.x;
        } else {
            return static_cast<long double>(x) * other.y - static_cast<long double>(y) * other.x;
        }
    }
    auto operator ^(const TPoint &other) const {
        if constexpr (std::is_integral_v<Type>) {
            return static_cast<long long>(x) * other.y - static_cast<long long>(y) * other.x;
        } else {
            return static_cast<long double>(x) * other.y - static_cast<long double>(y) * other.x;
        }
    }

    // the comparison of points is all follows:
    // the lowest point is considered the smallest, of which the leftmost is considered
    template<typename T1, typename T2> friend bool operator ==(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator !=(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator <(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator >(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator <=(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator >=(const TPoint<T1> &a, const TPoint<T2> &b);

    template<typename ST, typename T> friend ST& operator >>(ST& stream, TPoint<T> &p);
    template<typename ST, typename T> friend ST& operator <<(ST& stream, const TPoint<T> &p);
};

template<typename T1, typename T2> bool operator ==(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.x == b.x && a.y == b.y;
}
template<typename T1, typename T2> bool operator !=(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.x != b.x || a.y != b.y;
}
template<typename T1, typename T2> bool operator <(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.y < b.y || (a.y == b.y && a.x < b.x);
}
template<typename T1, typename T2> bool operator >(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.y > b.y || (a.y == b.y && a.x > b.x);
}
template<typename T1, typename T2> bool operator <=(const TPoint<T1> &a, const TPoint<T2> &b) { return !(a > b); }
template<typename T1, typename T2> bool operator >=(const TPoint<T1> &a, const TPoint<T2> &b) { return !(a < b); }

template<typename T> TPoint<T> operator +(const TPoint<T> &a, const TPoint<T> &b) {
    return TPoint<T>(a.x + b.x, a.y + b.y);
}
template<typename T> TPoint<T> operator -(const TPoint<T> &a, const TPoint<T> &b) {
    return TPoint<T>(a.x - b.x, a.y - b.y);
}

template<typename ST, typename T> ST& operator >>(ST& stream, TPoint<T> &p) { return stream >> p.x >> p.y; }
template<typename ST, typename T> ST& operator <<(ST& stream, const TPoint<T> &p) { return stream << p.x << ' ' << p.y; }

template<typename T> std::string to_string(const TPoint<T> &p) {
    return "(" + to_string(p.x) + ", " + to_string(p.y) + ")";
}
