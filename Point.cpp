#include <bits/stdc++.h>


template<typename Type>
class TPoint {
 public:
    Type x, y;
    int id;

    static constexpr Type eps = static_cast<Type>(1e-9);

    constexpr TPoint() : x(0), y(0), id(-1) {}
    template<typename T1, typename T2>
    TPoint(const T1 &x, const T2 &y) :
    x(static_cast<Type>(x)), y(static_cast<Type>(y)), id(-1) {}
    template<typename T1, typename T2, typename T3>
    TPoint(const T1 &x, const T2 &y, const T3 &id) :
    x(static_cast<Type>(x)), y(static_cast<Type>(y)), id(static_cast<int>(id)) {}

    void rotate_clockwise() { swap(x, y); y *= -1; }
    void rotate_counterclockwise() { swap(x, y); x *= -1; }

    [[nodiscard]] auto len2() const -> decltype(int64_t{}) {
        if constexpr (std::is_integral_v<Type>)
            return static_cast<long long>(x) * x + static_cast<long long>(y) * y;
        return static_cast<long double>(x) * x + static_cast<long double>(y) * y;
    }
    [[nodiscard]] long double len() const { return sqrt(len2()); }

    [[nodiscard]] bool is_upper() const { return y > eps || (abs(y) <= eps && x > eps); }

    [[nodiscard]] int cmp_polar(const TPoint &other) const {
        if (const bool a = is_upper(), b = other.is_upper(); a != b) { return a ? -1 : 1; }
        Type cross_prod = *this ^ other;
        return cross_prod > eps ? -1 : cross_prod < -eps ? 1 : 0;
    }

    tuple<Type, Type, int> operator()() const { return make_tuple(x, y, id); }

    template<typename T> TPoint &operator +=(const TPoint<T> &other) {
        x += static_cast<Type>(other.x); y += static_cast<Type>(other.y); return *this; }
    template<typename T> TPoint &operator -=(const TPoint<T> &other) {
        x -= static_cast<Type>(other.x); y -= static_cast<Type>(other.y); return *this; }
    template<typename T> TPoint &operator *=(const T &value) {
        x *= static_cast<Type>(value); y *= static_cast<Type>(value); return *this; }
    template<typename T> TPoint &operator /=(const T &value) {
        x /= static_cast<Type>(value); y /= static_cast<Type>(value); return *this; }
    TPoint operator -() { return TPoint(-x, -y, id); }

    template<typename T, typename VT> friend TPoint<T> operator *(const TPoint<T> &p, const VT &value);
    template<typename T, typename VT> friend TPoint<T> operator /(const TPoint<T> &p, const VT &value);

    template<typename T1, typename T2> friend bool operator ==(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator !=(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator <(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator >(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator <=(const TPoint<T1> &a, const TPoint<T2> &b);
    template<typename T1, typename T2> friend bool operator >=(const TPoint<T1> &a, const TPoint<T2> &b);

    template<typename T> friend TPoint<T> operator +(const TPoint<T> &a, const TPoint<T> &b);
    template<typename T> friend TPoint<T> operator -(const TPoint<T> &a, const TPoint<T> &b);

    template<typename T> friend T operator *(const TPoint<T> &a, const TPoint<T> &b);
    template<typename T> friend T operator %(const TPoint<T> &a, const TPoint<T> &b);
    template<typename T> friend T operator ^(const TPoint<T> &a, const TPoint<T> &b);

    template<typename ST, typename T> friend ST& operator >>(ST& stream, TPoint<T> &p);
    template<typename ST, typename T> friend ST& operator <<(ST& stream, TPoint<T> &p);
};

template<typename T, typename VT> TPoint<T> operator *(const TPoint<T> &p, const VT &value) {
    p.x *= value; p.y *= value; return p; }
template<typename T, typename VT> TPoint<T> operator /(const TPoint<T> &p, const VT &value) {
    p.x /= value; p.y /= value; return p; }

template<typename T1, typename T2> bool operator ==(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.x == b.x && a.y == b.y; }
template<typename T1, typename T2> bool operator !=(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.x != b.x || a.y != b.y; }
template<typename T1, typename T2> bool operator <(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.y < b.y || (a.y == b.y && a.x < b.x); }
template<typename T1, typename T2> bool operator >(const TPoint<T1> &a, const TPoint<T2> &b) {
    return a.y > b.y || (a.y == b.y && a.x > b.x); }
template<typename T1, typename T2> bool operator <=(const TPoint<T1> &a, const TPoint<T2> &b) { return !(a > b); }
template<typename T1, typename T2> bool operator >=(const TPoint<T1> &a, const TPoint<T2> &b) { return !(a < b); }

template<typename T> TPoint<T> operator +(const TPoint<T> &a, const TPoint<T> &b) {
    return TPoint<T>(a.x + b.x, a.y + b.y); }
template<typename T> TPoint<T> operator -(const TPoint<T> &a, const TPoint<T> &b) {
    return TPoint<T>(a.x - b.x, a.y - b.y); }

template<typename T> T operator *(const TPoint<T> &a, const TPoint<T> &b) { return a.x * b.x + a.y * b.y; }
template<typename T> T operator %(const TPoint<T> &a, const TPoint<T> &b) { return a.x * b.y - a.y * b.x; }
template<typename T> T operator ^(const TPoint<T> &a, const TPoint<T> &b) { return a.x * b.y - a.y * b.x; }

template<typename ST, typename T> ST& operator >>(ST& stream, TPoint<T> &p) { return stream >> p.x >> p.y; }
template<typename ST, typename T> ST& operator <<(ST& stream, TPoint<T> &p) { return stream << p.x << ' ' << p.y;}

template<typename T> string to_string(const TPoint<T> &p) { return "(" + to_string(p.x) + ", " + to_string(p.y) + ")"; }