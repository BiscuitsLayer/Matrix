#include "Circuit.hpp"

bool operator < (Edge lhs, Edge rhs) {
    if (lhs.first == rhs.first) {
        return lhs.second < rhs.second;
    }
    return lhs.first < rhs.first;
}

bool operator != (RV lhs, RV rhs) {
    return (lhs.data_ != rhs.data_);
}

std::ostream& operator << (std::ostream& stream, RV rv) {
    stream << std::setw (0);
    stream << "(" << rv.Resistance () << "R, " << rv.Voltage () << "V)";
    return stream;
}