#include "Circuit.hpp"

bool operator < (Edge lhs, Edge rhs) {
    if (lhs.first == rhs.first) {
        return lhs.second < rhs.second;
    }
    return lhs.first < rhs.first;
}