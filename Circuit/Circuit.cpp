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

void DFS::VertexEntry (int cur) {
    used_[cur] = true;
    curPath_.push_back (cur);
}

void DFS::VertexOutro (int cur) {
    curPath_.pop_back ();
    used_[cur] = false;
}

bool DFS::PushIfDifferent () {
    int temp = cyclesSets_.size ();
    cyclesSets_.insert ({ curPath_.begin (), curPath_.end () });
    if (cyclesSets_.size () != temp) {
        //  Size increased, => different cycle
        cycles_.push_back (curPath_);
        return true;
    }
    return false;
}

void DFS::Step (int cur, int prev) {
    if (used_[cur]) {
        //  it is a complete cycle only if the first vertex equals the last one
        if (cur != curPath_.front ()) {
            return;
        }
        curPath_.push_back (cur);
        PushIfDifferent ();
        curPath_.pop_back ();
        return;
    }
    else {
        VertexEntry (cur);
        for (int i = 0; i < table_->Shape ().first; ++i) {
            if (((*table_).At (cur, i) != RV {}) && (i != cur) && (i != prev)) {
                Step (i, cur);
            }
        }
        VertexOutro (cur);
    }
}

std::vector <std::vector <Vertex>> DFS::GetCycles () {
    for (int i = 0; i < used_.size (); ++i) {
        Step (i, -1);
    }
    return cycles_;
}

void Circuit::ComputeMaxIdx () {
    for (int i = 0; i < adjTable_.Shape ().first; ++i) {
        for (int j = 0; j < adjTable_.Shape ().second; ++j) {
            if (adjTable_.At (i, j) != RV {}) {
                Edge edge = { i, j };
                bool containsReversedEdge = false;
                int variableIdx = GetVariableIdx (edge, containsReversedEdge);
                maxIdx_ = std::max <int> (maxIdx_, variableIdx);
            }
        }                
    }
}

int Circuit::GetVariableIdx (Edge edge, bool& containsReversedEdge) {
    if (edge.first == edge.second) {
        throw std::invalid_argument ("Cyclic edge");
    }
    static int counter = 0;
    int ans = 0;
    try {
        ans = edgesToVariables_.at (edge);
        return ans;
    }
    catch (std::out_of_range &ex) {
        try {
            std::swap (edge.first, edge.second);
            ans = edgesToVariables_.at (edge);
            containsReversedEdge = true;
            return ans;
        }
        catch (std::out_of_range &ex) {
            ans = counter++;
            edgesToVariables_[edge] = ans;
        }
    }
    return ans;
}

PairMatrix Circuit::FirstKhLaw () {
    int adjTableSize = adjTable_.Shape ().first;   //  to avoid static_cast
    Linear::Matrix <double> lhs { adjTableSize, maxIdx_ + 1 };
    Linear::Matrix <double> rhs { adjTableSize, 1, 0 };
    for (int i = 0; i < adjTable_.Shape ().first; ++i) {
        for (int j = 0; j < adjTable_.Shape ().second; ++j) {
            if (adjTable_.At (i, j) != RV {}) {
                Edge edge = { i, j };
                bool containsReversedEdge = false;
                int variableIdx = GetVariableIdx (edge, containsReversedEdge);
                lhs.At (i, variableIdx) = 1 * (containsReversedEdge ? -1 : 1);
            }
        }
    }
    return { lhs, rhs };
}

PairMatrix Circuit::SecondKhLaw () {
    int cyclesSize = cycles_.size ();   //  to avoid static_cast
    Linear::Matrix <double> lhs { cyclesSize, maxIdx_ + 1 };
    Linear::Matrix <double> rhs { cyclesSize, 1, 0 };
    for (int i = 0; i < cycles_.size (); ++i) {
        for (int j = 0; j < cycles_[i].size () - 1; ++j) {
            Edge edge = { cycles_[i][j], cycles_[i][j+1] };
            bool containsReversedEdge = false;
            int variableIdx = GetVariableIdx (edge, containsReversedEdge);
            lhs.At (i, variableIdx) = adjTable_.At (edge.first, edge.second).Resistance () * (containsReversedEdge ? -1 : 1);
            rhs.At (i, 0) += adjTable_.At (edge.first, edge.second).Voltage ();
        }
    }
    return { lhs, rhs };
}

PairMatrix Circuit::Execute () {
    PairMatrix temp = FirstKhLaw ();
    lhs_ = temp.first;
    rhs_ = temp.second;
    temp = SecondKhLaw ();
    lhs_.AppendRows (temp.first);
    rhs_.AppendRows (temp.second);
    return { lhs_, rhs_ };
}