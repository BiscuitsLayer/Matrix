#pragma once

#include "../Matrix/Matrix.hpp"
#include <vector>
#include <map>
#include <set>

using Vertex = int;
using Edge = std::pair <Vertex, Vertex>;

//  So we can make them keys in std::map
bool operator < (Edge lhs, Edge rhs);

struct RV {  //  Resistance + voltage
    public:
        std::pair <double, double> data_ {};
        RV ():
            data_ ({ -1, -1 })
            {}
        RV (double resistance, double voltage):
            data_ ({ resistance, voltage })
            {}
        RV (std::pair <double, double> pair):
            data_ (pair)
            {}
            
        double& Resistance () {
            return data_.first;
        }
        double& Voltage () {
            return data_.second;
        }
};

bool operator != (RV lhs, RV rhs);
std::ostream& operator << (std::ostream& stream, RV rv);

using PairMatrix = std::pair <Linear::Matrix <double>, Linear::Matrix <double>>;

class DFS final {
    private:
        //  WORK
        Linear::Matrix <RV>* table_ {};
        std::vector <bool> used_ {};
        std::vector <int> curPath_ {};
        std::set <std::set <Vertex>> cyclesSets_ {};

        //  RESULT
        std::vector <std::vector <Vertex>> cycles_ {};

        void VertexEntry (int cur) {
            used_[cur] = true;
            curPath_.push_back (cur);
        }

        void VertexOutro (int cur) {
            curPath_.pop_back ();
            used_[cur] = false;
        }

        bool PushIfDifferent () {
            int temp = cyclesSets_.size ();
            cyclesSets_.insert ({ curPath_.begin (), curPath_.end () });
            if (cyclesSets_.size () != temp) {
                //  Size increased, => different cycle
                cycles_.push_back (curPath_);
                return true;
            }
            return false;
        }

        void Step (int cur, int prev) {
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

    public:
        DFS (Linear::Matrix <RV>* table):
            table_ (table),
            used_ (table->Shape ().first, false),
            curPath_ ({})
            {}

        std::vector <std::vector <Vertex>> GetCycles () {
            for (int i = 0; i < used_.size (); ++i) {
                Step (i, -1);
            }
            return cycles_;
        }
};

class Circuit final {
    private:
        //  GIVEN
        Linear::Matrix <RV> adjTable_ {};   //  TODO: Matrix можно

        //  TOOLS
        DFS dfs_;

        //  CALCULATED
        int maxIdx_ = 0;
        std::map <Edge, int> edgesToVariables_ {};
        std::vector <std::vector <Vertex>> cycles_ {};

        //  RESULT
        Linear::Matrix <double> lhs_ {};
        Linear::Matrix <double> rhs_ {};

        void ComputeMaxIdx () {
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

    public:
        Circuit (const Linear::Matrix <RV>& adjTable):
            adjTable_ (adjTable),
            dfs_ (&adjTable_),
            edgesToVariables_ ({}),
            cycles_ (dfs_.GetCycles ())
            {
                ComputeMaxIdx ();
            }

        int GetVariableIdx (Edge edge, bool& containsReversedEdge) {
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

        PairMatrix FirstKhLaw () {
            int adjTableSize = adjTable_.Shape ().first;   //  to avoid static_cast
            Linear::Matrix <double> lhs { adjTableSize, maxIdx_ + 1 };
            Linear::Matrix <double> rhs { adjTableSize, 1, 0 };
            for (int i = 0; i < adjTable_.Shape ().first; ++i) {
                for (int j = 0; j < adjTable_.Shape ().second; ++j) {
                    if (adjTable_.At (i, j) != RV {}) {
                        Edge edge = { i, j };
                        bool containsReversedEdge = false;
                        int variableIdx = GetVariableIdx (edge, containsReversedEdge);
                        lhs.At (i, variableIdx) = 1 * (containsReversedEdge ? 1 : -1);
                    }
                }
            }
            return { lhs, rhs };
        }

        PairMatrix SecondKhLaw () {
            int cyclesSize = cycles_.size ();   //  to avoid static_cast
            Linear::Matrix <double> lhs { cyclesSize, maxIdx_ + 1 };
            Linear::Matrix <double> rhs { cyclesSize, 1, 0 };
            for (int i = 0; i < cycles_.size (); ++i) {
                for (int j = 0; j < cycles_[i].size () - 1; ++j) {
                    Edge edge = { cycles_[i][j], cycles_[i][j+1] };
                    bool containsReversedEdge = false;
                    int variableIdx = GetVariableIdx (edge, containsReversedEdge);
                    lhs.At (i, variableIdx) = adjTable_.At (edge.first, edge.second).Resistance () * (containsReversedEdge ? 1 : -1);
                    rhs.At (i, 0) += adjTable_.At (edge.first, edge.second).Voltage ();
                }
            }
            return { lhs, rhs };
        }

        PairMatrix Execute () {
            PairMatrix temp = FirstKhLaw ();
            lhs_ = temp.first;
            rhs_ = temp.second;

            temp = SecondKhLaw ();
            lhs_.AppendRows (temp.first);
            rhs_.AppendRows (temp.second);

            return { lhs_, rhs_ };

        }

};