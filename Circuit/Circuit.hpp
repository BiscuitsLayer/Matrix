#pragma once

#include "../MatrixLib.hpp"
#include <vector>
#include <map>
#include <set>

using Vertex = int;
using Edge = std::pair <Vertex, Vertex>;

bool operator < (Edge lhs, Edge rhs) {
    if (lhs.first == rhs.first) {
        return lhs.second < rhs.second;
    }
    return lhs.first < rhs.first;
}

using RV = std::pair <double, double>; //   Resistance + voltage
const RV emptyRV = { -1, -1 };

class DFS final {
    private:
        //  WORK
        std::vector <std::vector <RV> >* table_ {};
        std::vector <bool> used_ {};
        std::vector <int> curPath_ {};
        std::set <std::set <Vertex> > cyclesSets_ {};

        //  RESULT
        std::vector <std::vector <Vertex> > cycles_ {};

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
                for (int i = 0; i < table_->size (); ++i) {
                    if (((*table_)[cur][i] != emptyRV) && (i != cur) && (i != prev)) {
                        Step (i, cur);
                    }
                }
                VertexOutro (cur);
            }
        }

    public:
        DFS (std::vector <std::vector <RV> >* table):
            table_ (table),
            used_ (table->size (), false),
            curPath_ ({})
            {}

        std::vector <std::vector <Vertex> > GetCycles () {
            for (int i = 0; i < used_.size (); ++i) {
                Step (i, -1);
            }
            return cycles_;
        }
};

class Circuit final {
    private:
        //  GIVEN
        std::vector <std::vector <RV> > adjTable_ {};   //  TODO: Matrix можно

        //  TOOLS
        DFS dfs_;

        //  CALCULATED
        int maxIdx_ = 0;
        std::map <Edge, int> edgesToVariables_ {};
        std::vector <std::vector <Vertex> > cycles_ {};

        //  RESULT
        Linear::Matrix <double> lhs_ {};
        Linear::Matrix <double> rhs_ {};

        void ComputeMaxIdx () {
            for (int i = 0; i < adjTable_.size (); ++i) {
                for (int j = 0; j < adjTable_[i].size (); ++j) {
                    if (adjTable_[i][j] != emptyRV) {
                        Edge edge = { i, j };
                        bool containsReversedEdge = false;
                        int variableIdx = GetVariableIdx (edge, containsReversedEdge);
                        maxIdx_ = std::max <int> (maxIdx_, variableIdx);
                    }
                }                
            }
        }

    public:
        Circuit (const std::vector <std::vector <RV> >& adjTable):
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

        Linear::Matrix <double> FirstKhLaw () {
            Linear::Matrix <double> lhs { adjTable_.size (), maxIdx_ + 1 };
            Linear::Matrix <double> rhs { adjTable_.size (), 1, 0 };
            for (int i = 0; i < adjTable_.size (); ++i) {
                for (int j = 0; j < adjTable_[i].size (); ++j) {
                    if (adjTable_[i][j] != emptyRV) {
                        Edge edge = { i, j };
                        bool containsReversedEdge = false;
                        int variableIdx = GetVariableIdx (edge, containsReversedEdge);
                        lhs.At (i, variableIdx) = 1 * (containsReversedEdge ? 1 : -1);
                    }
                }
            }
            DEBUG (lhs);
            DEBUG (rhs);
        }

        Linear::Matrix <double> SecondKhLaw () {
            Linear::Matrix <double> lhs { cycles_.size (), maxIdx_ + 1 };
            Linear::Matrix <double> rhs { cycles_.size (), 1, 0 };
            for (int i = 0; i < cycles_.size (); ++i) {
                for (int j = 0; j < cycles_[i].size () - 1; ++j) {
                    Edge edge = { cycles_[i][j], cycles_[i][j+1] };
                    bool containsReversedEdge = false;
                    int variableIdx = GetVariableIdx (edge, containsReversedEdge);
                    lhs.At (i, variableIdx) = adjTable_[edge.first][edge.second].first * (containsReversedEdge ? 1 : -1);
                    rhs.At (i, 0) += adjTable_[edge.first][edge.second].second;
                }
            }
            DEBUG (lhs);
            DEBUG (rhs);
        }

};