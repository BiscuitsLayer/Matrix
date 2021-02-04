#pragma once

//  SYSTEM
#include <vector>
#include <map>
#include <set>

//  MATRIX
#include "../Matrix/Matrix.hpp"

//  TYPEDEFS
using Vertex = int;
using Edge = std::pair <Vertex, Vertex>;
using PairMatrix = std::pair <Linear::Matrix <double>, Linear::Matrix <double>>;

//  So we can make them keys in std::map
bool operator < (Edge lhs, Edge rhs);

struct RV {  //  Resistance + voltage
    public:
        //  DATA
        std::pair <double, double> data_ {};

        //  CTORS
        RV ():
            data_ ({ -1, -1 })
            {}
        RV (double resistance, double voltage):
            data_ ({ resistance, voltage })
            {}
        RV (std::pair <double, double> pair):
            data_ (pair)
            {}

        //  GETTERS
        double& Resistance ()   { return data_.first; }
        double& Voltage ()      { return data_.second; }
};

//  OVERLOADED OPERATORS
bool operator != (RV lhs, RV rhs);
std::ostream& operator << (std::ostream& stream, RV rv);

class DFS final {
    private:
        //  COMPUTATIONS
        Linear::Matrix <RV>* table_ {};
        std::vector <bool> used_ {};
        std::vector <int> curPath_ {};
        std::set <std::set <Vertex>> cyclesSets_ {};

        //  RESULT
        std::vector <std::vector <Vertex>> cycles_ {};

        //  VERTEX AND PATH METHODS
        void VertexEntry (int cur);
        void VertexOutro (int cur);
        bool PushIfDifferent ();

        //  ALGORITHM
        void Step (int cur, int prev);
    public:
        //  CTOR
        DFS (Linear::Matrix <RV>* table):
            table_ (table),
            used_ (table->Shape ().first, false),
            curPath_ ({})
            {}

        //  CYCLES
        std::vector <std::vector <Vertex>> GetCycles ();
};

class Circuit final {
    private:
        //  GIVEN
        Linear::Matrix <RV> adjTable_ {};

        //  TOOLS
        DFS dfs_;

        //  COMPUTATIONS
        int maxIdx_ = 0;
        std::map <Edge, int> edgesToVariables_ {};
        std::vector <std::vector <Vertex>> cycles_ {};
        void ComputeMaxIdx ();

        //  RESULT
        Linear::Matrix <double> lhs_ {};
        Linear::Matrix <double> rhs_ {};
    public:
        //  CTOR
        Circuit (const Linear::Matrix <RV>& adjTable):
            adjTable_ (adjTable),
            dfs_ (&adjTable_),
            edgesToVariables_ ({}),
            cycles_ (dfs_.GetCycles ())
            {
                ComputeMaxIdx ();
            }

        //  VARIABLE INDEX
        int GetVariableIdx (Edge edge, bool& containsReversedEdge);

        //  KIRCHHOFF'S LAWS
        PairMatrix FirstKhLaw ();
        PairMatrix SecondKhLaw ();

        //  EXECUTE
        PairMatrix Execute ();

};