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

class Circuit {
    private:
        std::vector <std::vector <RV> > adjTable_ {};   //  TODO: Matrix можно
        std::set <std::set <Vertex> > cyclesSets_ {};   //  Auxiliary set to make cycles unique
        std::vector <std::vector <Vertex> > cycles_ {};
        std::map <Edge, int> edgesToVariables_;
        
    public:
        Circuit (const std::vector <std::vector <RV> >& adjTable):
            adjTable_ (adjTable),
            cyclesSets_ ({}),
            cycles_ ({}),
            edgesToVariables_ ({})
            {}

        void DeepFirstSearch (int cur, int prev, std::vector <bool>* used, std::vector <int>* curPath) {
            if ((*used)[cur]) {
                //  it is a complete cycle only if the first vertex equals the last one
                if (cur != (*curPath)[0]) {
                    return;
                }

                int temp = cyclesSets_.size ();
                curPath->push_back (cur);
                cyclesSets_.insert ({ curPath->begin (), curPath->end () });
                if (cyclesSets_.size () != temp) {
                    //  We have found a different cycle!
                    cycles_.push_back (*curPath);
                }
                
                curPath->pop_back ();
                return;
            }
            else {
                (*used)[cur] = true;
                curPath->push_back (cur);

                RV emptyRV = { -1, -1 };
                for (int i = 0; i < adjTable_.size (); ++i) {
                    if ((adjTable_[cur][i] != emptyRV) && (i != cur) && (i != prev)) {
                        DeepFirstSearch (i, cur, used, curPath);
                    }
                }

                curPath->pop_back ();
                (*used)[cur] = false;
            }
        }

        void GetCycles () {
            int size = adjTable_.size ();
            std::vector <bool> used (size, false);
            std::vector <int> curPath {};

            for (int i = 0; i < size; ++i) {
                DeepFirstSearch (i, -1, &used, &curPath);
            }
            //  show ans
            std::cerr << "Cycles:" << std::endl;
            for (auto& cycle : cycles_) {
                std::cout << cycle << std::endl;
            }
        }

        int GetVariable (Edge edge, bool& containsReversedEdge) {
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
            int adjTableSize = adjTable_.size ();
            std::vector <std::vector <double> > lhs (adjTableSize, std::vector <double> (adjTableSize));    //   main part of the Matrix
            std::vector <double> rhs (adjTableSize, 0);                                                     //  additional column
            int maxIdx = 0;
            RV emptyRV = { -1, -1 };
            for (int i = 0; i < adjTableSize; ++i) {
                //  calculate number of variables to resize
                for (int j = 0; j < adjTableSize; ++j) {
                    if (adjTable_[i][j] != emptyRV) {
                        Edge edge = { i, j };
                        bool containsReversedEdge = false;
                        int variableIdx = GetVariable (edge, containsReversedEdge);
                        maxIdx = std::max <int> (maxIdx, variableIdx);
                    }
                }
                //  resize
                lhs[i].resize (maxIdx + 1);
                //  write the law
                for (int j = 0; j < adjTableSize; ++j) {
                    if (adjTable_[i][j] != emptyRV) {
                        Edge edge = { i, j };
                        bool containsReversedEdge = false;
                        int variableIdx = GetVariable (edge, containsReversedEdge);
                        lhs[i][variableIdx] = 1 * (containsReversedEdge ? 1 : -1);
                    }
                }
            }
            //  show ans
            for (int i = 0; i < lhs.size (); ++i) {
                std::cerr << lhs[i] << "= " << rhs[i] << std::endl;
            }
        }

        Linear::Matrix <double> SecondKhLaw () {
            int cyclesSize = cycles_.size ();
            std::vector <std::vector <double> > lhs (cyclesSize, std::vector <double> {});       //   main part of the Matrix
            std::vector <double> rhs (cyclesSize, 0);                                            //  additional column
            int maxIdx = 0;
            for (int i = 0; i < cyclesSize; ++i) {
                //  calculate number of variables to resize
                for (int j = 0; j < cycles_[i].size () - 1; ++j) {
                    Edge edge = { cycles_[i][j], cycles_[i][j+1] };
                    bool containsReversedEdge = false;
                    int variableIdx = GetVariable (edge, containsReversedEdge);
                    maxIdx = std::max <int> (maxIdx, variableIdx);
                }
                //  resize
                lhs[i].resize (maxIdx + 1);
                //  get resistance coefficients
                for (int j = 0; j < cycles_[i].size () - 1; ++j) {
                    Edge edge = { cycles_[i][j], cycles_[i][j+1] };
                    bool containsReversedEdge = false;
                    int variableIdx = GetVariable (edge, containsReversedEdge);
                    lhs[i][variableIdx] = adjTable_[edge.first][edge.second].first * (containsReversedEdge ? 1 : -1);
                    rhs[i] += adjTable_[edge.first][edge.second].second;
                }
            }
            //  show ans
            for (int i = 0; i < lhs.size (); ++i) {
                lhs[i].resize (maxIdx + 1, 0);
                std::cerr << lhs[i] << "= " << rhs[i] << std::endl;
            }
        }

};