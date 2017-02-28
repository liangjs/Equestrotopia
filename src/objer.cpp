#include "objer.h"
#include <map>
#include <cctype>
#include <fstream>
#include <cstdlib>

namespace Equestria {
    extern std::map<std::string, int> mtlIndex; // from model.h

    typedef std::vector<std::string>::iterator vsi_t;

    Point __getObjPoint(vsi_t itb, vsi_t ite) {
        if (ite - itb >= 3) {
            double x = atof((itb++)->c_str());
            double y = atof((itb++)->c_str());
            double z = atof(itb->c_str());
            return Point(x, y, z);
        }
        double x = atof((itb++)->c_str());
        double y = atof(itb->c_str());
        return Point(x, y, 0);
    }

    void strSplit(const std::string &str, std::vector<std::string> &ans) {
        std::string now;
        ans.clear();
        for (auto c : str)
            if (isspace(c)) {
                if (!now.empty())
                    ans.push_back(now);
                now.clear();
            }
            else
                now += c;
        if (!now.empty())
            ans.push_back(now);
    }

    void objRead(const std::string &file) {
        std::ifstream fin(file);
        std::string line;
        std::vector<Point> vlist, vnlist, vtlist;
        int material = -1;
        vlist.push_back(Point());
        vnlist.push_back(Point());
        vtlist.push_back(Point());
        while (getline(fin, line)) {
            std::vector<std::string> split;
            strSplit(line, split);
            if (split.empty())
                continue;
            std::string op = split.front();
            if (op == "#") // comment
                ;
            else if (op == "v") // vertex
                vlist.push_back(__getObjPoint(++split.begin(), split.end()));
            else if (op == "vn") // vertex normal
                vnlist.push_back(__getObjPoint(++split.begin(), split.end()));
            else if (op == "vp") // Parameter space vertices
                ;
            else if (op == "vt") // texture coordinate
                vtlist.push_back(__getObjPoint(++split.begin(), split.end()));
            else if (op == "f") { // face
                std::vector<Point> pl, nl, tl;
                for (auto it = ++split.begin(); it != split.end(); ++it) {
                    std::string &tmp = *it;
                    std::string::size_type pos1 = tmp.find('/');
                    if (pos1 == std::string::npos) {
                        pl.push_back(vlist[atoi(tmp.c_str())]);
                        nl.push_back(Point());
                        tl.push_back(Point());
                    }
                    else {
                        pl.push_back(vlist[atoi(tmp.substr(0, pos1).c_str())]);
                        std::string::size_type pos2 = tmp.find('/', pos1 + 1);
                        if (pos2 == std::string::npos) {
                            tl.push_back(vtlist[atoi(tmp.substr(pos1 + 1).c_str())]);
                            nl.push_back(Point());
                        }
                        else {
                            tl.push_back(vtlist[atoi(tmp.substr(pos1 + 1, pos2 - pos1 - 1).c_str())]);
                            nl.push_back(vnlist[atoi(tmp.substr(pos2 + 1).c_str())]);
                        }
                    }
                }
                polygon.push_back(Polygon(pl, nl, tl, material));
            }
            else if (op == "mtllib")   // external .mtl file
                ;
            else if (op == "usemtl") // use material
                material = mtlIndex[*split.rbegin()];
            else if (op == "o") // object
                ;
            else if (op == "g") // group
                ;
            else if (op == "s") // smooth shading
                ;
        }
        fin.close();
    }

}
