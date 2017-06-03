#include "input.h"
#include "brdf.h"
#include <map>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace Equestria
{
    std::vector<Polygon *> polygon;
    std::map<std::string, int> mtlIndex;
    std::vector<Material> material;
    std::map<std::string, int> txIdx;
    std::vector<Texture> texture;

    void readModel(const std::string &file)
    {
        std::ifstream fin(file);
        int num;
        fin >> num; // read obj (and load mtl file according to "mtllib" in obj)
        std::string fname;
        for (int i = 0; i < num; ++i) {
            fin >> fname;
            objRead(fname);
        }
        fin.close();
        for (auto i : mtlIndex) { // read brdf
            bool ok = BRDF::read_brdf((i.first + ".brdf").c_str(), material[i.second].brdf);
            if (!ok) {
                std::cerr << "cannot read BRDF file \"" << i.first << ".brdf\"" << std::endl;
                material[i.second].brdf = NULL;
            }
        }
    }

    void rotateModel(const std::string &fname)
    {
        FILE *rotfile = fopen(fname.c_str(), "r");
        if (rotfile == NULL) {
            std::cerr << "no rotate file \"" << fname << "\"" << std::endl;
            return;
        }
        int flag;
        fscanf(rotfile, "%d", &flag);
        if (!flag)
            return;
        double dr, ax, ay, az;
        fscanf(rotfile, "%lf%lf%lf%lf", &ax, &ay, &az, &dr);
        Point axis(ax, ay, az);
        for (auto &p : polygon)
            p->rotate(dr, axis);
    }

    typedef std::vector<std::string>::iterator vsi_t;

    Point __getObjPoint(vsi_t itb, vsi_t ite)
    {
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

    void strSplit(const std::string &str, std::vector<std::string> &ans)
    {
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

    void objRead(const std::string &file)
    {
        //std::cout << "Loading " << file << " ..." << std::endl;
        std::ifstream fin(file);
        if (!fin) {
            std::cerr << "could not read object file \"" << file << "\"" << std::endl;
            return;
        }
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
            if (op[0] == '#') // comment
                ;
            else if (op == "v") // vertex
                vlist.push_back(__getObjPoint(++split.begin(), split.end()));
            else if (op == "vn") // vertex normal
                vnlist.push_back(__getObjPoint(++split.begin(), split.end()));
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
                polygon.push_back(new Polygon(pl, nl, tl, material));
            }
            else if (op == "mtllib") { // external .mtl file
                for (vsi_t i = split.begin() + 1; i != split.end(); ++i)
                    mtlRead(*i);
            }
            else if (op == "usemtl") // use material
                material = mtlIndex[split.back()];
            else if (op == "g" || op == "o" || op == "s")
                ;//std::cout << "ignore operator \"" << op << "\"" << std::endl;
            else
                std::cerr << "operator \"" << op << "\" not supported" << std::endl;
        }
        fin.close();
    }

    void mtlRead(const std::string &file)
    {
        //std::cout << "Loading " << file << " ..." << std::endl;
        std::ifstream fin(file);
        if (!fin) {
            std::cerr << "could not read material file \"" << file << "\"" << std::endl;
            return;
        }
        std::string line;
        Material *pm = NULL;
        while (getline(fin, line)) {
            std::vector<std::string> split;
            strSplit(line, split);
            if (split.empty())
                continue;
            std::string op = split.front();
            if (op == "newmtl") {
                int num = mtlIndex.size();
                mtlIndex.insert(std::make_pair(split.back(), num));
                material.push_back(Material());
                pm = &material.back();
            }
            else if (op[0] == '#') // comment
                ;
            else if (op == "Ka")
                pm->mtl.Ka = __getObjPoint(++split.begin(), split.end());
            else if (op == "Kd")
                pm->mtl.Kd = __getObjPoint(++split.begin(), split.end());
            else if (op == "Ks")
                pm->mtl.Ks = __getObjPoint(++split.begin(), split.end());
            else if (op == "Tf")
                pm->mtl.Tf = __getObjPoint(++split.begin(), split.end());
            else if (op == "Ns")
                pm->mtl.Ns = atof(split.back().c_str());
            else if (op == "Ni")
                pm->mtl.Ni = atof(split.back().c_str());
            else if (op == "d")
                pm->mtl.Tr = 1 - atof(split.back().c_str());
            else if (op == "Tr")
                pm->mtl.Tr = atof(split.back().c_str());
            else if (op == "illum")
                pm->mtl.illum = atoi(split.back().c_str());
            else if (op == "map_Ka") {
                readTexture(split.back());
                pm->mtl.mapKa = txIdx[split.back()];
            }
            else if (op == "map_Kd") {
                readTexture(split.back());
                pm->mtl.mapKd = txIdx[split.back()];
            }
            else if (op == "map_Ks" || op == "refl") {
                readTexture(split.back());
                pm->mtl.mapKs = txIdx[split.back()];
            }
            //else if (op == "map_bump" || op == "bump")
            //    pm->mtl.mapBump = split.back();
            else
                std::cerr << "operator \"" << op << "\" not supported" << std::endl;
        }
    }

    void readTexture(const std::string &file)
    {
        if (txIdx.find(file) != txIdx.end())
            return;
        txIdx.insert(std::make_pair(file, texture.size()));
        Texture tx;
        FILE *fin = fopen(file.c_str(), "rb");
        if (!fin) {
            std::cerr << "could not read texture file \"" << file << "\"" << std::endl;
            return;
        }
        char c1, c2;
        fscanf(fin, "%c%c\n", &c1, &c2);
        if (c1 != 'P' || c2 != '6') {
            std::cerr << "texture file \"" << file << "\" format error" << std::endl;
            fclose(fin);
            return;
        }
        fscanf(fin, "%d %d\n", &tx.w, &tx.h);
        int tmp;
        fscanf(fin, "%d\n", &tmp);
        if (tmp != 255) {
            std::cerr << "texture file \"" << file << "\" format error" << std::endl;
            fclose(fin);
            return;
        }
        tx.buf = new unsigned char[3 * tx.w * tx.h];
        fread(tx.buf, sizeof(unsigned char), tx.w * tx.h * 3, fin);
        fclose(fin);
        texture.push_back(tx);
    }

    Point Material::MTL::getKa(double u, double v)
    {
        return mapKa != -1 ? texture[mapKa].get(u, v) : Ka;
    }
    Point Material::MTL::getKd(double u, double v)
    {
        return mapKd != -1 ? texture[mapKd].get(u, v) : Kd;
    }
    Point Material::MTL::getKs(double u, double v)
    {
        return mapKs != -1 ? texture[mapKs].get(u, v) : Ks;
    }
}
