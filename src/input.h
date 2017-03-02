#ifndef INPUT_H
#define INPUT_H

#include <vector>
#include <map>
#include <string>
#include "geometry.h"

namespace Equestria {
    extern std::vector<Polygon> polygon;
    extern std::map<std::string, int> mtlIndex;

    struct Material {
        double *brdf;
        struct MTL {
            Point Ka; // ambient
            Point Kd; // diffuse
            Point Ks; // specular
            double Ns; // specular exponent, range between 0 and 1000
            double d, Tr; // d for dissolved, Tr for transparent, d + Tr = 1
            int illum; /* 0. Color on and Ambient off
                          1. Color on and Ambient on
                          2. Highlight on
                          3. Reflection on and Ray trace on
                          4. Transparency: Glass on, Reflection: Ray trace on
                          5. Reflection: Fresnel on and Ray trace on
                          6. Transparency: Refraction on, Reflection: Fresnel off and Ray trace on
                          7. Transparency: Refraction on, Reflection: Fresnel on and Ray trace on
                          8. Reflection on and Ray trace off
                          9. Transparency: Glass on, Reflection: Ray trace off
                          10. Casts shadows onto invisible surfaces */
            std::string mapKa; // ambient texture map
            std::string mapKd; // diffuse texture map, most of time be the same as mapKa
            std::string mapKs; // specular color texture map
            std::string mapNs; // specular highlight component
            std::string mapd; // alpha texture map
            std::string mapBump; // bump map (which by default uses luminance channel of the image)
            std::string mapDisp; // displacement map
            std::string mapDecal; // stencil decal texture (defaults to 'matte' channel of the image)
        } mtl;
    };
    extern std::vector<Material> material;

    void readModel(const std::string &file);
    void objRead(const std::string &file); // read obj file and save polygons to "polygon"
    void strSplit(const std::string &str, std::vector<std::string> &ans);
    void mtlRead(const std::string &file); // read mtl file and save to "material"
}

#endif
