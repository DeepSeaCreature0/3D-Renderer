#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <memory>
#include <limits>
#include <algorithm>
#include <cassert>
using namespace std;

# define PI 3.141592653589793
# define infinity 1e8
# define MAX_RAY_DEPTH 5

// 3D vector for coordinate and RGB color
template<typename T>
class Vec3
{
public:
    T x{}, y{}, z{}; // x,y,z coordinates

    Vec3() = default;
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    // convert to unit vector
    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);
            x *= invNor;
            y *= invNor;
            z *= invNor;
        }
        return *this;
    }

    // operation overloading //
    
    // scalar multiplication
    Vec3 operator*(T f) const { return {x * f, y * f, z * f}; }
    // (Hadamard) multiplication of vectors
    Vec3 operator*(const Vec3& v) const { return {x * v.x, y * v.y, z * v.z}; }
    // dot-product of vectors
    T dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    // vector subtraction
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    // vector addition
    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    // inplace addition
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    // inplace cross-product
    Vec3& operator*=(const Vec3& v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
    // reverse vector direction by 180 deg
    Vec3 operator-() const { return {-x, -y, -z}; }
    // vector length square
    T length2() const { return x * x + y * y + z * z; }
    // vector length
    T length() const { return sqrt(length2()); }
    /* 
        overloads the insertion operator << for the ostream class,
        allowing objects of the Vec3 class to be output to a stream
        (like std::cout) in a custom format.
    */
    friend ostream& operator<<(ostream& os, const Vec3& v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

using Vec3f = Vec3<float>;

// sphere object - can reflect and refract light
class Sphere
{
public:
    // using  "{}" for default initialization
    Vec3f center;
    float radius{}, radius2{};
    Vec3f surfaceColor, emissionColor{};
    float transparency{}, reflection{};

    
    Sphere(const Vec3f& c, float r, const Vec3f& sc, float refl = 0, float transp = 0, const Vec3f& ec = {})
        : center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
          transparency(transp), reflection(refl)
    {}

    // return true if ray intersect with any object in its path else false
    bool intersect(const Vec3f& rayorig, const Vec3f& raydir, float& t0, float& t1) const
    {
        Vec3f l = center - rayorig; // vector pointing from rayorig to center of sphere
        float tca = l.dot(raydir);  // projection of the vector l onto the ray direction raydir
        if (tca < 0) return false;  // the sphere is behind the ray's origin

        float d2 = l.dot(l) - tca * tca; // squared distance from the sphere's center to the closest point on the ray to the center
        if (d2 > radius2) return false; // ray misses the sphere

        float thc = sqrt(radius2 - d2); // Tangent Half-Chord: half the length of the chord made by the ray when passing through sphere
        t0 = tca - thc; // distance from rayorig to point of contact when entering sphere
        t1 = tca + thc; // distance from rayorig to point of contact to sphere when exiting sphere
        
        return true;
    }
};

float mix(float a, float b, float mix)
{
    return b * mix + a * (1 - mix);
}

// trace ray path
Vec3f trace(const Vec3f& rayorig, const Vec3f& raydir, const vector<Sphere>& spheres, int depth)
{
    float tnear = infinity; 
    const Sphere* sphere = nullptr;
    // Find Closest Intersection
    for (const auto& sph : spheres) {
        float t0 = infinity, t1 = infinity;
        if (sph.intersect(rayorig, raydir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &sph;
            }
        }
    }

    if (!sphere) return {2};

    Vec3f surfaceColor = {0};
    Vec3f phit = rayorig + raydir * tnear; // hit point
    Vec3f nhit = phit - sphere->center; // normal at the hit point
    nhit.normalize();

    float bias = 1e-4;  // A small offset to avoid self-intersection
    bool inside = false;
    // check if the ray is inside the sphere
    if (raydir.dot(nhit) > 0) { 
        nhit = -nhit;
        inside = true;
    }

    // check if the sphere is reflective or transparent and trace reflectionand refraction rays
    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
        float facingratio = -raydir.dot(nhit);  // cosine of the angle between the ray direction and the normal
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1); // change the mix value to tweak the effect
        Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);   // The reflection direction
        refldir.normalize();

        Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
        Vec3f refraction = {0};

        if (sphere->transparency) {
            float ior = 1.1;    // Index of refraction of the material.
            float eta = inside ? ior : 1 / ior;     //The ratio of refractive indices
            float cosi = -nhit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi); // Determines if refraction is possible (k >= 0). Calculated using Snell's law
            Vec3f refrdir = raydir * eta + nhit * (eta * cosi - sqrt(k));   // The refraction direction
            refrdir.normalize();
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
        }

        // mix of reflection and refraction colors
        surfaceColor = (
            reflection * fresneleffect +
            refraction * (1 - fresneleffect) * sphere->transparency
        ) * sphere->surfaceColor;
    }
    // If the sphere is not reflective or transparent, we compute direct illumination from light sources
    else {
        // check if its a light source
        for (const auto& sph : spheres) {
            if (sph.emissionColor.x > 0) {
                Vec3f transmission = {1};   // Represents whether the light reaches the surface
                Vec3f lightDirection = sph.center - phit;   //direction from the intersection point (phit) to the light source
                lightDirection.normalize();

                //  Check for shadows
                for (const auto& j : spheres) {
                    if (&sph != &j) {
                        float t0, t1;
                        if (j.intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                            transmission = {0};
                            break;
                        }
                    }
                }

                surfaceColor += sphere->surfaceColor * transmission *
                max(0.0f, nhit.dot(lightDirection)) * sph.emissionColor;    // nhit.dot(lightDirection)=cosine
            }
        }
    }
    
    return surfaceColor + sphere->emissionColor;
}

void render(const vector<Sphere>& spheres)
{
    int width = 640, height = 480;
    auto image = make_unique<Vec3f[]>(width * height);
    auto pixel = image.get();
    float invWidth = 1.0f / width, invHeight = 1.0f / height;
    float fov = 50, aspectratio = width / static_cast<float>(height);   // field of view angle and ratio of the image's width to its height
    float angle = tan(PI * 0.5f * fov / 180.0f);    // determine the extent of the camera's view in the image plane

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5f) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5f) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace({0}, raydir, spheres, 0);
        }
    }

    ofstream ofs("raytracer_output.ppm", ios::out | ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height; ++i) {
        ofs << static_cast<unsigned char>(min(1.0f, image[i].x) * 255) <<
               static_cast<unsigned char>(min(1.0f, image[i].y) * 255) <<
               static_cast<unsigned char>(min(1.0f, image[i].z) * 255);
    }
}

int main()
{
    vector<Sphere> spheres;

    spheres.emplace_back(Vec3f(0.0f, -10004, -20), 10000, Vec3f(0.0f, 1.0f, 0.0f)); // ground sphere
    spheres.emplace_back(Vec3f(0.0f, 0, -20), 4, Vec3f(1.00f, 0.32f, 0.36f), 1, 0.5f);
    spheres.emplace_back(Vec3f(5.0f, -1, -15), 2, Vec3f(0.90f, 0.76f, 0.46f));
    spheres.emplace_back(Vec3f(5.0f, 0, -25), 3, Vec3f(0.65f, 0.77f, 0.97f));
    spheres.emplace_back(Vec3f(-5.5f, 0, -15), 3, Vec3f(0.90f, 0.90f, 0.90f));
    spheres.emplace_back(Vec3f(0.0f, 10, 30), 3, Vec3f(0.00f, 0.00f, 0.00f), 0, 0.0f, Vec3f(1.5));  //light source

    render(spheres);
    
    return 0;
}

