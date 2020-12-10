/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Nikhil Johny Karuthedath
 * *************************
*/

#ifdef WIN32
    #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
    #include <GL/gl.h>
    #include <GL/glut.h>
#elif defined(__APPLE__)
    #include <OpenGL/gl.h>
    #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <imageIO.h>
#include <iostream>
#include <vector>
#include <cmath>

#include <string.h>
#ifdef WIN32
    #define strcasecmp _stricmp
#endif


using namespace std;

#define MAX_TRIANGLES 200000
#define MAX_SPHERES 1000
#define MAX_LIGHTS 1000

char * filename = NULL;

// different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

// the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

double aspect_ratio = ((double) WIDTH) / ((double)HEIGHT);

#define PI 3.14159265
#define LEAST_POSSIBLE_T 1e8

// Soft shadow 
#define EXTRA_LIGHTS 4.0

// Antialiasing
#define RAY_SAMPLES 1.0

bool enableSoftShadows = false;

// Color struct
struct Color {
  double r, g, b;

  Color() { r = g = b = 0; }

  Color(double r1, double g1, double b1): r(r1), g(g1), b(b1) {}

  Color addColor(Color color) {
    return Color(r + color.r, g + color.g, b + color.b);
  }

  Color multiplyColor(double scale) {
    return Color(r * scale, g * scale, b * scale);
  }

  Color multiplyColor(Color color) {
    return Color(r * color.r, g * color.g, b * color.b);
  }

  void clampColorTo1() {
    if (r > 1.0) r = 1.0;
    if (g > 1.0) g = 1.0;
    if (b > 1.0) b = 1.0;
  }
};

// Vector 3D struct
struct Vector3D {
  double x, y, z;

  Vector3D() { x = y = z = 0.0; }

  Vector3D(double x1, double y1, double z1): x(x1), y(y1), z(z1) {}
  
  Vector3D negate() { 
    return scalarMultiply(-1.0f); 
  }

  Vector3D addVector(Vector3D vector) {
    return Vector3D(x + vector.x, y + vector.y, z + vector.z);
  }

  Vector3D subtract(Vector3D vector) {
    return addVector(vector.negate());
  }

  Vector3D scalarMultiply(double scalar) {
    return Vector3D(scalar * x, scalar * y, scalar * z);
  }

  double length(){
    return sqrt(x * x + y * y + z * z);
  }

  Vector3D normalize() {
    return scalarMultiply(1.0f / length());
  }

  double dotProduct(Vector3D vector) {
    return x * vector.x + y * vector.y + z * vector.z;
  }

  Vector3D crossProduct(Vector3D vector) {
    return Vector3D(y * vector.z - z * vector.y, -(x * vector.z - z * vector.x), x * vector.y - y * vector.x);
  }
};

struct Vertex {
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle {
  Vertex v[3];
};

struct Sphere{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// Get triangle vertex
Vector3D getTriangleVertex(int triangleIndex, int vertex) {
  return Vector3D(triangles[triangleIndex].v[vertex].position[0],
            triangles[triangleIndex].v[vertex].position[1],
            triangles[triangleIndex].v[vertex].position[2]);
}

// Get triangle vertices
void getTriangleVertices(int triangleIndex, Vector3D &v0, Vector3D &v1, Vector3D &v2) {
  v0 = getTriangleVertex(triangleIndex, 0);
  v1 = getTriangleVertex(triangleIndex, 1);
  v2 = getTriangleVertex(triangleIndex, 2);
}

// Get light position
Vector3D getLightPosition(int lightIndex) {
  return Vector3D(lights[lightIndex].position[0],
          lights[lightIndex].position[1],
          lights[lightIndex].position[2]);
}

// Get light color
Color getLightColor(int lightIndex) {
  return Color(lights[lightIndex].color[0],
         lights[lightIndex].color[1],
         lights[lightIndex].color[2]);
}

// Calculate shadow ray = light position - intersection point
Vector3D getShadowRay(Vector3D intersectionPoint, int lightIndex) {
  return getLightPosition(lightIndex).subtract(intersectionPoint).normalize();
}

// Calculate the normal of a triangle
Vector3D getTriangleNormal(int triangleIndex) {
  Vector3D v0, v1, v2;
  getTriangleVertices(triangleIndex, v0, v1, v2);

  Vector3D v0v1 = v1.subtract(v0);
  Vector3D v0v2 = v2.subtract(v0);
  
  return v0v1.crossProduct(v0v2);
}

// Calculate reflection vector, r = 2 * (l . n) * n - l
Vector3D getReflectionVector(Vector3D l, Vector3D n) {
  return (n.scalarMultiply(2 * l.dotProduct(n)).subtract(l)).normalize();
}

// Calculate interpolated color = alpha * c0 + beta * c1 + gamma * c2
Color getInterpolatedColor (Color c0, Color c1, Color c2, double alpha, double beta, double gamma) {
  return c0.multiplyColor(alpha).addColor(c1.multiplyColor(beta)).addColor(c2.multiplyColor(gamma));
}

// Using barycentric coordinate find alpha, beta and gamma and check if intersection point is inside the triangle
bool isTriangleIntersected(int triangleIndex, Vector3D ray, double tValue, Vector3D origin, double &alpha, double &beta, double &gamma) {
  // Get triangle details
  Vector3D normal = getTriangleNormal(triangleIndex);
  Vector3D v0, v1, v2;
  getTriangleVertices(triangleIndex, v0, v1, v2);

  Vector3D intersectionPoint = origin.addVector(ray.scalarMultiply(tValue));

  // Total area
  Vector3D v2ToV0 = v2.subtract(v0);
  Vector3D v1ToV0 = v1.subtract(v0);
  Vector3D perpendicular = v2ToV0.crossProduct(v1ToV0);
  double area = normal.dotProduct(perpendicular);

  // Area of V2V1P
  Vector3D v2ToIntersectionPoint = intersectionPoint.subtract(v2);
  Vector3D v1ToIntersectionPoint = intersectionPoint.subtract(v1);
  perpendicular = v2ToIntersectionPoint.crossProduct(v1ToIntersectionPoint);
  double area1 = normal.dotProduct(perpendicular);
  alpha = area1 / area;
  if (alpha < 0) return false;

  // Area of V0V2P
  Vector3D v0ToIntersectionPoint = intersectionPoint.subtract(v0);
  v2ToIntersectionPoint = intersectionPoint.subtract(v2);
  perpendicular = v0ToIntersectionPoint.crossProduct(v2ToIntersectionPoint);
  double area2 = normal.dotProduct(perpendicular);
  beta = area2 / area;
  if (beta < 0) return false;

  gamma = 1 - alpha - beta;
  if (gamma < 0) return false;
 
  return true;
}

// Get sphere details based on index
void getSphereDetails (int sphereIndex, double &radius, double &xc, double &yc, double &zc) {
  radius = spheres[sphereIndex].radius;
  xc = spheres[sphereIndex].position[0];
  yc = spheres[sphereIndex].position[1];
  zc = spheres[sphereIndex].position[2];
}

// Check if the shadow ray is blocked by a sphere or triangle
bool isShadowRayBlocked(Vector3D shadowRay, Vector3D intersection, int lightIndex) {
  double epsilon = 1e-8;
  Vector3D light = getLightPosition(lightIndex);

  // Distance between the light and intersection
  double distanceIntersectionToLight = light.subtract(intersection).length();

  // Loop over all spheres
  for (int sphereIndex = 0; sphereIndex < num_spheres; sphereIndex++) {
    double xc, yc, zc, radius;
    getSphereDetails(sphereIndex, radius, xc, yc, zc);

    double b = 2 * (shadowRay.x * (intersection.x - xc) + shadowRay.y * (intersection.y - yc) + shadowRay.z * (intersection.z - zc));
    double c = pow(intersection.x - xc, 2.0) + pow(intersection.y - yc, 2.0) + pow(intersection.z - zc, 2.0) - pow(radius, 2.0);
    double root = b * b - 4 * c;

    // If there are no intersections, continue
    if (root < 0) continue;

    // If there are intersections, return true
    double t0 = (- b + sqrt(root)) / 2.0;
    double t1 = (- b - sqrt(root)) / 2.0;
    if (t0 > epsilon && shadowRay.scalarMultiply(t0).length() < distanceIntersectionToLight) return true;
    if (t1 > epsilon && shadowRay.scalarMultiply(t1).length() < distanceIntersectionToLight) return true;
  }

  // Loop over all triangles
  for (int triangleIndex = 0; triangleIndex < num_triangles; triangleIndex++) {
    Vector3D n = getTriangleNormal(triangleIndex).normalize();
    double nDotShadowRay = n.dotProduct(shadowRay);

    // Ray and the plane are parallel
    if (-epsilon < nDotShadowRay && nDotShadowRay < epsilon) continue;

    // get t for intersection with the plane
    Vector3D v0 = getTriangleVertex(triangleIndex, 0);
    double d = n.dotProduct(v0);
    double tValue = (d - n.dotProduct(intersection)) / nDotShadowRay;

    // The intersection is behind ray origin
    if (tValue <= epsilon) continue;
      
    // Check if the triangle in the plane is intersected
    double alphaTemp, betaTemp, gammaTemp;
    bool isIntersected = isTriangleIntersected(triangleIndex, shadowRay, tValue, intersection, alphaTemp, betaTemp, gammaTemp);
    if (tValue > epsilon &&
      shadowRay.scalarMultiply(tValue).length() < distanceIntersectionToLight && isIntersected) return true;
  }

  return false;
}

// Calculate sphere color using phong model
Color getSphereColor(Vector3D ray, Vector3D shadowRay, int lightIndex, double minT, int minTIndex) {
  double xc, yc, zc, radius;
  Vector3D normal, reflection;

  // Get the sphere and its details
  getSphereDetails(minTIndex, radius, xc, yc, zc);
  normal = ray.scalarMultiply(minT).subtract(Vector3D(xc, yc, zc)).scalarMultiply(1 / radius);
  reflection = getReflectionVector(shadowRay, normal);

  // Calculate the diffuse color
  double lDotN = shadowRay.dotProduct(normal);
  if (lDotN < 0) lDotN = 0.0;

  Color diffuse = Color(spheres[minTIndex].color_diffuse[0], spheres[minTIndex].color_diffuse[1], spheres[minTIndex].color_diffuse[2]);
  Color finalDiffuse = diffuse.multiplyColor(lDotN);

  // Calculate the specular color
  double rDotV = reflection.dotProduct(ray.negate());
  if (rDotV < 0) rDotV = 0.0;

  Color specular = Color(spheres[minTIndex].color_specular[0], spheres[minTIndex].color_specular[1], spheres[minTIndex].color_specular[2]);
  Color finalSpecular = specular.multiplyColor(pow(rDotV, spheres[minTIndex].shininess));

  // Calculate phong model
  Color color = getLightColor(lightIndex).multiplyColor(finalDiffuse.addColor(finalSpecular));
  return color;
}

// Calculate triangle color using phong model
Color getTriangleColor(Vector3D ray, Vector3D shadowRay, int lightIndex, int minTIndex, double alpha, double beta, double gamma) {
  Vector3D n0, n1, n2;
  n0 = Vector3D(triangles[minTIndex].v[0].normal[0], triangles[minTIndex].v[0].normal[1], triangles[minTIndex].v[0].normal[2]);
  n1 = Vector3D(triangles[minTIndex].v[1].normal[0],triangles[minTIndex].v[1].normal[1],triangles[minTIndex].v[1].normal[2]);
  n2 = Vector3D(triangles[minTIndex].v[2].normal[0],triangles[minTIndex].v[2].normal[1],triangles[minTIndex].v[2].normal[2]);

  // Get final normal = alpha * n0 + beta * n1 + gamma * n2
  Vector3D normal = n0.scalarMultiply(alpha).addVector(n1.scalarMultiply(beta)).addVector(n2.scalarMultiply(gamma)).normalize();

  // Get reflection vector
  Vector3D reflection = getReflectionVector(shadowRay, normal);

  // Calculate the diffuse color
  Color v0Kd, v1Kd, v2Kd;
  v0Kd = Color(triangles[minTIndex].v[0].color_diffuse[0], triangles[minTIndex].v[0].color_diffuse[1], triangles[minTIndex].v[0].color_diffuse[2]);
  v1Kd = Color(triangles[minTIndex].v[1].color_diffuse[0], triangles[minTIndex].v[1].color_diffuse[1], triangles[minTIndex].v[1].color_diffuse[2]);
  v2Kd = Color(triangles[minTIndex].v[2].color_diffuse[0], triangles[minTIndex].v[2].color_diffuse[1], triangles[minTIndex].v[2].color_diffuse[2]);

  double ln = shadowRay.dotProduct(normal);
  if (ln < 0) ln = 0;
  if (ln > 1) ln = 1.0;
  Color kd = getInterpolatedColor(v0Kd, v1Kd, v2Kd, alpha, beta, gamma);
  Color diffuse = kd.multiplyColor(ln);;

  // Calculate the specular color
  Color v0Ks, v1Ks, v2Ks;
  v0Ks = Color(triangles[minTIndex].v[0].color_specular[0], triangles[minTIndex].v[0].color_specular[1], triangles[minTIndex].v[0].color_specular[2]);
  v1Ks = Color(triangles[minTIndex].v[1].color_specular[0], triangles[minTIndex].v[1].color_specular[1], triangles[minTIndex].v[1].color_specular[2]);
  v2Ks = Color(triangles[minTIndex].v[2].color_specular[0], triangles[minTIndex].v[2].color_specular[1], triangles[minTIndex].v[2].color_specular[2]);
  double rv = reflection.dotProduct(ray.negate());
  if (rv < 0) rv = 0.0;
  if (rv > 1) rv = 1.0;
  Color ks = getInterpolatedColor(v0Ks, v1Ks, v2Ks, alpha, beta, gamma);

  double shininesses[3] = { triangles[minTIndex].v[0].shininess, triangles[minTIndex].v[1].shininess, triangles[minTIndex].v[2].shininess };
  double interpolateShininess = alpha * shininesses[0] + beta * shininesses[1] + (gamma) * shininesses[2];

  Color specular = ks.multiplyColor(pow(rv, interpolateShininess));

  // Calculate phong model
  Color light = getLightColor(lightIndex);
  Color color = light.multiplyColor(diffuse.addColor(specular));
  return color;
}

// Calculate ray color using the phong model
Color getRayColor(Vector3D ray, bool isMinTForSphere, double minT, int minTIndex, double alpha, double beta, double gamma) {
  Vector3D shadowRay;
  Color rayColor;

  // If ray intersects with some object in the scene
  if (minT != LEAST_POSSIBLE_T) {
    // Loop over all lights in the scene
    for (int lightIndex = 0; lightIndex < num_lights; lightIndex++) {
      // Shoot a shadow ray to the light source
      shadowRay = getShadowRay(ray.scalarMultiply(minT), lightIndex);

      // If shadow ray is not blocked, calculate the phong model
      bool isBlocked = isShadowRayBlocked(shadowRay, ray.scalarMultiply(minT), lightIndex);
      if (!isBlocked) {
        // If the first intersection is with a sphere, calculate phong model for the sphere
        if (isMinTForSphere) {
          rayColor = rayColor.addColor(getSphereColor(ray, shadowRay, lightIndex, minT, minTIndex));
        }
        // Else calculate phong model for the triangle
        else {
          rayColor = rayColor.addColor(getTriangleColor(ray, shadowRay, lightIndex, minTIndex, alpha, beta, gamma));
        }
      }
    }
  }
  // If ray missed everything in the scene, display a background color
  else {
    Color white = Color(1.0, 1.0, 1.0);
    rayColor = rayColor.addColor(white);
  }
  return rayColor;
}


// Check if the ray intersects the sphere
void checkSphereIntersection(Vector3D ray, double &leastT, int &minTIndex) {
  double xc, yc, zc, radius, root, b, c;
  double t0, t1, minT;

  // Loop over all spheres
  for (int sphereIndex = 0; sphereIndex < num_spheres; sphereIndex++) {
    // Get the sphere and its details
    getSphereDetails(sphereIndex, radius, xc, yc, zc);

    // Calculate roots
    b = -2.0 * (ray.x * xc + ray.y * yc + ray.z * zc);
    c = pow(xc, 2.0) + pow(yc, 2.0) + pow(zc, 2.0) - pow(radius, 2.0);
    root = b * b - 4 * c;

    // If there is an intersection, find the first intersection
    if (root >= 0) {
      t0 = (- b + sqrt(root)) / 2.0;
      t1 = (- b - sqrt(root)) / 2.0;

      // Find the minimun t value
      if (t0 > 0 && t1 > 0) {
        if (t0 < t1)
          minT = t0;
        else
          minT = t1;
      }
      else if (t0 > 0) minT = t0;
      else if (t1 > 0) minT = t1;
      else continue;

      if (minT < leastT) {
        leastT = minT;
        minTIndex = sphereIndex;
      }
    }
  }
}

// Check if the ray intersects the triangle
void checkTriangleIntersection(Vector3D ray, double &leastT, int &minTIndex, bool &isMinTForSphere, double &alpha, double &beta, double &gamma) {
  double minT;
  double alphaTemp, betaTemp, gammaTemp;
  double epsilon = 1e-8;

  // Loop over all triangles
  for (int triangleIndex = 0; triangleIndex < num_triangles; triangleIndex++) {
    Vector3D n = getTriangleNormal(triangleIndex).normalize();
    double nDotRay = n.dotProduct(ray);

    if (nDotRay < epsilon || epsilon < nDotRay) {
      Vector3D v0 = getTriangleVertex(triangleIndex, 0);
      minT = n.dotProduct(v0) / nDotRay;

      // Check for triangle intersections
      if (minT > 0 && minT < leastT) {
        if (isTriangleIntersected(triangleIndex, ray, minT, Vector3D(0, 0, 0), alpha, beta, gamma)) {
          leastT = minT;
          minTIndex = triangleIndex;
          isMinTForSphere = false;
          alphaTemp = alpha;
          betaTemp = beta;
          gammaTemp = gamma;
        }
      }
    }
  }
  alpha = alphaTemp;
  beta = betaTemp;
  gamma = gammaTemp;
}

// Generate random additional lights to get soft shadows
void simulateSoftShadows() {
  int numOfLights = num_lights;

  // Loop over all light sources
  for (int lightIndex = 0; lightIndex < num_lights; lightIndex++) {
    for (int i = 0; i < EXTRA_LIGHTS; i++) {
      // Get a random position
      double delta = 0.1;
      static double numbers[] = { -0.8, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.8} ;
      int index = rand() % (sizeof numbers / sizeof *numbers);
      double random = numbers[index];
      double randomPosition = random * delta;

      // Randomly distribute additional lights
      lights[numOfLights].position[0] = lights[lightIndex].position[0] + randomPosition;
      lights[numOfLights].position[1] = lights[lightIndex].position[1] + randomPosition;
      lights[numOfLights].position[2] = lights[lightIndex].position[2] + randomPosition;

      // Distribute light intensity
      lights[numOfLights].color[0] = lights[lightIndex].color[0] / EXTRA_LIGHTS;
      lights[numOfLights].color[1] = lights[lightIndex].color[1] / EXTRA_LIGHTS;
      lights[numOfLights].color[2] = lights[lightIndex].color[2] / EXTRA_LIGHTS;

      numOfLights++;
    }

    lights[lightIndex].color[0] = lights[lightIndex].color[0] / EXTRA_LIGHTS;
    lights[lightIndex].color[1] = lights[lightIndex].color[1] / EXTRA_LIGHTS;
    lights[lightIndex].color[2] = lights[lightIndex].color[2] / EXTRA_LIGHTS;
  }

  num_lights = numOfLights;
}

// Get ray color after checking for any intersections
Color checkRayIntersections(Vector3D ray) {
  double alpha, beta, gamma;
  double minT;
  double leastT = LEAST_POSSIBLE_T;
  int minTIndex;
  bool isMinTForSphere = true;

  // Calculate closest intersection of ray to objects in the scene
  checkSphereIntersection(ray, leastT, minTIndex);
  checkTriangleIntersection(ray, leastT, minTIndex, isMinTForSphere, alpha, beta, gamma);
  minT = leastT;

  // Calculate phong lighting to get final pixel color for the ray
  return getRayColor(ray, isMinTForSphere, minT, minTIndex, alpha, beta, gamma);
}

// Shoot a ray and calculate its color based on its intersections
Color shootRay(int x, int i, int y, int j) {
  double deltaRaySamples = 1 / RAY_SAMPLES;

  // Get image plane dimensions
  double fovBy2Degrees = fov / 2 * PI / 180.0;
  double displacementX = -aspect_ratio * tan(fovBy2Degrees);
  double displacementY = -tan(fovBy2Degrees);
  double imagePlaneSize = 2 * tan(fovBy2Degrees) / ((double)HEIGHT);

  // Calculate ray
  Vector3D ray;
  ray.x = imagePlaneSize * ((double) x + i * deltaRaySamples) + displacementX;
  ray.y = imagePlaneSize * ((double) y + j * deltaRaySamples) + displacementY;
  ray.z = -1.0f;
  ray = ray.normalize();

  return checkRayIntersections(ray);
}

// Draw scene
void draw_scene() {
  if(enableSoftShadows) {
    simulateSoftShadows();
  }

  // Loop over all pixels
  for(unsigned int x = 0; x < WIDTH; x++) {
    glPointSize(2.0);    
    glBegin(GL_POINTS);
    for(unsigned int y = 0; y < HEIGHT; y++) {
      Color pixelColor;
      // Loop over all sample rays
      for (int i = 0; i < RAY_SAMPLES; i++) {
        for (int j = 0; j < RAY_SAMPLES; j++) {
          // Get color for this ray
          Color color = shootRay(x, i, y, j);
          pixelColor = pixelColor.addColor(color);
        }
      }

      // Calculate final pixel color by adding up and clamping the color
      double redChannel = pixelColor.r / pow(RAY_SAMPLES, 2) + ambient_light[0];
      double greenChannel = pixelColor.g / pow(RAY_SAMPLES, 2) + ambient_light[1];
      double blueChannel = pixelColor.b / pow(RAY_SAMPLES, 2) + ambient_light[2];

      Color finalPixelColor = Color(redChannel, greenChannel, blueChannel);
      finalPixelColor.clampColorTo1();
      plot_pixel(x, y, finalPixelColor.r * 255, finalPixelColor.g * 255, finalPixelColor.b * 255);
    }

    glEnd();
    glFlush();
  }

  printf("Done!\n");
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg() {
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found) {
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3]) {
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r) {
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi) {
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv) {
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display() {}

void init() {
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle() {
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv) {
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
