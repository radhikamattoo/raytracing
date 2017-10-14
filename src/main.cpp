// Radhika Mattoo, rm3485@nyu.edu
// Computer Graphics Fall 2017
// Assignment 1

// C++ include
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;

// Returns discriminant for given data
double get_discriminant(Vector3d &ray_origin, Vector3d &ray_direction, Vector3d &sphere_center, double sphere_radius)
{
  // Find discriminant to determine if there's a solution
  double A_ = ray_direction.dot(ray_direction);
  double B_ = 2 * ray_direction.dot((ray_origin - sphere_center));
  double C_ = ((ray_origin - sphere_center).dot((ray_origin - sphere_center))) - std::pow(sphere_radius,2);
  double discriminant =  (std::pow(B_,2) - (4 * A_ * C_));

  return discriminant;
}
// Gets a ray's intersection point with sphere
Vector3d get_intersection(double discriminant, Vector3d &ray_direction, Vector3d &ray_origin, Vector3d &sphere_center)
{
  double A_ = ray_direction.dot(ray_direction);
  double B_ = 2 * ray_direction.dot((ray_origin - sphere_center));

  double T_pos = (-(B_) + std::sqrt(discriminant))/(2 * A_);
  double T_neg = (-(B_) - std::sqrt(discriminant))/(2 * A_);

  // Plug smaller but positive t into  p(t) = e + td and find the point of intersection
  // Smaller t is 'in front'
  Vector3d ray_intersection;

  if(T_pos < T_neg){
     ray_intersection = ray_origin + (T_pos * ray_direction);
  }else{
     ray_intersection = ray_origin + (T_neg * ray_direction);
  }
  return ray_intersection;
}

// Determines the color for a given pixel using its ray's intersection point
double get_pixel_color(bool diffuse, double discriminant, double sphere_radius, Vector3d &origin, Vector3d light_position, Vector3d &ray_direction, Vector3d &ray_origin, Vector3d &sphere_center)
{
  double pixel_value;

  Vector3d ray_intersection = get_intersection(discriminant, ray_direction, ray_origin, sphere_center);

  // Compute normal at the intersection point
  Vector3d ray_normal = (ray_intersection - sphere_center)/sphere_radius;

  // Get L (pointing towards light source)
  Vector3d L = (light_position - ray_intersection).normalized().transpose();

  if(diffuse){
    // Simple diffuse model
    pixel_value =  ray_normal.dot(L);

    // Clamp to zero
    pixel_value = max(pixel_value,0.);

  }else{
    // Ambient lighting is b/w 0 and 1
    double ambient = 0.01;

    // Phong exponent
    int phong = 10;

    // Get V (pointing towards the camera)
    Vector3d V = (origin - ray_intersection).normalized().transpose();

    // Get H (bisector of V and L)
    Vector3d H = (V + L).normalized().transpose();

    double diffuse = ray_normal.dot(L);
    diffuse = max(diffuse, 0.);

    // Get specular value
    double specular = ray_normal.dot(H);
    specular = max(specular, 0.);
    specular = std::pow(specular, phong);

    // Sum them
    pixel_value = ambient + diffuse + specular;
  }
  return pixel_value;
}

vector<float> split_line(string line, bool F)
{
  string extracted;
  vector<float> data;
  int z = 0;

  // # of faces is the first character in every line, start at 2 to skip
  if(F){
    z = 2;
  }
  for(int i = z; i <= line.length(); i++){

    char val = line[i];

    if(val == ' ' || i == line.length()){ // Finished building int
      // Convert to int and push to data vector
      data.push_back(atof(extracted.c_str()));
      extracted = "";
    }else{ // Still building int
      extracted.push_back(val);
    }
  }
  return data;
}

// Iterates through given OFF file and fills V & F matrices
pair<MatrixXd, MatrixXd> read_off_data(string filename, bool enlarge)
{
  // Load file
  string line;
  ifstream stream(filename.c_str());
  getline(stream, line); //first line is OFF

  // Get data from 2nd line
  getline(stream, line);
  vector<float> data = split_line(line, false);

  // Extract metadata into vars
  int vertices = data[0];
  int faces = data[1];
  MatrixXd V = MatrixXd::Zero(vertices, 3);
  MatrixXd F = MatrixXd::Zero(faces, 3);
  vector<float> line_data;

  // Fill V & F matrices from file
  for(int v = 0; v < vertices; v++){
    getline(stream, line);
    line_data = split_line(line, false);

    for(int j = 0; j < 3; j++){
      if(enlarge){
        V(v,j) = (line_data[j]*3);
      }else{
        V(v,j) = (line_data[j]/10);
      }
      if(!enlarge){
        V(v,0) += 0.1;
      }else{
        V(v,0) -= 0.15;
        V(v,1) -= 0.1;
      }
    }

  }

  for(int f = 0; f < faces; f++){
    getline(stream, line);
    line_data = split_line(line, true);
    for(int j = 0; j < 3; j++){
      F(f,j) = line_data[j];
    }
  }

  // Construct pair and return
  pair<MatrixXd, MatrixXd> matrices(V, F);
  return matrices;

}
// Solve for x, given Ax = b
vector<float> solver(Vector3d &a_coord, Vector3d &b_coord, Vector3d &c_coord, Vector3d &ray_direction, Vector3d &ray_origin)
{
  // Construct matrices/vectors to solve for
  Matrix3f A_;
  Vector3f b_;

  Vector3d a_minus_b = a_coord - b_coord;
  Vector3d a_minus_c = a_coord - c_coord;
  Vector3d a_minus_e = a_coord - ray_origin;

  A_ << a_minus_b[0], a_minus_c[0], ray_direction[0],   a_minus_b[1], a_minus_c[1], ray_direction[1], a_minus_b[2], a_minus_c[2], ray_direction[2];
  b_ << a_minus_e[0], a_minus_e[1], a_minus_e[2];

  // cout << "Here is the matrix A:\n" << A_ << endl;
  // cout << "Here is the vector b:\n" << b_ << endl;

  Vector3f sol = A_.colPivHouseholderQr().solve(b_);
  if(!(A_*sol).isApprox(b_, 0.003)){
    sol[0] = -1;
    sol[1] = -1;
    sol[2] = -1;
  }
  vector<float> solutions;

  // cout << "Here is the solution: " << sol << endl;
  solutions.push_back(sol[0]);
  solutions.push_back(sol[1]);
  solutions.push_back(sol[2]);
  return solutions;
}

// Translates indices from F matrix into 3D coordinates from V
vector<Vector3d> get_triangle_coordinates(MatrixXd &V, MatrixXd &F, unsigned row)
{
  Vector3f coordinates;
  for(unsigned y=0; y< F.cols();y++)
  {
    // Get F indices from row
    coordinates[y] = F(row,y);
  }
  // Now have indices to index into V with
  float a_component = coordinates[0];
  float b_component = coordinates[1];
  float c_component = coordinates[2];

  // cout << a_component << endl;
  // cout << b_component << endl;
  // cout <<"C index:" << c_component << endl;

  // Get triangle coordinates
  Vector3d a_coord = RowVector3d(V(a_component,0),V(a_component, 1),V(a_component, 2));
  Vector3d b_coord = RowVector3d(V(b_component,0),V(b_component, 1),V(b_component, 2));
  Vector3d c_coord = RowVector3d(V(c_component,0),V(c_component, 1),V(c_component, 2));

  // cout << "a_coord: " <<  a_coord << endl;
  // cout << "b_coord: " << b_coord << endl;
  // cout << "c_coord: " << c_coord << endl;

  vector<Vector3d> ret;
  ret.push_back(a_coord);
  ret.push_back(b_coord);
  ret.push_back(c_coord);
  return ret;
}
// Gets pixel value for intersection point
float compute_triangle_pixel(Vector3d &ray_intersection, Vector3d &ray_normal, Vector3d &light_position, Vector3d &origin)
{
  // Get L (pointing towards light source)
  Vector3d L = (light_position - ray_intersection).normalized().transpose();

  // Ambient lighting is b/w 0 and 1
  double ambient = 0.01;

  // Phong exponent
  int phong = 10;

  // Get V (pointing towards the camera)
  Vector3d V = (origin - ray_intersection).normalized().transpose();

  // Get H (bisector of V and L)
  Vector3d H = (V + L).normalized().transpose();

  double diffuse = ((light_position - ray_intersection).normalized().transpose()) * ray_intersection;
  diffuse = max(diffuse, 0.);

  // Get specular value
  double specular = ray_normal.dot(H);
  specular = max(specular, 0.);
  specular = std::pow(specular, phong);

  // Sum them & return
  return ambient + diffuse + specular;
}

// 1.4 Ray Tracing Triangle Meshes
bool does_intersect(float t, float u, float v)
{
  if(t > 0 && u >= 0 && v >=0 && (u+v) <= 1){ return true; }
  return false;
}

bool is_a_shadow(Vector3d &ray_intersection, Vector3d &light_position, MatrixXd &F_bumpy, MatrixXd &V_bumpy, MatrixXd &F_bunny, MatrixXd &V_bunny)
{
  float epsilon = 0.01;
  Vector3d shadow_origin = ray_intersection;
  Vector3d shadow_direction = (light_position - ray_intersection).normalized().transpose();
  shadow_origin = shadow_origin + (epsilon * shadow_direction);

  bool intersects = false;

  for(unsigned x=0; x<(F_bumpy.rows()); x++)
  {
    // Get 3D coordinates
    vector<Vector3d> coords_bumpy = get_triangle_coordinates(V_bumpy, F_bumpy, x);
    Vector3d a_coord_bumpy = coords_bumpy[0];
    Vector3d b_coord_bumpy = coords_bumpy[1];
    Vector3d c_coord_bumpy = coords_bumpy[2];

    vector<Vector3d> coords_bunny = get_triangle_coordinates(V_bunny, F_bunny, x);
    Vector3d a_coord_bunny = coords_bunny[0];
    Vector3d b_coord_bunny = coords_bunny[1];
    Vector3d c_coord_bunny = coords_bunny[2];

    // Get u, v, and t
    vector<float> solutions_bumpy = solver(a_coord_bumpy, b_coord_bumpy, c_coord_bumpy, shadow_direction, shadow_origin);
    float u_bumpy =  solutions_bumpy[0];
    float v_bumpy =  solutions_bumpy[1];
    float t_bumpy =  solutions_bumpy[2];

    vector<float> solutions_bunny = solver(a_coord_bunny, b_coord_bunny, c_coord_bunny, shadow_direction, shadow_origin);
    float u_bunny =  solutions_bunny[0];
    float v_bunny =  solutions_bunny[1];
    float t_bunny =  solutions_bunny[2];

    // Check for intersection
    bool bunny_intersects = does_intersect(t_bunny, u_bunny,v_bunny);
    bool bumpy_intersects = does_intersect(t_bumpy, u_bumpy, v_bumpy);
    if(bunny_intersects || bumpy_intersects){
      intersects = true;
      break;
    }
  }
  return intersects;
}

// 1.1 Ray Tracing Spheres
void part1()
{
  std::cout << "Part 1.1: Ray Tracing Spheres" << std::endl;
  MatrixXd C = MatrixXd::Zero(800,800);
  const std::string filename("part1.png");

  MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

  // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
  Vector3d origin(-1,1,1);
  Vector3d x_displacement(2.0/A.cols(),0,0);
  Vector3d y_displacement(0,-2.0/A.rows(),0);

  // CHANGE THIS LINE FOR NON-CENTERED SPHERE
  Vector3d sphere_center = RowVector3d(-0.5,0,0);
  Vector3d sphere_center_2 = RowVector3d(0.5,0,0);
  Vector3d sphere_center_3 = RowVector3d(0.0,0.5,0);
  const double sphere_radius = 0.4;

  const Vector3d light_position(-1,1,1);
  const bool diffuse = true;

  for (unsigned i=0;i<C.cols();i++)
  {
      for (unsigned j=0;j<C.rows();j++)
      {
          // Prepare the ray
          Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
          Vector3d ray_direction = RowVector3d(0,0,-1);

          // Find discriminant and intersection point
          double discriminant =  get_discriminant(ray_origin, ray_direction, sphere_center, sphere_radius);
          double discriminant_2 =  get_discriminant(ray_origin, ray_direction, sphere_center_2, sphere_radius);
          double discriminant_3 =  get_discriminant(ray_origin, ray_direction, sphere_center_3, sphere_radius);

          if(discriminant >= 0){
            C(i,j) = get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center);
          }else if(discriminant_2 >=0){
            C(i,j) = get_pixel_color(diffuse, discriminant_2, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center_2);

          }else if(discriminant_3 >=0){
            C(i,j) = get_pixel_color(diffuse, discriminant_3, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center_3);

          }else{
            C(i,j) = 1.0;
          }
          // Disable the alpha mask for this pixel
          A(i,j) = 1;
      } // inner loop
  } // outer loop
  std::cout << "Saving 1.1 sphere to png" << std::endl;
  write_matrix_to_png(C,C,C,A,filename);
}

// 1.2 Shading
void part2()
{
  std::cout << "Part 1.2: Shading" << std::endl;

    const string filename("part2.png");

    // Store the color
    MatrixXd R = MatrixXd::Zero(800,800);
    MatrixXd G = MatrixXd::Zero(800,800);
    MatrixXd B = MatrixXd::Zero(800,800);

    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    // One sphere should have specular, and one only diffuse lighting
    Vector3d sphere_center = RowVector3d(-0.5,0,0);
    Vector3d sphere_center_2 = RowVector3d(0.5,0,0);
    Vector3d sphere_center_3 = RowVector3d(0.0,0.5,0);
    const double sphere_radius = 0.4;

    const Vector3d light_position(-1,2,1);
    const Vector3d light_position_2(1, -2, -1);

    bool diffuse;
    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

            // Find discriminants to determine if there's a solution
            double discriminant =  get_discriminant(ray_origin, ray_direction, sphere_center, sphere_radius);
            double discriminant_2 = get_discriminant(ray_origin, ray_direction, sphere_center_2, sphere_radius);
            double discriminant_3 = get_discriminant(ray_origin, ray_direction, sphere_center_3, sphere_radius);

            if(discriminant >= 0){
              diffuse = true; // one sphere should have only diffuse
              R(i,j) = get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center);
              R(i,j) += get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center);
              R(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else if(discriminant_2 >= 0){
              diffuse = false; // the other sphere should be specular
              G(i,j) = get_pixel_color(diffuse, discriminant_2, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center_2);
              G(i,j) += get_pixel_color(diffuse, discriminant_2, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center_2);
              G(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else if(discriminant_3 >= 0){
              diffuse = false; // the other sphere should be specular
              B(i,j) = get_pixel_color(diffuse, discriminant_3, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center_3);
              B(i,j) += get_pixel_color(diffuse, discriminant_3, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center_3);
              B(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else{
              R(i,j) = 1.0;
              G(i,j) = 1.0;
              B(i,j) = 1.0;
            }
            // Disable the alpha mask for this pixel
            A(i,j) = 1;

        } // inner loop
    } // outer loop
    std::cout << "Saving sphere 1.2 to png" << std::endl;
    write_matrix_to_png(R,G,B,A,filename);
}
// 1.3 Perspective Projection
void part3()
{
  std::cout << "Part 1.3: Perspective Projection" << std::endl;

    const std::string filename("part3.png");

    // Store the color
    MatrixXd R = MatrixXd::Zero(800,800);
    MatrixXd G = MatrixXd::Zero(800,800);
    MatrixXd B = MatrixXd::Zero(800,800);

    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    // One sphere should have specular, and one only diffuse lighting
    Vector3d sphere_center = RowVector3d(-0.5,0,0);
    Vector3d sphere_center_2 = RowVector3d(0.5,0,0);
    Vector3d sphere_center_3 = RowVector3d(0.0,0.5,0);
    const double sphere_radius = 0.4;

    const Vector3d light_position(-1,1,1);
    const Vector3d light_position_2(1, -2, -1);

    bool diffuse;
    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin;

            // Formula taken from textbook
            double focal_length = 1.0;
            Vector3d w(0,0,-1);
            Vector3d ray_direction =  (focal_length * w) + (double(i)*x_displacement + double(j)*y_displacement);

            // Find discriminants to determine if there's a solution
            double discriminant =  get_discriminant(ray_origin, ray_direction, sphere_center, sphere_radius);
            double discriminant_2 = get_discriminant(ray_origin, ray_direction, sphere_center_2, sphere_radius);
            double discriminant_3 = get_discriminant(ray_origin, ray_direction, sphere_center_3, sphere_radius);

            if(discriminant >= 0){
              diffuse = true; // one sphere should have only diffuse
              R(i,j) = get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center);
              R(i,j) += get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center);
              R(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else if(discriminant_2 >= 0){
              diffuse = false; // the other sphere should be specular
              G(i,j) = get_pixel_color(diffuse, discriminant_2, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center_2);
              G(i,j) += get_pixel_color(diffuse, discriminant_2, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center_2);
              G(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else if(discriminant_3 >= 0){
              diffuse = false; // the other sphere should be specular
              B(i,j) = get_pixel_color(diffuse, discriminant_3, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center_3);
              B(i,j) += get_pixel_color(diffuse, discriminant_3, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center_3);
              B(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else{
              R(i,j) = 1.0;
              G(i,j) = 1.0;
              B(i,j) = 1.0;
            }
            // Disable the alpha mask for this pixel
            A(i,j) = 1;
        } // inner loop
    } // outer loop
    std::cout << "Saving sphere 1.3 to png" << std::endl;
    write_matrix_to_png(R,G,B,A,filename);
}
// Returns a <bool,float> describing which object it intersected with and the color to set
// True means it hit the BUNNY, False means it hit  BUMPY
pair<bool,float> reflection_intersection(Vector3d &r, Vector3d &intersection_point, Vector3d &light_position, MatrixXd &F_bumpy, MatrixXd &V_bumpy, MatrixXd &F_bunny, MatrixXd &V_bunny)
{
  // R is our ray direction, intersection point is our origin
  float epsilon = 0.01;
  Vector3d ray_origin = intersection_point;
  Vector3d ray_direction = r;
  ray_origin = ray_origin + (epsilon * ray_direction);


  Vector3d origin(-1,1,1);
  pair<bool, float> pixel_data(false,-100);
  // float specular_color =
  // Iterate through all faces and determine if intersection occurs
  for(unsigned x=0; x<(F_bumpy.rows()); x++)
  {
    // Get 3D coordinates
    vector<Vector3d> coords_bumpy = get_triangle_coordinates(V_bumpy, F_bumpy, x);
    Vector3d a_coord_bumpy = coords_bumpy[0];
    Vector3d b_coord_bumpy = coords_bumpy[1];
    Vector3d c_coord_bumpy = coords_bumpy[2];


    vector<Vector3d> coords_bunny = get_triangle_coordinates(V_bunny, F_bunny, x);
    Vector3d a_coord_bunny = coords_bunny[0];
    Vector3d b_coord_bunny = coords_bunny[1];
    Vector3d c_coord_bunny = coords_bunny[2];

    // Get u, v, and t
    vector<float> solutions_bumpy = solver(a_coord_bumpy, b_coord_bumpy, c_coord_bumpy, ray_direction, ray_origin);
    float u_bumpy =  solutions_bumpy[0];
    float v_bumpy =  solutions_bumpy[1];
    float t_bumpy =  solutions_bumpy[2];

    vector<float> solutions_bunny = solver(a_coord_bunny, b_coord_bunny, c_coord_bunny, ray_direction, ray_origin);
    float u_bunny =  solutions_bunny[0];
    float v_bunny =  solutions_bunny[1];
    float t_bunny =  solutions_bunny[2];

    // Check for intersection
    bool bunny_intersects = does_intersect(t_bunny, u_bunny,v_bunny);
    bool bumpy_intersects = does_intersect(t_bumpy, u_bumpy, v_bumpy);

    if(bunny_intersects){
      Vector3d ray_intersection = ray_origin + (t_bunny * ray_direction);
      Vector3d ray_normal = ((b_coord_bunny - a_coord_bunny).cross(c_coord_bunny - a_coord_bunny)).normalized();
      float pixel_value = compute_triangle_pixel(ray_intersection, ray_normal, light_position, origin);
      pixel_data = make_pair(true, pixel_value);
      break;
    }else if(bumpy_intersects){
      Vector3d ray_intersection = ray_origin + (t_bumpy * ray_direction);
      Vector3d ray_normal = ((b_coord_bumpy - a_coord_bumpy).cross(c_coord_bumpy - a_coord_bumpy)).normalized();
      float pixel_value = compute_triangle_pixel(ray_intersection, ray_normal, light_position, origin);
      pixel_data = make_pair(false, pixel_value);
      break;
    }
  } //end for loop
  return pixel_data;

}

void part4(bool shadows, bool reflections)
{
    if(shadows){
      cout << "Part 1.5: Shadows" << endl;
    }else if(reflections){
      cout << "Part 1.6: Reflections" << endl;
    }else{
      cout << "Part 1.4: Ray Tracing Triangle Meshes" << endl;
    }
    const std::string filename("part4-perspective-shifted-bunny-no-shadows.png");

    // Create data matrices from OFF files
    pair<MatrixXd, MatrixXd> bumpy = read_off_data("../data/bumpy_cube.off", false);
    MatrixXd V_bumpy = bumpy.first;
    MatrixXd F_bumpy = bumpy.second;

    pair<MatrixXd, MatrixXd> bunny = read_off_data("../data/bunny.off", true);
    MatrixXd V_bunny = bunny.first;
    MatrixXd F_bunny = bunny.second;

    // Store the color
    MatrixXd R = MatrixXd::Zero(800,800);
    MatrixXd G = MatrixXd::Zero(800,800);
    MatrixXd B = MatrixXd::Zero(800,800);

    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/A.cols(),0,0);
    Vector3d y_displacement(0,-2.0/A.rows(),0);

    // Vector3d light_position(-1,1,1);
    Vector3d light_position(1,-1,1);
    // mirror is the plane where y = -0.6
    float mirror = -0.6;
    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {
          // Prepare the ray
          Vector3d ray_origin;
          Vector3d ray_direction;
          if(reflections){
            ray_origin = origin;

            // Formula taken from textbook
            double focal_length = 1.0;
            Vector3d w(0,0,-1);
            ray_direction =  (focal_length * w) + (double(i)*x_displacement + double(j)*y_displacement);
          }else{
            ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            ray_direction = RowVector3d(0,0,-1);
          }

          // Iterate through all faces and determine if intersection occurs
          for(unsigned x=0; x<(F_bumpy.rows()); x++)
          {
            // Get 3D coordinates
            vector<Vector3d> coords_bumpy = get_triangle_coordinates(V_bumpy, F_bumpy, x);
            Vector3d a_coord_bumpy = coords_bumpy[0];
            Vector3d b_coord_bumpy = coords_bumpy[1];
            Vector3d c_coord_bumpy = coords_bumpy[2];


            vector<Vector3d> coords_bunny = get_triangle_coordinates(V_bunny, F_bunny, x);
            Vector3d a_coord_bunny = coords_bunny[0];
            Vector3d b_coord_bunny = coords_bunny[1];
            Vector3d c_coord_bunny = coords_bunny[2];

            // Get u, v, and t
            vector<float> solutions_bumpy = solver(a_coord_bumpy, b_coord_bumpy, c_coord_bumpy, ray_direction, ray_origin);
            float u_bumpy =  solutions_bumpy[0];
            float v_bumpy =  solutions_bumpy[1];
            float t_bumpy =  solutions_bumpy[2];

            vector<float> solutions_bunny = solver(a_coord_bunny, b_coord_bunny, c_coord_bunny, ray_direction, ray_origin);
            float u_bunny =  solutions_bunny[0];
            float v_bunny =  solutions_bunny[1];
            float t_bunny =  solutions_bunny[2];

            // Check for intersection
            bool bunny_intersects = does_intersect(t_bunny, u_bunny,v_bunny);
            bool bumpy_intersects = does_intersect(t_bumpy, u_bumpy, v_bumpy);

            if(bunny_intersects){
              Vector3d ray_intersection = ray_origin + (t_bunny * ray_direction);
              Vector3d ray_normal = ((b_coord_bunny - a_coord_bunny).cross(c_coord_bunny - a_coord_bunny)).normalized();
              if(shadows){ // Do we want to check for shadows?
                if( is_a_shadow(ray_intersection,light_position, F_bumpy, V_bumpy, F_bunny, V_bunny)){
                  R(i,j) = 0.0;
                }else{
                  R(i,j) = compute_triangle_pixel(ray_intersection, ray_normal, light_position, origin);
                }
              }else{
                R(i,j) = compute_triangle_pixel(ray_intersection, ray_normal, light_position, origin);
              }
              A(i,j) = 1.0;
              break;
            }else if(bumpy_intersects){
              Vector3d ray_intersection = ray_origin + (t_bumpy * ray_direction);
              Vector3d ray_normal = ((b_coord_bumpy - a_coord_bumpy).cross(c_coord_bumpy - a_coord_bumpy)).normalized();
              if(shadows){ // Do we want to check for shadows?
                if( is_a_shadow(ray_intersection,light_position, F_bumpy, V_bumpy, F_bunny, V_bunny) ){
                  G(i,j) = 0.0;
                }
                else{
                  G(i,j) = compute_triangle_pixel(ray_intersection, ray_normal, light_position, origin);
                }
              }else{
                G(i,j) = compute_triangle_pixel(ray_intersection, ray_normal, light_position, origin);
              }
              A(i,j) = 1.0;
              break;
            }
          } // F matrix for loop
          if(reflections){
            float t_reflection = (mirror - ray_origin[1])/ray_direction[1];
            if(t_reflection > 0){
              // If there's a reflection, find r and perform intersection again
              Vector3d ray_intersection_reflection = ray_origin + (t_reflection * ray_direction);
              // Vector3d ray_normal = ray_intersection_reflection.normalized().transpose();
              Vector3d ray_normal = RowVector3d(0,1,0); // Normal from the y-plane is a unit vector pointing directly upwards
              Vector3d d_direction = (ray_intersection_reflection - origin).normalized().transpose();
              Vector3d r_direction = d_direction - (2 *(d_direction.dot(ray_normal)) * ray_normal);

              // Perform ray intersection on r
              pair<bool,float> result = reflection_intersection(r_direction, ray_intersection_reflection, light_position, F_bumpy, V_bumpy, F_bunny, V_bunny);
              bool is_bunny = result.first;
              float pixel_value = result.second;
              if(pixel_value != -100){
                if(is_bunny){
                  R(i,j) = pixel_value;
                }else{
                  G(i,j) = pixel_value;
                }
                A(i,j) = 1.0;
              } // pixel_value if
            } // t > 0 if
          } // reflections if

        } // inner for loop
        if(i % 200 == 0){
          cout << "At outer index: " << i << endl;
        }
      } // outer for loop
      std::cout << "Saving sphere 1.4/1.5 to png" << std::endl;
      write_matrix_to_png(R,G,B,A,filename);
}


int main()
{
    // part1();             // 1.1
    // part2();             // 1.2
    // part3();             // 1.3
    // part4(false, false); // 1.4
    // part4(true, false); // 1.5
    part4(false, true); // 1.6


    return 0;
}
