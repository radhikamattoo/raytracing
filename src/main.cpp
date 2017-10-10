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
  int i = 0;

  // # of faces is the first character in every line, start at 2 to skip
  if(F){
    i = 2;
  }
  for(i; i <= line.length(); i++){
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
pair<MatrixXd, MatrixXd> read_off_data(string filename)
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
      V(v,j) = line_data[j];
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
vector<float> solver(Vector3d &a_coord, Vector3d &b_coord, Vector3d &c_coord, Vector3d &ray_direction, Vector3d &ray_origin)
{
  // Construct matrices/vectors to solve for
  Matrix3f A_;
  Vector3f b_;
  A_ << (a_coord - b_coord)[0], (a_coord - c_coord)[0], ray_direction[0],   (a_coord - b_coord)[1], (a_coord - c_coord)[1], ray_direction[1],  (a_coord - b_coord)[2], (a_coord - c_coord)[2], ray_direction[2];
  b_ << (a_coord - ray_origin)[0], (a_coord - ray_origin)[1], (a_coord - ray_origin)[2];

  Vector3f sol = A_.colPivHouseholderQr().solve(b_);
  vector<float> solutions;

  solutions.push_back(sol[0]);
  solutions.push_back(sol[1]);
  solutions.push_back(sol[2]);
  return solutions;
}

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

  // Get triangle coordinates
  Vector3d a_coord = RowVector3d(V(a_component,0),V(a_component, 1),V(a_component, 2));
  Vector3d b_coord = RowVector3d(V(b_component,0),V(b_component, 1),V(b_component, 2));
  Vector3d c_coord = RowVector3d(V(c_component,0),V(c_component, 1),V(c_component, 2));

  vector<Vector3d> ret;
  ret.push_back(a_coord);
  ret.push_back(b_coord);
  ret.push_back(c_coord);
  return ret;
}
float get_triangle_pixel(Vector3d &ray_intersection, Vector3d &ray_normal, Vector3d &light_position, Vector3d &origin)
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

  double diffuse = ray_normal.dot(L);
  diffuse = max(diffuse, 0.);

  // Get specular value
  double specular = ray_normal.dot(H);
  specular = max(specular, 0.);
  specular = std::pow(specular, phong);

  // Sum them & return
  return ambient + diffuse + specular;

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


// 1.4 Ray Tracing Triangle Meshes
void part4()
{
    cout << "Part 1.4: Ray Tracing Triangle Meshes" << endl;
    const std::string filename("part4.png");

    // Create data matrices from OFF files
    pair<MatrixXd, MatrixXd> bumpy = read_off_data("data/bumpy_cube.off");
    MatrixXd V_bumpy = bumpy.first;
    MatrixXd F_bumpy = bumpy.second;

    pair<MatrixXd, MatrixXd> bunny = read_off_data("data/bunny.off");
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

    Vector3d light_position(-1,1,1);
    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {
          // Prepare the ray (orthographic)
          // TODO: Make perspective?
          Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
          Vector3d ray_direction = RowVector3d(0,0,-1);

          // Iterate through all faces and determine if intersection occurs
          // TODO: Do for both bumpy and bunny matrices
          for(unsigned x=0; x<F_bumpy.rows(); x++)
          {
            // Get coordinates for given face
            vector<Vector3d> coords = get_triangle_coordinates(V_bumpy, F_bumpy, x);
            Vector3d a_coord = coords[0];
            Vector3d b_coord = coords[1];
            Vector3d c_coord = coords[2];

            // Get u, t, and v
            vector<float> solutions = solver(a_coord, b_coord, c_coord, ray_direction, ray_origin);
            float u =  solutions[0];
            float t =  solutions[1];
            float v =  solutions[2];

            // Check for intersection
            if(t > 0 && u >= 0 && v >=0 && (u+v) <= 1)
            {
              Vector3d ray_intersection = ray_origin + (t * ray_direction);
              Vector3d ray_normal = (b_coord - a_coord).cross(c_coord - a_coord)/(((b_coord - a_coord).cross(c_coord - a_coord)).norm());
              R(i,j) = get_triangle_pixel(ray_intersection, ray_normal, light_position, origin);
            }else{
              R(i,j) = 1.0;
              G(i,j) = 1.0;
              B(i,j) = 1.0;
            } // else
          } // outer for loop
          A(i,j) = 1.0;
        } // inner for loop
        if(i % 200 == 0){
          cout << "At outer index: " << i << endl;
        }
      } // outer for loop
      std::cout << "Saving sphere 1.4 to png" << std::endl;
      write_matrix_to_png(R,G,B,A,filename);

}
// 1.5 Shadows
void part5()
{

}

// 1.6 Reflections on the Floor
void part6()
{

}
int main()
{
    // part1();
    // part2();
    // part3();
    part4();

    return 0;
}
