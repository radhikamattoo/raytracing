// Radhika Mattoo, rm3485@nyu.edu
// Computer Graphics Fall 2017
// Assignment 1

// C++ include
#include <iostream>
#include <string>
#include <vector>

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
  // because the smaller t is 'in front'
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

//
void perspective_projection()
{

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
  Vector3d sphere_center = RowVector3d(0,0,0);
  const double sphere_radius = 0.9;

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
          if(discriminant >= 0){
            C(i,j) = get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center);
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

    const std::string filename("part2.png");

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
    const double sphere_radius = 0.4;

    const Vector3d light_position(-1,2,1);
    const Vector3d light_position_2(1, -2, -1);

    bool diffuse;
    for (unsigned i=0;i<A.cols();i++)
    {
        for (unsigned j=0;j<A.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin;

            // Formula taken from textbook
            Vector3d dir(0,0,-1);
            Vector3d ray_direction =  dir + (double(i)*x_displacement + double(j)*y_displacement);

            // Find discriminants to determine if there's a solution
            double discriminant =  get_discriminant(ray_origin, ray_direction, sphere_center, sphere_radius);
            double discriminant_2 = get_discriminant(ray_origin, ray_direction, sphere_center_2, sphere_radius);

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

int main()
{
    // part1();
    // part2();
    part3();

    return 0;
}
