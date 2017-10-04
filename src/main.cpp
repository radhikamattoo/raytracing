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

void part1()
{
    // std::cout << "Part 1: Writing a grid png image" << std::endl;

    const std::string filename("part1.png");
    Eigen::MatrixXd M(800,800);

    // Draw a grid, each square has a side of e pixels
    const int e = 50;
    const double black = 0;
    const double white = 1;

    for (unsigned wi = 0; wi<M.cols();++wi)
        for (unsigned hi = 0; hi < M.rows(); ++hi)
            M(hi,wi) = (lround(wi / e) % 2) == (lround(hi / e) % 2) ? black : white;

    // Write it in a png image. Note that the alpha channel is reversed to make the white (color = 1) pixels transparent (alhpa = 0)
    write_matrix_to_png(M,M,M,1.0-M.array(),filename);
}


void part2()
{
  std::cout << "Part 2: Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("part2.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/C.cols(),0,0);
    Vector3d y_displacement(0,-2.0/C.rows(),0);

    // TODO: Render multiple spheres with different colors
    // One should have specular  and one should have diffuse lighting
    Vector3d sphere_center = RowVector3d(-0.5,0,0);
    Vector3d sphere_center_2 = RowVector3d(0.5,0,0);
    const double sphere_radius = 0.4;

    // TODO: Add a second light source
    const Vector3d light_position(-1,2,1);
    const Vector3d light_position_2(1, -2, -1);

    bool diffuse;
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

            // Find discriminants to determine if there's a solution
            double discriminant =  get_discriminant(ray_origin, ray_direction, sphere_center, sphere_radius);
            double discriminant_2 = get_discriminant(ray_origin, ray_direction, sphere_center_2, sphere_radius);

            if(discriminant >= 0){
              diffuse = true; // one sphere should have only diffuse
              C(i,j) = get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center);
              C(i,j) += get_pixel_color(diffuse, discriminant, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center);
              C(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else if(discriminant_2 >= 0){
              diffuse = false; // the other sphere should be specular
              C(i,j) = get_pixel_color(diffuse, discriminant_2, sphere_radius, origin, light_position, ray_direction, ray_origin, sphere_center_2);
              C(i,j) += get_pixel_color(diffuse, discriminant_2, sphere_radius, origin, light_position_2, ray_direction, ray_origin, sphere_center_2);
              C(i,j) -= 0.1; //only supposed to count ambient lighting once
            }else{
              C(i,j) = 1.0;
            }

            // Disable the alpha mask for this pixel
            A(i,j) = 1;

        } // inner loop
    } // outer loop

    // Save to png
    std::cout << "Saving to png";
    // if(diffuse){
      // write_matrix_to_png(C,C,C,A,filename);
    // }else{
      // add color!
      write_matrix_to_png(C,C,C,A,filename);
    // }
}


int main()
{
    part1();
    part2();

    return 0;
}
