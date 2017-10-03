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

void part0()
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

void part1()
{
  std::cout << "Part 2: Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("part2.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    Vector3d origin(-1,1,1);
    Vector3d x_displacement(2.0/C.cols(),0,0);
    Vector3d y_displacement(0,-2.0/C.rows(),0);

    // Construct sphere & give non-origin center
    const double sphere_x = -.5;
    const double sphere_y = .3;
    const double sphere_z = -2;

    Vector3d sphere_center = RowVector3d(sphere_x,sphere_y,sphere_z);
    const double sphere_radius = 0.9;

    // Single light source
    const Vector3d light_position(-1,1,1);

    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            Vector3d ray_direction = RowVector3d(0,0,-1);

            // Find discriminant to determine if there's a solution
            double A_ = ray_direction.dot(ray_direction);
            double B_ = 2 * ray_direction.dot((ray_origin - sphere_center));
            double C_ = ((ray_origin - sphere_center).dot((ray_origin - sphere_center))) - std::pow(sphere_radius,2);
            double discriminant =  (std::pow(B_,2) - (4 * A_ * C_));

            // The ray hit the sphere, compute the exact intersection point
            if(discriminant >= 0){
              // std::cout << "Ray hit sphere\n";
              double T_pos = (-(B_) + std::sqrt(discriminant))/(2 * A_);
              double T_neg = (-(B_) - std::sqrt(discriminant))/(2 * A_);

              // Using this line to initialize ray_intersection
              Vector3d ray_intersection;

              // Plug smaller but positive t into  p(t) = e + td and find the point of intersection
              // because the smaller t is 'in front'
              //FIXME:
              if(T_pos < T_neg){
                 ray_intersection = ray_origin + (T_pos * ray_direction);
              }else{
                 ray_intersection = ray_origin + (T_neg * ray_direction);
              }

              // Compute normal at the intersection point
              Vector3d ray_normal = (ray_intersection - sphere_center)/sphere_radius;
              // Vector3d ray_normal = 2 * (ray_intersection - sphere_center);
              // Vector3d ray_normal = (ray_intersection - sphere_center)/sphere_radius;

              // Simple diffuse model
              C(i,j) =  ray_normal.dot((light_position-ray_intersection).normalized().transpose());

              // Clamp to zero
              C(i,j) = max(C(i,j),0.);


              // Disable the alpha mask for this pixel
              A(i,j) = 1;

            }
        }
    }
    // Save to png
    std::cout << "Saving to png";
    write_matrix_to_png(C,C,C,A,filename);

}

void get_intersection_point()
{

}

int main()
{
    part0();
    part1();

    return 0;
}
