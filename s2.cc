//
// <YOUR NAME>
// Computational Vision Homework 4
// Program 2
//

#include "image.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;
using namespace ComputerVisionProjects;


// @brief Implementation of Program 2 to compute light source directions
// @param input_params filename the name of the input file containing X0 Y0 radius
// @param input_sphere_filenames the names of the the three sphere files
// @param output_directions_filename the name of the output directions filename
void ComputeAndSaveDirections(const string &input_params_filename, const array<string,3> &input_sphere_filenames,
                              const string &output_directions_filename) {
  cout << "Compute And Save Directions" << endl;
  cout << "Input params filename: " << input_params_filename << endl;
  cout << "Sphere image filename 1: " << input_sphere_filenames[0] << endl;
  cout << "Sphere image filename 2: " << input_sphere_filenames[1] << endl;
  cout << "Sphere image filename 3: " << input_sphere_filenames[2] << endl;
  cout << "Output directions filename: " << output_directions_filename << endl;
}

int main(int argc, char **argv){  
  if (argc != 6) {
    printf("Usage: %s  {input parameters file} {image 1} {image 2} {image 3} {output directions file}\n", argv[0]);
    return 0;
  }

  const string input_params_filename(argv[1]);
  const array<string, 3> input_sphere_filenames{string{argv[2]}, string{argv[3]}, string{argv[4]}};
  const string output_directions_filename(argv[5]);

  ComputeAndSaveDirections(input_params_filename, input_sphere_filenames, output_directions_filename);

  return 0;
}
