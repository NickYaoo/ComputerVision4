//
// <YOUR NAME>
// Computational Vision Homework 4
// Program 1
//

#include "image.h"
#include "DisjSets.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;
using namespace ComputerVisionProjects;


 // @brief Implementation of Task Program 1 to compute sphre params
 //  
 // @param input_sphere_filename the name of the input sphere image
 // @param threshold threshold for binarization
 // @param output_params_filename the name of file to store X0 Y0 radius
void ComputeAndSaveSphereParametrs(const string &input_sphere_filename, int threshold, const string &output_params_filename) {
  cout << "Computing Sphere Parameters" << endl;
  cout << "Input filename: " << input_sphere_filename << endl;
  cout << "Threshold: " << threshold << endl;
  cout << "Output filename: " << output_params_filename << endl;
}

int main(int argc, char **argv){  
  if (argc != 4) {
    printf("Usage: %s {input original image} {input threshold value} {output parameters file}\n", argv[0]);
    return 0;
  }
  
  const string input_filename(argv[1]);
  const int threshold = atoi(argv[2]);
  const string output_params_filename(argv[3]);

  ComputeAndSaveSphereParametrs(input_filename, threshold, output_params_filename);

  return 0;
}
