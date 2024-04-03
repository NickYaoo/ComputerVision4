//
// <YOUR NAME>
// Computational Vision Homework 3
// Program 3
//

#include "image.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;
using namespace ComputerVisionProjects;

// @brief Computes Normals and Albedos
// @param input_directions_filename the name of the input directions filename (as computed from s2)
// @param input_object_filenames the names of the the three object files
// @param step the step parameter (for skipping in normal visualization)
// @param threshold the threshold parameter (for determining bright pixels)
// @param output_normals_filename the name of output normals image
// @param output_albedos_filename the name of output albedos image
void ComputeAndSaveNormalsAndAlbedoImages(const string &input_directions_filename, const array<string,3> &input_object_filenames, int step, int threshold,
                                          const string &output_normals_filename, const string &output_albedos_filename) {
  cout << "Compute Normals and Albedos using Photometric Stereo" << endl;
  cout << "Input directions filename: " << input_directions_filename << endl;
  cout << "Object image filename 1: " << input_object_filenames[0] << endl;
  cout << "Object image filename 2: " << input_object_filenames[1] << endl;
  cout << "Object image filename 3: " << input_object_filenames[2] << endl;
  cout << "Step parameter: " << step << endl;
  cout << "Threshold parameter: " << threshold << endl;
  cout << "Output normals filename: " << output_normals_filename << endl;
  cout << "Output albedos filename: " << output_albedos_filename << endl;
  // .. Code that calls other functions
}

int main(int argc, char **argv){
  if (argc != 9) {
    printf("Usage: %s  {input directions filename} {object 1 filename} {object 2 filename} {object 3 filename} {step} {threshold} {output normals filename} {output albedos filename}\n", argv[0]);
    return 0;
  }
  
  const string input_directions_filename{argv[1]};
  const array<string, 3> input_object_filenames{string{argv[2]}, string{argv[3]}, string{argv[4]}};
  const int step = atoi(argv[5]);
  const int threshold = atoi(argv[6]);
  const string output_normals_filename{argv[7]};
  const string output_albedos_filename{argv[8]};
  
  ComputeAndSaveNormalsAndAlbedoImages(input_directions_filename, input_object_filenames, step, threshold, output_normals_filename, output_albedos_filename);
  
  return 0;
}
