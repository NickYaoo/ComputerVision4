//
// Nicklaus Yao
// Computational Vision Homework 4
// Program 2
//

#include "image.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

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

  // Read in all 3 spheres
  // Read in sphere_1
  Image sphere_1;
  if (!ReadImage(input_sphere_filenames[0], &sphere_1)) {
    cout <<"Can't open file " << &input_sphere_filenames[0] << endl;
  }
    // Read in sphere_2
  Image sphere_2;
  if (!ReadImage(input_sphere_filenames[1], &sphere_2)) {
    cout <<"Can't open file " << &input_sphere_filenames[1] << endl;
  }
    // Read in sphere_3
  Image sphere_3;
  if (!ReadImage(input_sphere_filenames[2], &sphere_3)) {
    cout <<"Can't open file " << &input_sphere_filenames[2] << endl;
  }

  // Read in values from S1 - yields center_x, center_y, radius
  std::ifstream file(input_params_filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << input_params_filename << std::endl;
    }

  vector<vector<double>> s1_params;
  // typical sStream structure for reading in lines and each value
  std::string line;
  while (std::getline(file, line)) {
      std::vector<double> row;
      std::istringstream iss(line);
      int value;
      while (iss >> value) {
          row.push_back(value);
      }
      s1_params.push_back(row);
    }

  // Testing for Read In
  // cout << s1_params[0][0] << " " << s1_params[0][1] << " " << s1_params[0][2] << endl;

  // Create datastructure {sphere_1: x, y, z, b1}
  //                      {sphere_2: x, y, z, b2}
  //                      {sphere_2: x, y, z, b2}
  vector<vector<double>> normals(3, vector<double>(4,0));

  // Compute the Brightness Point of Each Sphere - stored in b1,b2,b3
  // Brute Force, iterate through each point of each Sphere - NOTE: Dimensions are identical
  for (int i = 0; i < sphere_1.num_rows(); i++){
    for(int j = 0; j < sphere_1.num_columns(); j++){
      // Sphere_1
      if (normals[0][3] < sphere_1.GetPixel(i,j)){
        // set brightness level to new and update our x and y values
        normals[0][3] = sphere_1.GetPixel(i,j);
        normals[0][0] = i;
        normals[0][1] = j;
      }
      // Sphere_2
      if (normals[1][3] < sphere_2.GetPixel(i,j)){
        // set brightness level to new and update our x and y values
        normals[1][3] = sphere_2.GetPixel(i,j);
        normals[1][0] = i;
        normals[1][1] = j;
      }
      // Sphere_3
      if (normals[2][3] < sphere_3.GetPixel(i,j)){
        // set brightness level to new and update our x and y values
        normals[2][3] = sphere_3.GetPixel(i,j);
        normals[2][0] = i;
        normals[2][1] = j;
      }
    }
  }

  // Testing for Brightest Pixel
  /**
  for(int i = 0; i < 3; i++){
    cout << "Sphere " << i << ": ";
    for(int j = 0; j < 4; j++){
      cout << normals[i][j] << " ";
    }
    cout << endl;
  }
  */

  // Calculate Non-unit normal for Each Sphere
  for(int i = 0; i < 3; i++){
    //cout << "Sphere " << i << ": ";
    // Compute (X-X0)
    normals[i][0] = normals[i][0] - s1_params[0][0];
    // Compute (Y-Y0)
    normals[i][1] = normals[i][1] - s1_params[0][1];
    // Compute for (Z-Z0)
    normals[i][2] = sqrt(pow(s1_params[0][2],2) - pow(normals[i][0],2) - pow(normals[i][1],2));
        //for(int j = 0; j < 4; j++){
      //cout << normals[i][j] << " ";
    //}
    //cout << endl;
  }


  // Compute the Unit Normal for Each Sphere
  for(int i = 0; i < 3; i++){
      // Computing Unit Normal
      double scalar = sqrt(pow(normals[i][0],2) + pow(normals[i][1],2) + pow(normals[i][2],2));
      normals[i][0] = normals[i][0] / scalar;
      normals[i][1] = normals[i][1] / scalar;
      normals[i][2] = normals[i][2] / scalar;

      // Computing Normal Scaled by Brightness
      for(int j = 0; j < 3; j++){
        normals[i][j] = normals[i][j] * normals[i][3];
      }
  }

  // Write to a textfile
  std::ofstream output_info(output_directions_filename);
  for(int sphere = 0; sphere < 3; sphere++){
    output_info << normals[sphere][0] <<  " " << normals[sphere][1] << " " << normals[sphere][2] << endl;
  }
  output_info.close();
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
