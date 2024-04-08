//
// Nicklaus Yao
// Computational Vision Homework 3
// Program 3
//

#include "image.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>

using namespace std;
using namespace ComputerVisionProjects;

// 2 by 2 determinant
double det2by2(vector<vector<double>>& matrix){
  // fix -0 edge case
  double det = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
  if (det == -0){
    return 0;
  }
  return det;
}
// Assume 3 by 3 Matrix
vector<vector<double>> InverseMatrix(const vector<vector<double>>& matrix){
  // Use First Row across * noninclusive 2 by 2 determinant method
  double det3by3 = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
                  matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
                  matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
  // Check Determinant Calculations
  //cout << "Determinant of 3 by 3 Matrix: " << det3by3 <<endl;

  // Out of Scope Check: We cannot invert a matrix of determinant 0
  if(det3by3 == 0){
    cout << "Cannot invert matrix with determinant 0" << endl;
    return {};
  }
  // Matrix of Minors + Cofactor Matrix
  vector<vector<double>> minor(3,vector<double>(3));
  // matrix to keep track of 2 by 2s
  vector<vector<double>> sub(2, vector<double>(2));
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      // keeps track of which row we populate for sub matrix
      int sub_row = 0;
      for(int row = 0; row < 3; row++){
        // Similar to above, we want noninclusive rows/ cols determinant
        if(row != i){
          // keeps track of which col we want to populate for sub matrix
          int sub_col = 0;
          for(int col = 0; col < 3; col++){
            if(col != j){
              // Populate the 2 by 2 matrix
              sub[sub_row][sub_col] = matrix[row][col];
              sub_col += 1;
            }
          }
          sub_row += 1;
        }
      }
      // Checking sub matrix:
      /**
      for(int r = 0; r < 2; r++){
        for(int c = 0; c<2; c++){
          cout << sub[r][c] << " ";
        }
        cout << endl;
      }
      cout << endl;
      */
      // Calculate determinant of 2 by 2 and then place into minor matrix
      if((i+j) % 2 == 0){
        minor[i][j] = det2by2(sub);
      }
      else{
        minor[i][j] = -1 * det2by2(sub);
      }
    }
  }
  // Testing for Minor matrix
  /**
  cout << "Testing for Cofactor Matrix:" << endl;
  for(int i =0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      cout << minor[i][j] << " ";
    }
    cout << endl;
  }
  */

  // Transpose Matrix ==> Adjont Matrix
  // cout << "Inverse:" << endl;
  vector<vector<double>> inverse(3,vector<double>(3));
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      // Transpose
      inverse[i][j] = minor[j][i];
      // Compute Inverse matrix by dividing Transposed Matrix by the 3 by 3 determinant
      // Fixing -0 Edge Case:
      if(inverse[i][j] / det3by3 == -0){
        inverse[i][j] = 0;
      }
      else{
        inverse[i][j] = inverse[i][j] / det3by3;
      }
      //cout << inverse[i][j] << " ";
    }
    // cout << endl;
  }
  
  // final return
  return inverse;
}
vector<vector<double>> compute_n(const vector<vector<double>>& matrix3by3, const vector<vector<double>> matrix3by1){
  // Handle Error Cases
  if(matrix3by3.size() != 3 || matrix3by3[0].size() != 3){
    cout << "Did not input 3 by 3 Matrix First" << endl;
    return {};
  }
  if(matrix3by1.size() != 3 || matrix3by1[0].size() != 1){
    cout << "Did not input 3 by 1 Matrix Second" << endl;
    return {};
  }
  // Dot Product Calculations
  vector<vector<double>> output(3, vector<double>(1));
  for(int row = 0; row < 3; row++){
    for(int col = 0; col < 3; col++){
      output[row][0] += (matrix3by3[row][col] * matrix3by1[col][0]);
    }
    //cout << output[row][0] << " ";
  }
  //cout << endl;
  return output;
}
double compute_albedo(vector<vector<double>>& n){
  // Should be a 3 by 1 matrix ==> matrix n
  // Handle Edge Case
  if(n.size() != 3 || n[0].size() != 1){
    cout << "Did not input 3 by 1 Matrix" << endl;
    return 0;
  }
  return sqrt(pow(n[0][0],2) + pow(n[1][0],2) + pow(n[2][0],2));
}
vector<vector<double>> compute_unitNormal(vector<vector<double>>& n){
  // Should be a 3 by 1 matrix ==> matrix n
  // Handle Edge Case
  if(n.size() != 3 || n[0].size() != 1){
    cout << "Did not input 3 by 1 Matrix" << endl;
    return {};
  }
  double length = compute_albedo(n);
  vector<vector<double>> unitNormal(3,vector<double>(1));
  //cout << endl << "UnitNormal:" << endl;
  for(int i = 0; i < 3; i++){
    unitNormal[i][0] = n[i][0]/length;
    //cout << unitNormal[i][0] << " ";
  }
  return unitNormal;
}
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


  // Read in S2 input and convert into 3 by 3 matrix called S2
  std::ifstream file(input_directions_filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << input_directions_filename << std::endl;
    }

  vector<vector<double>> s2_params;
  // typical sStream structure for reading in lines and each value
  std::string line;
  while (std::getline(file, line)) {
      std::vector<double> row;
      std::istringstream iss(line);
      double value;
      while (iss >> value) {
          row.push_back(value);
      }
      s2_params.push_back(row);
    }
  /**
  cout << "Testing Read in Values:" << endl;
  // Testing for Read In
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      cout << s2_params[i][j] << " ";
    }
    cout << endl;
  }
  */

  // Read in images for the Spheres
    // Read in sphere_1
  Image sphere_1;
  if (!ReadImage(input_object_filenames[0], &sphere_1)) {
    cout <<"Can't open file " << &input_object_filenames[0] << endl;
  }
    // Read in sphere_2
  Image sphere_2;
  if (!ReadImage(input_object_filenames[1], &sphere_2)) {
    cout <<"Can't open file " << &input_object_filenames[1] << endl;
  }
    // Read in sphere_3
  Image sphere_3;
  if (!ReadImage(input_object_filenames[2], &sphere_3)) {
    cout <<"Can't open file " << &input_object_filenames[2] << endl;
  }

  // Does not matter which sphere we take, all have the same dimensions
  int num_rows = sphere_3.num_rows();
  int num_cols = sphere_3.num_columns();

  // Albedo Image
  Image albedo_img;
  albedo_img.AllocateSpaceAndSetSize(num_rows, num_cols);
  albedo_img.SetNumberGrayLevels(255);

  // Needle Map Image - base is input Image 1
  Image needle_map = sphere_1;


  // greatest albedo - used to normalize albedo to be between 0-255
  double max_albedo = 0.0;

  // hold the albedo values because SetPixel needs int but albedo is a double before scaling
  vector<vector<double>> albedo_space(num_rows, vector<double>(num_cols));
  // thresholding
  for(int x = 0; x < num_rows; x++){
    for(int y = 0; y < num_cols; y++){
      // convert to double - GetPixel is int
      double i1 = sphere_1.GetPixel(x,y)/(1.0);
      double i2 = sphere_2.GetPixel(x,y)/(1.0);
      double i3 = sphere_3.GetPixel(x,y)/(1.0);

      if(i1 > threshold & i2 > threshold & i3 > threshold){
        // create intensity vector
        vector<vector<double>> intensity = {{i1}, {i2}, {i3}};

        // testing for intensity size:
        // cout << intensity.size() << " " << intensity[0].size() << endl;

        // solve for n
        vector<vector<double>> n = compute_n(InverseMatrix(s2_params), intensity);
        // calculate albedo
        double albedo = compute_albedo(n);

        // determine max
        if (max_albedo < albedo){
          max_albedo = albedo;
        }
        // populate the albedo space - all none changed is 0.0
        albedo_space[x][y] = albedo;

        // Needle Map Creation

        // Gridpoints are Black ==> 
        // Normals are White ==> Unit Normal Vector ==> X10 value to make visible
        // Normal ==> S^-1 [I1, I2, I3] = N ==> Unit Normal == N/Length(N)
        vector<vector<double>> unit_normal = compute_unitNormal(n);
        if(x % step == 0 && y % step == 0 ){
          // scale unit normal vector
          unit_normal[0][0] = unit_normal[0][0] * 10;
          unit_normal[1][0] = unit_normal[1][0] * 10;
          // Draw normal Line
          DrawLine(x, y, x + unit_normal[0][0], y + unit_normal[1][0], 255, &needle_map);
          // make gridPoint Black
          needle_map.SetPixel(x,y,0);
          // make white circle around the gridpoint - similar to sample image
          needle_map.SetPixel(x+1, y, 255);
          needle_map.SetPixel(x-1, y, 255);
          needle_map.SetPixel(x, y+1, 255);
          needle_map.SetPixel(x, y-1, 255);
        }
      }
      else{
        albedo_img.SetPixel(x,y,0); 
      }
    }
  }
  // Standardize albedo values to be between 0-255
  //cout << "Show albedo_space values:" << endl;
  for(int x = 0; x < num_rows; x++){
    for(int y = 0; y < num_cols; y++){
      if(albedo_space[x][y] != 0.0){
        // cout << albedo_space[x][y] << " ";
        albedo_img.SetPixel(x,y, (255 * albedo_space[x][y] / max_albedo));
      }
    }
  }

  // Writing Out for Needle Map
  if (!WriteImage(output_normals_filename, needle_map)){
    cout << "Can't write to file " << output_normals_filename << endl;
  }

  // Writing Out for Albedo Image
  if (!WriteImage(output_albedos_filename, albedo_img)){
    cout << "Can't write to file " << output_albedos_filename << endl;
  }
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
