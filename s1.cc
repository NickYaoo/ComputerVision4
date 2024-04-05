//
// Nicklaus Yao
// Computational Vision Homework 4
// Program 1
//

#include "image.h"
#include "DisjSets.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <vector>

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

  // Read in image
  Image image;
  if (!ReadImage(input_sphere_filename, &image)) {
    cout <<"Can't open file " << input_sphere_filename << endl;
  }

  // Assuming Orthographic Projection - projects circle onto image plane
  // Obtain Clean Binary Circle - Experiment with Thresholding
  vector<double> dict(3,0);

  // We do not need to produce an output of the Binary Circle, but we can for checking
    // Create output image
  Image output;
  output.AllocateSpaceAndSetSize(image.num_rows(), image.num_columns());
  output.SetNumberGrayLevels(255);

 for(int i = 0; i < image.num_rows(); i++){
    for(int j = 0; j < image.num_columns(); j++){
      int color = image.GetPixel(i,j);

      // we check if the pixel is above the threshold
      if(color > threshold){
        output.SetPixel(i,j,255);
        // we need to perform these operations regardless if we made a new object or not
          // add to Area
          dict[2] += 1;
          // add i to avg-x-cord ==> we take the sume of the x-cords first and then we will later find the average
          dict[0] += i;
          // add j to avg-y-cord ==> same as above
          dict[1] += j;
      }
      else{
        output.SetPixel(i,j,0);
      }
    }
 }
     // UPDATE STEP: COMPUTE THE CENTROID using count/area aka find the avg_x and avg_y
    // Use iterators to loop through all keys of the dictionary
      // Update x-center, y-center
      double area = dict[2];
      dict[0] = dict[0] / area; // avg_x
      dict[1] = dict[1] / area; //avg_y


  if (!WriteImage("Binary_Output.pgm", output)){
    cout << "Can't write to file " << "Binary_Output" << endl;
  }

  // Compute the Radius
  // (Rightmost - leftmost + Lowest-Highest) / 2 = Diameter
  // We will go in reverse for the number of rows because we want to find Lowest First
  int highest = image.num_rows()-1;
  int leftmost = image.num_columns()-1;
  int rightmost = 0;
  int lowest = 0;
  bool lowestFound = false;
  for(int i = image.num_rows()-1; i >= 0; i--){
    for(int j = 0; j < image.num_columns(); j++){
      if(output.GetPixel(i,j) == 255){
        // the first occurance of a valid pixel on row i is the lowest possible
        // after finding, we do not care, set lowestFound to true
        if (!lowestFound){
          lowest = i;
          lowestFound = true;
        }
        // each time we find a smaller highest or a smaller leftmost, we update
        // each time we find a greater than current rightmost, we update
        if(j < leftmost){
          leftmost = j;
        }
        if(j > rightmost){
          rightmost = j;
        }
        if(i < highest){
          highest = i;
        }
      }
    }
   }
   cout << rightmost << " " << leftmost << " " << lowest << " " << highest << " " << std::endl;
  double radius = ((rightmost - leftmost) + (lowest - highest)) / 4;

  
  // Write to textfile
  std::ofstream output_info(output_params_filename);
  output_info << int(dict[0]) << " " << int(dict[1]) << " " << int(radius);
  output_info.close();
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
