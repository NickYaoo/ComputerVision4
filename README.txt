Nicklaus Yao
Computational Vision Homework 4

**Thresholds should also be at the file "thresholds.txt" (two thresholds, one for s1 and one for s3).
1. All four parts of the assignment have been completed.
2. Have some problems when ...
3. All programs runs as directed in assignment details.
4. Input files same as the ones provided.


----------------------
To compile in Linux:
----------
 
   make all

----------

Program_1 s1
    - Threshold used is: 89
    - Threshold is stored in the first line of threshold.txt
    - Output: X and Y coordinates of Center of Sphere(Image), radius of the sphere
    - X and Y are found using by first thresholding the image and then computing the centroid
    - Radius is found by taking the Average distance between the highest and lowest pixel of the circle
        and the leftmost and rightmost pixel of the circle as the diameter. Then we divide by 2 for the radius.

To run:
---------
./s1 sphere_0.pgm {input int threshold} s1_output.txt

Program_2 s2
    - The formula used for computation of sphere normals was: 1/(N1^2+ N2^2+ N3^2) * {N1, N2, N3} = Unit normals
        - Note: This is the unit normals of the brightest pixel X, Y of each of the photos
        - The non unit normal is found by {X-X0, Y-Y0, sqrt(r^2 - (X-X0)^2 - (Y-Y0)^2)} where
            X0, Y0, and r was the output of Program s1 and the input file(s1_output.txt)
    - Output: X, Y, Z components of the unit normal vector scaled by the brightness of the brightest pixel

To run:
---------
./s2 s1_output.txt sphere_1.pgm sphere_2.pgm sphere_3.pgm s2_output.txt

Program_3 s3
    - Threshold used is: 87
    - Threshold is stored in the second line in thresholds.txt
    - We find the albedo and normal of each point that is above the threshold in all three sphere images
        by taking the product of the inverse of the matrix given by s2_output.txt and the vector {I1, I2, I3}
        - {I1,I2,I3} is the brightness of the pixel in question in every single image - all should be above threshold
    - Output: 
        - Normal Needle Map: Created a gridmap (spaced with the step parameter) of black gridpoints with a white line representing
            the normal of the pixel 
        - Albedo Image: Calcuate the albedo of each pixel above the threshold and scale it to fit 0-255 so that the pixel is this alebdo value scaled

To run:
---------
./s3 s2_output.txt sphere_1.pgm sphere_2.pgm sphere_3.pgm {input int step parameter } {input int threshold} output_normals.pgm output_albedo.pgm