/********************************************//**
* @mainpage AreaCon: A C++ Library for Area-Constrained Partitioning
* @details AreaCon is a light-weight, C++ library for carrying out area-constrained partitioning operations in floating-point precision. The library is primarily a numerical implementation of the area-constrained partitioning algorithm by Patel, Frasca, and Bullo, 2014 (see http://www.areacon.org/References ).
                                               
    Low-level source code documentation is provided here; however, additional discussion about how to use the code along with numerical examples can be found on the project's Wiki ( https://github.com/jrpeters/AreaCon/wiki ). 
                                               
    Please cite AreaCon whenever possible. A bibtex citation is provided below:
    @verbatim
    @Misc{AreaCon:16,
    author =        {J.R. Peters and Contributors},
    title =         {The {AreaCon} Library},
    howpublished =  {\texttt{http://www.areacon.org}},
    year =          {2016},
    note =          {v. 1},
    abstract =      {A C++ library for area-constrained partitioning.},
     }
    @endverbatim
                                               
    For more information on the code and information on how to contact me, please see http://www.areacon.org
* @author Jeffrey R. Peters ( http://www.jeffreyrpeters.com , https://github.com/jrpeters )
* @version 1.0
* @date Feb. 2016
* @copyright Copyright &copy; 2016. The Regents of the University of California. All rights reserved. Licensed pursuant to the terms and conditions available for viewing at: http://opensource.org/licenses/BSD-3-Clause .
                                               
* @copyright The Clipper library ( www.angusj.com/delphi/clipper.php ), as used by AreaCon, is also subject to the following:
                                               
* @copyright Copyright &copy; 2016. The Regents of the University of California. Distributed under the Boost Software License, Version 1.0. (See http://www.boost.org/LICENSE_1_0.txt ).
***********************************************/


#ifndef __AreaCon__areacon__
#define __AreaCon__areacon__

#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "clipper.hpp"




namespace AreaCon {
    // Point Class-------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    /**
     * The base class for defining a point in the two-dimensional plane.
     * @author Jeffrey R. Peters
     */
    class Point
    {
    public:
        //@{
        double x;
        double y;
        //@}
        static double Robustness_Constant;/**<A constant used to enhance numerical robustness. Loosely, when a Euclidean distance is less than the robustness constant, the distance is considered to be 0.*/
        //@{
        /**
         * Default Constructor
         * @param[in] x The x coordinate
         * @param[in] y The y coordinate
         */
        Point(const double x = INFINITY, const double y = INFINITY);
        //@}
        
        //@{
        /**
         * Calculates the Euclidean norm of the point
         * @return Euclidean norm
         */
        double Norm() const;
        /**
         * Adds the point Test to the current point component-wise
         * @return Sum = Point(Test.x+x, Test.y+y)
         */
        Point AddPoint(const Point Test) const;
        /**
         * Returns a point (vector) of length Norm that points in a direction perpendicular to the line between Test and the current point.
         * @param[in] Test The test point
         * @param[in] Norm The desired norm of the resulting point (vector)
         * @return Perpendicular Point
         */
        Point FindPerpDirection(const Point Test, const double Norm) const;
        /**
         * Reverses the orientation of the current point (vector), i.e., multiplies by -1
         */
        void FlipDirection(void);
        /**
         * Finds the perpendicular distance of the point to the line connecting Test1, Test2
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test Point
         * @return The perpendicular distance
         */
        double PerpDistanceToLine(const Point Test1, const Point Test2) const;
        /**
         * Tests to see if the points Test1, Test2 are numerically collinear with the point object.
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @return Indicator of collinearity
         */
        bool AreCollinear(const Point Test1, const Point Test2)const;
        /**
         * Tests to see if the point object lies (numerically) on the line between Test1, Test2
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @return Indicator stating if point object is between the test points
         */
        bool AreBetween(const Point Test1, const Point Test2) const;
        /**
         * Multiplies the point object (vector) by a constant factor
         * @param[in] factor The multiplication factor
         */
        void Mult(const double factor);
        //@}
        
        //@{
        /**
         * Tests to see if Test1, Test2 are equal.
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @return Equality indicator
         */
        static bool IsEqual(const Point Test1, const Point Test2);
        /**
         * Finds the Euclidean distance between the points Test1, Test2.
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @return Euclidean distance
         */
        static double Distance(const Point Test1, const Point Test2);
        /**
         * Finds a point along the line connecting Test1, Test2 that is a normalized distance away from the point Test1, e.g, distance = 0.5 corresponds to the point that is exactly half-way between Test1, Test2.
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @param[in] distance The normalized distance along the line connecting Test1, Test2
         * @return The resulting point along the line
         */
        static Point FindPointAlongLine(const Point Test1, const Point Test2, const double distance);
        /**
         * Finds the Euclidean norm of the point Test
         * @param[in] Test The test point
         * @return The Euclidean distance
         */
        static double Norm(const Point Test);
        /**
         * Adds the points Test1, Test2 component-wise
         * @param[in] Test1 The first summand
         * @param[in] Test2 The second summand
         * @return Sum of the points Test1, Test2
         */
        static Point AddPoints(const Point Test1, const Point Test2);
        /**
         * Returns a point (vector) of length Norm that points in a direction perpendicular to the line between Test1 and Test2
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @param[in] Norm The desired norm of the resulting point (vector)
         * @return Perpendicular Point
         */
        static Point FindPerpDirection(const Point Test1, const Point Test2, const double Norm);
        /**
         * Finds the perpendicular distance of Test3 to the line connecting Test1, Test2
         * @param[in] Test1 The first point defining the line
         * @param[in] Test2 The second point defining the line
         * @param[in] Test3 The test point
         * @return The perpendicular distance
         */
        static double PerpDistanceToLine(const Point Test1, const Point Test2, const Point Test3);
        /**
         * Tests to see if Test3 lies (numerically) on the line between Test1, Test2
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @param[in] Test3 The third test point
         * @return Indicator stating if point object is between the test points
         */
        static bool AreCollinear(const Point Test1, const Point Test2, const Point Test3);
        /**
         * Tests to see if Test3 lies (numerically) on the line between Test1, Test2
         * @param[in] Test1 The first test point
         * @param[in] Test2 The second test point
         * @param[in] Test3 The third test point
         * @return Indicator stating if point object is between the test points
         */
        static bool AreBetween(const Point Test1, const Point Test2, const Point Test3);
        /**
         * Returns the end-points of the line which forms the intersection of the lines connecting p1, p2 and p3, p4, respectively. Note that p1,p2,p3,p4 must be collinear. If the lines do not intersect, or only intersect at a single point, the return vector will contain less than 2 entries
         * @param[in] p1 The first end-point of the first line
         * @param[in] p2 The second end-point of the first line
         * @param[in] p3 The first end-point of the second line
         * @param[in] p4 The second end-point of the second line
         * @return A vector containing the end-points of the intersection line;
         */
        static std::vector<Point> FindCollinearIntersection(const Point p1, const Point p2, const Point p3, const Point p4);
        
    };
    // Poly Class-------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    /**
     * The base class for defining (convex) polygons.
     * @author Jeffrey R. Peters
     */
    class Poly
    {
    public:
        //@{
        /**
         * Default Constructor
         * @param[in] Vertices The points defining the vertices of the polygon in counter-clockwise order (the first vertex is not repeated).
         */
        Poly(std::vector<Point> Vertices = {});
        //@}
        /**
         * A function used for setting the vertices of the polygon.
         * @param[in] Vertices The points defining the vertices of the polygon in counter-clockwise order (the first vertex is not repeated).
         * @param[in] GetExtrema Flag that determines whether the extrema should be calculated and re-set
         */
        void SetVertices(const std::vector<Point> Vertices, const bool GetExtrema = true);
        /**
         @return The list of Vertices
         */
        std::vector<Point> GetVertices(void) const;
        /**
         * @return The number of vertices.
         */
        int GetNVertices(void) const;
        /**
         *  Determines if the point Test lies within the polygon.
         * @param[in] Test The test point
         * @return Indicator of whether or not Test is inside the polygon
         */
        bool pnpoly(const Point Test) const;
        /**
         * Returns the extreme x and y values.
         * @param[out] minx, miny, maxx, maxy
         */
        void GetExtrema(double &minx, double &miny, double &maxx, double &maxy) const;
        
    private:
        double minx;/**< The minimum x coordinate of the polygon*/
        double miny;/**< The minimum y coordinate of the polygon*/
        double maxx;/**< The maximum x coordinate of the polygon*/
        double maxy;/**< The maximum y coordinate of the polygon*/
        std::vector<Point> Vertices;/**< The Vertices of the polygon*/
        int NPoly;/**< The number of vertices that the polygon has*/
        /**
         * A function used to initialize the polygon by calculating extrema and running checks for dimensional consistency.
         */
        void InitializePoly(void);
    };
    // Mult_Array Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    /**
     * Container for storing two-dimensional arrays of size NxN.
     * @author Jeffrey R. Peters
     */
    class Mult_Array
    {
    public:
        //@{
        /**
         * Constructor
         * @param[in] N The size of the array
         */
        Mult_Array(const int N);
        /** Copy Constructor
         * @param[in] obj Element to by copied
         */
        Mult_Array(const Mult_Array &obj);
        /** Copy Assignment Operator
         * @param[in] obj Element to by copied
         */
        Mult_Array& operator=(const Mult_Array &obj);
        /**
         * Destructor. Cleans up dynamically allocated memory.
         */
        ~Mult_Array(void);
        //@}
        double **Array;/**<A two-dimensional array*/
        const int N;/**<The size of the array.*/
    };
    // Delaunay Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    /**
     * Container for storing Delaunay (duel) graphs.
     * @author Jeffrey R. Peters
     */
    class DelaunayGraph
    {
    public:
        //@{
        /**
         * Constructor.
         * @param[in] NRegions The number of regions in the partition or the number of nodes in the Delaunay graph
         */
        DelaunayGraph(const int NRegions);
        /** Copy Constructor
         * @param[in] obj Element to by copied
         */
        DelaunayGraph(const DelaunayGraph &obj);
        /** Copy Assignment Operator
         * @param[in] obj Element to by copied
         */
        DelaunayGraph& operator=(const DelaunayGraph &obj);
        /**
         * Destructor. Cleans up dynamically allocated memory.
         */
        ~DelaunayGraph();
        //@}
        
        Point ***Graph;/**<A multi-dimensional array whose (i,j)-th entry holds the endpoints of the line segment that is shared by region i and region j. If the two regions do not share a common edge, then at least 1 index of the (i,j)-th entry will equal INFINITY.*/
        const int NRegions;/**<The number of regions under consideration.*/
    };
    // Int_Params Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    /**
     * A container for values used in quickly calculating area integrals over polygonal regions of interest.
     * @author Jeffrey R. Peters
     */
    class Int_Params
    {
    public:
        //@{
        std::vector<double> Coefficient_a, Coefficient_b, Coefficient_c, Coefficient_d;/**<Coefficients used in quickly calculating area integrals. Usually populated as a part of the function Density.FindIntegralCoefficients.*/
        std::vector<double> Int, Intx, Inty;/**<Parameters representing area integrals over grid squares. Int represents the total integral, Intx represents the integral of x*f(x,y), and Inty represents the integral of y*f(x,y). Usually populated as a part of the function Density.FindIntegralVector.*/
        double Unweighted_Area; /**<The overall area of some polygonal region of interest*/
        //@}
        
        //@{
        /*
         * Default Constructor.
         * @param[in] Coefficient_a A coefficient used in quickly calculating area integrals
         * @param[in] Coefficient_b A coefficient used in quickly calculating area integrals
         * @param[in] Coefficient_c A coefficient used in quickly calculating area integrals
         * @param[in] Coefficient_d A coefficient used in quickly calculating area integrals
         * @param[in] Int Vector representing the total integral over individual grid squares
         * @param[in] Intx Vector representing the integral of the function x*f(x,y) over individual grid squares
         * @param[in] Inty Vector representing the integral of the function y*f(x,y) over individual grid squares
         * @param[in] Unweighted_Area The overall area of some polygonal region
         */
        Int_Params(std::vector<double> Coefficient_a = {}, std::vector<double> Coefficient_b= {}, std::vector<double> Coefficient_c= {}, std::vector<double> Coefficient_d= {}, std::vector<double> Int = {}, std::vector<double> Intx = {}, std::vector<double> Inty = {}, double UnweightedArea = 0);
        //@}
        /**
         * Checks to make sure that the user-input parameters satisfy the required bounds
         */
        void CheckParameters(void) const;
    };
    // Parameters Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    class Parameters
    {
    public:
        //@{
        /**
         * Default Constructor. Default values are recommended to produce reasonable solutions with reasonable efficiency in most scenarios.
         * @param[in] line_int_step Spacing parameter used for calculating line integrals
         * @param[in] weights_step Step-size used in updating the weighting parameters
         * @param[in] centers_step Step-size used in updating the center locations
         * @param[in] volume_tolerance Parameter used in determining whether desired volumes have been acheived
         * @param[in] convergence_criterion Parameter used as an algorithmic stopping criterion
         * @param[in] max_iterations_volume Upper bound on the number of volumetric iterations
         * @param[in] max_iterations_centers Upper bound on the number of centroidal movement iterations
         * @param[in] Volume_Lower_Bound A lower bound on the weighted area of each region
         * @param[in] Robustness_Constant A constant used to enhace numerical robustness (see Point class)
         */
        Parameters(const double line_int_step = 0.1,const double weights_step = 0.1, const double centers_step = 1,const double volume_tolerance = 0.002, const double convergence_criterion = 0.02, const int max_iterations_volume = 200, const int max_iterations_centers = 500, const double Volume_Lower_Bound = 10e-6, const double Robustness_Constant = 10e-8);
        //@}
        //@{
        const double line_int_step;/**<Spacing parameter used for calculating line integrals*/
        const double weights_step; /**<Step-size used in updating the weighting parameters*/
        const double centers_step; /**<Step-size used in updating the center locations*/
        const double volume_tolerance; /**<Parameter used in determining whether desired volumes have been acheived*/
        const double convergence_criterion;/**<Parameter used as an algorithmic stopping criterion*/
        const int max_iterations_volume;/**<Upper bound on the number of volumetric iterations*/
        const int max_iterations_centers;/**<Upper bound on the number of centroidal movement iterations*/
        
        const double Volume_Lower_Bound;
        const double Robustness_Constant;
        //@}
    private:
        /**
         * Checks to see if user-input parameters satisfy the required bounds
         */
        void CheckParameters(void);
    };
    
    // Density Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    /**
     * The base class for defining probability density functions.
     * @author Jeffrey R. Peters
     */
    class Density
    {
    public:
        //@{
        /**
         * Default Constructor.
         */
        Density();
        /**
         * @param[in] Region The (convex) polygonal region of interest.
         * @param[in] Nx The number of grid points in the x direction.
         * @param[in] Ny The number of grid points in the y direction.
         * @param[in] Values A vector containing the value of the density function at the grid-point locations (the value at the (i,j)-th grid point is stored in the (Ny*i+j)-th entry of Values.
         */
        Density(const Poly Region, const int Nx = 0, const int Ny = 0, const std::vector<double> Values = {});
        //@}
        /** Function used to set a new polygonal region of interest.
         * @param[in] Region The (convex) polygonal region of interest.
         * @param[in] Nx The number of grid points in the x direction.
         * @param[in] Ny The number of grid points in the y direction.
         * @param[in] Values A vector containing the value of the density function at the grid-point locations (the value at the (i,j)-th grid point is stored in the (Ny*i+j)-th entry of Values.
         */
        void SetNewRegion(const Poly Region, const int Nx = 0, const int Ny = 0, const std::vector<double> Values = {});
        /**
         * Function used to re-set the values of the grid-parameters.
         * @param[in] Nx The number of grid points in the x direction.
         * @param[in] Ny The number of grid points in the y direction.
         * @param[in] Values A vector containing the value of the density function at the grid-point locations (the value at the (i,j)-th grid point is stored in the (Ny*i+j)-th entry of Values.
         */
        void SetParameters(const int Nx,const int Ny,const std::vector<double> Values);
        /**
         * @return Nx
         */
        int GetNx(void) const {return Nx;};
        /**
         * @return Ny
         */
        int GetNy(void) const {return Ny;};;
        /**
         * @return Region
         */
        Poly GetRegion(void) const {return Region;};;
        /**
         * @return GridInRegion
         */
        std::vector<bool> GetGridInRegion(void) const {return GridInRegion;};
        /**
         * @return Integral
         */
        Int_Params GetIntegral(void) const {return Integral;};
        /**
         * Returns the extreme x and y values.
         * @param[out] minx, miny, maxx, maxy
         */
        void GetExtrema (double &minx, double &miny, double &maxx, double &maxy) const {minx = this->minx; miny = this->miny; maxx = this->maxx; maxy = this->maxy;}
        /**
         * Calculates the line integral of the density function over the straight line connecting points p1, p2
         * @param[in] spacing The spacing between evaluation points
         * @param[in] p1, p2 The endpoints of the line in question
         * @return The value of the line integral
         */
        double LineIntegral(double spacing, const Point &p1, const Point &p2) const;
        /**
         * Evaluates the integral of the density over the polygon Region
         * @param[in] Region The polygon over which the integral is evaluated
         * @return The weighted area of the region
         */
        double CalculateWeightedArea(const Poly Region) const;
        /**
         * Calculates the centroid of the polygon Region with respect to the density
         * @param[in] Region The polygon of interest
         * @param[in] Volume The total volume of the region in question
         * @return The location of the centroid
         */
        Point CalculateCentroid(const Poly Region, const double &Volume) const;
        /**
         * Sets the lower volume bound (default = 0). This bound is used to avoid numerical instability in partition calculations.
         * @param[in] VolumeLowerBound The new bound value;
         */
        void SetVolumeLowerBound(const double VolumeLowerBound);
        /**
         * Returns the lower volume bound.
         * @return Volume_Lower_Bound;
         */
        double GetVolumeLowerBound(void);
        
        void WriteToFile(const std::string filename)const;

    private:
        Poly Region;/**<The region of interest*/
        int Nx;/**<The number of grid points in the x direction*/
        int Ny;/**<The number of grid points in the y direction*/
        double dx;/**<The grid spacing: dx = maxx-minx/(Nx-1)*/
        double dy;/**<The grid spacing: dy = maxy-miny/(Ny-1)*/
        double minx;/**< The minimum x coordinate of the polygon*/
        double miny;/**< The minimum y coordinate of the polygon*/
        double maxx;/**< The maximum x coordinate of the polygon*/
        double maxy;/**< The maximum y coordinate of the polygon*/
        double Volume_Lower_Bound;/**<A lower bound on any calculated volume (default = 0). This parameter is used to avoid numerical instability in partition calculations.*/
        std::vector<double> Values;/**<A vector containing the value of the density function at the grid-point locations (the value at the (i,j)-th grid point is stored in the (Ny*i+j)-th entry of Values.*/
        std::vector<bool> GridInRegion;/**<A vector whose entries indicate whether or not grid points lie within Region. If the (i,j)-th grid point lies within the polygonal region of interest, then the (Ny*i+j)-th entry of GridInRegion is true, otherwise it is false.*/
        Int_Params Integral;/**<Container that holds parameters relevant to quickly calculating area integrals.*/
        
        /**
         * A function that checks consistency of parameter sizes.
         */
        void CheckParameterSizes(void);
        /**
         * Finds and sets the extrema for the polygon Region, i.e., sets minx, maxx, miny, maxy.
         */
        void SetExtrema(void);
        /**
         * Calculates the value of dx and dy.
         */
        void Setdxy(void);
        /**
         * Performs pre-processing steps that are necessary for fast integration.
         */
        void PreprocessIntegral(void);
        /**
         * Creates the coefficients that are used in numerically evaluating integrals. Results are stored in the associated Int_Params container Integral
         */
        void CreateIntegralCoefficients(void);
        /**
         * Pre-calculates and stores the value of relevant integrals over individual grid squares. Results are stored in the associated Int_Params container Integral
         * @return The value of the total integral of the density under the region of interest.
         */
        double CreateIntegralVector(void);
        /**
         * Makes sure that the values stored in the integral vector are normalized so that their sum is equal to 1 over the region of interest.
         * @param[in] Total The value of the total integral of the density under the region of interest.
         */
        void NormalizeIntegralVector(const double &Total);
        
        /**
         * Uses interpolation to find the value of the density at the point Test, which is not necessarily a grid point.
         * @param[in] Test The point of interest. 
         * @return The value of the density function at Test
         */
        double InterpolateValue(const Point &Test) const;
        
        /**
         * Returns the world coordinates of the grid-point associated with the ii-th entry of the vector Values.
         * @return The world coordinates of the associated grid point.
         */
        Point ConvertIndextoWorld(const int ii) const;
    };
    // Partition Class--------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------
    /**
     * The base class for defining a Partition.
     * @author Jeffrey R. Peters
     */
    class Partition
    {
    public:
        //@{
        /**
         * Constructors.
         */
        /**
         * @param[in] NRegions The number of regions desired.
         * @param[in] Prior The density function goverining partition creation
         * @param[in] desired_area A vector containing the desired areas of the regions in the resulting configuration. Defaults to equal area.
         * @param[in] Alg_Params Various algorithmic parameters
         */
        Partition(int NRegions = 0, Density Prior = Density(), std::vector<double> desired_area = {}, Parameters Alg_Params = Parameters());
        //@}
        //@{
        /**
         * Sets the partition variables.
         * @param[in] NRegions The number of regions desired.
         * @param[in] Prior The density function goverining partition creation
         * @param[in] desired_area A vector containing the desired areas of the regions in the resulting configuration. Defaults to equal area.
         */
        void SetPartitionVariables(int NRegions = 0, Density Prior = Density(), std::vector<double> desired_area = {});
        /**
         * Used to initialize algorithmic process variables Centers and Weights. If no input is given, default centers and weights are created.
         * @param[in] Centers The initial center locations.
         * @param[in] Weights The initial weight values.
         */
        void InitializePartition(std::vector<Point> Centers = {},std::vector<double> Weights = {});
        /**
         * @return The current value of Covering
         */
        std::vector<Poly> GetCovering(void){return Covering;}
        /**
         * @return The current value of Centers
         */
        std::vector<Point> GetCenters(void){return Centers;}
        /**
         * @return The current value of Weights
         */
        std::vector<double> GetWeights(void){return Weights;}
        /**
         * The main function used for calculating partitions. Partitions are calculated and the resultant configuration is stored in the containers Centers and Covering. If WriteToFile = true, then the evolution of the centers and partitions will be written to the files filename_centers and filename_partitions, respectively.
         * @param[in] WriteToFile Flag indicating if result should be written to file
         * @param[in] filename_partition Output filename to be written with partition data
         * @param[in] filename_centers Output filename to be written with center data
         */
        void CalculatePartition(bool WriteToFile, std::string filename_partition = "", std::string filename_centers = "" );
        //@}
    private:
        //@{
        std::vector<Point> Centers; /**<The vector of center locations.*/
        std::vector<Poly> Covering; /**<The covering of the area of interest.*/
        std::vector<double> Weights; /**<The vector of weights associated with each area.*/
        //@}
        //@{
        const Parameters Alg_Params;/**<Algorithmic parameters.*/
        std::vector<double> desired_area;/**<A vector specifying the desired areas of the resultant configurations.*/
        Density Prior;/**<The prior probability density function.*/
        int NRegions;/**<The number of regions desired.*/
        //@}
        //@{
        /**
         * Checks Parameters to ensure consistent sizes
         */
        void CheckParams(void);
        /**
         * Creates default centers within the region defined by the polygon Region. The parameter multiplier is used to create center points which are ensured to lie within the polygon.
         * @param[in] Region The region of interest
         * @param[in] multiplier A parameter used in center creation
         * @return A flag indicating whether centers were successfully created
         */
        bool CreateDefaultCenters(const Poly Region, const double multiplier);
        /**
         * Creates default centers within the region defined by the polygon Region. The parameter multiplier is used to create center points which are ensured to lie within the polygon.
         * @param[in] Region The region of interest
         * @param[in] initial_multiplier A parameter used in center creation
         * @param[in] max_steps An upper bound on the number of center creation attempts
         */
        void CreateDefaultCenters(const Poly Region, const double initial_multiplier, const int max_steps);
        /**
         * Creates the power diagram generated from the current values of Centers and Weights.
         */
        bool CreatePowerDiagram(void);
        /**
         * A function used to eliminate redundancies in the covering that may arise due to numerical error.
         * @param[in] tolerance A tolerance for determining whether distinct vertices should be combined into a single vertex
         * @param[in] mult A multiplier that affects the degree of numerical accuracy
         */
        void CleanCovering(const double tolerance, const long int &mult);
        /**
         * Creates the Delaunay graph (Duel Graph) based on the partition stored in Covering.
         * @param[out] Delaunay The Delaunay graph.
         */
        void CreateDelaunayGraph(DelaunayGraph &Delaunay) const;
        /**
         * Update the center locations based on the current configuration
         * @param[in] volumes The current volumes of the regions in Covering
         * @return The sum of the squares of the Euclidean distance moved by each of the centers
         */
        double GradientStepCenter(const std::vector<double> &volumes);
        /**
         * Update the center locations based on the current configuration (This function is used for the initial center step, and can speed convergence in cases where the center stepsize parameter is not chosen equal to 1).
         * @param[in] temp_step A step-size that overwrites the value indicated in Alg_Params.
         * @param[in] volumes The current volumes of the regions in Covering
         * @return The sum of the squares of the Euclidean distance moved by each of the centers
         */
        void GradientStepCenter(const double &temp_step, const std::vector<double> &volumes);
        /**
         * Update the weights based on the current configuration.
         * @param[in] volumes The current volumes of the regions in Covering
         * @param[in] SharedEdges The current Delaunay graph.
         */
        void GradientStepWeights(const std::vector<double> &volumes, const DelaunayGraph &SharedEdges);
        /**
         * Calculates a measure of volumetric error
         * @param[in] volumes The volumes of the regions in Covering.
         * @return A measure of the volumetric error.
         */
        double CalculateError(const std::vector<double> &volumes);
        /**
         * A function that calculates the volumes of the regions in Covering
         * @return A vector containing the volumes of the regions in Covering.
         */
        std::vector<double> CalculateVolumes(void);
    };
    
}//the AreaCon namespace
#endif /* defined(__AreaCon__areacon__) */
