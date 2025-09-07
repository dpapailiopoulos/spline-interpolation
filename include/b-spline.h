#ifndef B_SPLINE_H
#define B_SPLINE_H

#include <iostream>
#include <array>
#include <cmath>
#include <eigen3/Eigen/Dense>

#include "solve_tridiagonal.h"

namespace spline_interpolation
{
using PointsArray = Eigen::MatrixX2d;
using Point = Eigen::RowVector2d;

// All boundary conditions that have been implemented for the b-spline
enum BoundaryCondition {
    OpenLoop,
    ClosedLoop,
    BoundaryConditionSize,
};

// Contains the current spline segment index and the segment parameter value for a specific global
// parameter value. Used for sampling points along the spline.
struct LocalSplineIndex {
    int spline_segment;
    double local_parameter;
};

class BSpline {
private:
    // A vector of all the points we want to interpolate the spline through
    PointsArray target_points_;

    // Number of points to be interpolated
    int spline_size_;

    // Check boundary condition struct
    BoundaryCondition boundary_condition_;

    /* The matrix representing the system of equations used to calculate the spline.
     * Most coefficient will be zero so the matrix is sparse and we can proccess it to make it banded as well,
     * lowering solving time by using the Thomas algorithm implemented in the solve_tridiagonal.h file.
     */
    TriDiagonal spline_tridiagonal_matrix_;

    /* These are the values that should be passed to the Thomas algorithm. It mostly contains the interpolated points.
     * For the OpenLoop boundary condition the extra values are the derivatives at the end points that can be defined
     * within the solveSpline() method.
     * For the ClosedLoop boundary condition no extra values must be determined since the last three control points
     * are wrapped around the first three in order to satisfy the equality of the second derivatives at the end points.
     */
    PointsArray values_matrix_;

    // The B-Spline control points
    PointsArray control_points_;

    // Creates the tridiagonal matrix that represents the system of equations for solving the spline interpolation
    void createTridiagonalMatrix();

    // Solves the system of equations and calculates the values of the control points defining the interpolation curve
    void solveSpline();

    // Checks if the first and last points in the target_points_ matrix have the same coordinates
    void checkClosedLoop();

    /* This method returns the current spline segment and the local parameter of said segment based on
     * the global parameter given as argument. Note that both the local and global parameter should
     * have values in the range [0,1], with 0 indicating the first interpolated point and 1 the last.
     * If an out of range parameter value is given a runtime error will be thrown.
     */
    LocalSplineIndex getLocalIndex(double parameter);

public:
    // Constructor
    BSpline(const PointsArray& targets, const BoundaryCondition condition);

    // Deconstructor
    ~BSpline();

    // --- Setters ---

    // Sets new target points to be interpolated and recalculates the entire spline. It is possible to change the
    // boundary condition as well. If not specified it will not change.
    void setTargetPoints(const PointsArray& targets,
        const BoundaryCondition condition = BoundaryCondition::BoundaryConditionSize);

    // --- Getters ---

    // Returns the second derivative of the spline based on the given spline parameter
    Point getSecondDerivative(double parameter);
    // Returns the first derivative of the spline based on the given spline parameter
    Point getFirstDerivative(double parameter);
    // Returns the spline coordinates of the point represented by the given parameter
    Point getPoint(double parameter);
    // Returns the spline curvature at the point represented by the given parameter
    double getCurvature(double parameter);
    // Returns the spline tangent angle in radians at the point represented by the given parameter
    double getTangent(double parameter);
    // Returns the calculated control points
    PointsArray getControlPoints();
    // Samples through the range of the curve returning as many sample points as defined by the resolution
    PointsArray getSplineCurve(const long int resolution = 2000);
};
} // namespace spline_interpolation
#endif