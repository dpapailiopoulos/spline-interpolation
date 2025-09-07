#include <fstream>
#include "cubic_spline.h"
#include "read_path.h"

namespace spline_interpolation
{
CubicSpline::CubicSpline(const PointsArray& targets, const BoundaryCondition condition) {
    target_points_ = targets;
    spline_size_ = target_points_.rows();
    boundary_condition_ = condition;

    createTridiagonalMatrix();
    solveSpline();
}

CubicSpline::~CubicSpline() {}

void CubicSpline::createTridiagonalMatrix() {
    // All the in between points have the same condition of equal first derivative
    spline_tridiagonal_matrix_.resize(spline_size_, spline_size_);
    spline_tridiagonal_matrix_.diagonal(0).setConstant(4.0);
    spline_tridiagonal_matrix_.diagonal(1).setConstant(1.0);
    spline_tridiagonal_matrix_.diagonal(-1).setConstant(1.0);

    defineBoundaryCondition();
}

void CubicSpline::defineBoundaryCondition() {
    // The boundary conditions are defined by the first and last row of the tridiagonal matrix
    switch (boundary_condition_)
    {
    case BoundaryCondition::ClosedLoop:
        checkClosedLoop();
        // Set first row
        spline_tridiagonal_matrix_.diagonal(0)(0) = 1.0;
        spline_tridiagonal_matrix_.diagonal(1)(0) = 1.0;
        // Set last row
        spline_tridiagonal_matrix_.diagonal(0)(spline_size_ - 1) = -3.0;
        spline_tridiagonal_matrix_.diagonal(-1)(spline_size_ - 2) = -1.0;
        break;
    case BoundaryCondition::Anchored:
        // Set first row
        spline_tridiagonal_matrix_.diagonal(0)(0) = 1.0;
        spline_tridiagonal_matrix_.diagonal(1)(0) = -1.0;
        // Set last row
        spline_tridiagonal_matrix_.diagonal(0)(spline_size_ - 1) = 1.0;
        spline_tridiagonal_matrix_.diagonal(-1)(spline_size_ - 2) = -1.0;
        break;
    case BoundaryCondition::NaturalSpline:
        // Set first row
        spline_tridiagonal_matrix_.diagonal(0)(0) = 1.0;
        // Set last row
        spline_tridiagonal_matrix_.diagonal(0)(spline_size_ - 1) = -2.0;
        break;
    default:
        throw std::runtime_error("CubicSpine -> Unkown boundary condition given");
        break;
    }
}

void CubicSpline::solveSpline() {
    values_matrix_.resize(spline_size_, 2);
    values_matrix_.setZero();
    for (int i{ 1 }; i < spline_size_ - 1; i++)
    {
        // Results from in between conditions
        values_matrix_.row(i) = target_points_.row(i - 1) - 2 * target_points_.row(i) + target_points_.row(i + 1);
    }

    // For closed-loop we have equal first and second derivatives at both end points.
    // The system is no longer tridiagonal but we can solve it using the Sherman-Morrison formula.
    // For more info: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variants
    if (boundary_condition_ == BoundaryCondition::ClosedLoop)
    {
        values_matrix_.row(0) = target_points_.row(1) - target_points_.row(0) - target_points_.row(spline_size_ - 1) + target_points_.row(spline_size_ - 2);

        // Temporary Vectors for Sherman-Morrison formula
        Eigen::VectorXd u{spline_size_};
        u.setZero();
        u(0) = 1.0;
        u(spline_size_ - 1) = 1.0;

        Eigen::VectorXd v{spline_size_};
        v.setZero();
        v(0) = 1.0;
        v(spline_size_ - 2) = 1.0;
        v(spline_size_ - 1) = 2.0;

        Eigen::MatrixXd y{ solveTridiagonal(spline_tridiagonal_matrix_, 6 * values_matrix_) };
        Eigen::VectorXd q{ solveTridiagonal(spline_tridiagonal_matrix_, u) };

        second_derivatives_ = y - (q * v.transpose() * y) / (1 + v.transpose() * q);
    }
    else
    {
        // For other cases we just solve the tridiagonal system
        second_derivatives_ = solveTridiagonal(spline_tridiagonal_matrix_, 6 * values_matrix_);
    }
}

void CubicSpline::checkClosedLoop() {
    if (target_points_.row(0) != target_points_.row(spline_size_ - 1))
    {
        throw std::runtime_error("CubicSpline -> Closed-Loop condition specified while first and last points do not have the same coordinates");
    }
}

LocalSplineIndex CubicSpline::getLocalIndex(double parameter) {
    if (parameter < 0 or parameter > 1)
    {
        throw std::runtime_error("CubicSpline -> Parameter value is out of range");
    }
    LocalSplineIndex local_index{};

    // Defines the interval of the spline parameter on which the spline transitions from one curve
    // segment to the next.
    double parameter_interval{ 1.0 / (spline_size_ - 1) };
    // Current spline segment based on parameter value given
    local_index.spline_segment = parameter / parameter_interval;
    // Current segment parameter
    local_index.local_parameter = std::fmod(parameter, parameter_interval) / parameter_interval;

    /* Due to floating point error, at the last point we may transition to spline segment equal to
     * the number of interpolation point, which is actually out of range since the number of
     * interpolation points is N+1 and of spline segments  is N. So we check for very small values
     * of the local parameter in which case we set that to 1 and transition to the previous spline
     * segment.
     */
    if (local_index.local_parameter < 1e-10 and local_index.spline_segment != 0)
    {
        local_index.local_parameter = 1.0;
        local_index.spline_segment--;
    }
    return local_index;
}

// --- Setters ---

void CubicSpline::setTargetPoints(const PointsArray& targets, const BoundaryCondition condition) {
    target_points_ = targets;
    spline_size_ = target_points_.rows();
    if (condition != BoundaryCondition::BoundaryConditionSize) { boundary_condition_ = condition; }

    createTridiagonalMatrix();
    solveSpline();
}

// --- Getters ---

Point CubicSpline::getSecondDerivative(double parameter) {
    LocalSplineIndex index{ getLocalIndex(parameter) };
    int i{ index.spline_segment };
    double u{ index.local_parameter };
    Point derivative_value{};

    // Equation of second derivative as a function of the already calculated ones
    derivative_value = second_derivatives_.row(i) + (second_derivatives_.row(i + 1) - second_derivatives_.row(i)) * u;

    return derivative_value;
}

Point CubicSpline::getFirstDerivative(double parameter) {
    LocalSplineIndex index{ getLocalIndex(parameter) };
    int i{ index.spline_segment };
    double u{ index.local_parameter };
    Point derivative_value{};

    // Equation of first derivative as a function of second derivatives
    derivative_value = target_points_.row(i + 1) - target_points_.row(i) - 1.0 / 6.0 * second_derivatives_.row(i + 1) - 1.0 / 3.0 * second_derivatives_.row(i) \
        + second_derivatives_.row(i) * u \
        + 1.0 / 2.0 * (second_derivatives_.row(i + 1) - second_derivatives_.row(i)) * std::pow(u, 2);

    return derivative_value;
}

Point CubicSpline::getPoint(double parameter) {
    LocalSplineIndex index{ getLocalIndex(parameter) };
    int i{ index.spline_segment };
    double u{ index.local_parameter };
    Point spline_point{};

    // Equation of point coordinates as a function of second derivatives
    spline_point = target_points_.row(i) \
        + (target_points_.row(i + 1) - target_points_.row(i) - 1.0 / 6.0 * second_derivatives_.row(i + 1) - 1.0 / 3.0 * second_derivatives_.row(i)) * u \
        + 1.0 / 2.0 * second_derivatives_.row(i) * std::pow(u, 2) \
        + 1.0 / 6.0 * (second_derivatives_.row(i + 1) - second_derivatives_.row(i)) * std::pow(u, 3);

    return spline_point;
}

PointsArray CubicSpline::getSecondDerivativeAtTargets() {
    return second_derivatives_;
}

PointsArray CubicSpline::getSplineCurve(const long int resolution) {
    PointsArray spline_curve{ resolution, 2 };
    // Iterate through all points and get spline coordinates
    for (long int i{ 0 }; i <= resolution - 1; i++)
    {
        // i is in the range [0,resolution-1], need to transfer that to [0,1]
        double param{ static_cast<double>(i) / static_cast<double>(resolution - 1) };
        spline_curve.row(i) = getPoint(param);
    }
    return spline_curve;
}

double CubicSpline::getCurvature(double parameter) {
    Point first_der{ getFirstDerivative(parameter) };
    Point second_der{ getSecondDerivative(parameter) };
    double curvature{
        (first_der(0) * second_der(1) - first_der(1) * second_der(0)) / (std::pow(std::pow(first_der(0), 2) + std::pow(first_der(1), 2), 3.0 / 2.0))
    };
    return curvature;
}

double CubicSpline::getTangent(double parameter) {
    Point first_der{ getFirstDerivative(parameter) };
    return std::atan2(first_der(1), first_der(0));
}
} // namespace spline_interpolation


// Run for test path
int main() {
    spline_interpolation::PointsArray midpoints{ readTrack("../point_paths/test_path.txt") };
    midpoints.conservativeResize(midpoints.rows() + 1, midpoints.cols());
    midpoints.row(midpoints.rows() - 1) = midpoints.row(0);

    spline_interpolation::CubicSpline spline{midpoints, spline_interpolation::BoundaryCondition::ClosedLoop};

    std::cout << "At s = 0.67" << '\n';
    std::cout << "Spline coordinates: " << spline.getPoint(0.67) << '\n';
    std::cout << "Spline curvature: " << spline.getCurvature(0.67) << '\n';
    std::cout << "Spline tangent: " << spline.getTangent(0.67) << "\n\n";

    std::cout << "Start: x' = " << spline.getFirstDerivative(0) << " , x'' = " << spline.getSecondDerivative(0) << '\n';
    std::cout << "End: x' = " << spline.getFirstDerivative(1) << " , x'' = " << spline.getSecondDerivative(1) << "\n\n";

    std::ofstream spline_curve("cubic_spline_points.txt");
    spline_curve << spline.getSplineCurve(2000);
    spline_curve.close();

    return 0;
}