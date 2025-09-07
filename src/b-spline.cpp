#include <fstream>
#include "b-spline.h"
#include "read_path.h"

namespace spline_interpolation
{
BSpline::BSpline(const PointsArray& targets, const BoundaryCondition condition) {
    target_points_ = targets;
    spline_size_ = target_points_.rows();
    boundary_condition_ = condition;

    createTridiagonalMatrix();
    solveSpline();
}

BSpline::~BSpline() {}

void BSpline::createTridiagonalMatrix() {
    switch (boundary_condition_)
    {
    case BoundaryCondition::ClosedLoop:
        checkClosedLoop();
        // Closed loop b-spline requires two extra control points
        spline_tridiagonal_matrix_.resize(spline_size_ + 2, spline_size_ + 2);
        spline_tridiagonal_matrix_.diagonal(0).setConstant(1.0);
        spline_tridiagonal_matrix_.diagonal(1).setConstant(4.0);
        spline_tridiagonal_matrix_.diagonal(2).setConstant(1.0);

        // Set last three rows
        spline_tridiagonal_matrix_.row(spline_size_ - 1).setZero();
        spline_tridiagonal_matrix_.row(spline_size_).setZero();
        spline_tridiagonal_matrix_.row(spline_size_ + 1).setZero();

        spline_tridiagonal_matrix_(spline_size_ - 1, 0) = -1.0;
        spline_tridiagonal_matrix_(spline_size_ - 1, spline_size_ - 1) = 1.0;

        spline_tridiagonal_matrix_(spline_size_, 1) = -1.0;
        spline_tridiagonal_matrix_(spline_size_, spline_size_) = 1.0;

        spline_tridiagonal_matrix_(spline_size_ + 1, 2) = -1.0;
        spline_tridiagonal_matrix_(spline_size_ + 1, spline_size_ + 1) = 1.0;

        break;
    case BoundaryCondition::OpenLoop:
        spline_tridiagonal_matrix_.resize(spline_size_, spline_size_);
        spline_tridiagonal_matrix_.diagonal(-1).setConstant(1.0);
        spline_tridiagonal_matrix_.diagonal(0).setConstant(4.0);
        spline_tridiagonal_matrix_.diagonal(1).setConstant(1.0);

        // Set first and last row
        spline_tridiagonal_matrix_(0, 1) = 2.0;
        spline_tridiagonal_matrix_(spline_size_ - 1, spline_size_ - 2) = 2.0;

        break;
    default:
        throw std::runtime_error("B-spline -> Unkown boundary condition given");
        break;
    }
}

void BSpline::solveSpline() {
    control_points_.resize(spline_size_ + 2, 2);
    if (boundary_condition_ == BoundaryCondition::ClosedLoop)
    {
        values_matrix_.resize(spline_size_ + 2, 2);
        values_matrix_.setZero();
        values_matrix_.block(0, 0, spline_size_ - 1, 2) = target_points_.block(0, 0, spline_size_ - 1, 2);

        control_points_ = spline_tridiagonal_matrix_.colPivHouseholderQr().solve(6 * values_matrix_);
    }
    else
    {
        Point first_point_derivative{};
        first_point_derivative << 0.0, 0.0;

        Point last_point_derivative{};
        first_point_derivative << 0.0, 0.0;

        values_matrix_.resize(spline_size_, 2);
        values_matrix_ = target_points_;
        values_matrix_.row(0) += first_point_derivative / (3 * (spline_size_ - 1));
        values_matrix_.row(spline_size_ - 1) -= last_point_derivative / (3 * (spline_size_ - 1));

        control_points_.block(1, 0, spline_size_, 2) = solveTridiagonal(spline_tridiagonal_matrix_, 6 * values_matrix_);
        control_points_.row(0) = control_points_.row(2) - 2 * first_point_derivative / (spline_size_ - 1);
        control_points_.row(spline_size_ + 1) = control_points_.row(spline_size_ - 1) - 2 * last_point_derivative / (spline_size_ - 1);
    }
}

void BSpline::checkClosedLoop() {
    if (target_points_.row(0) != target_points_.row(spline_size_ - 1))
    {
        throw std::runtime_error("B-spline -> Closed-Loop condition specified while first and last points do not have the same coordinates");
    }
}

LocalSplineIndex BSpline::getLocalIndex(double parameter) {
    if (parameter < 0 or parameter > 1)
    {
        throw std::runtime_error("B-spline -> Parameter value is out of range");
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

void BSpline::setTargetPoints(const PointsArray& targets, const BoundaryCondition condition) {
    target_points_ = targets;
    spline_size_ = target_points_.rows();
    if (condition != BoundaryCondition::BoundaryConditionSize) { boundary_condition_ = condition; }

    createTridiagonalMatrix();
    solveSpline();
}

// --- Getters ---

Point BSpline::getSecondDerivative(double parameter) {
    LocalSplineIndex index{ getLocalIndex(parameter) };
    int i{ index.spline_segment };
    double u{ index.local_parameter };
    Point derivative_value{};

    double b_m2{ u };
    double b_m1{ 1 - 3 * u };
    double b_0{ -2 + 3 * u };
    double b_1{ 1 - u };

    // Equation of second derivative as a function of the already calculated ones
    derivative_value = b_m2 * control_points_.row(i + 3) + b_m1 * control_points_.row(i + 2) \
        + b_0 * control_points_.row(i + 1) + b_1 * control_points_.row(i);

    return derivative_value;
}

Point BSpline::getFirstDerivative(double parameter) {
    LocalSplineIndex index{ getLocalIndex(parameter) };
    int i{ index.spline_segment };
    double u{ index.local_parameter };
    Point derivative_value{};

    double b_m2{ std::pow(u,2) / 2.0 };
    double b_m1{ (1 + 2 * u - 3 * std::pow(u,2)) / 2.0 };
    double b_0{ (-4 * u + 3 * std::pow(u,2)) / 2.0 };
    double b_1{ (-1 + 2 * u - std::pow(u,2)) / 2.0 };

    // Equation of first derivative as a function of second derivatives
    derivative_value = b_m2 * control_points_.row(i + 3) + b_m1 * control_points_.row(i + 2) \
        + b_0 * control_points_.row(i + 1) + b_1 * control_points_.row(i);

    return derivative_value;
}

Point BSpline::getPoint(double parameter) {
    LocalSplineIndex index{ getLocalIndex(parameter) };
    int i{ index.spline_segment };
    double u{ index.local_parameter };
    Point spline_point{};

    double b_m2{ std::pow(u,3) / 6.0 };
    double b_m1{ (1 + 3 * u + 3 * std::pow(u,2) - 3 * std::pow(u,3)) / 6.0 };
    double b_0{ (4 - 6 * std::pow(u,2) + 3 * std::pow(u,3)) / 6.0 };
    double b_1{ (1 - 3 * u + 3 * std::pow(u,2) - std::pow(u,3)) / 6.0 };

    // Equation of point coordinates as a function of second derivatives
    spline_point = b_m2 * control_points_.row(i + 3) + b_m1 * control_points_.row(i + 2) \
        + b_0 * control_points_.row(i + 1) + b_1 * control_points_.row(i);

    return spline_point;
}

PointsArray BSpline::getControlPoints() {
    return control_points_;
}

PointsArray BSpline::getSplineCurve(const long int resolution) {
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

double BSpline::getCurvature(double parameter) {
    Point first_der{ getFirstDerivative(parameter) };
    Point second_der{ getSecondDerivative(parameter) };
    double curvature{
        (first_der(0) * second_der(1) - first_der(1) * second_der(0)) / (std::pow(std::pow(first_der(0), 2) + std::pow(first_der(1), 2), 3.0 / 2.0))
    };
    return curvature;
}

double BSpline::getTangent(double parameter) {
    Point first_der{ getFirstDerivative(parameter) };
    return std::atan2(first_der(1), first_der(0));
}
} // namespace spline_interpolation


// Random stuff for testing
int main() {
    spline_interpolation::PointsArray midpoints{ readTrack("../point_paths/test_path.txt") };
    midpoints.conservativeResize(midpoints.rows() + 1, midpoints.cols());
    midpoints.row(midpoints.rows() - 1) = midpoints.row(0);

    spline_interpolation::BSpline spline{midpoints, spline_interpolation::BoundaryCondition::ClosedLoop};

    std::cout << "At s = 0.67" << '\n';
    std::cout << "Spline coordinates: " << spline.getPoint(0.67) << '\n';
    std::cout << "Spline curvature: " << spline.getCurvature(0.67) << '\n';
    std::cout << "Spline tangent: " << spline.getTangent(0.67) << "\n\n";

    std::cout << "Start: x' = " << spline.getFirstDerivative(0) << " , x'' = " << spline.getSecondDerivative(0) << '\n';
    std::cout << "End: x' = " << spline.getFirstDerivative(1) << " , x'' = " << spline.getSecondDerivative(1) << "\n\n";

    std::ofstream spline_curve("b_spline_points.txt");
    spline_curve << spline.getSplineCurve(2000);
    spline_curve.close();

    return 0;
}