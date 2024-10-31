#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

class Line3D
{
public:
    /// Constructor
    /// Initialise the number of intervals and step size for finite difference
    Line3D() : N_interval_trapezoidal(100), FD_step_for_tangent(1.0e-8) {};

    /// Position vector as a function of the curve parameter
    /// sigma: r(sigma)
    /// Pure virtual function
    virtual std::vector<double> position(const double& sigma) const = 0;

    /// Tangent vector as a function of the curve parameter
    /// sigma: dr/dsigma
    /// Virtual
    /// Can be numerically approximated
    virtual std::vector<double> tangent(const double& sigma) const
    {
        /// Calculate position + step and position - step
        std::vector<double> pos_for = position(sigma + FD_step_for_tangent);
        std::vector<double> pos_back = position(sigma - FD_step_for_tangent);

        // Create a vector for x, y and z
        std::vector<double> tangent(3);

        /// Loop to approximate each component
        for (unsigned i = 0; i < 3; i++)
        {
            tangent[i] = (pos_for[i] - pos_back[i]) / (2 * FD_step_for_tangent);
        }

        return tangent;
    };

    /// Output n_point points along the curve,
    /// at equally spaced values of sigma between sigma=a and sigma=b.
    /// filename specifies the name of the file that will contain n_point lines,
    /// each containing the three coordinates (x, y, z) of the point.
    void output(const std::string& filename, const double& a, const double& b, const unsigned& n_point) const
    {
        /// Open the file
        std::ofstream file(filename);

        /// Throw an error in case file could not be opened
        if (!file.is_open())
        {
            std::cerr << "Error: could not open the file" << std::endl;
            return;
        }

        /// Calculate the step size
        double delta_sigma = (b - a) / (n_point - 1);
        
        /// Loop to compute positions along the curve
        for (unsigned i = 0; i < n_point; ++i) {
            
            /// Calculate current value of sigma
            double sigma = a + i * delta_sigma;
            
            /// Calculate the position for the current sigma
            std::vector<double> pos = position(sigma);

            /// Write coordinates to the file
            file << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }

        /// Close the file
        file.close();
    }

    /// Virtual
    /// Can be numerically approximated
    /// Return the arclength along the curve between sigma=a and sigma=b.
    virtual double arclength(const double& a, const double& b) const
    {
        // Step size h definition  
        double h = (b - a) / N_interval_trapezoidal;

        // Calculation of the contribution of the first and last pointsl
        double integral = 0.5 * tangent_norm(a);
        integral += 0.5 * tangent_norm(b);

        // Perform trapezoidal rule iteration step
        // to calculate the interior points and the sum of the integral
        for (unsigned i = 1; i < N_interval_trapezoidal; ++i)
        {
            // Calculate the x value
            double sigma = a + i * h;

            // Sum up the functions values
            integral += tangent_norm(sigma);
        }

        /// Return the arclenght
        return h * integral;
    }

protected:
    unsigned N_interval_trapezoidal;
    double FD_step_for_tangent;

    /// Initialize auxilary function to calculate tangent norm 
    double tangent_norm(const double& sigma) const
    {
        /// Create a vector for x, y and z
        /// Calls tangent function using 'this' pointer as used within the class 
        std::vector<double> t = this->tangent(sigma);

        /// Return normalized tangent
        return std::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    }
};

class Helix : public Line3D
{
public:
    // Constructor for helix with given radius and pitch
    Helix(const double& radius, const double& h) : R(radius), H(h) {}

    // Override position function for helix
    std::vector<double> position(const double& sigma) const override
    {
        /// Create a vector for x, y and z
        std::vector<double> pos(3);
        pos[0] = R * std::cos(sigma); /// x
        pos[1] = R * std::sin(sigma); /// y
        pos[2] = H * sigma; /// z
        return pos;
    }

private:
    double R;
    double H;
};