#include "analysis_exercise4.h" 

int main() {
    Helix h(1.37, 1.12); // Create a helix with R = 1.37 and H = 1.12

    double L = h.arclength(-0.3, 2.2); // Calculate arclength between a = -0.3 and b = 2.2
    std::cout << "Arclength of the helix: " << L << std::endl; // Should now reflect the new N_interval_trapezoidal

    h.output("result.dat", -0.3, 2.2, 100);

    return 0;
}
