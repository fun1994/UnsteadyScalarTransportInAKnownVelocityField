#include "UnsteadyScalarTransportInAKnownVelocityField.h"

void test(std::string temporal, std::string index) {
	UnsteadyScalarTransportInAKnownVelocityField USTIAKVF(1.2, 0.1, 20, 20, 0.0001, temporal, 1e-8, 1.0);
	Data data;
	USTIAKVF.solve(data);
	std::cout << "temporal=" << temporal << " time=" << data.t[data.t.size() - 1] << std::endl;
	data.save(index);
}

void test() {
	test("explicit Euler", "1");
	test("implicit Euler", "2");
	test("Crank-Nicolson", "3");
	test("three-time-level", "4");
}

int main() {
	test();
	return 0;
}
