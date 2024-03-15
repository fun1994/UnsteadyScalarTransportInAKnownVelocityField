#pragma once
#include "Data.h"

class UnsteadyScalarTransportInAKnownVelocityField {
	double rho;
	double Gamma;
	int Nx;
	int Ny;
	double dx;
	double dy;
	double dt;
	std::string temporal;
	double tol;
	double omega;
	void initialize(Data& data) {
		data.x.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			data.x[i] = (i + 0.5) * dx;
		}
		data.y.resize(Ny);
		for (int i = 0; i < Ny; i++) {
			data.y[i] = (i + 0.5) * dy;
		}
		data.t.push_back(0.0);
		data.phi.resize(1);
		data.phi[0].resize(Nx);
		for (int i = 0; i < Nx; i++) {
			data.phi[0][i].resize(Ny);
		}
	}
	void initialize(std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P, std::vector<std::vector<double>>& phi) {
		A_W.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_W[i].resize(Ny);
		}
		A_S.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_S[i].resize(Ny);
		}
		A_P.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_P[i].resize(Ny);
		}
		A_N.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_N[i].resize(Ny);
		}
		A_E.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_E[i].resize(Ny);
		}
		Q_P.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			Q_P[i].resize(Ny);
		}
		phi.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			phi[i].resize(Ny);
		}
	}
	void discretize(Data& data, std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P, bool flag) {
		if (temporal == "explicit Euler") {
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					A_P[i][j] += rho * dx * dy;
					Q_P[i][j] += rho * dx * dy * data.phi[data.phi.size() - 1][i][j];
					if (i == Nx - 1) {
						Q_P[i][j] += -rho * dy * dt * data.phi[data.phi.size() - 1][i][j];
					}
					else {
						Q_P[i][j] += -rho * (data.x[i] + data.x[i + 1]) * dy * dt / 4 * (data.phi[data.phi.size() - 1][i][j] + data.phi[data.phi.size() - 1][i + 1][j]);
					}
					if (i > 0) {
						Q_P[i][j] += rho * (data.x[i - 1] + data.x[i]) * dy * dt / 4 * (data.phi[data.phi.size() - 1][i - 1][j] + data.phi[data.phi.size() - 1][i][j]);
					}
					if (j < Ny - 1) {
						Q_P[i][j] += rho * (data.y[j] + data.y[j + 1]) * dx * dt / 4 * (data.phi[data.phi.size() - 1][i][j] + data.phi[data.phi.size() - 1][i][j + 1]);
					}
					if (j > 0) {
						Q_P[i][j] += -rho * (data.y[j - 1] + data.y[j]) * dx * dt / 4 * (data.phi[data.phi.size() - 1][i][j - 1] + data.phi[data.phi.size() - 1][i][j]);
					}
					if (i < Nx - 1) {
						Q_P[i][j] += Gamma * dy * dt / dx * (data.phi[data.phi.size() - 1][i + 1][j] - data.phi[data.phi.size() - 1][i][j]);
					}
					if (i == 0) {
						Q_P[i][j] += -2 * Gamma * dy * dt / dx * (data.phi[data.phi.size() - 1][i][j] - (1 - data.y[j]));
					}
					else {
						Q_P[i][j] += -Gamma * dy * dt / dx * (data.phi[data.phi.size() - 1][i][j] - data.phi[data.phi.size() - 1][i - 1][j]);
					}
					if (j == Ny - 1) {
						Q_P[i][j] += -2 * Gamma * dx * dt / dy * data.phi[data.phi.size() - 1][i][j];
					}
					else {
						Q_P[i][j] += Gamma * dx * dt / dy * (data.phi[data.phi.size() - 1][i][j + 1] - data.phi[data.phi.size() - 1][i][j]);
					}
					if (j > 0) {
						Q_P[i][j] += -Gamma * dx * dt / dy * (data.phi[data.phi.size() - 1][i][j] - data.phi[data.phi.size() - 1][i][j - 1]);
					}
				}
			}
		}
		else if (temporal == "implicit Euler" || (temporal == "three-time-level" && flag)) {
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					A_P[i][j] += rho * dx * dy;
					Q_P[i][j] += rho * dx * dy * data.phi[data.phi.size() - 1][i][j];
					if (i == Nx - 1) {
						A_P[i][j] += rho * dy * dt;
					}
					else {
						A_P[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy * dt / 4;
						A_E[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy * dt / 4;
					}
					if (i > 0) {
						A_W[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy * dt / 4;
						A_P[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy * dt / 4;
					}
					if (j < Ny - 1) {
						A_P[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx * dt / 4;
						A_N[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx * dt / 4;
					}
					if (j > 0) {
						A_S[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx * dt / 4;
						A_P[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx * dt / 4;
					}
					if (i < Nx - 1) {
						A_P[i][j] += Gamma * dy * dt / dx;
						A_E[i][j] += -Gamma * dy * dt / dx;
					}
					if (i == 0) {
						A_P[i][j] += 2 * Gamma * dy * dt / dx;
						Q_P[i][j] += 2 * Gamma * dy * dt / dx * (1 - data.y[j]);
					}
					else {
						A_W[i][j] += -Gamma * dy * dt / dx;
						A_P[i][j] += Gamma * dy * dt / dx;
					}
					if (j == Ny - 1) {
						A_P[i][j] += 2 * Gamma * dx * dt / dy;
					}
					else {
						A_P[i][j] += Gamma * dx * dt / dy;
						A_N[i][j] += -Gamma * dx * dt / dy;
					}
					if (j > 0) {
						A_S[i][j] += -Gamma * dx * dt / dy;
						A_P[i][j] += Gamma * dx * dt / dy;
					}
				}
			}
		}
		else if (temporal == "Crank-Nicolson") {
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					A_P[i][j] += rho * dx * dy;
					Q_P[i][j] += rho * dx * dy * data.phi[data.phi.size() - 1][i][j];
					if (i == Nx - 1) {
						A_P[i][j] += rho * dy * dt / 2;
						Q_P[i][j] += -rho * dy * dt / 2 * data.phi[data.phi.size() - 1][i][j];
					}
					else {
						A_P[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy * dt / 8;
						A_E[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy * dt / 8;
						Q_P[i][j] += -rho * (data.x[i] + data.x[i + 1]) * dy * dt / 8 * (data.phi[data.phi.size() - 1][i][j] + data.phi[data.phi.size() - 1][i + 1][j]);
					}
					if (i > 0) {
						A_W[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy * dt / 8;
						A_P[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy * dt / 8;
						Q_P[i][j] += rho * (data.x[i - 1] + data.x[i]) * dy * dt / 8 * (data.phi[data.phi.size() - 1][i - 1][j] + data.phi[data.phi.size() - 1][i][j]);
					}
					if (j < Ny - 1) {
						A_P[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx * dt / 8;
						A_N[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx * dt / 8;
						Q_P[i][j] += rho * (data.y[j] + data.y[j + 1]) * dx * dt / 8 * (data.phi[data.phi.size() - 1][i][j] + data.phi[data.phi.size() - 1][i][j + 1]);
					}
					if (j > 0) {
						A_S[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx * dt / 8;
						A_P[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx * dt / 8;
						Q_P[i][j] += -rho * (data.y[j - 1] + data.y[j]) * dx * dt / 8 * (data.phi[data.phi.size() - 1][i][j - 1] + data.phi[data.phi.size() - 1][i][j]);
					}
					if (i < Nx - 1) {
						A_P[i][j] += Gamma * dy * dt / (2 * dx);
						A_E[i][j] += -Gamma * dy * dt / (2 * dx);
						Q_P[i][j] += Gamma * dy * dt / (2 * dx) * (data.phi[data.phi.size() - 1][i + 1][j] - data.phi[data.phi.size() - 1][i][j]);
					}
					if (i == 0) {
						A_P[i][j] += Gamma * dy * dt / dx;
						Q_P[i][j] += Gamma * dy * dt / dx * (1 - data.y[j]);
						Q_P[i][j] += -Gamma * dy * dt / dx * (data.phi[data.phi.size() - 1][i][j] - (1 - data.y[j]));
					}
					else {
						A_W[i][j] += -Gamma * dy * dt / (2 * dx);
						A_P[i][j] += Gamma * dy * dt / (2 * dx);
						Q_P[i][j] += -Gamma * dy * dt / (2 * dx) * (data.phi[data.phi.size() - 1][i][j] - data.phi[data.phi.size() - 1][i - 1][j]);
					}
					if (j == Ny - 1) {
						A_P[i][j] += Gamma * dx * dt / dy;
						Q_P[i][j] += -Gamma * dx * dt / dy * data.phi[data.phi.size() - 1][i][j];
					}
					else {
						A_P[i][j] += Gamma * dx * dt / (2 * dy);
						A_N[i][j] += -Gamma * dx * dt / (2 * dy);
						Q_P[i][j] += Gamma * dx * dt / (2 * dy) * (data.phi[data.phi.size() - 1][i][j + 1] - data.phi[data.phi.size() - 1][i][j]);
					}
					if (j > 0) {
						A_S[i][j] += -Gamma * dx * dt / (2 * dy);
						A_P[i][j] += Gamma * dx * dt / (2 * dy);
						Q_P[i][j] += -Gamma * dx * dt / (2 * dy) * (data.phi[data.phi.size() - 1][i][j] - data.phi[data.phi.size() - 1][i][j - 1]);
					}
				}
			}
		}
		else if (temporal == "three-time-level" && !flag) {
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					A_P[i][j] += rho * dx * dy;
					Q_P[i][j] += 4 * rho * dx * dy / 3 * data.phi[data.phi.size() - 1][i][j] - rho * dx * dy / 3 * data.phi[data.phi.size() - 2][i][j];
					if (i == Nx - 1) {
						A_P[i][j] += 2 * rho * dy * dt / 3;
					}
					else {
						A_P[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy * dt / 6;
						A_E[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy * dt / 6;
					}
					if (i > 0) {
						A_W[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy * dt / 6;
						A_P[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy * dt / 6;
					}
					if (j < Ny - 1) {
						A_P[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx * dt / 6;
						A_N[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx * dt / 6;
					}
					if (j > 0) {
						A_S[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx * dt / 6;
						A_P[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx * dt / 6;
					}
					if (i < Nx - 1) {
						A_P[i][j] += 2 * Gamma * dy * dt / (3 * dx);
						A_E[i][j] += -2 * Gamma * dy * dt / (3 * dx);
					}
					if (i == 0) {
						A_P[i][j] += 4 * Gamma * dy * dt / (3 * dx);
						Q_P[i][j] += 4 * Gamma * dy * dt / (3 * dx) * (1 - data.y[j]);
					}
					else {
						A_W[i][j] += -2 * Gamma * dy * dt / (3 * dx);
						A_P[i][j] += 2 * Gamma * dy * dt / (3 * dx);
					}
					if (j == Ny - 1) {
						A_P[i][j] += 4 * Gamma * dx * dt / (3 * dy);
					}
					else {
						A_P[i][j] += 2 * Gamma * dx * dt / (3 * dy);
						A_N[i][j] += -2 * Gamma * dx * dt / (3 * dy);
					}
					if (j > 0) {
						A_S[i][j] += -2 * Gamma * dx * dt / (3 * dy);
						A_P[i][j] += 2 * Gamma * dx * dt / (3 * dy);
					}
				}
			}
		}
	}
	void SOR(Data& data, std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P, std::vector<std::vector<double>>& phi) {
		while (true) {
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					phi[i][j] = Q_P[i][j] - A_P[i][j] * data.phi[data.phi.size() - 1][i][j];
					if (i < Nx - 1) {
						phi[i][j] -= A_E[i][j] * data.phi[data.phi.size() - 1][i + 1][j];
					}
					if (i > 0) {
						phi[i][j] -= A_W[i][j] * phi[i - 1][j];
					}
					if (j < Ny - 1) {
						phi[i][j] -= A_N[i][j] * data.phi[data.phi.size() - 1][i][j + 1];
					}
					if (j > 0) {
						phi[i][j] -= A_S[i][j] * phi[i][j - 1];
					}
					phi[i][j] *= omega / A_P[i][j];
					phi[i][j] += data.phi[data.phi.size() - 1][i][j];
				}
			}
			double res = residual(data, phi);
			data.phi[data.phi.size() - 1] = phi;
			if (res < tol) {
				break;
			}
		}
	}
	double residual(Data& data, std::vector<std::vector<double>>& phi) {
		double res = 0.0;
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				res = res < abs(data.phi[data.phi.size() - 1][i][j] - phi[i][j]) ? abs(data.phi[data.phi.size() - 1][i][j] - phi[i][j]) : res;
			}
		}
		return res;
	}
	double residual(Data& data) {
		double res = 0.0;
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				res = res < abs(data.phi[data.phi.size() - 1][i][j] - data.phi[data.phi.size() - 2][i][j]) ? abs(data.phi[data.phi.size() - 1][i][j] - data.phi[data.phi.size() - 2][i][j]) : res;
			}
		}
		return res;
	}
public:
	UnsteadyScalarTransportInAKnownVelocityField(double rho, double Gamma, int Nx, int Ny, double dt, std::string temporal, double tol, double omega) :rho(rho), Gamma(Gamma), Nx(Nx), Ny(Ny), dx(1.0 / Nx), dy(1.0 / Ny), dt(dt), temporal(temporal), tol(tol), omega(omega) {}
	void solve(Data& data) {
		initialize(data);
		bool flag = true;
		while (true) {
			std::vector<std::vector<double>> A_W, A_S, A_P, A_N, A_E, Q_P, phi;
			initialize(A_W, A_S, A_P, A_N, A_E, Q_P, phi);
			discretize(data, A_W, A_S, A_P, A_N, A_E, Q_P, flag);
			flag = false;
			data.t.push_back(data.t[data.t.size() - 1] + dt);
			data.phi.push_back(data.phi[data.phi.size() - 1]);
			SOR(data, A_W, A_S, A_P, A_N, A_E, Q_P, phi);
			double res = residual(data);
			if (res < tol) {
				break;
			}
		}
	}
};
