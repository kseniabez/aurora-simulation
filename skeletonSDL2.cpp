#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include "SDL2Auxiliary.h"

using namespace std;
using glm::vec3;
// ---------------------------------------------------------
// STRUCTURES

struct FluidCube {
	int size;
	float dt;
	float diff;
	float visc;

	std::vector<float> s;
	std::vector<float> density;

	std::vector<float> Vx;
	std::vector<float> Vy;
	std::vector<float> Vz;

	std::vector<float> Vx0;
	std::vector<float> Vy0;
	std::vector<float> Vz0;

	FluidCube(int N, float diffusion, float viscosity, float deltaTime)
		: size(N), diff(diffusion), visc(viscosity), dt(deltaTime) {

		// Resizing the vectors
		s.resize(size * size * size);
		density.resize(size * size * size);
		Vx.resize(size * size * size);
		Vy.resize(size * size * size);
		Vz.resize(size * size * size);
		Vx0.resize(size * size * size);
		Vy0.resize(size * size * size);
		Vz0.resize(size * size * size);

		// Debugging outputs to check vector sizes
		std::cout << "FluidCube initialized with size " << size << "x" << size << "x" << size << std::endl;
		std::cout << "Vectors resized to size: " << s.size() << " (s), "
			<< density.size() << " (density), "
			<< Vx.size() << " (Vx), "
			<< Vy.size() << " (Vy), "
			<< Vz.size() << " (Vz), "
			<< Vx0.size() << " (Vx0), "
			<< Vy0.size() << " (Vy0), "
			<< Vz0.size() << " (Vz0)" << std::endl;
	}
};

// ---------------------------------------------------------
// GLOBAL VARIABLES
const int SCREEN_WIDTH = 100;
const int SCREEN_HEIGHT = 100;
SDL2Aux* sdlAux;
int t;
int N = 100;
FluidCube* cube;

#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)

// ---------------------------------------------------------
// FUNCTION DECLARATIONS
void Draw();
void Update();
void FluidCubeAddDensity(FluidCube* cube, int x, int y, int z, float amount);
void set_bnd(int b, vector<float>& x);
void lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter);
void diffuse(int b, vector<float>& x, vector<float>& x0, float diff, float dt, int iter);
void project(vector<float>& velocX, vector<float>& velocY, vector<float>& velocZ, vector<float>& p, vector<float>& div, int iter);
static void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY, vector<float>& velocZ, float dt);
void FluidCubeAddVelocity(FluidCube* cube, int x, int y, int z, float amountX, float amountY, float amountZ);
void FluidCubeStep(FluidCube* cube);

// ---------------------------------------------------------
// FUNCTION DEFINITIONS
int main(int argc, char* argv[])
{
	t = SDL_GetTicks();

	sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);

	cube = new FluidCube(N, 0.0001f, 0.0001f, 0.1f);

	if (cube == nullptr) {
		std::cerr << "Failed to allocate FluidCube!" << std::endl;
	}
	else {
		std::cout << "FluidCube successfully created!" << std::endl;
	}


	while (!sdlAux->quitEvent()) {
		Draw();
		Update();
	}
	sdlAux->saveBMP("screenshot.bmp");
	return 0;
}

void Update()
{
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t) / 1000.0f;
	t = t2;  
	std::cout << "Render time: " << dt << " ms." << endl;

	FluidCubeAddDensity(cube, N / 2, N / 2, N / 2, 1.0f);
	FluidCubeAddVelocity(cube, N / 2, N / 2, N / 2, 1.0f, 0.0f, 0.0f); 

	FluidCubeStep(cube);
}

void Draw()
{
	sdlAux->clearPixels();
	for (int y = 0; y < SCREEN_HEIGHT; ++y) {
		for (int x = 0; x < SCREEN_WIDTH; ++x) {

			int index = IX(x, y, 50);  

			float densityVal = cube->density[index];

			sdlAux->putPixel(x, y, densityVal * vec3(0.0f, 1.0f, 0.0f));
		}
	}
	sdlAux->render();
}
void FluidCubeAddDensity(FluidCube* cube, int x, int y, int z, float amount)
{
	cube->density[IX(x, y, z)] += amount;
}
void FluidCubeAddVelocity(FluidCube* cube, int x, int y, int z, float amountX, float amountY, float amountZ)
{

	int index = IX(x, y, z);

	cube->Vx[index] += amountX;
	cube->Vy[index] += amountY;
	cube->Vz[index] += amountZ;
}

void FluidCubeStep(FluidCube* cube)
{
	int N = cube->size;
	float visc = cube->visc;
	float diff = cube->diff;
	float dt = cube->dt;

	diffuse(1, cube->Vx0, cube->Vx, visc, dt, 4);
	diffuse(2, cube->Vy0, cube->Vy, visc, dt, 4);
	diffuse(3, cube->Vz0, cube->Vz, visc, dt, 4);

	project(cube->Vx0, cube->Vy0, cube->Vz0, cube->Vx, cube->Vy, 4);

	advect(1, cube->Vx, cube->Vx0, cube->Vx0, cube->Vy0, cube->Vz0, dt);
	advect(2, cube->Vy, cube->Vy0, cube->Vx0, cube->Vy0, cube->Vz0, dt);
	advect(3, cube->Vz, cube->Vz0, cube->Vx0, cube->Vy0, cube->Vz0, dt);

	project(cube->Vx, cube->Vy, cube->Vz, cube->Vx0, cube->Vy0, 4);

	diffuse(0, cube->s, cube->density, diff, dt, 4);
	advect(0, cube->density, cube->s, cube->Vx, cube->Vy, cube->Vz, dt);
}

void set_bnd(int b, vector<float>& x)
{
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
			x[IX(i, j, N - 1)] = b == 3 ? -x[IX(i, j, N - 2)] : x[IX(i, j, N - 2)];
		}
	}
	for (int k = 1; k < N - 1; k++) {
		for (int i = 1; i < N - 1; i++) {
			x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
			x[IX(i, N - 1, k)] = b == 2 ? -x[IX(i, N - 2, k)] : x[IX(i, N - 2, k)];
		}
	}
	for (int k = 1; k < N - 1; k++) {
		for (int j = 1; j < N - 1; j++) {
			x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
			x[IX(N - 1, j, k)] = b == 1 ? -x[IX(N - 2, j, k)] : x[IX(N - 2, j, k)];
		}
	}

	x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)]
		+ x[IX(0, 1, 0)]
		+ x[IX(0, 0, 1)]);

	x[IX(0, N - 1, 0)] = 0.33f * (x[IX(1, N - 1, 0)]
		+ x[IX(0, N - 2, 0)]
		+ x[IX(0, N - 1, 1)]);

	x[IX(0, 0, N - 1)] = 0.33f * (x[IX(1, 0, N - 1)]
		+ x[IX(0, 1, N - 1)]
		+ x[IX(0, 0, N - 1)]);

	x[IX(0, N - 1, N - 1)] = 0.33f * (x[IX(1, N - 1, N - 1)]
		+ x[IX(0, N - 2, N - 1)]
		+ x[IX(0, N - 1, N - 2)]);
	x[IX(N - 1, 0, 0)] = 0.33f * (x[IX(N - 2, 0, 0)]
		+ x[IX(N - 1, 1, 0)]
		+ x[IX(N - 1, 0, 1)]);
	x[IX(N - 1, N - 1, 0)] = 0.33f * (x[IX(N - 2, N - 1, 0)]
		+ x[IX(N - 1, N - 2, 0)]
		+ x[IX(N - 1, N - 1, 1)]);
	x[IX(N - 1, 0, N - 1)] = 0.33f * (x[IX(N - 2, 0, N - 1)]
		+ x[IX(N - 1, 1, N - 1)]
		+ x[IX(N - 1, 0, N - 2)]);
	x[IX(N - 1, N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1, N - 1)]
		+ x[IX(N - 1, N - 2, N - 1)]
		+ x[IX(N - 1, N - 1, N - 2)]);
}
void lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter)
{
	float cRecip = 1.0 / c;
	for (int k = 0; k < iter; k++) {
		for (int m = 1; m < N - 1; m++) {
			for (int j = 1; j < N - 1; j++) {
				for (int i = 1; i < N - 1; i++) {
					x[IX(i, j, m)] =
						(x0[IX(i, j, m)]
							+ a * (x[IX(i + 1, j, m)]
								+ x[IX(i - 1, j, m)]
								+ x[IX(i, j + 1, m)]
								+ x[IX(i, j - 1, m)]
								+ x[IX(i, j, m + 1)]
								+ x[IX(i, j, m - 1)]
								)) * cRecip;
				}
			}
		}
		set_bnd(b, x);
	}
}
void diffuse(int b, vector<float>& x, vector<float>& x0, float diff, float dt, int iter)
{
	float a = dt * diff * (N - 2) * (N - 2);
	lin_solve(b, x, x0, a, 1 + 6 * a, iter);
}
static void project(vector<float>& velocX, vector<float>& velocY, vector<float>& velocZ, vector<float>& p, vector<float>& div, int iter)
{
	for (int k = 1; k < N - 1; k++) {
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
				div[IX(i, j, k)] = -0.5f * (
					velocX[IX(i + 1, j, k)]
					- velocX[IX(i - 1, j, k)]
					+ velocY[IX(i, j + 1, k)]
					- velocY[IX(i, j - 1, k)]
					+ velocZ[IX(i, j, k + 1)]
					- velocZ[IX(i, j, k - 1)]
					) / N;
				p[IX(i, j, k)] = 0;
			}
		}
	}
	set_bnd(0, div);
	set_bnd(0, p);
	lin_solve(0, p, div, 1, 6, iter);

	for (int k = 1; k < N - 1; k++) {
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
				velocX[IX(i, j, k)] -= 0.5f * (p[IX(i + 1, j, k)]
					- p[IX(i - 1, j, k)]) * N;
				velocY[IX(i, j, k)] -= 0.5f * (p[IX(i, j + 1, k)]
					- p[IX(i, j - 1, k)]) * N;
				velocZ[IX(i, j, k)] -= 0.5f * (p[IX(i, j, k + 1)]
					- p[IX(i, j, k - 1)]) * N;
			}
		}
	}
	set_bnd(1, velocX);
	set_bnd(2, velocY);
	set_bnd(3, velocZ);
}
static void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY, vector<float>& velocZ, float dt)
{
	float i0, i1, j0, j1, k0, k1;

	float dtx = dt * (N - 2);
	float dty = dt * (N - 2);
	float dtz = dt * (N - 2);

	float s0, s1, t0, t1, u0, u1;
	float tmp1, tmp2, tmp3, x, y, z;

	float Nfloat = N;
	float ifloat, jfloat, kfloat;
	int i, j, k;

	for (k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
		for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
			for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
				tmp1 = dtx * velocX[IX(i, j, k)];
				tmp2 = dty * velocY[IX(i, j, k)];
				tmp3 = dtz * velocZ[IX(i, j, k)];
				x = ifloat - tmp1;
				y = jfloat - tmp2;
				z = kfloat - tmp3;

				if (x < 0.5f) x = 0.5f;
				if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
				i0 = floorf(x);
				i1 = i0 + 1.0f;
				if (y < 0.5f) y = 0.5f;
				if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
				j0 = floorf(y);
				j1 = j0 + 1.0f;
				if (z < 0.5f) z = 0.5f;
				if (z > Nfloat + 0.5f) z = Nfloat + 0.5f;
				k0 = floorf(z);
				k1 = k0 + 1.0f;

				s1 = x - i0;
				s0 = 1.0f - s1;
				t1 = y - j0;
				t0 = 1.0f - t1;
				u1 = z - k0;
				u0 = 1.0f - u1;

				int i0i = i0;
				int i1i = i1;
				int j0i = j0;
				int j1i = j1;
				int k0i = k0;
				int k1i = k1;

				d[IX(i, j, k)] =

					s0 * (t0 * (u0 * d0[IX(i0i, j0i, k0i)]
						+ u1 * d0[IX(i0i, j0i, k1i)])
						+ (t1 * (u0 * d0[IX(i0i, j1i, k0i)]
							+ u1 * d0[IX(i0i, j1i, k1i)])))
					+ s1 * (t0 * (u0 * d0[IX(i1i, j0i, k0i)]
						+ u1 * d0[IX(i1i, j0i, k1i)])
						+ (t1 * (u0 * d0[IX(i1i, j1i, k0i)]
							+ u1 * d0[IX(i1i, j1i, k1i)])));
			}
		}
	}
	set_bnd(b, d);
}