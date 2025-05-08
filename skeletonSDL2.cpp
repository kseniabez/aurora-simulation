#include <iostream>
#include <sstream>
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

	std::vector<float> Vx0;
	std::vector<float> Vy0;

	FluidCube(int N, float diffusion, float viscosity, float deltaTime)
		: size(N), diff(diffusion), visc(viscosity), dt(deltaTime) {

		// Resizing the vectors
		s.resize(size * size);
		density.resize(size * size);
		Vx.resize(size * size);
		Vy.resize(size * size);
		Vx0.resize(size * size);
		Vy0.resize(size * size);

		// Debugging outputs to check vector sizes
		std::cout << "FluidCube initialized with size " << size << "x" << size << "x" << size << std::endl;
		std::cout << "Vectors resized to size: " << s.size() << " (s), "
			<< density.size() << " (density), "
			<< Vx.size() << " (Vx), "
			<< Vy.size() << " (Vy), "
			<< Vx0.size() << " (Vx0), "
			<< Vy0.size() << " (Vy0), " << std::endl;
	}
};

// ---------------------------------------------------------
// GLOBAL VARIABLES
const int SCREEN_WIDTH = 200;
const int SCREEN_HEIGHT = 200;
SDL2Aux* sdlAux;
int t;
int N = 100;
FluidCube* cube;

int frameCount = 0;
vector<bool> is_occupied;
#define IX(x, y) ((x) + (y) * N)

int stageDuration = 100;

// ---------------------------------------------------------
// FUNCTION DECLARATIONS
void Draw();
void Update();
void FluidCubeAddDensity(FluidCube* cube, int x, int y, float amount);
void set_bnd(int b, vector<float>& x);
void lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter);
void diffuse(int b, vector<float>& x, vector<float>& x0, float diff, float dt, int iter);
void project(vector<float>& velocX, vector<float>& velocY, vector<float>& p, vector<float>& div, int iter);
static void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY, float dt);
void FluidCubeAddVelocity(FluidCube* cube, int x, int y, float amountX, float amountY);
void FluidCubeStep(FluidCube* cube);
void AuroraPhase(FluidCube* cube, int phase, float centerX, float centerY, float radius, float falloff, float strength);

// ---------------------------------------------------------
// FUNCTION DEFINITIONS
int main(int argc, char* argv[])
{
	t = SDL_GetTicks();

	sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);

	cube = new FluidCube(N, 0.000001f, 0.0001f, 0.1f);

	if (cube == nullptr) {
		std::cerr << "Failed to allocate FluidCube!" << std::endl;
	}
	else {
		std::cout << "FluidCube successfully created!" << std::endl;
	}

	is_occupied = vector<bool>(N * N, false);

	while (!sdlAux->quitEvent()) {
		Draw();
		Update();

		std::ostringstream oss;
		oss << frameCount++;
		std::string str = oss.str();

		std::string filename = "video/screenshot_" + str + ".bmp";
		sdlAux->saveBMP(filename.c_str());
	}
	sdlAux->saveBMP("screenshot.bmp");
	return 0;
}

void Update()
{
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t) / 1000.0f;
	t = t2;
	std::cout << "Render time: " << dt * 1000.0f << " ms." << std::endl;

	int phase = (frameCount / stageDuration) % 4;

	for (int j = 0; j < N && frameCount == 0; ++j) {
		for (int i = N / 2 - 1; i < N / 2 + 1; i++)
		{
			FluidCubeAddDensity(cube, j, i, 2.0f);
		}
	}

	AuroraPhase(cube, phase, N / 2, N / 2, 10.0f, 20.0f, 0.002f);
	FluidCubeStep(cube);
	frameCount++;
}

void AuroraPhase(FluidCube* cube, int phase, float centerX, float centerY, float radius, float falloff, float strength)
{
	for (int j = 1; j < N - 1; ++j) {
		for (int i = 1; i < N - 1; ++i) {
			float dx = i - centerX;
			float dy = j - centerY;
			float dist = sqrtf(dx * dx + dy * dy);

			float weight = 0.0f;
			if (dist < radius)
				weight = 1.0f;
			else if (dist < radius + falloff)
				weight = 1.0f - (dist - radius) / falloff;

			if (weight <= 0.0f) continue;

			int col = frameCount % (stageDuration * 4);

		switch (phase) {
			case 0:
				/*for (int k = N / 2 - 1; k < N / 2 + 1; k++)
				{
					FluidCubeAddDensity(cube, col, k, 0.1f);
				}*/
				break;

			case 1: {
				// Simple vertical oscillation (init wiggle)
				float amplitude = strength * 1.0f;
				float frequency = 2.0f * M_PI / N;
				float omega = 0.1f;
				float vx = 0.0f;
				float vy = amplitude * sinf(frequency * i - omega * frameCount) * weight;
				FluidCubeAddVelocity(cube, i, j, vx, vy);
				break;
			}

			case 2: {
				// Twisting with radial pull
				float swirl_strength = strength * expf(-frameCount / 400.0f);
				float radial_pull = strength * 0.2f;
				float dist_safe = dist + 1e-5f;

				float velX = (-dy / dist_safe * swirl_strength - dx / dist_safe * radial_pull) * weight;
				float velY = (+dx / dist_safe * swirl_strength - dy / dist_safe * radial_pull) * weight;

				FluidCubeAddVelocity(cube, i, j, velX, velY);
				break;
			}

			case 3: {
				// Full spiral vortex
				float spiral_strength = strength / 2 * expf(-frameCount / 500.0f);
				float dist_safe = dist + 1e-5f;

				float angle = atan2f(dy, dx);
				float velX = -sinf(angle) * spiral_strength * weight;
				float velY = cosf(angle) * spiral_strength * weight;

				FluidCubeAddVelocity(cube, i, j, velX, velY);
				break;
			}
			}
		}
	}
}


void Draw()
{
	sdlAux->clearPixels();
	for (int y = 0; y < SCREEN_HEIGHT; ++y) {
		for (int x = 0; x < SCREEN_WIDTH; ++x) {
			int indX = x * N / SCREEN_WIDTH;
			int indY = y * N / SCREEN_HEIGHT;
			int index = IX(indX, indY);  

			float densityVal = cube->density[index];

			sdlAux->putPixel(x, y, densityVal * vec3(0.0f, 1.0f, 0.0f));
		}
	}
	sdlAux->render();
}
void FluidCubeAddDensity(FluidCube* cube, int x, int y, float amount)
{
	cube->density[IX(x, y)] += amount;
}
void FluidCubeAddVelocity(FluidCube* cube, int x, int y, float amountX, float amountY)
{

	int index = IX(x, y);

	cube->Vx[index] += amountX;
	cube->Vy[index] += amountY;
}

void FluidCubeStep(FluidCube* cube)
{
	int N = cube->size;
	float visc = cube->visc;
	float diff = cube->diff;
	float dt = cube->dt;

	diffuse(1, cube->Vx0, cube->Vx, visc, dt, 4);
	diffuse(2, cube->Vy0, cube->Vy, visc, dt, 4);

	project(cube->Vx0, cube->Vy0, cube->Vx, cube->Vy, 4);

	advect(1, cube->Vx, cube->Vx0, cube->Vx0, cube->Vy0, dt);
	advect(2, cube->Vy, cube->Vy0, cube->Vx0, cube->Vy0, dt);

	project(cube->Vx, cube->Vy, cube->Vx0, cube->Vy0, 4);

	diffuse(0, cube->s, cube->density, diff, dt, 4);
	advect(0, cube->density, cube->s, cube->Vx, cube->Vy, dt);
}

void set_bnd(int b, vector<float>& x)
{
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			int idx = IX(i, j);
			if (is_occupied[idx]) {
				// Average values of non-occupied neighbors
				float sum = 0.0f;
				int count = 0;
				if (!is_occupied[IX(i - 1, j)]) { sum += x[IX(i - 1, j)]; count++; }
				if (!is_occupied[IX(i + 1, j)]) { sum += x[IX(i + 1, j)]; count++; }
				if (!is_occupied[IX(i, j - 1)]) { sum += x[IX(i, j - 1)]; count++; }
				if (!is_occupied[IX(i, j + 1)]) { sum += x[IX(i, j + 1)]; count++; }
				if (count > 0)
					x[idx] = sum / count;
				else
					x[idx] = 0.0f;
			}
		}
	}

	// Original boundary code for outer edges
	for (int i = 1; i < N - 1; i++) {
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N - 1, i)] = b == 1 ? -x[IX(N - 2, i)] : x[IX(N - 2, i)];
	}

	// Corner handling (you can keep as-is or improve)
	x[IX(0, 0)] = 0.33f * (x[IX(1, 0)] + x[IX(0, 1)] + x[IX(0, 0)]);
	x[IX(0, N - 1)] = 0.33f * (x[IX(1, N - 1)] + x[IX(0, N - 2)] + x[IX(0, N - 1)]);
	x[IX(N - 1, 0)] = 0.33f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)] + x[IX(N - 1, 0)]);
	x[IX(N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)] + x[IX(N - 1, N - 1)]);
}
void lin_solve(int b, vector<float>& x, vector<float>& x0, float a, float c, int iter)
{
	float cRecip = 1.0 / c;
	for (int k = 0; k < iter; k++) {
			for (int j = 1; j < N - 1; j++) {
				for (int i = 1; i < N - 1; i++) {
					x[IX(i, j, m)] =
						(x0[IX(i, j)]
							+ a * (x[IX(i + 1, j)]
								+ x[IX(i - 1, j)]
								+ x[IX(i, j + 1)]
								+ x[IX(i, j - 1)]
								+ x[IX(i, j)]
								+ x[IX(i, j)]
								)) * cRecip;
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
static void project(vector<float>& velocX, vector<float>& velocY, vector<float>& p, vector<float>& div, int iter)
{
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
				div[IX(i, j, k)] = -0.5f * (
					velocX[IX(i + 1, j)]
					- velocX[IX(i - 1, j)]
					+ velocY[IX(i, j + 1)]
					- velocY[IX(i, j - 1)]
					) / N;
				p[IX(i, j, k)] = 0;
			}
		}
	set_bnd(0, div);
	set_bnd(0, p);
	lin_solve(0, p, div, 1, 6, iter);

	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
				- p[IX(i - 1, j)]) * N;
			velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
				- p[IX(i, j - 1)]) * N;
		}
	}
	set_bnd(1, velocX);
	set_bnd(2, velocY);
}
static void advect(int b, vector<float>& d, vector<float>& d0, vector<float>& velocX, vector<float>& velocY, float dt)
{
	float i0, i1, j0, j1;

	float dtx = dt * (N - 2);
	float dty = dt * (N - 2);

	float s0, s1, t0, t1;
	float tmp1, tmp2, x, y;

	float Nfloat = N;
	float ifloat, jfloat;
	int i, j;

		for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
			for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
				tmp1 = dtx * velocX[IX(i, j)];
				tmp2 = dty * velocY[IX(i, j)];
				x = ifloat - tmp1;
				y = jfloat - tmp2;

				if (x < 0.5f) x = 0.5f;
				if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
				i0 = floorf(x);
				i1 = i0 + 1.0f;
				if (y < 0.5f) y = 0.5f;
				if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
				j0 = floorf(y);
				j1 = j0 + 1.0f;

				s1 = x - i0;
				s0 = 1.0f - s1;
				t1 = y - j0;
				t0 = 1.0f - t1;

				int i0i = i0;
				int i1i = i1;
				int j0i = j0;
				int j1i = j1;

				d[IX(i, j)] =
					s0 * (t0 * d0[IX(i0i, j0i)]
						+ (t1 * d0[IX(i0i, j1i)]))
					+ s1 * (t0 * d0[IX(i1i, j0i)]
						+ (t1 * d0[IX(i1i, j1i)]));

				
			}
		}
	set_bnd(b, d);
}