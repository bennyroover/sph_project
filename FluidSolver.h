#ifndef FLUIDSOLVER
#define FLUIDSOLVER

#include <glm/glm.hpp>
#include <glm/gtc/random.hpp>
#include <vector>
#include <forward_list>

#define PI 3.14159
#define EPS 1e-6
#define BOUND_DAMPING -0.5f
//#define BOUND_DAMPING -1.0f

#define MAX_SPEED 1000.0f

struct Particle {
	glm::vec3 pos;
	glm::vec3 vel;
	glm::vec3 force;
	glm::vec3 force_prev;
	glm::vec3 pos_prev;
	//glm::vec3 vel_prev;
	float p;
	float rho;
	float rho_prev;
	int cell; // index of cell where particle is located
	std::forward_list<int> neighbors;

	Particle() {
		pos = glm::vec3(0.0f, 0.0f, 0.0f);
		pos_prev = glm::vec3(0.0f, 0.0f, 0.0f);
		vel = glm::vec3(0.0f, 0.0f, 0.0f);
		force = glm::vec3(0.0f, 0.0f, 0.0f);
	}
	Particle(glm::vec3 _pos) {
		pos = _pos;
		pos_prev = glm::vec3(0.0f, 0.0f, 0.0f);
		vel = glm::vec3(0.0f, 0.0f, 0.0f);
		force = glm::vec3(0.0f, 0.0f, 0.0f);
		force_prev = glm::vec3(0.0f, 0.0f, 0.0f);
		rho = 1.0f;
		rho_prev = 1.0f;
		p = 0.0f;
		cell = -1;
	}
};

struct Cell {
	int size; // number of members
	std::forward_list<int> members;
	glm::vec3 center;
};


class FluidSolver {
public:

	FluidSolver(int N, float width, float length, float height);
	~FluidSolver();

	void initParticles();
	void reset();
	void step();

	void setN(int _N);
	int getN();

	void addParticle();
	void deleteParticle();

	void setR(float R);
	void setmu(float mu);
	void setrho(float rho0);
	void setk(float k);
	void setG(float G);
	void setDimensions(float l, float w, float h);
	void setConfig(int c);

	float getR();
	float getmu();
	float getrho();
	float getk();
	float getG();
	float getTimeStep();

	glm::vec3 getPos(int i);
	glm::vec3 getVel(int i);

	void equilibrate();
	void printParticles();

	std::vector<Particle> particles;

private:
	std::vector<Cell> cells;

	int N;    // number of particles
	int config; // initial particle configuration
	float dt; // time step

	float m; // particle mass
	float R; // cutoff radius
	float Rsq; // cutoff radius squared
	float cs;
	float rho0; // equilibrium density
	float k; // bulk modulus
	float mu; // viscosity
	float G; // gravity

	glm::vec3 grav; // [m/s^2]

	float width;
	float length;
	float height;

	float xb;
	float yb;
	float zb;

	void calcForces();
	void calcDensityPressure();
	void integrate();

	float pressure(float rho);
	glm::vec3 force_pressure(int i, int j, glm::vec3 rn, float rm);
	glm::vec3 force_viscous(int i, int j, float rm);
	glm::vec3 force_surface(int i, int j);
	glm::vec3 force_lj(int i, int j, glm::vec3 rn, float rm);

	void euler(int i);
	void leapfrog(int i);
	void boundary(int i);
	void checkNan(int i);

	int Xcells; // 
	int Ycells; //
	int Zcells; //
	int Ncells; // 
	float cell_width;
	void initCells(); //
	int mapCell(glm::vec3 p); //
	glm::ivec3 mapPosInt(int i, int nx, int ny, int nz);
	int gridMap(int x, int y, int z, int nx, int nz);
	glm::vec3 mapPos(int i, int nx, int ny, int nz, float d);  //
	glm::vec3 mapCellPos(int i);
	void setCell(int i); //

	float eq_steps;

	float POLY6;
	float SPIKY;
	float SPIKY_GRAD;
	float VISC;
	float VISC_LAP;

	float LJ;
	float ep;
	float sig;

	float poly6(float rsq);
	float spiky(float r);
	float spiky_grad(float r);
	float visc_lap(float r);

	float lenSqr(glm::vec3 v);
	void printVec3(glm::vec3 v);
	void printCells();
};

#endif