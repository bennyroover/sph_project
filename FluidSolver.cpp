#include "FluidSolver.h"
#include <cmath>

FluidSolver::FluidSolver(int _N, float _width, float _length, float _height) {
	N = _N;
	width = _width;
	length = _length;
	height = _height;

	xb = 0.5f*width;
	yb = 0.5f*height;
	zb = 0.5f*length;

	dt = 0.0005f; // time

	R = 0.06f; // cutoff distance
	Rsq = R * R;
	rho0 = 100.0f;  // rest density
	k = 1.0f;    // bulk modulus
	m = rho0 * PI*Rsq*R; // mass
	mu = 1000.0f; // viscosity

	POLY6 = 315.0f / (64.0f*PI*glm::pow(R, 9.0f));
	SPIKY = 15.0f / (PI*glm::pow(R, 6.0f));
	SPIKY_GRAD = -45.0f / (PI*glm::pow(R, 6.0f));
	VISC = 15.0f / (2.0f*PI*glm::pow(R, 3.0f));
	VISC_LAP = 45.0f / (PI*glm::pow(R, 6.0f));

	ep = 0.001f;
	sig = R;
	LJ = 48.0f*ep;

	eq_steps = 10;

	G = 1000.0f;
	grav = G * glm::vec3(0.0f, -1.0f, 0.0f);

	cell_width = R;
	config = 0;

	initCells();
	initParticles();
}

FluidSolver::~FluidSolver() {
}

void FluidSolver::reset() {
	particles.clear();
	//cells.clear();
	initParticles();
	initCells();
}

int FluidSolver::getN() {
	return N;
}

void FluidSolver::setN(int _N) {

	if (_N > N) {
		while (N < _N) {
			addParticle();
		}
	} else {
		while (N > _N) {
			deleteParticle();
		}
	}
}

void FluidSolver::setConfig(int _config) {
	config = _config;
}

void FluidSolver::addParticle() {
	float x = glm::linearRand(-1.0f*xb, xb);
	float y = glm::linearRand(-1.0f*yb, yb);
	float z = glm::linearRand(-1.0f*zb, zb);
	glm::vec3 pos = glm::vec3(x, y, z);

	Particle particle(glm::vec3(x, y, z));
	particle.cell = mapCell(pos);
	particles.push_back(particle);
	N++;
	cells[particle.cell].members.push_front(N-1);
	cells[particle.cell].size++;
}

void FluidSolver::deleteParticle() {
	cells[particles[N-1].cell].members.remove(N-1);
	cells[particles[N-1].cell].size--;
	particles.pop_back();
	N--;
}

glm::vec3 FluidSolver::getPos(int i) {
	return particles[i].pos;
}

glm::vec3 FluidSolver::getVel(int i) {
	return particles[i].vel;
}

void FluidSolver::setR(float _R) {
	R = _R;
	Rsq = R * R;
	cell_width = R;
	initCells();
}

void FluidSolver::setmu(float _mu) {
	mu = _mu;
}

void FluidSolver::setrho(float _rho0) {
	rho0 = _rho0;
	m = rho0 * PI*Rsq*R;
}

void FluidSolver::setk(float _k) {
	k = _k;
}

void FluidSolver::setG(float _G) {
	G = _G;
	grav = G * glm::vec3(0.0f, -1.0f, 0.0f);
}

void FluidSolver::setDimensions(float l, float w, float h) {
	length = l;
	width = w;
	height = h;

	xb = 0.5f*width;
	yb = 0.5f*height;
	zb = 0.5f*length;

	initCells();
}

float FluidSolver::getR() {
	return R;
}

float FluidSolver::getmu() {
	return mu;
}

float FluidSolver::getrho() {
	return rho0;
}

float FluidSolver::getk() {
	return k;
}

float FluidSolver::getTimeStep() {
	return dt;
}

void FluidSolver::initParticles() {
	float mag = 0.01f;
	glm::vec3 pos;
	glm::vec3 jiggle;
	for (int i = 0; i < N; i++) {
		if (config == 0) {
			float x = glm::linearRand(-1.0f*xb, xb);
			float y = glm::linearRand(-1.0f*yb, yb);
			float z = glm::linearRand(-1.0f*zb, zb);
			pos = glm::vec3(x,y,z);
		} else if (config == 1) {
			float x = glm::linearRand(0.5f*xb, xb);
			float y = glm::linearRand(-1.0f*yb, 0.5f*yb);
			float z = glm::linearRand(-0.8f*zb, 0.8f*zb);
			pos = glm::vec3(x, y, z);
		} else if (config == 2) {
			pos = glm::ballRand(0.2f*width);
		}
		Particle particle(pos);
		particle.rho_prev = m * poly6(0.0f);
		particle.cell = mapCell(pos);
		particles.push_back(particle);
		cells[particle.cell].members.push_front(i);
		cells[particle.cell].size++;
	}
	for (int i = 0; i < eq_steps; i++) {
		equilibrate();
	}
}

void FluidSolver::initCells() {
	cells.clear();
	Xcells = glm::ceil(width / cell_width);
	Ycells = glm::ceil(height / cell_width);
	Zcells = glm::ceil(length / cell_width);
	Ncells = Xcells * Ycells * Zcells;

	for (int i = 0; i < Ncells; i++) {
		Cell cell;
		cell.center = mapCellPos(i);
		cell.size = 0;
		cells.push_back(cell);
	}
}


// 
void FluidSolver::equilibrate() {
	glm::vec3 delta;
	float len;

	step();

	for (int i = 0; i < particles.size(); i++) {
		delta = particles[i].pos - particles[i].pos_prev;
		len = glm::length(delta);
		if (len > R) {
			particles[i].pos = particles[i].pos_prev + 2.0f*R*delta/len;
		}
		if (glm::length(particles[i].vel) > MAX_SPEED) {
			particles[i].vel = glm::vec3(0.0f, 0.0f, 0.0f);
		}
	}
}

void FluidSolver::step() {
	calcDensityPressure();
	calcForces();
	//printf("FORCE\n");
	//printParticles();
	integrate();
	//printf("INTEGRATE\n");
	//printParticles();
}

//
void FluidSolver::calcDensityPressure() {
	float rsq, rm;

	for (int i = 0; i < particles.size(); i++) {
		particles[i].force_prev = particles[i].force;
		particles[i].force = glm::vec3(0.0f, 0.0f, 0.0f);
		particles[i].neighbors.clear();

		particles[i].rho_prev = particles[i].rho;
		particles[i].rho = poly6(0.0f);

		int ncell;
		glm::ivec3 id = mapPosInt(particles[i].cell, Xcells, Ycells, Zcells);
		std::forward_list<int>::iterator it;

		for (int yi = -1; yi <= 1; yi++) {
			for (int zi = -1; zi <= 1; zi++) {
				for (int xi = -1; xi <= 1; xi++) {
					ncell = gridMap(id.x+xi, id.y+yi, id.z+zi, Xcells, Zcells);
					if (ncell < 0 or ncell >= Ncells) continue;
					if (cells[ncell].size == 0) continue;
					
					for (it = cells[ncell].members.begin(); it != cells[ncell].members.end(); ++it) {
						if (*it == i) continue;
						rm = glm::length(particles[*it].pos - particles[i].pos);
						rsq = rm * rm;
						if (rsq < Rsq) {
							particles[i].neighbors.push_front(*it);
							particles[i].rho += poly6(rsq);
							//particles[i].rho += spiky(rm);
						}
					}
				}
			}
		}
		particles[i].rho *= m;
		particles[i].p = pressure(particles[i].rho);
	}
}

float FluidSolver::pressure(float rho) {
	return k * (rho - rho0); // ideal gas
	//return k*(glm::pow(rho/rho0,7.0f) - 1.0f); // tait equation
}

// F_i = Fp + Fv + Fg
void FluidSolver::calcForces() {
	glm::vec3 rn; // normalize vector between particles
	float rm;	  // distance between particles
	std::forward_list<int>::iterator it;

	for (int i = 0; i < particles.size(); i++) {
		for (it = particles[i].neighbors.begin(); it != particles[i].neighbors.end(); ++it) {
			rn = particles[*it].pos - particles[i].pos;
			rm = glm::length(rn);
			rn = glm::normalize(rn);
			particles[i].force = particles[i].force + force_pressure(i, *it, rn, rm) + force_viscous(i, *it, rm);
			particles[*it].force = particles[*it].force + force_pressure(*it, i, -1.0f*rn, rm) + force_viscous(*it, i, rm);
			particles[*it].neighbors.remove(i);
		}
		particles[i].force = particles[i].force + particles[i].rho*grav;
	}
}

// Fp = -m * sum_j [(p_i+p_j)/(2*rho_j)*kernel_grad(r_i-r_j)
glm::vec3 FluidSolver::force_pressure(int i, int j, glm::vec3 rn, float rm) {
	return -1.0f*rn*m*(particles[i].p + particles[j].p) / (2.0f*particles[j].rho) * spiky_grad(rm);
}

// Fv = mu*m * sum_j [(v_j - v_i)/rho_j]*kernel_lap(r_i-r_j)
glm::vec3 FluidSolver::force_viscous(int i, int j, float rm) {
	return (mu*m/particles[j].rho)*(particles[j].vel - particles[i].vel)*visc_lap(rm);
}

glm::vec3 FluidSolver::force_surface(int i, int j) {
	return glm::vec3(0.0f,0.0f,0.0f);
}

glm::vec3 FluidSolver::force_lj(int i, int j, glm::vec3 rn, float r) {
	return (LJ*rn / r)*(glm::pow(sig / r, 12.0f) - 0.5f*glm::pow(sig / r, 6.0f));
}

void FluidSolver::integrate() {
	glm::vec3 temp;
	for (int i = 0; i < particles.size(); i++) {
		temp = particles[i].pos;
		//euler(i);
		leapfrog(i);
		boundary(i);
		//checkNan(i);
		particles[i].pos_prev = temp;
		setCell(i);
	}
}

// 3d point ot 1d index
int FluidSolver::mapCell(glm::vec3 p) {
	glm::ivec3 comps = glm::floor((p + glm::vec3(xb, yb, zb)) / cell_width);
	if (comps.x >= Xcells) {
		comps.x = Xcells - 1;
	}
	else if (comps.x < 0) {
		comps.x = 0;
	}
	if (comps.y >= Ycells) {
		comps.y = Ycells - 1;
	}
	else if (comps.y < 0) {
		comps.y = 0;
	}
	if (comps.z >= Zcells) {
		comps.z = Zcells - 1;
	}
	else if (comps.z < 0) {
		comps.z = 0;
	}
	return gridMap(comps.x, comps.y, comps.z, Xcells, Zcells);
}

// 3d index to 1d index
int FluidSolver::gridMap(int x, int y, int z, int nx, int nz) {
	return x + z * nx + y * nx*nz;
}

//
glm::ivec3 FluidSolver::mapPosInt(int i, int nx, int ny, int nz) {
	//return glm::ivec3(i % nx, i / (nx*nz), (i / nz) % nz);
	return glm::ivec3(i % nx, i / (nx*ny), (i / nz) % nz);
}

glm::vec3 FluidSolver::mapPos(int i, int nx, int ny, int nz, float d) {
	return d * glm::vec3(mapPosInt(i, nx, ny, nz)) - glm::vec3(xb, yb, zb);
}

glm::vec3 FluidSolver::mapCellPos(int i) {
	return mapPos(i, Xcells, Ycells, Zcells, cell_width);
}

void FluidSolver::setCell(int i) {

	if (mapCell(particles[i].pos) != particles[i].cell) {
		cells[particles[i].cell].members.remove(i);
		cells[particles[i].cell].size--;
		particles[i].cell = mapCell(particles[i].pos);
		cells[particles[i].cell].members.push_front(i);
		cells[particles[i].cell].size++;
	}
}

void FluidSolver::boundary(int i) {
	if (particles[i].pos.x < -1.0f*xb) {
		particles[i].pos.x = -1.0f*xb + EPS;
		particles[i].vel.x *= BOUND_DAMPING;
	}
	else if (particles[i].pos.x > xb) {
		particles[i].pos.x = xb - EPS;
		particles[i].vel.x *= BOUND_DAMPING;
	}
	if (particles[i].pos.y < -1.0f*yb) {
		particles[i].pos.y = -1.0f*yb + EPS;
		particles[i].vel.y *= BOUND_DAMPING;
	}
	else if (particles[i].pos.y > yb) {
		particles[i].pos.y = yb - EPS;
		particles[i].vel.y *= BOUND_DAMPING;
	}
	if (particles[i].pos.z < -1.0f*zb) {
		particles[i].pos.z = -1.0f*zb + EPS;
		particles[i].vel.z *= BOUND_DAMPING;
	}
	else if (particles[i].pos.z > zb) {
		particles[i].pos.z = zb - EPS;
		particles[i].vel.z *= BOUND_DAMPING;
	}
}

void FluidSolver::checkNan(int i) {
	if (isnan(particles[i].pos.x)) {
		particles[i].pos.x = 0.0f;
	}
	if (isnan(particles[i].pos.y)) {
		particles[i].pos.y = 0.0f;
	}
	if (isnan(particles[i].pos.z)) {
		particles[i].pos.z = 0.0f;
	}
	if (isnan(particles[i].vel.x)) {
		particles[i].vel.x = 0.0f;
	}
	if (isnan(particles[i].vel.y)) {
		particles[i].vel.y = 0.0f;
	}
	if (isnan(particles[i].vel.z)) {
		particles[i].vel.z = 0.0f;
	}
	if (isnan(particles[i].force.x)) {
		particles[i].force.x = 0.0f;
	}
	if (isnan(particles[i].force.y)) {
		particles[i].force.y = 0.0f;
	}
	if (isnan(particles[i].force.z)) {
		particles[i].force.z = 0.0f;
	}
}


// implicit euler
void FluidSolver::euler(int i) {
	particles[i].vel = particles[i].vel + dt * particles[i].force / particles[i].rho;
	particles[i].pos = particles[i].pos + dt * particles[i].vel;
}

//
void FluidSolver::leapfrog(int i) {
	particles[i].pos = particles[i].pos + dt * (particles[i].vel + 0.5f*dt*particles[i].force_prev / particles[i].rho_prev);
	particles[i].vel = particles[i].vel + 0.5f*dt*(particles[i].force_prev / particles[i].rho_prev + particles[i].force / particles[i].rho);
}

float FluidSolver::poly6(float rsq) {
	if (rsq < Rsq) {
		return POLY6 * glm::pow((Rsq - rsq), 3);
	}
	return 0.0f;
}

float FluidSolver::spiky(float r) {
	if (r < R) {
		return SPIKY * glm::pow((R - r), 3);
	}
	return 0.0f;
}

float FluidSolver::spiky_grad(float r) {
	if (r <= R) {
		return SPIKY_GRAD * glm::pow(R - r, 2.0f);
	}
	return 0.0f;
}

float FluidSolver::visc_lap(float r) {
	if (r <= R) {
		return VISC_LAP * (R - r);
	}
	return 0.0f;
}


float FluidSolver::lenSqr(glm::vec3 v) {
	return v.x*v.x + v.y*v.y + v.z*v.z;
}

void FluidSolver::printParticles() {
	std::forward_list<int>::iterator it;

	printf("PARTICLES\n");
	for (int i = 0; i < N; i++) {
		printf("%d cell=%d ", i, particles[i].cell);
		printf("neighbors: ");
		for (it = particles[i].neighbors.begin(); it != particles[i].neighbors.end(); ++it) {
			printf("%d ", *it);
		}
		//printf("rho=%.1f p=%.1f\n", particles[i].rho, particles[i].p);
		//printVec3(particles[i].pos);
		//printVec3(particles[i].vel);
		//printVec3(particles[i].force);
		//printVec3(mapCellPos(particles[i].cell));
		printf("\n");
	}
}

void FluidSolver::printCells() {
	std::forward_list<int>::iterator it;
	int num = 0;
	printf("\nCELLS ncells=%d\n", Ncells);
	for (int i = 0; i < Ncells; i++) {
		if (cells[i].size != 0) {
			printf("Cell %d: ", i);
			for (it = cells[i].members.begin(); it != cells[i].members.end(); ++it) {
				printf("%d ", *it);
				num++;
			}
			printf("\n");
		}
	}
	printf("N=%d n=%d\n", N, num);
}

void FluidSolver::printVec3(glm::vec3 v) {
	printf("%.2f %.2f %.2f\n", v.x, v.y, v.z);
}


