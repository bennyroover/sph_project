#pragma once

// from lab

#ifndef MYGLCANVAS_H
#define MYGLCANVAS_H

#define SEGMENTS 10
#define RADIUS 0.02
#define TIME_RATIO 60000

#if defined(__APPLE__)
#  include <OpenGL/gl3.h> // defines OpenGL 3.0+ functions
#else
#  if defined(WIN32)
#    define GLEW_STATIC 1
#  endif
#  include <GL/glew.h>
#endif
#include <FL/glut.h>
#include <FL/glu.h>
#include <glm/glm.hpp>
#include <time.h>
#include <iostream>

#include "ShaderManager.h"

#include "FluidSolver.h"
#include <Windows.h>

#define TO_RADIANS(a) (a*PI/180.0f)

class MyGLCanvas : public Fl_Gl_Window {
public:
	int wireframe, orbits, grid;
	glm::vec3 rotVec;
	glm::vec3 eyePosition;
	glm::vec3 lookatPoint;
	glm::vec3 lightPos;
	float scaleFactor;

	FluidSolver* fluid;

	int viewAngle;
	float clipNear;
	float clipFar;

	MyGLCanvas(int x, int y, int w, int h, const char *l = 0);
	~MyGLCanvas();

	void start();
	void reset();
	void setTimeScaling();

	void setParams();
	void updateBox();

	float N; // number of particles
	float R; // cutoff radius
	float rho0; // equilibrium density
	float k; // bulk modulus
	float mu; // viscosity
	float G; // gravity
	float timeScaling; // time scale
	float width;
	float length;
	float height;

	int configuration;

	void reloadShaders();


private:
	void draw();
	void drawScene();
	int handle(int);
	void resize(int x, int y, int w, int h);
	void updateCamera(int width, int height);

	unsigned int VBO_box;
	unsigned int VAO_box;

	unsigned int VBO_fluid;
	unsigned int VAO_fluid;

	void initShaders();
	ShaderManager* myShaderManager;

	void draw_particles();
	void draw_points();
	void draw_spheres();
	void draw_sphere(glm::vec3 p);


	void init_box();
	void draw_box();
	
	void init_particles();
	void draw_particles_shader();

	void startCounter();
	float getCurrentTime();

	int running;

	double PCFreq = 0.0;
	__int64 CounterStart = 0;

	float currentTime;
	float t;
	float dtSim;
	float accumulator;

	bool firstTime;
};

#endif // !MYGLCANVAS_H