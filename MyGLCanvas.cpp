#include "MyGLCanvas.h"
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

MyGLCanvas::MyGLCanvas(int x, int y, int w, int h, const char *l) : Fl_Gl_Window(x, y, w, h, l) {
	lookatPoint = glm::vec3(0.0f, 0.0f, 0.0f);
	scaleFactor = 1.0f;
	mode(FL_OPENGL3 | FL_RGB | FL_ALPHA | FL_DEPTH | FL_DOUBLE);
	rotVec = glm::vec3(0.0f, 0.0f, 0.0f);
	lightPos = glm::vec3(0.0, 2.0, 0.0);
	wireframe = 0;
	grid = 0;
	rotVec.x = rotVec.y = rotVec.z = 0.0f;
	eyePosition.x = 0.0f;
	eyePosition.y = 0.5f;
	eyePosition.z = -1.5f;

	viewAngle = 60;
	clipNear = 0.01f;
	clipFar = 10.0f;

	N = 2000;
	configuration = 0; // set to random
	width = 1.0f;
	length = 1.0f;
	height = 0.5f;

	fluid = new FluidSolver(N, width, length, height);
	running = 0;

	R = 0.05f; // cutoff distance
	rho0 = 8000.0f;  // rest density
	k = 300.0f;    // bulk modulus
	mu = 300.0f; // viscosity
	G = 5000.0f; // gravity
	timeScaling = 1;
	setTimeScaling();

	setParams();

	firstTime = true;

	myShaderManager = new ShaderManager();
}

MyGLCanvas::~MyGLCanvas() {
	delete myShaderManager;
	delete fluid;
}


void MyGLCanvas::start() {
	running = true;
	startCounter();
}

void MyGLCanvas::reset() {
	fluid->reset();
	running = false;
}

void MyGLCanvas::setParams() {
	fluid->setN(N);
	fluid->setR(R);
	fluid->setmu(mu);
	fluid->setk(k);
	fluid->setrho(rho0);
	fluid->setG(G);
	fluid->setConfig(configuration);
}

void MyGLCanvas::updateBox() {
	fluid->setDimensions(length, width, height);
	init_box();
}

void MyGLCanvas::setTimeScaling() {
	dtSim = fluid->getTimeStep() / timeScaling * TIME_RATIO;
}

void MyGLCanvas::draw_particles() {
	//draw_points();
	draw_spheres();
}

void  MyGLCanvas::draw_points() {
	glm::vec3 p;

	glPointSize(3.0);

	glBegin(GL_POINTS);
	for (int i = 0; i < fluid->getN(); i++) {
		p = fluid->getPos(i);
		glColor3f(0.0f, 0.0f, 1.0f);
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
}

void MyGLCanvas::draw_spheres() {

	for (int i = 0; i < fluid->getN(); i++) {
		draw_sphere(fluid->getPos(i));
	}
}

void MyGLCanvas::draw_sphere(glm::vec3 p) {
	glPushMatrix();
	glTranslatef(p.x, p.y, p.z);
	//glColor3f(0.0f, 0.0f, 1.0f);
	glutSolidSphere(RADIUS, SEGMENTS, SEGMENTS);
	glPopMatrix();
}

void MyGLCanvas::draw() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (!valid()) {  //this is called when the GL canvas is set up for the first time...
		puts("establishing GL context");

		glViewport(0, 0, w(), h());
		updateCamera(w(), h());
		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
		glEnable(GL_DEPTH_TEST);
		glPolygonOffset(1, 1);
		if (firstTime == true) {
			firstTime = false;
			initShaders();
		}

		init_box();
		init_particles();
	}

	// Clear the buffer of colors in each bit plane.
	// bit plane - A set of bits that are on or off (Think of a black and white image)

	updateCamera(w(), h());
	draw_box();

	if (running) {
		double newTime = getCurrentTime();
		double frameTime = newTime - currentTime;

		if (frameTime > 250) {
			frameTime = 250;
		}

		currentTime = newTime;
		accumulator += frameTime;

		while (accumulator >= dtSim) {
			fluid->step();
			accumulator -= dtSim;
			t += dtSim;
		}	
	}
	draw_particles_shader();
}

void MyGLCanvas::init_particles() {
	float maxPointSize = 1.5;
	glUseProgram(myShaderManager->program);
	glUniform1f(glGetUniformLocation(myShaderManager->program, "maxPointSize"), maxPointSize);
	glPointSize(10);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_POINT_SPRITE);
	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	glGenVertexArrays(1, &VAO_fluid);
	glGenBuffers(1, &VBO_fluid);
}

void MyGLCanvas::draw_particles_shader() {
	int box = 0;
	int isBoxLocation = glGetUniformLocation(myShaderManager->program, "isBox");

	Particle *particle = &((fluid->particles)[0]);

	glUseProgram(myShaderManager->program);
	glUniform1i(isBoxLocation, box);

	glBindVertexArray(VAO_fluid);
	glBindBuffer(GL_ARRAY_BUFFER, VBO_fluid);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Particle) * fluid->getN(), &((fluid->particles)[0]), GL_DYNAMIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)0);
	glEnableVertexAttribArray(0);

	/*glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)16);
	glEnableVertexAttribArray(1);*/

	glDrawArrays(GL_POINTS, 0, fluid->getN());
}

void MyGLCanvas::init_box() {
	float x1 = -0.5f*width - RADIUS;
	float x2 = 0.5f*width + RADIUS;
	float y1 = -0.5f*height - RADIUS;
	float y2 = 0.5f*height + RADIUS;
	float z1 = -0.5f*length - RADIUS;
	float z2 = 0.5f*length + RADIUS;

	glm::vec3 vertices[24]{
		glm::vec3(x1, y1, z2), glm::vec3(x2, y1, z2),
		glm::vec3(x1, y1, z1), glm::vec3(x2, y1, z1),
		glm::vec3(x1, y1, z2), glm::vec3(x1, y1, z1),
		glm::vec3(x2, y1, z2), glm::vec3(x2, y1, z1),

		glm::vec3(x1, y2, z2), glm::vec3(x2, y2, z2),
		glm::vec3(x1, y2, z1), glm::vec3(x2, y2, z1),
		glm::vec3(x1, y2, z2), glm::vec3(x1, y2, z1),
		glm::vec3(x2, y2, z2), glm::vec3(x2, y2, z1),

		glm::vec3(x1, y1, z2), glm::vec3(x1, y2, z2),
		glm::vec3(x1, y1, z1), glm::vec3(x1, y2, z1),

		glm::vec3(x2, y1, z2), glm::vec3(x2, y2, z2),
		glm::vec3(x2, y1, z1), glm::vec3(x2, y2, z1),
	};

	glGenVertexArrays(1, &VAO_box);
	glGenBuffers(1, &VBO_box);

	glBindVertexArray(VAO_box);

	glBindBuffer(GL_ARRAY_BUFFER, VBO_box);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
	glEnableVertexAttribArray(0);
}

void MyGLCanvas::draw_box() {
	int box = 1;
	int isBoxLocation = glGetUniformLocation(myShaderManager->program, "isBox");

	glUseProgram(myShaderManager->program);
	glUniform1i(isBoxLocation, box);
	glBindVertexArray(VAO_box);
	glDrawArrays(GL_LINES, 0, 24);
}

void MyGLCanvas::initShaders() {

	myShaderManager->initShader("shaders/330/test.vert", "shaders/330/test.frag");
}

void MyGLCanvas::reloadShaders() {
	myShaderManager->resetShaders();

	myShaderManager->initShader("shaders/330/test.vert", "shaders/330/test.frag");

	invalidate();
}

void MyGLCanvas::startCounter() {
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		printf("QueryPerformanceFrequency failed!\n");

	PCFreq = double(li.QuadPart) / 1000.0;

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
	currentTime = getCurrentTime();
	t = 0.0;
	accumulator = 0.0;
}

float MyGLCanvas::getCurrentTime() {
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return float(li.QuadPart - CounterStart) / PCFreq;
}

void MyGLCanvas::updateCamera(int width, int height) {
	float xy_aspect;
	xy_aspect = (float)width / (float)height;
	float ratio = 1;

	glm::mat4 modelViewMatrix = glm::lookAt(eyePosition, lookatPoint, glm::vec3(0.0f, 1.0f, 0.0f));
	modelViewMatrix = glm::rotate(modelViewMatrix, (float)TO_RADIANS(rotVec.x), glm::vec3(1.0f, 0.0f, 0.0f));
	modelViewMatrix = glm::rotate(modelViewMatrix, (float)TO_RADIANS(rotVec.y), glm::vec3(0.0f, 1.0f, 0.0f));
	modelViewMatrix = glm::rotate(modelViewMatrix, (float)TO_RADIANS(rotVec.z), glm::vec3(0.0f, 0.0f, 1.0f));
	modelViewMatrix = glm::scale(modelViewMatrix, glm::vec3(scaleFactor, scaleFactor, scaleFactor));

	//SHADER: passing the projection matrix to the shader... the projection matrix will be called myProjectionMatrix in shader
	glUniformMatrix4fv(glGetUniformLocation(myShaderManager->program, "myModelviewMatrix"), 1, false, glm::value_ptr(modelViewMatrix));

	//passing the light position to the shader
	glUniform3fv(glGetUniformLocation(myShaderManager->program, "myLightPosition"), 1, glm::value_ptr(lightPos));

	// pass far plane
	glUniform1f(glGetUniformLocation(myShaderManager->program, "farPlane"), clipFar);
	
	//SHADER: passing the projection matrix to the shader... the projection matrix will be called myProjectionMatrix in shader
	glm::mat4 perspectiveMatrix;
	perspectiveMatrix = glm::perspective((float)TO_RADIANS(viewAngle), xy_aspect, clipNear, clipFar);
	glUniformMatrix4fv(glGetUniformLocation(myShaderManager->program, "myProjectionMatrix"), 1, false, glm::value_ptr(perspectiveMatrix));
}


int MyGLCanvas::handle(int e) {
	//printf("Event was %s (%d)\n", fl_eventnames[e], e);
#ifndef __APPLE__
	if (firstTime && e == FL_SHOW && shown()) {
		firstTime = 0;
		make_current();
		GLenum err = glewInit(); // defines pters to functions of OpenGL V 1.2 and above
		if (GLEW_OK != err) {
			/* Problem: glewInit failed, something is seriously wrong. */
			fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		}
		else {
			//SHADER: initialize the shader manager and loads the two shader programs
			initShaders();
		}
	}
#endif
	switch (e) {
	case FL_ENTER: cursor(FL_CURSOR_HAND); break;
	case FL_LEAVE: cursor(FL_CURSOR_DEFAULT); break;
	case FL_KEYUP:
		printf("keyboard event: key pressed: %c\n", Fl::event_key());
		switch (Fl::event_key()) {
		case 'w': eyePosition.y += 0.05f;  break;
		case 'a': eyePosition.x += 0.05f; break;
		case 's': eyePosition.y -= 0.05f;  break;
		case 'd': eyePosition.x -= 0.05f; break;
		}
		updateCamera(w(), h());
		break;
	case FL_MOUSEWHEEL:
		printf("mousewheel: dx: %d, dy: %d\n", Fl::event_dx(), Fl::event_dy());
		eyePosition.z += Fl::event_dy() * -0.05f;
		updateCamera(w(), h());
		break;
	}

	return Fl_Gl_Window::handle(e);
}

void MyGLCanvas::resize(int x, int y, int w, int h) {
	Fl_Gl_Window::resize(x, y, w, h);
	puts("resize called");
}