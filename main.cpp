/*  =================== File Information =================
	File Name: main.cpp
	Description:
	Author: Michael Shah

	Purpose: Driver for 3D program to load .ply models
	Usage:
	===================================================== */

// from lab

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/names.h>

#include "MyGLCanvas.h"

using namespace std;

#define MAX_N 10000
#define MAX_RHO 10000
#define MAX_MU 1000
#define MAX_K 1000
#define K_STEP 10
#define MAX_CUT 0.2
#define MIN_CUT 0.01
#define MAX_G 10000


class MyAppWindow;
MyAppWindow *win;

class MyAppWindow : public Fl_Window {
public:
	Fl_Button* startButton;
	Fl_Button* resetButton;
	Fl_Slider* rotXSlider;
	Fl_Slider* rotYSlider;
	Fl_Slider* rotZSlider;


	Fl_Slider* NSlider; // particle number
	Fl_Slider* muSlider;  // viscosity
	Fl_Slider* rhoSlider; // density
	Fl_Slider* kSlider; // bulk modulus
	Fl_Slider* GSlider; // gravity
	Fl_Slider* RSlider; // cutoff distance
	Fl_Slider* TimeSlider; // time scaling slider

	Fl_Slider* lengthSlider;
	Fl_Slider* widthSlider;
	Fl_Slider* heightSlider;

	Fl_Button* reloadButton;

	MyGLCanvas* canvas;

public:
	// APP WINDOW CONSTRUCTOR
	MyAppWindow(int W, int H, const char*L = 0);

	static void idleCB(void* userdata) {
		win->canvas->redraw();
	}

private:
	// Someone changed one of the sliders
	static void rotateCB(Fl_Widget* w, void* userdata) {
		float value = ((Fl_Slider*)w)->value();
		*((float*)userdata) = value;
	}

	static void colorCB(Fl_Widget* w, void* userdata) {
		float value = ((Fl_Slider*)w)->value();
		*((float*)userdata) = value;
	}

	static void toggleCB(Fl_Widget* w, void* userdata) {
		int value = ((Fl_Button*)w)->value();
		printf("value: %d\n", value);
		*((int*)userdata) = value;
	}

	static void startCB(Fl_Widget*w, void*data) {
		cout << "start" << endl;
		win->canvas->start();
	}

	static void resetCB(Fl_Widget*w, void*data) {
		cout << "reset" << endl;
		win->canvas->reset();
	}

	static void timeScaleCB(Fl_Widget*w, void*data) {
		float value = ((Fl_Slider*)w)->value();
		*((float*)data) = value;
		win->canvas->setTimeScaling();
	}

	static void paramsCB(Fl_Widget* w, void* userdata) {
		float value = ((Fl_Slider*)w)->value();
		*((float*)userdata) = value;
		win->canvas->setParams();
	}

	static void boxCB(Fl_Widget* w, void* userdata) {
		float value = ((Fl_Slider*)w)->value();
		*((float*)userdata) = value;
		win->canvas->updateBox();
	}

	static void configCB(Fl_Widget* w, void* userdata) {
		const char* value = ((Fl_Button*)w)->label();
		if (strcmp("Random", value) == 0) {
			win->canvas->configuration = 0;
		}
		else if (strcmp("Wave", value) == 0) {
			win->canvas->configuration = 1;
		}
		else if (strcmp("Ball", value) == 0) {
			win->canvas->configuration = 2;
		}
		win->canvas->setParams();
	}

	//static void reloadCB(Fl_Widget* w, void* userdata) {
	//	win->canvas->reloadShaders();
	//}

};


MyAppWindow::MyAppWindow(int W, int H, const char*L) : Fl_Window(W, H, L) {
	begin();
	// OpenGL window

	canvas = new MyGLCanvas(10, 10, w() - 310, h() - 20);

	Fl_Pack* pack = new Fl_Pack(w() - 300, 30, 150, h(), "");
	pack->box(FL_DOWN_FRAME);
	pack->labelfont(1);
	pack->type(Fl_Pack::VERTICAL);
	pack->spacing(0);
	pack->begin();

		Fl_Pack* buttonsPack = new Fl_Pack(w() - 100, 30, 150, h(), "");
		buttonsPack->box(FL_DOWN_FRAME);
		buttonsPack->labelfont(1);
		buttonsPack->type(Fl_Pack::VERTICAL);
		buttonsPack->spacing(0);
		buttonsPack->begin();

		startButton = new Fl_Button(0, 0, pack->w() - 20, 20, "Start");
		startButton->callback(startCB, (void*)this);

		resetButton = new Fl_Button(0, 0, pack->w() - 20, 20, "Reset");
		resetButton->callback(resetCB, (void*)this);

		buttonsPack->end();

		Fl_Pack* configPack = new Fl_Pack(w() - 100, 130, 150, h(), "");
		configPack->box(FL_DOWN_FRAME);
		configPack->labelfont(1);
		configPack->type(Fl_Pack::VERTICAL);
		configPack->spacing(0);
		configPack->begin();

		{ Fl_Round_Button* tmpButton = new Fl_Round_Button(0, 0, pack->w() - 20, 20, "Random");
			tmpButton->type(102);
			tmpButton->down_box(FL_ROUND_DOWN_BOX);
			tmpButton->value(1);
			tmpButton->callback((Fl_Callback*)configCB);
		}
		{ Fl_Round_Button* tmpButton = new Fl_Round_Button(0, 0, pack->w() - 20, 20, "Wave");
			tmpButton->type(102);
			tmpButton->down_box(FL_ROUND_DOWN_BOX);
			tmpButton->callback((Fl_Callback*)configCB);
		}
		{ Fl_Round_Button* tmpButton = new Fl_Round_Button(0, 0, pack->w() - 20, 20, "Ball");
			tmpButton->type(102);
			tmpButton->down_box(FL_ROUND_DOWN_BOX);
			tmpButton->callback((Fl_Callback*)configCB);
		}

		configPack->end();

	Fl_Box *NTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "# Particles");
	NSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
	NSlider->align(FL_ALIGN_TOP);
	NSlider->type(FL_HOR_SLIDER);
	NSlider->bounds(0, MAX_N);
	NSlider->step(1);
	NSlider->value(canvas->N);
	NSlider->callback(paramsCB, (void*)(&(canvas->N)));

	Fl_Box *TimeTextBox = new Fl_Box(0, 0, pack->w() - 20, 20, "Time Scaling");
	TimeSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
	TimeSlider->align(FL_ALIGN_TOP);
	TimeSlider->type(FL_HOR_SLIDER);
	TimeSlider->bounds(0.1, 2);
	TimeSlider->step(0.1);
	TimeSlider->value(canvas->timeScaling);
	TimeSlider->callback(timeScaleCB, (void*)(&(canvas->timeScaling)));

	Fl_Box *muTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Viscosity");
	muSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
	muSlider->align(FL_ALIGN_TOP);
	muSlider->type(FL_HOR_SLIDER);
	muSlider->bounds(0, MAX_MU);
	muSlider->step(10);
	muSlider->value(canvas->mu);
	muSlider->callback(paramsCB, (void*)(&(canvas->mu)));

	Fl_Box *rhoTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Density");
	muSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
	muSlider->align(FL_ALIGN_TOP);
	muSlider->type(FL_HOR_SLIDER);
	muSlider->bounds(1, MAX_RHO);
	muSlider->step(100);
	muSlider->value(canvas->rho0);
	muSlider->callback(paramsCB, (void*)(&(canvas->rho0)));

	Fl_Box *kTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Bulk Modulus");
	kSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
	kSlider->align(FL_ALIGN_TOP);
	kSlider->type(FL_HOR_SLIDER);
	kSlider->bounds(0, MAX_K);
	kSlider->step(K_STEP);
	kSlider->value(canvas->k);
	kSlider->callback(paramsCB, (void*)(&(canvas->k)));

	Fl_Box *GTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Gravity");
	GSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
	GSlider->align(FL_ALIGN_TOP);
	GSlider->type(FL_HOR_SLIDER);
	GSlider->bounds(0, MAX_G);
	GSlider->step(100);
	GSlider->value(canvas->G);
	GSlider->callback(paramsCB, (void*)(&(canvas->G)));

	Fl_Box *RTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Cutoff Distance");
	RSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
	RSlider->align(FL_ALIGN_TOP);
	RSlider->type(FL_HOR_SLIDER);
	RSlider->bounds(MIN_CUT, MAX_CUT);
	RSlider->step(0.001);
	RSlider->value(canvas->R);
	RSlider->callback(paramsCB, (void*)(&(canvas->R)));

	pack->end();

	Fl_Pack* packCol2 = new Fl_Pack(w() - 150, 30, 150, h(), "");
	packCol2->box(FL_DOWN_FRAME);
	packCol2->type(Fl_Pack::VERTICAL);
	packCol2->spacing(0);
	packCol2->begin();


		//slider for controlling rotation
		Fl_Box *rotXTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "RotateX");
		rotXSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
		rotXSlider->align(FL_ALIGN_TOP);
		rotXSlider->type(FL_HOR_SLIDER);
		rotXSlider->bounds(-359, 359);
		rotXSlider->step(1);
		rotXSlider->value(canvas->rotVec.x);
		rotXSlider->callback(rotateCB, (void*)(&(canvas->rotVec.x)));

		Fl_Box *rotYTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "RotateY");
		rotYSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
		rotYSlider->align(FL_ALIGN_TOP);
		rotYSlider->type(FL_HOR_SLIDER);
		rotYSlider->bounds(-359, 359);
		rotYSlider->step(1);
		rotYSlider->value(canvas->rotVec.y);
		rotYSlider->callback(rotateCB, (void*)(&(canvas->rotVec.y)));

		Fl_Box *rotZTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "RotateZ");
		rotZSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
		rotZSlider->align(FL_ALIGN_TOP);
		rotZSlider->type(FL_HOR_SLIDER);
		rotZSlider->bounds(-359, 359);
		rotZSlider->step(1);
		rotZSlider->value(canvas->rotVec.z);
		rotZSlider->callback(rotateCB, (void*)(&(canvas->rotVec.z)));

		Fl_Box *lengthTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Length");
		lengthSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
		lengthSlider->align(FL_ALIGN_TOP);
		lengthSlider->type(FL_HOR_SLIDER);
		lengthSlider->bounds(0.2, 2);
		lengthSlider->step(0.1);
		lengthSlider->value(canvas->length);
		lengthSlider->callback(boxCB, (void*)(&(canvas->length)));

		Fl_Box *widthTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Width");
		widthSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
		widthSlider->align(FL_ALIGN_TOP);
		widthSlider->type(FL_HOR_SLIDER);
		widthSlider->bounds(0.2, 2);
		widthSlider->step(0.1);
		widthSlider->value(canvas->length);
		widthSlider->callback(boxCB, (void*)(&(canvas->width)));

		Fl_Box *heightTextbox = new Fl_Box(0, 0, pack->w() - 20, 20, "Height");
		widthSlider = new Fl_Value_Slider(0, 0, pack->w() - 20, 20, "");
		widthSlider->align(FL_ALIGN_TOP);
		widthSlider->type(FL_HOR_SLIDER);
		widthSlider->bounds(0.2, 5);
		widthSlider->step(0.1);
		widthSlider->value(canvas->height);
		widthSlider->callback(boxCB, (void*)(&(canvas->height)));

	

	/*Fl_Pack* packShaders = new Fl_Pack(w() - 100, 230, 100, h(), "Shaders");
	packShaders->box(FL_DOWN_FRAME);
	packShaders->labelfont(1);
	packShaders->type(Fl_Pack::VERTICAL);
	packShaders->spacing(0);
	packShaders->begin();

	reloadButton = new Fl_Button(0, 0, pack->w() - 20, 20, "Reload");
	reloadButton->callback(reloadCB, (void*)this);

	packShaders->end();*/

	packCol2->end();

	end();
}


/**************************************** main() ********************/
int main(int argc, char **argv) {
	win = new MyAppWindow(1000, 600, "Fluid Animation");
	win->resizable(win);
	Fl::add_idle(MyAppWindow::idleCB);
	win->show();
	return(Fl::run());
}