

#ifndef _NBODY_RENDERER_
#define _NBODY_RENDERER_

#include <iostream>
#include <thread>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>
static GLFWwindow* window;

// Include GLM
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace glm;

const int OGL_SCREEN_W = 1400;
const int OGL_SCREEN_H = 1200;
const int FBO_MARGIN = 100;

// Helper Objects
#include "Axes.hpp"
#include "Octree.hpp"
#include "Plane.hpp"
#include "Points.hpp"
#include "shader.hpp"

#define CAMERA_CONTROL 1


/**
 * A class encapsualting all OpenGL aspects
 * required to render an N-body simulation.
 *
 * It assumes the body positions are centered about the origin with
 * radius on the scale of ~10 units.
 *
 * Visuals are representative of stars for gravitational N-body.
 *
 * General usage pattern to to create the renderer and then
 * call startRenderThread to start rendering using another thread.
 * The (main) simultation thread should then use updatePoints()
 * to pass new positions to the renderer, as well as
 * shouldClose() to determine if the user has interupted the
 * rendering window thus likely the entire simulation.
 */
class NBodyRenderer {

private :

	Points points;
	Axes ax;
	Plane plane;
	OGLOctree* tree;

	glm::mat4 Projection;
	glm::mat4 V;

	GLuint forwardFBO;
	GLuint forwardRenderTex;

	GLuint blitFBO;
	GLuint blitTex;

	int WhichBlurShader = 1;
	float blurDownScale = 2.0;
	GLuint blurFBOs[4];
	GLuint blurRenderTexs[4];
	GLuint blur1ProgramID;
	GLuint blur2ProgramID;
	GLuint BlurTexID = 0;
	GLuint Blur1HorizLoc; //BlurMultLoc;
	GLuint Blur2HorizLoc; //BlurMultLoc;

	GLuint luminanceFBO;
	GLuint luminanceTex;
	GLuint luminanceProgramID;
	GLuint luminanceTexID = 0;

	GLuint deferredVAOID;
	GLuint deferredVBO;
	GLuint deferredProgramID;
	GLuint forwardTexID = 0;
	GLuint bloomTexID = 1;
	GLuint bloom2TexID = 2;
	GLuint avgLumTexID = 3;

 	double lastTime = glfwGetTime();
	long framecount = 0;
	bool vsync;
	bool init;
	bool close;

	std::thread renderThread;

	glm::vec3 eye = {10.0f, 5.0f, 10.0f};
	glm::vec3 targ = {0.0f, 0.0f, 0.0f};
	glm::vec3 up = {0.0f, 1.0f, 0.0f};
	float radiusFromOrigin = glm::length(eye);


	/**
	 * Control camera with mouse/keyboard.
	 * "Move the camera" while world is stationary.
	 */
	void cameraControlsFreeRoam(glm::mat4& V);


	/**
	 * Control camera with mouse/keyboard.
	 * "Move the world" while camera is stationary.
	 */
	void cameraControlsGlobe(glm::mat4& V);


	/**
	 * Initialize and setup ogl stuff: FBOs, textures, shaders.
	 */
	void setupRenderStages();


	/**
	 * Initialize OpenGL context and everything else.
	 * @return 0 iff successfully setup.
	 */
	int OGL_setup();


public :

	/**
	 * Create the NBodyRenderer along with an OpenGL context.
	 *
	 * @param do_vsync, should the OpenGL window use vsync?
	 */
	NBodyRenderer(bool do_vsync = true) :
		ax(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(2.0f, 2.0f, 2.0f), true),
		plane(30, true, true),
		tree(NULL),
		vsync(do_vsync),
		init(false),
		close(false)
	{

	}


	/**
	 * Destructor. Closes renderer window and ends OpenGL context.
	 */
	~NBodyRenderer() {
		// Close OpenGL window and terminate GLFW
		if (tree != NULL) {
			delete tree;
		}
		if (renderThread.joinable()) {
			renderThread.join();
		}
	}


	/**
	 * Upate the positions of the points.
	 * @see Axes.
	 *
	 * @param N, the number of points.
	 * @param positions, an array of 3*n positional data.
	 * @param idx, a mapping of indices to sort Points colors by.
	 */
	void updatePoints(long N, const double* positions, const float* colors) {
		if (positions == NULL) {
			return;
		}
		if (points.numPoints() != N) {
			points = Points(N, positions, colors);
		} else {
			points.setVertexData(N, positions, colors);
		}
	}


	/**
	 * Upate the positions of the points.
	 * @see Axes.
	 *
	 * @param N, the number of points.
	 * @param positions, an array of 3*N positional data.
	 * @param m, an array of N masses.
	 * @param domainSize, the half-size of the bounding cube containing all points.
	 * @param idx, a mapping of indices to sort Points colors by.
	 */
	void updatePoints(long N, const double* positions, const double* m, double domainSize, float* colors) {
		if (positions == NULL) {
			return;
		}
		if (points.numPoints() != N) {
			points = Points(N, positions, colors);
		} else {
			points.setVertexData(N, positions, colors);
		}

#if defined(DEBUG_OCTREE_OGL) && DEBUG_OCTREE_OGL
		if (tree == NULL) {
			tree = new OGLOctree();
		}
		tree->updateOctree(N, positions, m, domainSize);
#endif
	}

	void setOctree(NBOctree_t* tree) {
		if (this->tree == NULL) {
			this->tree = new OGLOctree(tree);
		} else {
			this->tree->setOctree(tree);
		}
	}


	/**
	 * Render the NBody scene.
	 *
	 * @returns true iff the rendering window should be closed.
	 */
	int draw();


	/**
	 * Start the thread with the render loop,
	 * also initializes the OpenGL context on that thread.
	 */
	void startRenderThread() {
		renderThread = std::thread( [this](){
			while(this->draw() == 0) {}
			glfwTerminate();
		});
	}


	/**
	 * Interrupt the render loop and stop rendering.
	 * Terminates the render thread.
	 */
	void stopRenderThread() {
		close = true;
		if (renderThread.joinable()) {
			renderThread.join();
		}
	}


	/**
	 * Wait for rendering to finish without interrupting.
	 * That is, the user has terminated the rendering window.
	 */
	void joinRenderThread() {
		if (renderThread.joinable()) {
			renderThread.join();
		}
	}


	/**
	 * Determine if the rendering window was asked to be closed.
	 *
	 * @return true iff the render window was asked to be closed.
	 */
	bool shouldClose() {
		return close;
	}


	/**
	 * Determine if the last points update has been rendered yet.
	 *
	 * @return true iff the last positional update has been rendered.
	 */
	bool needsUpdate() {
		return points.pointsStale();
	}

};



#endif
