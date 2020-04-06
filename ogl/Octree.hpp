
#ifndef _OGL_OCTREE_
#define _OGL_OCTREE_

#include <mutex>
#include "parallel/include/Synchronized.hpp"
#include "NBodyOctree.h"


/**
 * A class to render the edges of
 * an Octree in OpenGL.
 */
class OGLOctree {

	NBOctree_t* tree;

	std::mutex m_mutex;

	glm::vec4 color = glm::vec4(0.3f, 0.3f, 0.3f, 0.5f);

	/**
	 * Recursively render an octree node and its children.
	 */
	void renderNode(NBOctreeNode_t* node) {
		if (node == NULL) {
			return;
		}

		float xmax = node->center[0] + node->size;
		float xmin = node->center[0] - node->size;
		float ymax = node->center[1] + node->size;
		float ymin = node->center[1] - node->size;
		float zmax = node->center[2] + node->size;
		float zmin = node->center[2] - node->size;

		//12 lines per node = 24 vertices

		//top face
		glVertex3f(xmax, ymax, zmax);
		glVertex3f(xmax, ymax, zmin);
		glVertex3f(xmax, ymax, zmin);
		glVertex3f(xmin, ymax, zmin);
		glVertex3f(xmin, ymax, zmin);
		glVertex3f(xmin, ymax, zmax);
		glVertex3f(xmin, ymax, zmax);
		glVertex3f(xmax, ymax, zmax);

		//bottom face;
		glVertex3f(xmax, ymin, zmax);
		glVertex3f(xmax, ymin, zmin);
		glVertex3f(xmax, ymin, zmin);
		glVertex3f(xmin, ymin, zmin);
		glVertex3f(xmin, ymin, zmin);
		glVertex3f(xmin, ymin, zmax);
		glVertex3f(xmin, ymin, zmax);
		glVertex3f(xmax, ymin, zmax);

		//the four verticles lines
		glVertex3f(xmax, ymin, zmax);
		glVertex3f(xmax, ymax, zmax);
		glVertex3f(xmax, ymin, zmin);
		glVertex3f(xmax, ymax, zmin);
		glVertex3f(xmin, ymin, zmin);
		glVertex3f(xmin, ymax, zmin);
		glVertex3f(xmin, ymin, zmax);
		glVertex3f(xmin, ymax, zmax);

		if (node->children != NULL) {
			for (int i = 0; i < 8; ++i) {
				renderNode(node->children[i]);
			}
		}

	}

public:

	/**
	 * OGLOctree default constructor. Sets underlying tree to empty.
	 */
	OGLOctree() : tree(NULL) {}


	/**
	 * Construct an OpenGL octree from an NBodyOctree.
	 * The OGLOctree takes ownership of the NBOctee
	 * and frees it on destruction.
	 * @param tree, the NBodyOctree.
	 */
	OGLOctree(NBOctree_t* tree_in) : tree(tree_in) {}


	/**
	 * Destructor of the OGLOctree. Frees the underlying NBOctree.
	 */
	~OGLOctree() {
		if (tree != NULL) {
			freeOctree_NB(tree);
		}
	}


	/**
	 * Update the underlying octree given new positions.
	 * This is synchronized with drawing so that other threads
	 * may meanwhile the render thread renders.
	 *
	 * @param N, the number of bodies.
	 * @param r, an array of 3*N values representing the bodies' positions.
	 * @param m, an array of N values representing the bodies' masses.
	 * @param domainSize, the size of the entire octree.
	 */
	void updateOctree(long N, const double* r, const double* m, double domainSize) {
		synchronized_nonrecursive(m_mutex) {
        	buildOctreeInPlace_NB(N, r, m, domainSize, &tree);
		}
	}


	/**
	 * Update the underlying octree structure.
	 * Reserves the octree here to it continues to render it
	 * until a new octree is set.
	 *
	 * @param tree_in, the new octree.
	 */
	void setOctree(NBOctree_t* tree_in) {
		freeOctree_NB(tree);
		tree = tree_in;
	}


	/**
	 * Draw the octree.
	 */
	void draw() {
		if (tree == NULL) {
			return;
		}

		glMatrixMode( GL_MODELVIEW );
		glPushMatrix();


		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glEnable(GL_DEPTH_TEST);
		// Accept fragment if it closer to the camera than the former one
		glDepthFunc(GL_LESS);

		glLineWidth(1.0f);
		glColor4f(color.x, color.y, color.z, color.a);

		synchronized_nonrecursive(m_mutex) {
			glBegin(GL_LINES);
			if (tree != NULL) {
				renderNode(tree->root);
			}
			glEnd();
		}

		glPopMatrix();

		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
	}

};


#endif
