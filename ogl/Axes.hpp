
#ifndef _OGL_AXES_
#define _OGL_AXES_

/**
 * An OpenGL set of axes.
 */
class Axes {

	glm::vec3 origin;
	glm::vec3 extents;
	bool includeNeg;

	float endcapSize = 0.075;
	glm::vec4 xcol = glm::vec4(1.0f, 0.0f, 0.0f, 0.15f);
	glm::vec4 ycol = glm::vec4(0.0f, 1.0f, 0.0f, 0.15f);
	glm::vec4 zcol = glm::vec4(0.25f, 0.25f, 1.0f, 0.4f);
	glm::vec4 negxcol = glm::vec4(0.5f, 0.0f, 0.0f, 0.15f);
	glm::vec4 negycol = glm::vec4(0.0f, 0.5f, 0.0f, 0.15f);
	glm::vec4 negzcol = glm::vec4(0.15f, 0.15f, 0.5f, 0.4f);

public:

	/**
	 * Construct a set of positive axes placed with their origin at orig in world space,
	 * each axis having a length specified by their respective dimension in ex.
	 *
	 * @param orig, the world coordinates of the axes's origin.
	 * @param ex, the extents of the three axes.
	 */
	Axes(glm::vec3 orig, glm::vec3 ex) : origin(orig), extents(ex), includeNeg(false) {}


	/**
	 * Construct a set of axes placed with their origin at orig in world space,
	 * each axis having a length specified by their respective dimension in ex.
	 *
	 * @param orig, the world coordinates of the axes's origin.
	 * @param ex, the extents of the three axes.
	 * @param posAndNeg, if true, draws negative axes as well.
	 */
	Axes(glm::vec3 orig, glm::vec3 ex, bool posAndNeg) : origin(orig), extents(ex), includeNeg(posAndNeg) {}


	/**
	 * Draw the axes.
	 */
	void draw() {

		glMatrixMode( GL_MODELVIEW );
		glPushMatrix();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		// glEnable(GL_DEPTH_TEST);
		// Accept fragment if it closer to the camera than the former one
		// glDepthFunc(GL_LESS);

		glLineWidth(3.0f);
		glBegin(GL_LINES);
		glColor4f(xcol.x, xcol.y, xcol.z, xcol.a);
		glVertex3f(origin.x, origin.y, origin.z);
		glVertex3f(origin.x + extents.x, origin.y, origin.z);
		glVertex3f(origin.x + extents.x, origin.y, origin.z);
		glVertex3f(origin.x + extents.x, origin.y, origin.z+endcapSize);
		glVertex3f(origin.x + extents.x, origin.y, origin.z);
		glVertex3f(origin.x + extents.x, origin.y, origin.z-endcapSize);

		if (includeNeg) {
			glColor4f(negxcol.x, negxcol.y, negxcol.z, negxcol.a);
			glVertex3f(origin.x, origin.y, origin.z);
			glVertex3f(origin.x - extents.x, origin.y, origin.z);
			glVertex3f(origin.x - extents.x, origin.y, origin.z);
			glVertex3f(origin.x - extents.x, origin.y, origin.z+endcapSize);
			glVertex3f(origin.x - extents.x, origin.y, origin.z);
			glVertex3f(origin.x - extents.x, origin.y, origin.z-endcapSize);
		}

		glColor4f(ycol.x, ycol.y, ycol.z, ycol.a);
		glVertex3f(origin.x, origin.y, origin.z);
		glVertex3f(origin.x, origin.y + extents.y, origin.z);
		glVertex3f(origin.x, origin.y + extents.y, origin.z);
		glVertex3f(origin.x, origin.y + extents.y, origin.z+endcapSize);
		glVertex3f(origin.x, origin.y + extents.y, origin.z);
		glVertex3f(origin.x, origin.y + extents.y, origin.z-endcapSize);

		if (includeNeg) {
			glColor4f(negycol.x, negycol.y, negycol.z, negycol.a);
			glVertex3f(origin.x, origin.y, origin.z);
			glVertex3f(origin.x, origin.y - extents.y, origin.z);
			glVertex3f(origin.x, origin.y - extents.y, origin.z);
			glVertex3f(origin.x, origin.y - extents.y, origin.z+endcapSize);
			glVertex3f(origin.x, origin.y - extents.y, origin.z);
			glVertex3f(origin.x, origin.y - extents.y, origin.z-endcapSize);
		}

		glColor4f(zcol.x, zcol.y, zcol.z, zcol.a);
		glVertex3f(origin.x, origin.y, origin.z);
		glVertex3f(origin.x, origin.y, origin.z + extents.z);
		glVertex3f(origin.x, origin.y, origin.z + extents.z);
		glVertex3f(origin.x+endcapSize, origin.y, origin.z + extents.z);
		glVertex3f(origin.x, origin.y, origin.z + extents.z);
		glVertex3f(origin.x-endcapSize, origin.y, origin.z + extents.z);

		if (includeNeg) {
			glColor4f(negzcol.x, negzcol.y, negzcol.z, negzcol.a);
			glVertex3f(origin.x, origin.y, origin.z);
			glVertex3f(origin.x, origin.y, origin.z - extents.z);
			glVertex3f(origin.x, origin.y, origin.z - extents.z);
			glVertex3f(origin.x+endcapSize, origin.y, origin.z - extents.z);
			glVertex3f(origin.x, origin.y, origin.z - extents.z);
			glVertex3f(origin.x-endcapSize, origin.y, origin.z - extents.z);
		}

		glEnd();

		glPopMatrix();

		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
	}

};


#endif
