
#ifndef _OGL_PLANE_
#define _OGL_PLANE_


/**
 * An OpenGL square plane.
 * Renders a flat quad of a particular size,
 * with grid lines along unit intervals.
 * The plane is always centered at the origin in world cooridnates.
 */
class Plane {

public:
	enum PLANE_WHICH {
		x,
		y,
		z
	};

private:

	PLANE_WHICH plane = PLANE_WHICH::x;

	glm::vec4 color = glm::vec4(0.20f, 0.20f, 0.20f, 1.0f);

	GLfloat size;

	bool softEdge;
	bool gridOnly;


public:

	/**
	 * Construct a default plane.
	 */
	Plane() : size(30.0), softEdge(false), gridOnly(false) {}


	/**
	 * Construct a plane of a particular size.
	 * @param sz, the side length of the plane.
	 */
	Plane(GLfloat sz) : size(sz), softEdge(false) {}


	/**
	 * Construct a plane of a particular size, with the option
	 * for a "soft edge" where the edge of the plane fades to 0 opacity.
	 *
	 * @param sz, the side length of the plane.
	 * @param softEdge_in, if true, edges of the plane fade to 0 opacity.
	 */
	Plane(GLfloat sz, bool softEdge_in) : size(sz), softEdge(softEdge_in) {}


	/**
	 * Construct a plane of a particular size, with the option
	 * for a "soft edge" where the edge of the plane fades to 0 opacity
	 * as well as the option to draw only grid lines.
	 *
	 * @param sz, the side length of the plane.
	 * @param softEdge_in, if true, edges of the plane fade to 0 opacity.
	 * @param gridOnly_in, if true, only the grid lines are drawn and not the plane itself.
	 */
	Plane(GLfloat sz, bool softEdge_in, bool gridOnly_in) :
		size(sz),
		softEdge(softEdge_in),
		gridOnly(gridOnly_in)
	{

	}


	/**
	 * Draw the plane.
	 */
	void draw() {

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_POLYGON_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		if (plane == PLANE_WHICH::x) {
			if (!gridOnly) {
				glBegin(GL_QUADS);
				glColor4f(color.x, color.y, color.z, color.w);
				glVertex3f(-size, 0.0f, -size);
				glVertex3f(size, 0.0f, -size);
				glVertex3f(size, 0.0f, size);
				glVertex3f(-size, 0.0f, size);
				glEnd();
			}

			glLineWidth(0.5f);
			glBegin(GL_LINES);
			if (softEdge) {
				for (int i = -size; i < size; ++i) {
					glColor4f(color.x, color.y, color.z, (size-fabs(i))/size);
					glVertex3f(1.0f*i, 0.0f, 0.0f);
					glColor4f(color.x, color.y, color.z, 0.0f);
					glVertex3f(1.0f*i, 0.0f, size);

					glColor4f(color.x, color.y, color.z, (size-fabs(i))/size);
					glVertex3f(1.0f*i, 0.0f, 0.0f);
					glColor4f(color.x, color.y, color.z, 0.0f);
					glVertex3f(1.0f*i, 0.0f, -size);


					glColor4f(color.x, color.y, color.z, (size-fabs(i))/size);
					glVertex3f(0.0f, 0.0f, 1.0f*i);
					glColor4f(color.x, color.y, color.z, 0.0f);
					glVertex3f(size, 0.0f, 1.0f*i);

					glColor4f(color.x, color.y, color.z, (size-fabs(i))/size);
					glVertex3f(0.0f, 0.0f, 1.0f*i);
					glColor4f(color.x, color.y, color.z, 0.0f);
					glVertex3f(-size, 0.0f, 1.0f*i);
				}
			} else {
				for (int i = -size; i < -size/2; ++i) {
					glVertex3f(1.0f, 0.0f, -size);
					glVertex3f(1.0f, 0.0f, size);
					glVertex3f(size, 0.0f, 1.0f*i);
					glVertex3f(-size, 0.0f, 1.0f*i);
				}
				for (int i = size/2; i < size; ++i) {
					glVertex3f(1.0f*i, 0.0f, -size);
					glVertex3f(1.0f*i, 0.0f, size);
					glVertex3f(size, 0.0f, 1.0f*i);
					glVertex3f(-size, 0.0f, 1.0f*i);
				}
			}
			glEnd();

		}

		glPopMatrix();
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_POLYGON_SMOOTH);
	}

};

#endif
