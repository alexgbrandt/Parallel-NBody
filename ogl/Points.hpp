
#ifndef _OGL_POINTS_
#define _OGL_POINTS_

#include <math.h>
#include "shader.hpp"

/**
 * A set of (continuously updating) points,
 * rendered in OpenGL, where each is represented
 * in fact by a square with a 2D Gaussian alpha component.
 *
 * Requires constants to be setup already: OGL_SCREEN_W, OGL_SCREEN_H
 */
class Points {

private:

	long nPoints;
	GLdouble* vertex_data;
	GLfloat* color_data;
	bool init;
	bool needsUpdate;

	GLuint VertexArrayID;
	GLuint programID;
	GLuint MatrixID;
	GLuint FlareSizeID;

	GLuint colorbuffer;
	GLuint vertexbuffer;

	GLuint pointTex;
	GLuint pointTexID = 0;
	int textureSize = 32;


	void copyToBuffer() {
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLdouble)*3*nPoints, vertex_data);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat)*4*nPoints, color_data);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		needsUpdate = false;
	}

	void setup(const GLdouble* vertex_data) {
		glGenVertexArrays(1, &VertexArrayID);
		glBindVertexArray(VertexArrayID);

		programID = LoadShaders( "shaders/PointShader.vertexshader", "shaders/PointShader.geometryshader", "shaders/PointShader.fragmentshader" );

		glUseProgram(programID);

		MatrixID = glGetUniformLocation(programID, "MVP");
		FlareSizeID = glGetUniformLocation(programID, "flare_size");
		GLfloat flareSize[2] = { ((float) textureSize) / (2.0f*OGL_SCREEN_W), ((float) textureSize) / (2.0f*OGL_SCREEN_H)};
		glUniform2f(FlareSizeID, flareSize[0], flareSize[1]);

	 	GLuint TextureID  = glGetUniformLocation(programID, "tex");
		glUniform1i(TextureID, pointTexID);
		glUseProgram(0);

		glGenBuffers(1, &vertexbuffer);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLdouble)*3*nPoints, vertex_data, GL_DYNAMIC_DRAW);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(
			0,                  // position attrib
			3,                  // size
			GL_DOUBLE,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// createColors();
		glGenBuffers(1, &colorbuffer);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*4*nPoints, color_data, GL_DYNAMIC_DRAW);

		// 2nd attribute buffer : colors
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(
			1,                  // color attrib
			4,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);


		//Create the point texture now;
		glGenTextures(1, &pointTex);
		glBindTexture(GL_TEXTURE_2D, pointTex);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, textureSize, textureSize, 0, GL_RED, GL_FLOAT, 0);
		// glCreateTextures(GL_TEXTURE_2D, 1, &pointTex);
		// glTextureStorage2D(pointTex, 1, GL_R32F, textureSize, textureSize);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		GLfloat texturePixels[textureSize*textureSize];
		//generate a discrete 2D Gaussian
		float sigma2 = textureSize * 0.5f;
		float A = 1.0f;
		for (int i = 0; i < textureSize; ++i) {
			float i1 = i-textureSize / 2.0f;
			for (int j = 0; j < textureSize; ++j) {
				float j1 = j - textureSize / 2.0f;
				texturePixels[i*textureSize + j] =
					pow(A*exp(-1.0f*((i1*i1)/(2*sigma2) + (j1*j1)/(2*sigma2))), 2.2);
			}
		}

		//fill server side texture data
		glTexSubImage2D(
			GL_TEXTURE_2D,		//target tex
			0, 					//level
			0,					//xoffset
			0,					//yoffset
			textureSize,		//width
			textureSize,		//height
			GL_RED, 			//texture format
			GL_FLOAT, 			//pixel data type
			texturePixels		//the data
		);
		glBindTexture(GL_TEXTURE_2D, 0);

		init = true;
	}

	// GLfloat emission[4] = {0.8f, 0.8f, 0.8f, 1.0f};
	// GLfloat specular[4] = {1.0f, 0.2f, 0.2f, 1.0f};


public:

	/**
	 * Default constructor for no points.
	 */
	Points() : nPoints(0), vertex_data(NULL), color_data(NULL), init(false), needsUpdate(false) {

	}

	/**
	 * Create a set of N points where the i'th point has position
	 * (pointData[3*i], pointData[3*i + 1], pointData[3*i+2]) in world coordinates.
	 *
	 * @param N, the number of points.
	 * @param pointData, A 3*N array of positional data.
	 */
	Points(long N, const GLdouble* pointData, const GLfloat* colorData) : nPoints(N),/* color_data(NULL),*/ init(false), needsUpdate(true) {
		vertex_data = (GLdouble*) malloc(sizeof(GLdouble)*N*3);
		memcpy(vertex_data, pointData, sizeof(GLdouble)*3*N);
		color_data = (GLfloat*) malloc(sizeof(GLfloat)*N*4);
		memcpy(color_data, colorData, sizeof(GLfloat)*4*N);
	}

	/**
	 * Move assignment.
	 * @param other, the other Points to move from.
	 */
	Points& operator=(Points&& other) {
		nPoints = other.nPoints;
		vertex_data = other.vertex_data;
		color_data = other.color_data;

		other.nPoints = 0;
		other.vertex_data = NULL;
		other.color_data = NULL;
		return *this;
	}

	/**
	 * Destructor.
	 */
	~Points() {
		free(vertex_data);
		free(color_data);
	}


	/**
	 * Get the number of points.
	 *
	 * @return the number of points.
	 */
	long numPoints() {
		return nPoints;
	}


	/**
	 * Determine if the current positional data has been used yet to render.
	 * This helps control copying where the generation of positions
	 * is faster than rendering.
	 *
	 * @returns true iff the current positional data has been rendered.
	 */
	bool pointsStale() {
		return !needsUpdate;
	}


	/**
	 * Set the vertex positions to pointData.
	 * pointData has same meaning as in the constructor.
	 *
	 * @param N, the number of points
	 * @param pointData, a 3*N array of positional data.
	 * @param colorData, a 4*N array of color data.
	 */
	void setVertexData(long N, const GLdouble* pointData, const GLfloat* colorData) {
		if (nPoints == 0 || N != nPoints) {
			return;
		}

		memcpy(vertex_data, pointData, sizeof(GLdouble)*nPoints*3);
		memcpy(color_data, colorData, sizeof(GLfloat)*nPoints*4);
		needsUpdate = true;
	}


	/**
	 * Draw the points.
	 * @param VP, the current Projection * View matrix.
	 */
	void draw(const glm::mat4& VP) {
		if (nPoints == 0) {
			return;
		}
		if (!init) {
			setup(vertex_data);
		}
		if (needsUpdate) {
			copyToBuffer();
		}

		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE);

		glBindVertexArray(VertexArrayID);

		glUseProgram(programID);
		// we assume identity model matrix here
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &VP[0][0]);

		glActiveTexture(GL_TEXTURE0);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, pointTex);
		// glBindTextureUnit(pointTexID, pointTex);

		glDrawArrays(GL_POINTS, 0, nPoints);

		glBindTexture(GL_TEXTURE_2D, 0);
		glUseProgram(0);
		glBindVertexArray(0);
		glDisable(GL_BLEND);
	}

};


#endif
