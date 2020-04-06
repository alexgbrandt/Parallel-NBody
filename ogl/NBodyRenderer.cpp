
#include "NBodyRenderer.hpp"

// "Move the camera" while world is stationary.
void NBodyRenderer::cameraControlsFreeRoam(glm::mat4& V) {

	static glm::vec3 position = eye;
	//exactly the angles to look at origin from (x,y,z) = (+v, +v/2, +v);
	static GLfloat theta = 1.25f*3.14159f;
	static GLfloat phi = 0.608f*3.14159f;

	static double mouseDownX;
	static double mouseDownY;
	static bool firstPress = true;
	static double lastTime = glfwGetTime();


	double dx = 0.0, dy = 0.0;
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if (state == GLFW_PRESS)
	{
		if (firstPress) {
			glfwGetCursorPos(window, &mouseDownX, &mouseDownY);
			firstPress = false;
			dx = 0.0;
			dy = 0.0;
		} else {
			double xpos, ypos;
			glfwGetCursorPos(window, &xpos, &ypos);

			dx = xpos - mouseDownX;
			dy = ypos - mouseDownY;

			mouseDownX = xpos;
			mouseDownY = ypos;
		}
	}
	if (state == GLFW_RELEASE) {
		firstPress = true;
	}

	// Avoid super jumpy motion
	// if (fabs(dx) < 10 && fabs(dy) < 10) {
		theta += 0.002f * dx;
		phi   += -0.002f * dy;
	// }

	glm::vec3 direction(
	    sin(phi) * sin(theta),
	    cos(phi),
	    sin(phi) * cos(theta)
	);
	glm::normalize(direction);

	glm::vec3 right = glm::normalize(glm::cross(direction, up));

	double currentTime = glfwGetTime();
	float speed = 1.5f;
	float deltaTime = (currentTime - lastTime);
	lastTime = currentTime;

	// Move forward
	if (glfwGetKey( window, GLFW_KEY_UP ) == GLFW_PRESS){
		position += direction * deltaTime * speed;
	}
	// Move backward
	if (glfwGetKey( window, GLFW_KEY_DOWN ) == GLFW_PRESS){
		position -= direction * deltaTime * speed;
	}
	// Strafe right
	if (glfwGetKey( window, GLFW_KEY_RIGHT ) == GLFW_PRESS){
		position += right * deltaTime * speed;
	}
	// Strafe left
	if (glfwGetKey( window, GLFW_KEY_LEFT ) == GLFW_PRESS){
		position -= right * deltaTime * speed;
	}

	eye = position;
	glm::vec3 curTarg = eye + direction;

	glm::mat4 Vtrans = glm::mat4(1.0f);
	Vtrans[3][0] = -eye.x;
	Vtrans[3][1] = -eye.y;
	Vtrans[3][2] = -eye.z;

	V = glm::mat4(1.0f);
	glm::vec3 D = glm::normalize(eye - curTarg);
	glm::vec3 R = glm::normalize(glm::cross(up, D));
	up = glm::normalize(glm::cross(D, R));

	V[0][0] = R.x;
	V[1][0] = R.y;
	V[2][0] = R.z;

	V[0][1] = up.x;
	V[1][1] = up.y;
	V[2][1] = up.z;

	V[0][2] = D.x;
	V[1][2] = D.y;
	V[2][2] = D.z;

	V = V * Vtrans;
}

// Rotate the world as if it is a globe instead of rotating oneself.
// "Move the world" while camera is stationary.
void NBodyRenderer::cameraControlsGlobe(glm::mat4& V) {

	static glm::vec3 position = eye;
	//exactly the angles to look at origin from (x,y,z) = (+v, +v/2, +v);
	static GLfloat theta = 1.25f*3.14159f;
	static GLfloat phi = 0.608f*3.14159f;

	static double mouseDownX;
	static double mouseDownY;
	static bool firstPress = true;
	static double lastTime = glfwGetTime();

	double dx = 0.0, dy = 0.0;
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if (state == GLFW_PRESS)
	{
		if (firstPress) {
			glfwGetCursorPos(window, &mouseDownX, &mouseDownY);
			firstPress = false;
			dx = 0.0;
			dy = 0.0;
		} else {
			double xpos, ypos;
			glfwGetCursorPos(window, &xpos, &ypos);

			dx = xpos - mouseDownX;
			dy = ypos - mouseDownY;

			dx *= -1; dy *= -1; //for globe movement

			mouseDownX = xpos;
			mouseDownY = ypos;
		}
	}
	if (state == GLFW_RELEASE) {
		firstPress = true;
	}

	//Avoid super jumpy motion
	// if (fabs(dx) < 10 && fabs(dy) < 10) {
		theta += 0.002f * dx;
		phi   += -0.002f * dy;
		if (theta > 2*_PI) {
			theta -= 2*_PI;
		}
		if (phi >= _PI) {
			phi = 0.9999999 * _PI;
		}
		if (phi <= -_PI) {
			phi = -0.9999999 * _PI;
		}
	// }

	glm::vec3 direction(
	    sin(phi) * sin(theta),
	    cos(phi),
	    sin(phi) * cos(theta)
	);
	glm::normalize(direction);

	float speed = 0.25f * radiusFromOrigin; //move faster further away
	double currentTime = glfwGetTime();
	float deltaTime = (currentTime - lastTime);
	lastTime = currentTime;

	// Move forward
	if (glfwGetKey( window, GLFW_KEY_UP ) == GLFW_PRESS) {
		radiusFromOrigin -= deltaTime * speed;
	}
	// Move backward
	if (glfwGetKey( window, GLFW_KEY_DOWN ) == GLFW_PRESS) {
		radiusFromOrigin += deltaTime * speed;
	}

	position = -1.0f * direction * radiusFromOrigin;
	V = glm::lookAt(position, targ, up);

}

void NBodyRenderer::setupRenderStages() {

	//FBOs and textures for blur shaders.
	glGenFramebuffers(4, blurFBOs);
	glGenTextures(4, blurRenderTexs);

	for (int k = 0; k < 2; ++k) {
		glBindTexture(GL_TEXTURE_2D, blurRenderTexs[k]);
		glTexImage2D(
			GL_TEXTURE_2D,
			0, //mipmap level
			GL_RGBA16F, //image format
			(OGL_SCREEN_W+FBO_MARGIN) / blurDownScale,  //width
			(OGL_SCREEN_H+FBO_MARGIN) / blurDownScale, //height
			0, //border always 09
			GL_RGBA, //texture type
			GL_FLOAT, //data type
			0 //pointer
		);
	}
	for (int k = 2; k < 4; ++k) {
		glBindTexture(GL_TEXTURE_2D, blurRenderTexs[k]);
		glTexImage2D(
			GL_TEXTURE_2D,
			0, //mipmap level
			GL_RGBA16F, //image format
			(OGL_SCREEN_W+FBO_MARGIN),  //width
			(OGL_SCREEN_H+FBO_MARGIN),  //height
			0, //border always 09
			GL_RGBA, //texture type
			GL_FLOAT, //data type
			0 //pointer
		);
	}

	for (int k = 0; k < 4; ++k) {
		glBindTexture(GL_TEXTURE_2D, blurRenderTexs[k]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindFramebuffer(GL_FRAMEBUFFER, blurFBOs[k]);
		glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, blurRenderTexs[k], 0);
	}
	glBindTexture(GL_TEXTURE_2D, 0);


	//FBO and texture for forward pass.
	glGenFramebuffers(1, &forwardFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, forwardFBO);

	glGenTextures(1, &forwardRenderTex);
	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, forwardRenderTex);
	glTexImage2DMultisample(
		GL_TEXTURE_2D_MULTISAMPLE,
		16, //numsamples
		GL_RGBA16F, //image format
		OGL_SCREEN_W+FBO_MARGIN, //width
		OGL_SCREEN_H+FBO_MARGIN, //height
		GL_FALSE //fixed sampling
	);

	// glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	// glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	// glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, forwardRenderTex, 0);
	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(status != GL_FRAMEBUFFER_COMPLETE) {
		fprintf(stderr, "NBodyRenderer: failed to create FBO\n");
		exit(ALLOC_ERROR);
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);


	//FBO and texture for blit then render to texture.
	glGenFramebuffers(1, &blitFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, blitFBO);

	glGenTextures(1, &blitTex);
	glBindTexture(GL_TEXTURE_2D, blitTex);
	glTexImage2D(
		GL_TEXTURE_2D,
		0, //mipmap level
		GL_RGBA16F, //image format
		(OGL_SCREEN_W+FBO_MARGIN),  //width
		(OGL_SCREEN_H+FBO_MARGIN),  //height
		0, //border always 0
		GL_RGBA, //texture type
		GL_FLOAT, //data type
		0 //pointer
	);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, blitTex, 0);

	status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(status != GL_FRAMEBUFFER_COMPLETE) {
		fprintf(stderr, "NBodyRenderer: failed to create FBO\n");
		exit(ALLOC_ERROR);
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);


	//Deferred blur shader
	blur1ProgramID = LoadShaders( "shaders/DeferredShader.vertexshader", "shaders/BlurShader.fragmentshader" );
	blur2ProgramID = LoadShaders( "shaders/DeferredShader.vertexshader", "shaders/BlurShader2.fragmentshader" );
	glUseProgram(blur1ProgramID);
	GLuint blurTexLoc = glGetUniformLocation(blur1ProgramID, "tex");
	glUniform1i(blurTexLoc, BlurTexID);
	GLuint BlurSizeLoc = glGetUniformLocation(blur1ProgramID, "size");
	glUniform2f(BlurSizeLoc, blurDownScale/(float) OGL_SCREEN_W, blurDownScale/(float) OGL_SCREEN_H);
	Blur1HorizLoc = glGetUniformLocation(blur1ProgramID, "mult");

	glUseProgram(blur2ProgramID);
	GLuint blur2TexLoc = glGetUniformLocation(blur2ProgramID, "tex");
	glUniform1i(blur2TexLoc, BlurTexID);
	GLuint blurOffsetScaleLoc = glGetUniformLocation(blur2ProgramID, "texOffset");
	glUniform2f(blurOffsetScaleLoc, 1.0f/(float) OGL_SCREEN_W, 1.0f/(float) OGL_SCREEN_H);
	GLuint Blur2SizeLoc = glGetUniformLocation(blur2ProgramID, "size");
	glUniform1i(Blur2SizeLoc, 15);
	Blur2HorizLoc = glGetUniformLocation(blur2ProgramID, "horizontalPass");
	GLuint BlurSigmaLoc = glGetUniformLocation(blur2ProgramID, "sigma");
	glUniform1f(BlurSigmaLoc, 2.f);

	glUseProgram(0);

	//FBO and texture for luminance
	// float base_width = OGL_SCREEN_W + 2.0;
	// float base_height = OGL_SCREEN_H + 2.0;
	// glGenFramebuffers(1, &luminanceFBO);
	// glBindFramebuffer(GL_FRAMEBUFFER, luminanceFBO);
	// glCreateTextures(GL_TEXTURE_2D, 1, &luminanceTex);
	// int lumLod = (int)floor(log2(max(base_width,base_height)/2));
	// glTextureStorage2D(
	// 	luminanceTex,
	// 	lumLod+1,
	// 	GL_R16F,
	// 	OGL_SCREEN_W/2.0,
	// 	OGL_SCREEN_H/2.0
	// );
	// glTextureParameteri(luminanceTex, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	// glTextureParameteri(luminanceTex, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	// glTextureParameteri(luminanceTex, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	// glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, luminanceTex, 0);
	// glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// // Luminance shader
	// luminanceProgramID = LoadShaders( "shaders/DeferredShader.vertexshader", "shaders/LuminanceShader.fragmentshader" );
	// glUseProgram(luminanceProgramID);
	// GLuint luminanceTexLoc = glGetUniformLocation(luminanceProgramID, "tex");
	// glUniform1i(luminanceTexLoc, luminanceTexID);
	// glUseProgram(0);

	glGenVertexArrays(1, &deferredVAOID);
	glBindVertexArray(deferredVAOID);
	glGenBuffers(1, &deferredVBO);
	glBindBuffer(GL_ARRAY_BUFFER, deferredVBO);
	GLfloat deferred_quad_vertex_data[] = {
        -1.0f, -1.0f, 0.0f,
         1.0f, -1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
         1.0f, -1.0f, 0.0f,
         1.0f,  1.0f, 0.0f
	};
	glBufferData(GL_ARRAY_BUFFER, sizeof(deferred_quad_vertex_data), deferred_quad_vertex_data, GL_STATIC_DRAW);

	//Deferred shader
	deferredProgramID = LoadShaders("shaders/DeferredShader.vertexshader", "shaders/DeferredShader.fragmentshader");
	glUseProgram(deferredProgramID);
	GLuint deferredTexLoc = glGetUniformLocation(deferredProgramID, "renderedTexture");
	glUniform1i(deferredTexLoc, forwardTexID);
	GLuint bloomTexLoc = glGetUniformLocation(deferredProgramID, "bloom");
	glUniform1i(bloomTexLoc, bloomTexID);
	GLuint bloom2TexLoc = glGetUniformLocation(deferredProgramID, "bloom2");
	glUniform1i(bloom2TexLoc, bloom2TexID);
	// GLuint lodID = glGetUniformLocation(deferredProgramID, "lum_lod");
	// glUniform1i(lodID, lumLod);
	// GLuint lumTexLoc = glGetUniformLocation(deferredProgramID, "lum");
	// glUniform1i(lumTexLoc, avgLumTexID);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(
		0,                  // position attrib
		3,                  // size per vertex
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
	);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

}


int NBodyRenderer::OGL_setup() {
	// Initialise GLFW
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	// glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow( OGL_SCREEN_W, OGL_SCREEN_H, "N-Body", NULL, NULL);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window.\n");
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	if (vsync) {
		glfwSwapInterval(2);
	}

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	// glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark blue background
	// glClearColor(0.1f, 0.1f, 0.2f, 1.0f);
	// Black background
	// If bg not fully black, deferred fragment shader will brighten bg color greatly.
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	Projection = glm::perspective(glm::radians(45.0f), ((float) OGL_SCREEN_W)/OGL_SCREEN_H, 0.001f, 10000.0f);
	V = glm::lookAt(eye, targ, up);

	GLenum err;
	while ((err= glGetError()) != GL_NO_ERROR)
	{
		fprintf(stderr, "before setup: %d\n", err);
	}
	setupRenderStages();
	while((err = glGetError()) != GL_NO_ERROR)
	{
		fprintf(stderr, "after setup: %d\n", err);
	}

	this->init = true;
	return 0;
}


int NBodyRenderer::draw() {
	if (close) {
		return close;
	}

	if (!init) {
		OGL_setup();
	}

	if (glfwGetKey(window, GLFW_KEY_ESCAPE ) == GLFW_PRESS || glfwWindowShouldClose(window)) {
		close = true;
		return close;
	}

	if ( CAMERA_CONTROL ) {
		cameraControlsGlobe(V);
	}

	////////////////////////////////
	//
	// Draw process is multipass:
	//   1/ Draw the points to a multi-samples texture, then blit
	//   2/ First blur pass computes a wide aura of light about each point
	//   3/ Second blur pass computes a gaussian/bloom effect about each point
	//   4/ Second forward pass draws in the axes and grid (multisampled), then blit
	//   5/ Final deferred pass blends the final forward blit and the two blurs
	//
	////////////////////////////////

	// Measure speed
	// double currentTime = glfwGetTime();
	// if ( currentTime - lastTime >= 1.0 ) {
	// 	fprintf(stderr, "%g fps\n", double(framecount) / (currentTime - lastTime));
	// 	framecount = 0;
	// 	lastTime = currentTime;
	// }

	GLenum err;

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(glm::value_ptr(Projection));

	glMatrixMode( GL_MODELVIEW );
	glPushMatrix();
	glLoadMatrixf(glm::value_ptr(V)); //load view only so that axes are at origin.


	/////////////////////////////////////////////////////////////////
	// Forward pass 1
	/////////////////////////////////////////////////////////////////

	glViewport(0,0, OGL_SCREEN_W+FBO_MARGIN, OGL_SCREEN_H+FBO_MARGIN);
	glBindFramebuffer(GL_FRAMEBUFFER, forwardFBO);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	points.draw(Projection * V);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	while((err = glGetError()) != GL_NO_ERROR)
	{
		fprintf(stderr, "Forward Pass 1 Error: %d\n", err);
	}


	/////////////////////////////////////////////////////////////////
	// Blit 1
	/////////////////////////////////////////////////////////////////
	glBindFramebuffer(GL_READ_FRAMEBUFFER, forwardFBO);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, blitFBO);
	glBlitFramebuffer(
		0, 0, OGL_SCREEN_W+FBO_MARGIN, OGL_SCREEN_H+FBO_MARGIN,
		0, 0, OGL_SCREEN_W+FBO_MARGIN, OGL_SCREEN_H+FBO_MARGIN,
		GL_COLOR_BUFFER_BIT, GL_NEAREST);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);


	/////////////////////////////////////////////////////////////////
	// Blur pass 1
	/////////////////////////////////////////////////////////////////
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);
	// glBindTexture(GL_TEXTURE_2D, sampleTex);

	glBindVertexArray(deferredVAOID);
	glViewport(0, 0, (OGL_SCREEN_W+FBO_MARGIN)/blurDownScale, (OGL_SCREEN_H+FBO_MARGIN)/blurDownScale);
	glUseProgram(blur1ProgramID);
	int flipflop = 0;
	int maxJ = 50;
	GLuint fbo;
	GLuint blurredTex1 = blitTex;
	for (int i = 0; i < 2; ++i) {
		glUniform2f(Blur1HorizLoc, (float) ((i+1)%2), (float) i);
		for (int j = 0; j < maxJ; ++j) {
			fbo = blurFBOs[flipflop];
			glBindFramebuffer(GL_FRAMEBUFFER, fbo);
			glBindTexture(GL_TEXTURE_2D, blurredTex1);
			// glBindTextureUnit(BlurTexID, sampleTex);
			glDrawArrays(GL_TRIANGLES, 0, 6);
			blurredTex1 = blurRenderTexs[flipflop]; //the one which was just written to
			flipflop = 1 - flipflop;
		}
	}
	glBindVertexArray(0);
	glUseProgram(0);
	while((err = glGetError()) != GL_NO_ERROR)
	{
		fprintf(stderr, "Blur Pass 1 Error: %d\n", err);
	}


	/////////////////////////////////////////////////////////////////
	// Blur 2 pass
	/////////////////////////////////////////////////////////////////
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);
	glBindVertexArray(deferredVAOID);
	glViewport(0, 0, (OGL_SCREEN_W+FBO_MARGIN), (OGL_SCREEN_H+FBO_MARGIN));
	glUseProgram(blur2ProgramID);
	GLuint blurredTex2 = blitTex;
	for (int i = 0; i < 2; ++i) {
		glUniform1i(Blur2HorizLoc, i);
		GLuint fbo = blurFBOs[2 + i];
		glBindFramebuffer(GL_FRAMEBUFFER, fbo);
		glBindTexture(GL_TEXTURE_2D, blurredTex2);
		// glBindTextureUnit(BlurTexID, sampleTex);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		blurredTex2 = blurRenderTexs[2 + i]; //the one which was just written to
	}
	glDisable(GL_TEXTURE_2D);
	glBindVertexArray(0);
	glUseProgram(0);
	while((err = glGetError()) != GL_NO_ERROR)
	{
		fprintf(stderr, "Blur Pass 2 Error: %d\n", err);
	}


	/////////////////////////////////////////////////////////////////
	// luminance averaging
	/////////////////////////////////////////////////////////////////
	// glBindVertexArray(deferredVAOID);
	// glViewport(0, 0, OGL_SCREEN_W+2, OGL_SCREEN_H+2);
	// glBindFramebuffer(GL_FRAMEBUFFER, luminanceFBO);
	// glUseProgram(luminanceProgramID);

	// glActiveTexture(GL_TEXTURE0);
	// glEnable(GL_TEXTURE_2D);
	// glBindTexture(GL_TEXTURE_2D, blitTex);
	// // glBindTextureUnit(luminanceTexID, blitTex);
	// glDrawArrays(GL_TRIANGLES, 0, 6);
	// glGenerateTextureMipmap(luminanceTex);
	// glBindFramebuffer(GL_FRAMEBUFFER, 0);
	// glUseProgram(0);
	// glBindVertexArray(0);


	/////////////////////////////////////////////////////////////////
	// forward pass number 2
	/////////////////////////////////////////////////////////////////
	glViewport(0, 0, OGL_SCREEN_W+FBO_MARGIN, OGL_SCREEN_H+FBO_MARGIN);
	glBindFramebuffer(GL_FRAMEBUFFER, forwardFBO);
	if (tree != NULL) {
		tree->draw();
	}
	plane.draw();
	ax.draw();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	while((err = glGetError()) != GL_NO_ERROR)
	{
		fprintf(stderr, "Forward Pass 2 Error: %d\n", err);
	}


	/////////////////////////////////////////////////////////////////
	// Blit 2
	/////////////////////////////////////////////////////////////////
	glBindFramebuffer(GL_READ_FRAMEBUFFER, forwardFBO);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, blitFBO);
	glBlitFramebuffer(
		0, 0, OGL_SCREEN_W+FBO_MARGIN, OGL_SCREEN_H+FBO_MARGIN,
		0, 0, OGL_SCREEN_W+FBO_MARGIN, OGL_SCREEN_H+FBO_MARGIN,
		GL_COLOR_BUFFER_BIT, GL_NEAREST);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);


	/////////////////////////////////////////////////////////////////
	// render from texture
	/////////////////////////////////////////////////////////////////
	glBindVertexArray(deferredVAOID);
	glViewport(0, 0, OGL_SCREEN_W, OGL_SCREEN_H);
	glUseProgram(deferredProgramID);

	glActiveTexture(GL_TEXTURE0 + forwardTexID);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, blitTex);
	glActiveTexture(GL_TEXTURE1);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, blurredTex1); //sampleTex is one of blurRenderTexs
	glActiveTexture(GL_TEXTURE0 + bloom2TexID);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, blurredTex2);
	// glActiveTexture(GL_TEXTURE0 + avgLumTexID);
	// glEnable(GL_TEXTURE_2D);
	// glBindTexture(GL_TEXTURE_2D, luminanceTex);
	// // glBindTextureUnit(deferredTexID, blitTex);
	// // glBindTextureUnit(bloomTexID, blurRenderTexs[0]);
	// // glBindTextureUnit(avgLumTexID, luminanceTex);

	glDrawArrays(GL_TRIANGLES, 0, 6);

	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE0);
	glDisable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	glUseProgram(0);
	glBindVertexArray(0);
	while((err = glGetError()) != GL_NO_ERROR)
	{
		fprintf(stderr, "Render From Texture Error: %d\n", err);
	}

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();


	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();


	++framecount;
	return close;
}
