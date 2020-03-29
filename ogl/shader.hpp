
/*************************
 *
 * Load and link OpenGL shaders.
 * Adapted from opengl-tutorial.org.
 *
 *************************/


#ifndef _SHADER_HPP_
#define _SHADER_HPP_

/**
 * Load an link vertex and fragment shaders, returning the compiled program id.
 */
GLuint LoadShaders(const char* vertex_file_path, const char* fragment_file_path);

/**
 * Load an link vertex, geometry, and fragment shaders, returning the compiled program id.
 */
GLuint LoadShaders(const char* vertex_file_path, const char* geometry_file_path, const char* fragment_file_path);

#endif
