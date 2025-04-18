#ifndef fvw_GL_H
#define fvw_GL_H

/* fvw_GL.h - model-independent graphics routines for fni_view(1). */
/* Last edited on 2025-01-20 14:25:08 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <assert.h>

#define __APPLE__ 0
#include <GL/glu.h>

#include <bool.h>

#define fvw_GL_default_window_corner_H  64
#define fvw_GL_default_window_corner_V  64
  /* Default screen coordinates of window's top left corner. */

#define fvw_GL_default_window_HSize  800
#define fvw_GL_default_window_VSize  600
  /* Default window size. */

void fvw_GL_initialize_libraries(int *argc, char** argv);
  /* Intializes the GL library.  Also parses any X11/glut command line 
    arguments in {argv}, and deletes them.  Therefore, it should be 
    called BEFORE the program's own option parsing. */

void fvw_GL_initialize_window
  ( char *title,
    void display(void),
    void reshape(int width, int height),
    void keyboard(unsigned char key, int x, int y),
    void passivemouse(int x, int y),
    void activemouse(int x, int y),
    void special(int key, int x, int y)
  );
  /* Intializes the GL window attributes with the given title
    and methods.  Assumes that {glutInit} has been called. */

void fvw_GL_start_viewing(void);
  /* Sets some graphics state parameters and starts the GL main loop. */


#endif
