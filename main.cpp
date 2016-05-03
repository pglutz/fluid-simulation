#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include<sys/stat.h>
#include<string>

#include "Particles.h"

#include <opencv2/opencv.hpp>


inline float clip(const float& n, const float& lower, const float& upper) 
{
    return glm::max(lower, glm::min(n, upper));
}

float theta = M_PI/8;
float phi = -M_PI/8+M_PI_2;
float dist = 2.5;
int width = 800;
int height = 800;
int frame = 0;
const int render_step = 3;
int mx, my;
char cur_time[512];
bool output_animation = false;

Particles particles;

void display(void);

void reshape(int width, int height);

void idle(void)
{
    particles.step();
    glutPostRedisplay();
    if(frame/render_step >= 300)
        return;
    if(frame%render_step == 0)
    {
      if (output_animation) {
        cv::Mat3b image(height, width);
        glReadPixels(0, 0, width, height, GL_BGR, GL_UNSIGNED_BYTE, image.data);
        cv::flip(image, image, 0);
        char fn[512];
        sprintf(fn, "%s%04d.png", cur_time, frame/render_step);
        cv::imwrite(fn, image);
      }
    }
    frame++;
}

void mouse(int button, int state, int x, int y);

void motion(int x, int y);

void keyboard(unsigned char c, int x, int y)
{
    switch(c)
    {
    case 'o' :
        break;
    }
}

int main(int argc, char** argv)
{

  //Options: -a for animation, -g for gravity, -iters for solver_iterations
  //-dt for time step, -h for smoothing radius, -r for rest density 
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-a") == 0) {
      output_animation = true;
    } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
      particles.gravity = std::stod(argv[i+1]);
    } else if (strcmp(argv[i], "-iters") == 0 && i + 1 < argc) {
      particles.solver_iterations = std::stoi(argv[i+1]);
    } else if (strcmp(argv[i], "-dt") == 0 && i + 1 < argc) {
      particles.dt = std::stod(argv[i+1]);
    } else if (strcmp(argv[i], "-h") == 0 && i + 1 < argc) {
      particles.set_h(std::stod(argv[i+1]));
    } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
      particles.rest = std::stod(argv[i+1]);
    } else if (strcmp(argv[i], "-grid") == 0) {
      particles.use_grid = true;
    }
  }

  if (output_animation) {
    //Create a directory for the images based on the current time
    time_t t = time(0);
    struct tm * now = localtime( & t );
    sprintf(cur_time, "result%d-%d-%d-%d:%d:%d/", now->tm_year + 1900,
	    now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min,
	    now->tm_sec);
    mkdir(cur_time, 0755); 
  }
  
    glutInit(&argc, argv);
  
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(width, height);

    (void)glutCreateWindow("Fluid Simulation");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutMainLoop();
    glutKeyboardFunc(keyboard);

    return EXIT_SUCCESS;
}

void reshape(int w, int h)
{
    width = w;
    height = h;
    glViewport(0, 0, w, h);
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // your drawing code goes here
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90, 1, 0.01, 100);
    gluLookAt(dist*sin(phi)*cos(theta), dist*cos(phi), dist*sin(phi)*sin(theta),
            0, 0, 0, 
            0, 1, 0);
    
    particles.render();

    glutSwapBuffers();
}

void mouse(int button, int state, int x, int y)
{
    if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        mx = x;
        my = y;
    }
}

void motion(int x, int y)
{
    int dx = x - mx;
    int dy = y - my;
    mx = x;
    my = y;
    if(abs(dx) > abs(dy))
        theta += dx*0.005;
    else
        phi -= dy*0.005;
    if(theta > 2*M_PI)
        theta -= 2*M_PI;
    if(theta < -2*M_PI)
        theta += 2*M_PI;
    phi = clip(phi, M_PI/12, M_PI*11/12);
    glutPostRedisplay();
}
