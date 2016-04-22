/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Particles.cpp
 * Author: swl
 * 
 * Created on April 15, 2016, 12:16 PM
 */

#include "Particles.h"

Particles::Particles() 
{
    poly6_h9 = 315.0/(64.0*M_PI*pow(h, 9));
    h2 = h*h;
    spiky_h6 = 45.0/(M_PI*pow(h,6));
    int nx = 10;
    int ny = 10;
    int nz = 10;
    float d = 0.1;
    for(int x=0; x<nx; x++)
    {
        for(int y=0; y<ny; y++)
        {
            for(int z=0; z<nz; z++)
            {
                Particle par;
                par.p = glm::dvec3((x+0.5-nx*0.5)*d, (y+0.5)*d-1.0, (z+0.5-nz*0.5)*d);
		par.v = glm::dvec3(0,0,0);
                particles.push_back(par);
            }
        }
    }
}

void Particles::step() {
  //Apply gravity to the particles
  for(int i = 0; i < particles.size(); i++)
  {
    particles[i].v += dt*glm::dvec3(0, -gravity, 0);
    particles[i].q = particles[i].p + dt*particles[i].v;
  }
  //Constraint solving loop
  for (int i = 0; i < solver_iterations; i++) {

    //Calculate lambda_i
    for (int j = 0; j < particles.size(); j++) {
      


    }

    //Calculate new positions and do collision detection
    for (int j = 0; j < particles.size(); j++) {
      

    }

    //Update positions
    for (int j = 0; j < particles.size(); j++) {
      particles[j].q += particles[i].delta_p;
    }
  }

  for (int i = 0; i < particles.size(); i++) {
    particles[i].v = (particles[i].q - particles[i].p)/dt;
    //TODO: Vorticity confinement and viscosity
    
    particles[i].p = particles[i].q;
  }


}

void Particles::render() const
{
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 10.0, 10.0, 10.0, 0.0 };
    glShadeModel (GL_SMOOTH);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glColor3f(0.2, 0.5, 0.8);
    glColorMaterial(GL_FRONT, GL_SPECULAR);
    glColor3f(0.9, 0.9, 0.9);
    glColorMaterial(GL_FRONT, GL_AMBIENT);
    glColor3f(0.2, 0.5, 0.8);
    
    for(const Particle &par : particles)
    {    
        
        glPushMatrix();
        glTranslatef(par.p.x, par.p.y, par.p.z);
        glutSolidSphere(0.05, 10, 10);
        glPopMatrix();
    }
    
    glPopAttrib();
}

