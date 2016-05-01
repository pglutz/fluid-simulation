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
#include <math.h>

Particles::Particles() 
{
    gravity =10.0;
    solver_iterations = 5;
    dt = .001;
    h = .4;
    rest = 400.0;
    epsilon = 0.01;
    k = .1;
    n = 4;
    q = .2*h;

    poly6_h9 = 315.0/(64.0*M_PI*pow(h, 9));
    h2 = h*h;
    spiky_h6 = 45.0/(M_PI*pow(h,6));
    W_dq = W_poly6(glm::dvec3(q));
    

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
                par.q = par.p; // initialize the tentative positions
		            par.v = glm::dvec3(0.0);
                particles.push_back(par);
            }
        }
    }
    initialize_walls(nx,ny,nz);
}

void Particles::initialize_walls(int nx, int ny, int nz) { // assumed to start at origin and constrained to move along one direction
    Walls bottom_wall;
    bottom_wall.n = dvec3(0.0,1.0,0.0);
    bottom_wall.p = dvec3(0.0,0.0,0.0);
    bottom_wall.v = dvec3(0.0,0.0,0.0);
    bottom_wall.p_prime = dvec3(nx*0.5,0,nz*0.5);
    walls.push_back(bottom_wall);

    // Walls left_wall;
    // bottom_wall.n = dvec3(1.0,0.0,0.0);
    // bottom_wall.p = dvec3(0.0,0.0,0.0);
    // bottom_wall.v = dvec3(0.0,0.0,0.0);
    // bottom_wall.p_prime = dvec3()
    // walls.push_back(left_wall);

    // Walls right_wall;
    // bottom_wall.n = dvec3(-1.0,0.0,0.0);
    // bottom_wall.p = dvec3(nx,0.0,0.0);
    // bottom_wall.v = dvec3(0.0,0.0,0.0);
    // walls.push_back(bottom_wall);

    // Walls front_wall;
    // bottom_wall.n = dvec3(0.0,0.0,1.0);
    // bottom_wall.p = dvec3(0.0,0.0,0.0);
    // bottom_wall.v = dvec3(0.0,0.0,0.0);
    // walls.push_back(front_wall);

    // Walls back_wall;
    // bottom_wall.n = dvec3(0.0,0.0,-1.0);
    // bottom_wall.p = dvec3(0.0,0.0,nz);
    // bottom_wall.v = dvec3(0.0,0.0,0.0);
    // walls.push_back(back_wall);

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
      particles[j].lambda = find_lambda(j);
    }

    //Calculate new positions and do collision detection
    for (int j = 0; j < particles.size(); j++) {
      particles[j].delta_p = find_delta_p(j);
      Particles particle_j = particles[j];
      //TODO: Collision detection
      for (int i = 0; i < walls.size(); ++i)
      {
        Walls wall_i = walls[i];
        if (ifCollision(wall_i,particle_j,delta_p)
        {
          
        }
      }

    }

    //Update positions
    for (int j = 0; j < particles.size(); j++) {
      particles[j].q += particles[j].delta_p;
    }
  }

  for (int i = 0; i < particles.size(); i++) {
    //Update velocity
    particles[i].v = (particles[i].q - particles[i].p)/dt;
    //TODO: Vorticity confinement and viscosity
    


    //Set position to predicted position
    particles[i].p = particles[i].q;
  }


}

// std::tuple<Particle, Particle> Particles::collisionResponse(Particle one, glm::dvec3 delta_p , Particle two) {
//   glm::dvec3 normal_collision = 
// }

bool Particles::ifCollision(Walls wall, Particle second, glm::dvec3 delta_p) {
  dvec3 tentative_posn = second.q + (delta_p * dt);
  double td_ref = delta_p;
  dvec3 d = glm::normalize(delta_p);
  double t_actual = delta_p/d ;
  double t = dot((wall.p_prime - second.q),wall.n)/(dot(wall.n, d) );
  if (t < 0)
    return false;

  if (t_actual > t )
  {
    
  }

}


double Particles::find_lambda(int i) {
    double C = 0.0;
    double nabla_C = 0.0;
    glm::dvec3 nabla_C_i = glm::dvec3(0.0);
    glm::dvec3 q_i = particles[i].q;

    //TODO: only sum over neighbors
    for (int j = 0; j < particles.size(); j++) {
        if (i != j) {
	    glm::dvec3 q_ij = q_i - particles[j].q;
	    C += W_poly6(q_ij);
	    glm::dvec3 delta_qj = W_spiky(q_ij);
	    nabla_C += dot(delta_qj, delta_qj);
	    nabla_C_i += delta_qj;
	  }
    }
    
    //if (C > 1000) {
    //std::cout << C << std::endl;
    //int a;
    //std::cin >> a;
    //}

    nabla_C += dot(nabla_C_i, nabla_C_i);
    nabla_C = nabla_C/(rest*rest);
    return -((C/rest) - 1.0)/(nabla_C + epsilon);
}

glm::dvec3 Particles::find_delta_p(int i) {
  double lambda_i = particles[i].lambda;
  glm::dvec3 q_i = particles[i].q;
  glm::dvec3 delta = glm::dvec3(0.0);
  double m = 1.0;

  //TODO: only sum over neighbors
  for (int j = 0; j < particles.size(); j++) {
    if (i != j) {
      //TODO: Add artificial pressure term
      glm::dvec3 q_ij = q_i - particles[j].q;
      double r2 = dot(q_ij, q_ij);
      if (r2 < h2) {
	m += 1;
	double s = -k*pow(W_poly6(q_ij)/W_dq, n);
	delta += (lambda_i + particles[j].lambda + s)*W_spiky(q_ij);
      }
    }
  }
  return (1.0/(m * rest))*delta;
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

