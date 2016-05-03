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
#include <glm/glm.hpp>

#define PI 3.14159265
#define E 2.718281828

Particles::Particles() 
{
    //Default values
    use_grid = false;
    gravity = 10;
    solver_iterations = 10;
    dt = .01;
    h = .4;
    rest = 810.0;
    epsilon = .01;
    k = .1;
    n = 4;
    q = .2*h;
    curr_t = 0.0;
    c = .01;
    
    poly6_h9 = 315.0/(64.0*M_PI*pow(h, 9));
    h2 = h*h;
    spiky_h6 = 45.0/(M_PI*pow(h,6));
    W_dq = W_poly6(glm::dvec3(q));

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
	            par.v = glm::dvec3(0.0);
              particles.push_back(par);
            }
        }
    }
    double x_bound_lower = (0.5-nx*0.5)*d -0.5;
    double z_bound = (0.5 - nz*0.5)*d -0.5;
    initialize_walls(x_bound_lower,ny,z_bound,d);
}

void Particles::initialize_walls(int nx, int ny, int nz, double d) {// constrained to move along one direction
    
    Walls bottom_wall;
    bottom_wall.n = glm::dvec3(0.0,1.0,0.0);
    // bottom_wall.p = glm::dvec3(0.0,-1.0,0.0);
    bottom_wall.p = -1.0;
    // bottom_wall.v = glm::dvec3(0.0,0.0,0.0);
    bottom_wall.ind = 1;
    bottom_wall.degrees = 0.0;
    walls.push_back(bottom_wall);


    Walls left_wall;
    left_wall.n = glm::dvec3(0.0,0.0,1.0);
    // left_wall.p = glm::dvec3(-1.0,0.0,0.0);
    left_wall.p = nz - 1.0;
    // left_wall.v = glm::dvec3(0.0,0.0,0.0);
    left_wall.degrees = 0.0;
    left_wall.ind = 2;
    walls.push_back(left_wall);


    Walls right_wall;
    right_wall.n = glm::dvec3(0.0,0.0,-1.0);
    // right_wall.p = glm::dvec3(nx*2.0,0.0,0.0);
    right_wall.p = nz +1.0;
    // right_wall.v = glm::dvec3(0.0,0.0,0.0);
    right_wall.ind = 2;
    right_wall.degrees = 0.0;
    walls.push_back(right_wall);


    Walls front_wall;
    front_wall.n = glm::dvec3(1.0,0.0,0.0);
    // front_wall.p = glm::dvec3(0.0,0.0,0.0);
    front_wall.p =nx - 1.0;
    // front_wall.v = glm::dvec3(0.0,0.0,0.0);
    front_wall.ind = 0;
    walls.push_back(front_wall);


    Walls back_wall;
    back_wall.n = glm::dvec3(-1.0,0.0,0.0);
    // back_wall.p = glm::dvec3(0.0,0.0,nz*2);
    back_wall.p = nx + 1.0;
    // back_wall.v = glm::dvec3(0.0,0.0,0.0);
    back_wall.ind = 0;
    walls.push_back(back_wall);

}

void Particles::set_h(double new_h) {
    h = new_h;
    q = .2*h;
    
    poly6_h9 = 315.0/(64.0*M_PI*pow(h, 9));
    h2 = h*h;
    spiky_h6 = 45.0/(M_PI*pow(h,6));
    
    W_dq = W_poly6(glm::dvec3(q));
}

void Particles::step() {

  // Moving the second wall (Simple Harmonic Motion)
  if (curr_t > 7) {
  double curr_deg = walls[1].degrees;
  walls[1].degrees = curr_deg + 20;
  walls[1].p = walls[1].p + (walls[1].p*cos( curr_deg * PI / 180.0 ))*0.1;//(exp(-0.01*curr_t)); // add something that damps the oscillation with time
  }

  //Apply gravity to the particles
  for (int i = 0; i < particles.size(); i++)
  {
    particles[i].v += dt*glm::dvec3(0.0, -gravity, 0.0);
    particles[i].q = particles[i].p + dt*particles[i].v;
  }

  //Hash particles to grid
  /*
  if (use_grid) {
    grid.clear();
    for (int i = 0; i < particles.size(); i++) {
      std::string key = round_to_str(particles[i].q);
      grid[key].push_back(particles[i]);
    }
  }
  */
  
  //Constraint solving loop
  for (int i = 0; i < solver_iterations; i++) {

    //Calculate lambda_i
    for (int j = 0; j < particles.size(); j++) {
      particles[j].lambda = find_lambda(j);
    }


    //Calculate new positions and do collision detection
    for (int j = 0; j < particles.size(); j++) {
      particles[j].delta_p = find_delta_p(j);

      int index;
      //TODO: Collision detection
      for (int i = 0; i < walls.size(); ++i)
      {
        Walls wall_i = walls[i];
        index = wall_i.ind;
        if ( wall_i.n[index] > 0 )
        {
          if (particles[j].q[index] < wall_i.p) {
           particles[j].q[index] = wall_i.p;
           reflect(particles[j],wall_i);
          }
        }
        else {    //negative normal
          if (particles[j].q[index] > wall_i.p) {
            particles[j].q[index] = wall_i.p;
            reflect(particles[j],wall_i);
          }
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
    //TODO: Vorticity confinement
    

    particles[i].delta_v = viscosity(i);
    //Set position to predicted position
    particles[i].p = particles[i].q;
  }
  for (int i = 0; i < particles.size(); i++) {
    particles[i].v += particles[i].delta_v;
  }

  curr_t += dt;

}

void Particles::reflect(Particle part, Walls wall) {
  glm::dvec3 w_i = part.p - part.q;
  glm::dvec3 w_o = 2.0*dot(w_i,wall.n)*wall.n - w_i;
  part.q = part.q + w_o;
}

bool Particles::ifCollision(Walls wall, Particle second, glm::dvec3 delta_p) {
  glm::dvec3 tentative_posn = second.q + (delta_p * dt);
  glm::dvec3 d = glm::normalize(delta_p);
  double t_actual = dot(d,delta_p);
  double t = 1.0; // dot((wall.p_prime - second.q),wall.n)/(dot(wall.n, d) );
  if (t < 0)
    return false;

  if (t_actual > t )
  {
    second.q = second.q + (t*d);  // return with corrected q value
  }
  else {
    second.q = second.q + (t_actual*d);
  }
  return true;

}


double Particles::find_lambda(int i) {
    double C = 0.0;
    double nabla_C = 0.0;
    glm::dvec3 nabla_C_i = glm::dvec3(0.0);
    glm::dvec3 q_i = particles[i].q;
    particles[i].d_W.clear();
    particles[i].W.clear();

    for (int j = 0; j < particles.size(); j++) {
      glm::dvec3 q_j = particles[j].q;
      if (q_i != q_j) {
	glm::dvec3 q_ij = q_i - q_j;
	double r2 = dot(q_ij, q_ij);
	if (r2 < h2) {
	  //Calculating density constraint
	  double W = W_poly6(q_ij);
	  C += W;

	  //Store computed value of W_poly6(q_i - q_j)
	  particles[i].W.push_back(W);

	  //Calculating gradient
	  glm::dvec3 delta_qj = W_spiky(q_ij);
	  nabla_C += dot(delta_qj, delta_qj);
	  nabla_C_i += delta_qj;

	  //Store computed value of W_spiky(q_i - q_j)
	  particles[i].d_W.push_back(delta_qj);

	} else {
	  particles[i].W.push_back(0);
	}
      } else {
	particles[i].W.push_back(0);
      }
    }  

      /*  
    } else {
      for (int t = -2; t <= 2; t++) {
        for (int u = -2; u <= 2; u++) {
	  for (int v = -2; v <= 2; v++) {
	    std::vector<Particle> list = grid[round_to_str(q_i + glm::dvec3(t*h, u*h, v*h))];
	    for (int j = 0; j < list.size(); j++) {
	      glm::dvec3 q_j = list[j].q;
	      if (q_i != q_j) {
		glm::dvec3 q_ij = q_i - q_j;
		double r2 = dot(q_ij, q_ij);
		if (r2 < h2) {
		  C += W_poly6(q_ij);
		  glm::dvec3 delta_qj = W_spiky(q_ij);
		  nabla_C += dot(delta_qj, delta_qj);
		  nabla_C_i += delta_qj;
		}
	      }
	    }
	  }
	}
      }
    }
      */

    nabla_C += dot(nabla_C_i, nabla_C_i);
    nabla_C = nabla_C/(rest*rest);
    return -((C/rest) - 1.0)/(nabla_C + epsilon);
}

glm::dvec3 Particles::find_delta_p(int i) {
  double lambda_i = particles[i].lambda;
  glm::dvec3 delta = glm::dvec3(0.0);
  double m = 1.0;

  for (int j = 0; j < particles.size(); j++) {
    double W = particles[i].W[j];
    if (W > 0) {
      //Artificial pressure term
      double s = -k*pow(W/W_dq, n);
      delta += (lambda_i + particles[j].lambda + s)*particles[i].d_W[(int) m - 1];
      m += 1;
    }
  }
  /*
  } else {
    for (int t = -2; t <= 2; t++) {
      for (int u = -2; u <= 2; u++) {
	for (int v = -2; v <= 2; v++) {
	  std::vector<Particle> list = grid[round_to_str(q_i + glm::dvec3(t*h, u*h, v*h))];
	  for (int j = 0; j < list.size(); j++) {
	    glm::dvec3 q_j = list[j].q;
	    if (q_i != q_j) {
	      glm::dvec3 q_ij = q_i - q_j;
	      double r2 = dot(q_ij, q_ij);
	      if (r2 < h2) {
		m += 1;
		double s = -k*pow(W_poly6(q_ij)/W_dq, n);
		delta += (lambda_i + list[j].lambda + s)*W_spiky(q_ij);
	      }
	    }
	  }
	}
      }
    }
  }
  */    

  return (1.0/(m*rest))*delta;
}

glm::dvec3 Particles::viscosity(int i) {
  glm::dvec3 v_i = particles[i].v;
  glm::dvec3 delta_v = glm::dvec3(0.0);
  int m = 1.0;
  for (int j = 0; j < particles.size(); j++) {
    double W = particles[i].W[j];
    if (W > 0) {
      m += 1;
      glm::dvec3 v_ij = particles[j].v - v_i;
      delta_v += W*v_ij;
    }
  }
  return (c/m)*delta_v;
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

