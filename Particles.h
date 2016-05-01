/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Particles.h
 * Author: swl
 *
 * Created on April 15, 2016, 12:16 PM
 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <glm/glm.hpp>
#include <vector>
#if defined(__APPLE_CC__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <math.h>
#endif

#include <iostream>
using namespace std;

class Particles {
public:
    Particles();
    void render() const;
    void step(); // simulate one frame
    bool ifCollision(glm::dvec3 first_posn, Particle second);
    void initialize_walls(int nx,int ny, int nz);

private:

    double W_poly6(glm::dvec3 r) {
      double r2 = dot(r, r);
      if (r2 > h2) {
	return 0.0;
      }
      return poly6_h9*pow((h2 - r2), 3);
    }

    glm::dvec3 W_spiky(glm::dvec3 r) {
      double norm_r = sqrt(dot(r, r));
      if (norm_r > h) {
	return glm::dvec3(0.0);
      }
      return -(spiky_h6*(h - norm_r)*(h - norm_r)/norm_r)*r;
    }

    double find_lambda(int i);
    glm::dvec3 find_delta_p(int i);

    struct Particle
    {
      //Particle position
      glm::dvec3 p;
      //Particle velocity
      glm::dvec3 v;
      //Temporary position
      glm::dvec3 q;
      //Lambda for constraint solving loop
      double lambda;
      //Position change
      glm::dvec3 delta_p;
    };

    struct Walls
    {
      //Wall normal
      glm::dvec3 n;
      //Wall velocity
      glm::dvec3 v;
      //Wall position
      glm::dvec3 p;
      // Wall center point 
      glm::dvec3 p_prime;
    };
    
    std::vector<Particle> particles;
    std::vector<Walls> walls; // boundaries
    double gravity;
    int solver_iterations;
    double dt;
    double h;
    double rest;
    double epsilon;
    double k;
    double q;
    double W_dq;
    double n;

    double poly6_h9;
    double h2;
    double spiky_h6;
};

#endif /* PARTICLES_H */

