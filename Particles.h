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

class Particles {
public:
    Particles();
    void render() const;
    void step(); // simulate one frame
private:

    double W_poly6(glm::dvec3 r) {
      double r2 = dot(r, r);
      return poly6_h9*pow((h2 - r2), 3);
    }

    glm::dvec3 W_spiky(glm::dvec3 r) {
      double norm_r = sqrt(dot(r, r));
      return (spiky_h6*(h - norm_r)*(h - norm_r)/norm_r)*r;
    }

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
    
    std::vector<Particle> particles;
    double gravity = .05;
    int solver_iterations = 10;
    double dt = .01;
    double h = .01;
    double poly6_h9;
    double h2;
    double spiky_h6;
};

#endif /* PARTICLES_H */

