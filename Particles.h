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

#include<sstream>
#include <iostream>
using namespace std;

#include <tr1/unordered_map>
using namespace std::tr1;

class Particles {
public:
    Particles();
    void render() const;
    void step(); // simulate one frame
    struct Walls
    {
      //Wall normal
      glm::dvec3 n;
      //Wall position
      double p;
      // dimension
      int ind;
      // for moving wall
      double degrees;
    };

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
    bool ifCollision(Walls wall, Particle second, glm::dvec3 delta_p);
    // bool ifCollision_1(Walls wall, Particle second, glm::dvec3 delta_p);
    void initialize_walls(int nx,int ny, int nz,double d);


    void set_h(double new_h);
    std::vector<Walls> walls; // boundaries
    double gravity;
    int solver_iterations;
    double dt;
    double h;
    double rest;
    bool use_grid;
    unordered_map<std::string, std::vector<Particle>> grid;
    std::vector<Particle> particles;   
    double curr_t;
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

    std::string round_to_str(glm::dvec3 v) {
      std::ostringstream strs;
      strs << floor(v.x/h) << floor(v.y/h) << floor(v.z/h);
      std::string str = strs.str();
      return str;
    }

    double find_lambda(int i);
    glm::dvec3 find_delta_p(int i);


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

