/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
using namespace std;

void ParticleFilter::init( double x,  double y,  double theta,  double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 20;
  default_random_engine gen;

  // Create a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);

  // Create normal distributions for y and psi
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_psi(theta, std[2]);

  // Clear all particles and weights
  particles.clear();
  weights.clear();
  for (unsigned int i = 0; i < num_particles; ++i) {
    Particle new_particle;
    
    // Sample  and from these normal distributions like this:
    //	 sample_x = dist_x(gen);
    //	 where "gen" is the random engine initialized earlier.

    new_particle.x      = dist_x(gen);
    new_particle.y      = dist_y(gen);
    new_particle.theta  = dist_psi(gen);
    new_particle.weight = 1.0;
    new_particle.id     = i;

    // Add this new particle to set of particles.
    particles.push_back(new_particle);
    weights.push_back(1.0/num_particles);
  }
  is_initialized = true;
}

void ParticleFilter::prediction( double delta_t,  double std_pos[],  double velocity,  double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;

  // Create normal distributions for x, y and psi
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for (unsigned int i = 0; i < num_particles; ++i) {
    Particle& current_particle = particles[i];
    
      if (fabs(yaw_rate)>0.001)
      {
        double theta_new   = current_particle.theta + yaw_rate*delta_t;
        double theta_old   = current_particle.theta;
        
        current_particle.x      += velocity/yaw_rate*(sin(theta_new)-sin(theta_old));
        current_particle.y      += velocity/yaw_rate*(cos(theta_old)-cos(theta_new));
        current_particle.theta   = theta_new;
      }
      else
      {
        current_particle.x     += velocity*cos(current_particle.theta)*delta_t;
        current_particle.y     += velocity*sin(current_particle.theta)*delta_t;
      }

      // Add noise:
      current_particle.x      += dist_x(gen);
      current_particle.y      += dist_y(gen);
      current_particle.theta  += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation( std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights( double sensor_range,  double std_landmark[],
		 std::vector<LandmarkObs> observations,  Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
  double sum_prob = 0.0;
  // loop thru all particles
  for (unsigned int i = 0; i < num_particles; ++i)
  {
    double xParticle = particles[i].x;
    double yParticle = particles[i].y;
    double theta     = particles[i].theta;
    double prob            = 1.0;
    
    //Predicted Landmarks (<=sensor_range)
    vector<Map::single_landmark_s> predicted;
    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); ++j) {
      double xLand = map_landmarks.landmark_list[j].x_f;
      double yLand = map_landmarks.landmark_list[j].y_f;
      
      if (dist(xLand, yLand, xParticle, yParticle) <= sensor_range) {
        predicted.push_back(map_landmarks.landmark_list[j]);
      }
    }

    //Transform measurements from vehicle to map
    for (unsigned int j = 0; j < observations.size(); ++j)
    {
      const double xVehicle = observations[j].x;
      const double yVehicle = observations[j].y;
      double min_dist       = sensor_range;
      double x_diff         = sensor_range;
      double y_diff         = sensor_range;
      
      double xMap = xParticle + xVehicle*cos(theta) - yVehicle*sin(theta);
      double yMap = yParticle + xVehicle*sin(theta) + yVehicle*cos(theta);
      
      for (unsigned int k = 0; k < predicted.size(); ++k)
      {
        double xLandmark = predicted[k].x_f;
        double yLandmark = predicted[k].y_f;
        double cur_dist = dist(xLandmark, yLandmark, xMap, yMap);
        
        if(cur_dist < min_dist)
        {
          min_dist = cur_dist;
          x_diff = xLandmark - xMap;
          y_diff = yLandmark - yMap;
        }
      }

      // update probability wrt each observation
      prob *= (1 / (2 * M_PI * std_landmark[0] * std_landmark[1])) *
      exp(-0.5 * ((x_diff * x_diff / (std_landmark[0] * std_landmark[0])) +
                   (y_diff * y_diff / (std_landmark[1] * std_landmark[1]))));
    }
    
    particles[i].weight = prob;
    sum_prob += prob;
    }
  // update the weights vector with normalized weights
  for (unsigned int i = 0; i < num_particles; ++i) {
    weights[i] = particles[i].weight / sum_prob;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<Particle> new_particles;
  default_random_engine gen;
  discrete_distribution<int> distribution(weights.begin(), weights.end());
  for (unsigned int i = 0; i < num_particles; ++i) {
    int weighted_index = distribution(gen);
    new_particles.push_back(particles[weighted_index]);
  }
  particles = new_particles;
  
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (unsigned int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
