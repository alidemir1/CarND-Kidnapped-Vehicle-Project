/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 20;  // TODO: Set the number of particles
  
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  for(int i = 0; i < num_particles; i++)
  {
    std::default_random_engine gen;
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles.push_back(p);
    weights.push_back(p.weight);
  }
  is_initialized = true;
  

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  double pred_x;
  double pred_y;
  double pred_theta;
  std::default_random_engine gen;
  
  for(unsigned int i = 0; i < particles.size(); i++)
  {
    Particle p = particles[i];
    
    
    
    //if yaw_rate too small, remove the factorization by making it equal to 1
    if(fabs(yaw_rate) < 0.0001)
    {
      pred_x = p.x + velocity * cos(p.theta) * delta_t;
      pred_y = p.y + velocity * sin(p.theta) * delta_t;
    }
    else
    {
      pred_x = p.x + (velocity / yaw_rate) * (sin(p.theta + (yaw_rate * delta_t)) - sin(p.theta));
      pred_y = p.y + (velocity / yaw_rate) * (cos(p.theta) - cos(p.theta + (yaw_rate * delta_t)));
    }
    pred_theta = p.theta + (yaw_rate * delta_t);
    
    
    std::normal_distribution<double> dist_x(pred_x, std_pos[0]);
  	std::normal_distribution<double> dist_y(pred_y, std_pos[1]);
  	std::normal_distribution<double> dist_theta(pred_theta, std_pos[2]);
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);  
   
    particles[i] = p;
   
    
  }
  

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations, double sensor_range, Particle p) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for(unsigned int i = 0; i < observations.size(); i++)
  {
    double min = sensor_range * sqrt(2);
    int idd = -1;
    for(unsigned int j = 0; j < predicted.size(); j++)
    {      
      double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      
      if(distance < min)
      {
        min = distance;
        idd = predicted[j].id;
      } 
    }
    
    observations[i].id = idd;
  }
    

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  total_weights = 0.0;
  
  for(unsigned int i = 0; i < particles.size(); i++)
  {
    vector<LandmarkObs> landmarks_on_map;
    for(unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++)
    {
      if(fabs(map_landmarks.landmark_list[k].x_f - particles[i].x) < sensor_range && fabs(map_landmarks.landmark_list[k].y_f - particles[i].y) < sensor_range)
      {
        landmarks_on_map.push_back(LandmarkObs {map_landmarks.landmark_list[k].id_i, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f});
      }
    }

    vector<LandmarkObs> obs_in_map_coordinates;
    
    for(unsigned int j = 0; j < observations.size(); j++)
    {    
      LandmarkObs L;
      L.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
      L.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
      L.id = j;
      obs_in_map_coordinates.push_back(L);
    }
    

    
    ParticleFilter::dataAssociation(landmarks_on_map, obs_in_map_coordinates, sensor_range, particles[i]);
    
    double prob = 1.0;
    
    for(unsigned int o = 0; o < obs_in_map_coordinates.size(); o++)
    {
      for(unsigned int n = 0; n < landmarks_on_map.size(); n++)
      {
        if(obs_in_map_coordinates[o].id == landmarks_on_map[n].id)
        {
          prob = prob * ((1  /  (2*M_PI*std_landmark[0] * std_landmark[1])) * exp(-1*((pow(obs_in_map_coordinates[o].x - landmarks_on_map[n].x, 2) / (2*std_landmark[0]*std_landmark[0])) + (pow(obs_in_map_coordinates[o].y - landmarks_on_map[n].y, 2) / (2*std_landmark[1]*std_landmark[1])))));
        }
      }  
    }
    
    particles[i].weight = prob;
 	total_weights += prob;
    
  }
  
  for (unsigned int q = 0; q < particles.size(); q++) {
    particles[q].weight /= total_weights;
    weights[q] = particles[q].weight;
  }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  	
  	
  
	vector<Particle> resampled_particles;

	// Create a generator to be used for generating random particle index and beta value
	std::default_random_engine gen;
	
	//Generate random particle index
	std::uniform_int_distribution<int> particle_index(0, particles.size() - 1);
	
  
	int current_index = particle_index(gen);
	
	double beta = 0.0;
	
	double max_weight_2 = 2.0 * *std::max_element(weights.begin(), weights.end());

	for (unsigned int i = 0; i < particles.size(); i++) {
		std::uniform_real_distribution<double> random_weight(0.0, max_weight_2);
		beta += random_weight(gen);
      
	  while (beta > weights[current_index]) {
	    beta -= weights[current_index];
	    current_index = (current_index + 1) % particles.size();
        
	  }
	  resampled_particles.push_back(particles[current_index]);
	}
	particles = resampled_particles;
  		
  	  
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}