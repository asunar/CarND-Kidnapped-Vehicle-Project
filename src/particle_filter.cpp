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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 10;
	std::default_random_engine gen;

	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);

	for(int i=0; i < num_particles; i++){
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
		weights.push_back(1);

	std::cout << particle.id << endl;
	std::cout << particle.x << endl;
	std::cout << particle.y << endl;
	std::cout << particle.theta << endl;
	std::cout << particle.weight << endl;		
	}



	is_initialized = true;
	std::cout << "Initialized" << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	for(int i = 0; i < num_particles; i++){
		Particle p = particles[i];

		double new_x = p.x;
		double new_y = p.y;
		double new_theta = p.theta;	

		if(fabs(yaw_rate) < 0.00001){
			new_x = p.x + velocity * delta_t * cos(p.theta); 
			new_y = p.y + velocity * delta_t * sin(p.theta); 
			new_theta = p.theta;
		} else {
			new_x = p.x + velocity/yaw_rate * (sin(p.theta + (yaw_rate * delta_t)) - sin(p.theta));
			new_y = p.y + velocity/yaw_rate * (cos(p.theta) - cos(p.theta + (yaw_rate * delta_t)));
			new_theta = p.theta + yaw_rate * delta_t;
		}

	    std::normal_distribution<double> noise_x(0, std_pos[0]);
	    std::normal_distribution<double> noise_y(0, std_pos[1]);
	    std::normal_distribution<double> noise_psi(0, std_pos[2]);
		
		particles[i].x = new_x + noise_x(gen);
		particles[i].y = new_y + noise_y(gen);
		particles[i].theta  = new_theta + noise_psi(gen);	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	double nearest_dist,o_x,o_y,pred_x,pred_y;
	for (int i = 0;i < observations.size();i++)
	{
		nearest_dist = std::numeric_limits<double>::max();
		o_x = observations[i].x;
		o_y = observations[i].y;
		for (int j = 0; j < predicted.size();j++)
		{
			pred_x = predicted[j].x;
			pred_y = predicted[j].y;
			if (dist(o_x,o_y,pred_x,pred_y) < nearest_dist)
			{
				nearest_dist = dist(o_x,o_y,pred_x,pred_y);
				observations[i].id = predicted[j].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	std::vector<LandmarkObs> transformed_observations;
	std::vector<LandmarkObs> particle_observations = observations;
	std::vector<LandmarkObs> inrange_landmarks;

	weights.clear();
	for (int i = 0; i < num_particles; i++)
	{
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;
		double global_x, global_y,obs_dist,landmark_dist;
		int l_id;

		transformed_observations.clear();
		LandmarkObs transformed_observation,g_landmark;

		particles[i].weight = 1; 

		for (int j = 0; j < particle_observations.size(); j++)
		{
			global_x = (particle_observations[j].x*cos(p_theta) - particle_observations[j].y*sin(p_theta)) + p_x; // Transformation into MAP coordinate system
			global_y = (particle_observations[j].x*sin(p_theta) + particle_observations[j].y*cos(p_theta)) + p_y;

			transformed_observation.x = global_x;
			transformed_observation.y = global_y;
			transformed_observation.id = particle_observations[j].id;
			transformed_observations.push_back(transformed_observation);

			std::cout << "Transformed Observation" << endl;
			std::cout << transformed_observation.id << endl;
			std::cout << transformed_observation.x << endl;
			std::cout << transformed_observation.y << endl;
		}
		inrange_landmarks.clear();

		for (int k = 0; k < map_landmarks.landmark_list.size();k++)
		{
			landmark_dist = dist(p_x,p_y,map_landmarks.landmark_list[k].x_f,map_landmarks.landmark_list[k].y_f);
			if (landmark_dist < sensor_range) // Considering only the landmarks with in sensor range
			{
				g_landmark.x = map_landmarks.landmark_list[k].x_f;
				g_landmark.y = map_landmarks.landmark_list[k].y_f;
				g_landmark.id = map_landmarks.landmark_list[k].id_i;
				inrange_landmarks.push_back(g_landmark);
			}
		}

		dataAssociation(inrange_landmarks,transformed_observations); 

		double o_x, o_y, landmark_x, landmark_y,dx,dy;
		int o_id, landmark_id;

		for (int m = 0; m < transformed_observations.size(); m++)
		{
			o_x = transformed_observations[m].x;
			o_y = transformed_observations[m].y;
			o_id =transformed_observations[m].id;

			for (int n = 0; n < inrange_landmarks.size(); n++)
			{
				if (inrange_landmarks[n].id == o_id)
				{
					landmark_x = inrange_landmarks[n].x;
					landmark_y = inrange_landmarks[n].y;

					break; 
				}

			}
			dx = o_x - landmark_x;
			dy = o_y - landmark_y;
			
			double d1 = (dx*dx)/(2.0*std_x*std_x);
			double d2 = (dy*dy)/(2.0*std_y*std_y);
			double p = (0.5/(M_PI*std_x*std_y))*exp(-(d1 + d2)); 

			particles[i].weight = particles[i].weight * p; 
			
		}

		weights.push_back(particles[i].weight);

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;

	std::vector<Particle> newParticles;
	std::discrete_distribution<> d(weights.begin(), weights.end());

	for (int i = 0; i < num_particles; i++) {
		int idx = d(gen);
		newParticles.push_back(particles[idx]);
	}
	particles = newParticles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
