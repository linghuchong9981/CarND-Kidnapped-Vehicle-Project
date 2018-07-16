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
    std::cout << "ParticleFilter init start" << std::endl;
    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i = 0; i < num_particles; ++i) {
        Particle particle;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0;
        particles.push_back(particle);
    }

    std::cout << "ParticleFilter init finish" << std::endl;

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    std::cout << "ParticleFilter update start" << particles[0].x << particles[0].y <<std::endl;
    default_random_engine gen;
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);


    //TODO check是引用还是传值
    for (int i = 0; i < num_particles; ++i) {
        Particle particle = particles[i];
        if (yaw_rate > 0.0001) {
            particle.x = particle.x +
                         (velocity / yaw_rate) * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta));
            particle.y = particle.y +
                         (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t));
            particle.theta += yaw_rate * delta_t;
        } else {
            particle.x = particle.x + velocity * delta_t * cos(particle.theta);
            particle.y = particle.y + velocity * delta_t * sin(particle.theta);
        }

        double x_noise = dist_x(gen);
        double y_noise = dist_y(gen);
        double theta_noise = dist_theta(gen);

        particle.x += x_noise;
        particle.y += y_noise;
        particle.theta += theta_noise;
        particles[i] = particle;
    }

    std::cout << "ParticleFilter update finish" << particles[0].x << particles[0].y << std::endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    if (predicted.size() == 0) {
        std::cout << "dataAssociation predicted size is 0" << std::endl;
        return;
    }

    std::vector<double> distance_vec;
    for (int i = 0; i < observations.size(); ++i) {
        if (predicted.size() == 0) {
            std::cout << "dataAssociation, observations is large than predicted" << std::endl;
            break;
        }
        //TODO 没有关联上怎么处理
        for (int j = 0; j < predicted.size(); ++j) {
            double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            distance_vec.push_back(distance);
        }
        int min_pos = distance(distance_vec.begin(), min_element(distance_vec.begin(), distance_vec.end()));
//        std::cout << "dataAssociation, min pos"<< min_pos << std::endl;
        observations[i].id = predicted[min_pos].id;
//        std::cout << "dataAssociation, erase" << std::endl;

    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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


    std::cout << "updateWeights start" << std::endl;
    if (observations.size() == 0) {
        std::cout << "observations size is zero" << std::endl;
        return;
    }

    for (int i = 0; i < num_particles; ++i) {
        Particle particle = particles[i];

        //transform
        std::vector<LandmarkObs> transformed_observations;
        for (int j = 0; j < observations.size(); ++j) {
            LandmarkObs obs = observations[j];
            obs.x = particle.x + cos(particle.theta) * obs.x - sin(particle.theta) * obs.y;
            obs.y = particle.y + sin(particle.theta) * obs.x + cos(particle.theta) * obs.y;
            transformed_observations.push_back(obs);
        }

//        std::cout << "transformed_observations finish" << std::endl;

        // get predict
        std::vector<LandmarkObs> predicted;

        for (int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
            Map::single_landmark_s landmark = map_landmarks.landmark_list[k];
            double distance = dist(particle.x, particle.y, landmark.x_f, landmark.y_f);
            if (distance <= sensor_range) {
                LandmarkObs obs;
                obs.id = landmark.id_i;
                obs.x = landmark.x_f;
                obs.y = landmark.y_f;
                predicted.push_back(obs);
            }
        }

//        std::cout << "predicted finish" << std::endl;

        dataAssociation(predicted, transformed_observations);

//        std::cout << "dataAssociation finish" << std::endl;
        //update weight
        for (int l = 0; l < transformed_observations.size(); ++l) {
            LandmarkObs observation = transformed_observations[l];
            Map::single_landmark_s landmark = find_landmark_data(map_landmarks, transformed_observations[l].id);
            double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
            double exponent = pow(observation.x - landmark.x_f, 2) / (2 * pow(std_landmark[0], 2))
                              + pow(observation.y - landmark.y_f, 2) / (2 * pow(std_landmark[1], 2));
            particle.weight *= gauss_norm * exp(exponent);
        }
//        std::cout << "update weight finish" << std::endl;

        particles[i] = particle;
    }

    std::cout << "ParticleFilter update" << std::endl;

}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution



    weights.clear();
    for (int i = 0; i <particles.size() ; ++i) {
        weights.push_back(particles[i].weight);
    }

    discrete_distribution<> distribution(weights.begin(),weights.end());
    default_random_engine gen;

    std::vector<Particle> updated;
    for (int j = 0; j < num_particles; ++j) {
        updated.push_back(particles[distribution(gen)]);
    }

    particles = move(updated);
    std::cout << "ParticleFilter resame" << std::endl;
}

Particle ParticleFilter::SetAssociations(Particle &particle, const std::vector<int> &associations,
                                         const std::vector<double> &sense_x, const std::vector<double> &sense_y) {
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
