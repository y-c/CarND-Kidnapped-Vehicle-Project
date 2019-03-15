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
      num_particles = 100;  // TODO: Set the number of particles
      std::default_random_engine generator;
    
    // Create normal (Gaussian) distributions for x, y, theta
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);

    for (int i = 0; i < num_particles; i++) {
          Particle p;
        p.id = i;
        p.x = dist_x(generator);
        p.y = dist_y(generator);
        p.theta = dist_theta(generator);
        p.weight = 1.0;
          
          // Add particle to list of particles
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
      std::default_random_engine generator;
    
    // Create normal (Gaussian) distributions for x, y, theta
    std::normal_distribution<double> noise_x(0.0, std_pos[0]);
    std::normal_distribution<double> noise_y(0.0, std_pos[1]);
    std::normal_distribution<double> noise_theta(0.0, std_pos[2]);

    if (fabs(yaw_rate) < 0.0001) {
      yaw_rate = 0.0001;
    }

    for (int i = 0; i < num_particles; i++) {
          particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) 
                            - sin(particles[i].theta)) + noise_x(generator);
          particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) 
                            - cos(particles[i].theta + yaw_rate * delta_t)) + noise_y(generator);
          particles[i].theta += yaw_rate * delta_t + noise_theta(generator);
      }
  }

  void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                       vector<LandmarkObs>& observations) {
    /**
     * TODO: Find the predicted measurement that is closest to each 
     *   observed measurement and assign the observed measurement to this 
     *   particular landmark.
     * NOTE: this method will NOT be called by the grading code. But you will 
     *   probably find it useful to implement this method and use it as a helper 
     *   during the updateWeights phase.
     */

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
      double var_range = std_landmark[0] * std_landmark[0];
    double var_bearing = std_landmark[1] * std_landmark[1];
    double det_std = std_landmark[0] * std_landmark[1];
    
    for (int i = 0; i < num_particles; i++) {
        particles[i].weight = 1;
        for (auto obs : observations) {
          double obs_x_trans = obs.x * cos(particles[i].theta) 
                                   - obs.y * sin(particles[i].theta) + particles[i].x;
          double obs_y_trans = obs.x * sin(particles[i].theta) 
                                   + obs.y * cos(particles[i].theta) + particles[i].y;
            double dist_min = sensor_range;
            double landmark_x = sensor_range;
            double landmark_y = sensor_range;    
            for (auto landmark : map_landmarks.landmark_list) {
          double dist = sqrt((obs_x_trans - landmark.x_f) * (obs_x_trans - landmark.x_f) 
                           + (obs_y_trans - landmark.y_f) * (obs_y_trans - landmark.y_f));
          if (dist < dist_min) {
            dist_min = dist;
            landmark_x = landmark.x_f;
            landmark_y = landmark.y_f;
          }
        } 
        double x_diff = obs_x_trans - landmark_x;
      double y_diff = obs_y_trans - landmark_y;
      double nom = exp(-0.5*((x_diff * x_diff) / var_range + (y_diff * y_diff) / var_bearing));
      double denom = 2 * M_PI * det_std;
      particles[i].weight *= nom / denom;
      } 
      weights[i] = particles[i].weight;
    }
  }

  void ParticleFilter::resample() {
    /**
     * TODO: Resample particles with replacement with probability proportional 
     *   to their weight. 
     * NOTE: You may find std::discrete_distribution helpful here.
     *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
     */
      std::default_random_engine generator;

      std::discrete_distribution<> weights_dis(weights.begin(), weights.end());
      // initialise new particle array
      vector<Particle> newParticles;
      // resample particles
      for (int i = 0; i < num_particles; i++) {
          newParticles.push_back(particles[weights_dis(generator)]);
      }
      particles = newParticles;
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
