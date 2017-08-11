#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

// Function which returns the closest way point within the given map
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{
  // Large number
	double closestLen = 100000;
	int closestWaypoint = 0;

  // Go through the way points and find the closest way point to the current
  // location of the car
	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x, y, map_x, map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}

  // Return the closest way point index
	return closestWaypoint;
}

// Function which returns the next way point
int NextWaypoint(double x, double y, double theta,
                 vector<double> maps_x, vector<double> maps_y,
                 vector<double> maps_dx, vector<double> maps_dy)
{
  // Get the closest way point
	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  // Get the corresponding co-ordinates
	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	// Heading vector
  double hx = map_x-x;
  double hy = map_y-y;

  // Normal vector:
  double nx = maps_dx[closestWaypoint];
  double ny = maps_dy[closestWaypoint];

  // Vector into the direction of the road (perpendicular to the normal vector)
  double vx = -ny;
  double vy = nx;

  // If the inner product of v and h is positive then we are behind the waypoint
  // so we do not need to increment closestWaypoint.
  // Otherwise we are beyond the waypoint and we need to increment closestWaypoint.
  double inner = hx*vx+hy*vy;
  if (inner < 0.0)
  {
    closestWaypoint++;
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y,vector<double> maps_dx, vector<double> maps_dy)
{
  // Get the next way points
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y, maps_dx, maps_dy);

  // Previous way point
  int prev_wp;
  prev_wp = next_wp - 1;

  // If the next way point is 0(circular), get the last way point in the vector
  if(next_wp == 0)
  {
    prev_wp  = maps_x.size() - 1;
  }

  // Distance in x, y
	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
  // Difference between current x,y and previous way point
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// Find the projection of x onto n
	double proj_norm = (x_x * n_x + x_y * n_y)/(n_x * n_x + n_y * n_y);
	double proj_x = proj_norm * n_x;
	double proj_y = proj_norm * n_y;

  // Get frenet d
	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  // See if d value is positive or negative by comparing it to a center point
	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  // Adjust 'd' based on center reference
	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// Calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

  // Return the 's' and 'd' co-ordinates
	return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
  // Initialize previous way point to be -1
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

  // Way point based on the previous way point
	int wp2 = (prev_wp + 1) % maps_x.size();
	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));

	// The x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);
	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading - pi()/2;
	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};
}

int main()
{
  uWS::Hub h;

  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  string map_file_ = "../data/highway_map.csv";
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line))
  {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // Start in lane 1(Middle lane)
  int lane = 1;
  int lane_change_wp = 0;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&lane_change_wp]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode)
  {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {
      auto s = hasData(data);

      if (s != "")
      {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry")
        {
          // Store the JSON data into local variables
          double car_x = j[1]["x"];
        	double car_y = j[1]["y"];
        	double car_s = j[1]["s"];
        	double car_d = j[1]["d"];
        	double car_yaw = j[1]["yaw"];
        	double car_speed = j[1]["speed"];
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];
        	vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

          // Reference velocity of 49.5 mph(below the 50mph speed limit)
        	double ref_vel = 49.5;

          // Get the previous path size
        	int prev_size = previous_path_x.size();

          // Initialize variables
        	int next_wp = -1;
        	double ref_x = car_x;
        	double ref_y = car_y;
        	double ref_yaw = deg2rad(car_yaw);

          // Check if the the size of the previous path is below 2
          // NOTE: This is the first time the path planner is going to run
        	if(prev_size < 2)
        	{
            // Get the next way point using the ego car's (x, y, yaw)
        		next_wp = NextWaypoint(ref_x, ref_y, ref_yaw,
                                   map_waypoints_x, map_waypoints_y,
                                   map_waypoints_dx,map_waypoints_dy);
        	}
          // Not the first time th path planner is running
        	else
        	{
            // Get the last 'x' from the previous path
    				ref_x = previous_path_x[prev_size - 1];
            // Get the previous to last 'x' from the previous path
    				double ref_x_prev = previous_path_x[prev_size - 2];
            // Get the last 'y' from the previous path
    				ref_y = previous_path_y[prev_size - 1];
            // Get the previous to last 'y' from the previous path
    				double ref_y_prev = previous_path_y[prev_size - 2];
            // Get the last heading direction of the previous path
    				ref_yaw = atan2((ref_y-ref_y_prev), (ref_x - ref_x_prev));

            // Get the next way point using the ego car's (x, y, yaw)
    				next_wp = NextWaypoint(ref_x, ref_y, ref_yaw,
                                   map_waypoints_x, map_waypoints_y,
                                   map_waypoints_dx, map_waypoints_dy);

            // Get the previous path end 's'
    				car_s = end_path_s;

            // Get the car speed from the previous path(1m/s = 2.23694mph)
    				car_speed = (sqrt((ref_x - ref_x_prev) * (ref_x - ref_x_prev) +
                              (ref_y - ref_y_prev) * (ref_y - ref_y_prev))/.02) * 2.237;
        	}

        	// Find ref_v to use for the respective lane
        	double closestDist_s = 100000;

          // Flag to indicate whether a lane change is required
        	bool change_lanes = false;

          // Go through the sensor fusion data of the ego car
        	for(int i = 0; i < sensor_fusion.size(); i++)
        	{
        		// Check if the car is in the same lane as the ego car
        		float d = sensor_fusion[i][6];
        		if(d < (2 + 4*lane + 2) && d > (2 + 4*lane - 2))
        		{
              // Get the (vx, vy) of that car
        			double vx = sensor_fusion[i][3];
        			double vy = sensor_fusion[i][4];

              // Speed of the car
        			double check_speed = sqrt(vx*vx + vy*vy);

              // 's' position of the car
        			double check_car_s = sensor_fusion[i][5];

              // Car's 's' position based on the speed
        			check_car_s += ((double)prev_size*.02 * check_speed);

              // Check if 's' value of the car greater than the ego car
              // but less than 's' gap required(30)
        			if((check_car_s > car_s) &&
                 ((check_car_s - car_s) < 30) &&
                 ((check_car_s - car_s) < closestDist_s))
        			{
                // Closest car found
        				closestDist_s = (check_car_s - car_s);

                // The gap between the car is greater than 20
        				if( (check_car_s - car_s) > 20)
        				{
        					// Match that cars speed
        					ref_vel = check_speed * 2.237;
                  // Change lanes
        					change_lanes = true;
        				}
        				else
        				{
        					// Go slightly slower than the cars speed(mph)
        					ref_vel = check_speed*2.237 - 5;
                  // Change lanes
        					change_lanes = true;
        				}
        			}
        		}
        	}

        	// Try to change lanes ego car is too close to car in front
        	if(change_lanes && ((next_wp - lane_change_wp) % map_waypoints_x.size() > 2))
        	{
            // Flag to tracked changed lanes
        		bool changed_lanes = false;

        		// First try to change to left lane(if not in left lane already)
        		if(lane != 0 && !changed_lanes)
        		{
              // Default to true for safe lane change
        			bool lane_safe = true;

        			for(int i = 0; i < sensor_fusion.size(); i++)
        			{
        				// Check for car in left lane
        				float d = sensor_fusion[i][6];
        				if(d < (2 + 4*(lane - 1) + 2) && d > (2 + 4*(lane - 1) - 2))
        				{
                  // Get the (vx, vy) of that car
        					double vx = sensor_fusion[i][3];
        					double vy = sensor_fusion[i][4];

                  // Speed of the car
        					double check_speed = sqrt(vx*vx + vy*vy);

                  // Current 's' position of the car
        					double check_car_s = sensor_fusion[i][5];

                  // Car's s position based on the speed
        					check_car_s+=((double)prev_size * .02 * check_speed);

                  // Distance between the cars
        					double dist_s = check_car_s - car_s;

                  // If the distance is less than the buffer of 20 then
                  // lane change is not safe
        					if(dist_s < 20 && dist_s > -20)
        					{
        						lane_safe = false;
        					}
        				}
        			}
              // Check if the flag for safe lane change is true
        			if(lane_safe)
        			{
                // Safe to change to left lane!
        				changed_lanes = true;
                // Change lane value
        				lane -= 1;
                // Way points for lane change
        				lane_change_wp = next_wp;
        			}
        		}

        		// Try to change to right lane(if not in right lane already)
        		if(lane != 2 && !changed_lanes)
        		{
              // Default to true for safe for lane change
        			bool lane_safe = true;
        			for(int i = 0; i < sensor_fusion.size(); i++)
        			{
        				// Check for car in right lane
        				float d = sensor_fusion[i][6];
        				if(d < (2 + 4*(lane + 1) + 2) && d > (2 + 4*(lane + 1) - 2))
        				{
                  // Get the (vx, vy) of that car
        					double vx = sensor_fusion[i][3];
        					double vy = sensor_fusion[i][4];

                  // Speed of the car
        					double check_speed = sqrt(vx*vx + vy*vy);

                  // Current 's' position of the car
        					double check_car_s = sensor_fusion[i][5];

                  // Car's s position based on the speed
        					check_car_s+=((double)prev_size * .02 * check_speed);

                  // Distance between the cars
        					double dist_s = check_car_s - car_s;

                  // If the distance is less than the buffer of 20 then
                  // lane change is not safe
        					if(dist_s < 20 && dist_s > -10)
        					{
        						lane_safe = false;
        					}
        				}
        			}
              // Check if the flag for safe lane change is true
        			if(lane_safe)
        			{
                // Safe to change to right lane!
        				changed_lanes = true;
                // Change lane value
        				lane += 1;
                // Way points for lane change
        				lane_change_wp = next_wp;
        			}
        		}
        	}

          // Vectors to store car's x and y
        	vector<double> ptsx;
        	vector<double> ptsy;

          // This is the first time the path planning algorithm is running
        	if(prev_size < 2)
        	{
            // Get the car's previous 'x' & 'y'
        		double prev_car_x = car_x - cos(car_yaw);
        		double prev_car_y = car_y - sin(car_yaw);

            // Add the previous and the current car positions
        		ptsx.push_back(prev_car_x);
        		ptsx.push_back(car_x);
        		ptsy.push_back(prev_car_y);
        		ptsy.push_back(car_y);
        	}
          // Not the first time the car is running
        	else
        	{
            // Add the previous 2 car positions
        		ptsx.push_back(previous_path_x[prev_size-2]);
        		ptsx.push_back(previous_path_x[prev_size-1]);
        		ptsy.push_back(previous_path_y[prev_size-2]);
        		ptsy.push_back(previous_path_y[prev_size-1]);
        	}

          // Get next 3 way points using the car's current 's' and lane position('d') at the middle of the lane
        	vector<double> next_wp0 = getXY(car_s+30, (2 + 4*lane),
                                          map_waypoints_s, map_waypoints_x, map_waypoints_y);
        	vector<double> next_wp1 = getXY(car_s+60, (2 + 4*lane),
                                          map_waypoints_s, map_waypoints_x, map_waypoints_y);
        	vector<double> next_wp2 = getXY(car_s+90, (2 + 4*lane),
                                          map_waypoints_s, map_waypoints_x, map_waypoints_y);

          // Push 3 more way points
        	ptsx.push_back(next_wp0[0]);
        	ptsx.push_back(next_wp1[0]);
        	ptsx.push_back(next_wp2[0]);
        	ptsy.push_back(next_wp0[1]);
        	ptsy.push_back(next_wp1[1]);
        	ptsy.push_back(next_wp2[1]);

          // Shift car reference angle to 0 degrees
        	for (int i = 0; i < ptsx.size(); i++)
        	{
        		double shift_x = ptsx[i] - ref_x;
        		double shift_y = ptsy[i] - ref_y;
				    ptsx[i] = ((shift_x * cos(0 - ref_yaw)) - (shift_y * sin(0 - ref_yaw)));
				    ptsy[i] = ((shift_x * sin(0 - ref_yaw)) + (shift_y * cos(0 - ref_yaw)));
          }

          // Create spline for car to follow based on the points generated(5)
        	tk::spline s;
        	s.set_points(ptsx, ptsy);

          // Vector of co-ordinates for the car to follow
        	vector<double> next_x_vals;
        	vector<double> next_y_vals;

          // Start from the previous path to allow for a smoother transition
        	for(int i = 0; i < previous_path_x.size(); i++)
        	{
        		next_x_vals.push_back(previous_path_x[i]);
        		next_y_vals.push_back(previous_path_y[i]);
        	}

          // Target 'x' position
        	double target_x = 30.0;
          // Get 'y' using the spline equation
        	double target_y = s(target_x);

          // Target distance to reach
        	double target_dist = sqrt((target_x) * (target_x) + (target_y) * (target_y));
        	double x_add_on = 0;

          // Populate 50 points for the car to follow using the spline
			    for (int i = 1; i <= 50 - previous_path_x.size(); i++)
          {
            // The car speed if below reference velocity
    				if(car_speed < ref_vel)
    				{
    					car_speed+=.224;
    				}
            // Car speed is above reference velocity
    				else if(car_speed > ref_vel)
    				{
    					car_speed-=.224;
    				}

            // Calculate the increments for 'x' based on the speed
    				double N = (target_dist / (.02*car_speed / 2.24));
            // New 'x' point
    				double x_point = x_add_on + (target_x)/N;
            // New 'y' point using the spline
    				double y_point = s(x_point);

            // Store for next iteration
    				x_add_on = x_point;

            // For converting from car co-ordinates to global co-ordinates
    				double x_ref = x_point;
    				double y_ref = y_point;

            // Global co-ordinates after conversion
    				x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
    				y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

            // New global 'x' & 'y' points with increments added
    				x_point += ref_x;
    				y_point += ref_y;

            // Add the points to the vector
    				next_x_vals.push_back(x_point);
    				next_y_vals.push_back(y_point);
    			}

        	json msgJson;
        	msgJson["next_x"] = next_x_vals;
        	msgJson["next_y"] = next_y_vals;

        	auto msg = "42[\"control\","+ msgJson.dump()+"]";

        	// This_thread::sleep_for(chrono::milliseconds(1000));
        	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      }
      else
      {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req,
                     char *data, size_t, size_t)
  {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req)
  {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length)
  {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  // Listen to port 4567 for communication
  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  // Run server
  h.run();
}
