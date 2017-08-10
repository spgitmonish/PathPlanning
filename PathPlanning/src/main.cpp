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

using namespace std;

#define EXPERIMENT 0

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{
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

// Function for calculating the distance based on (x1,y1) & (x2,y2)
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
int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
  // Get the closest way point
	int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

  // Get the corresponding co-ordinates
	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

  // Get the heading direction based on the next way point and
  // the car's current location
	double heading = atan2((map_y - y),(map_x - x));

  // Calculate the angle based on the current theta and the heading direction
	double angle = abs(theta - heading);

  // Check if the angle is greater than pi/4
	if(angle > pi()/4)
	{
    // Change the index to be the next one
		closestWaypoint++;
	}

  // Return the index
	return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
  // Get the next way points
	int next_wp = NextWaypoint(x, y, theta, maps_x,maps_y);

  // Previous way point
	int prev_wp;
	prev_wp = next_wp - 1;

  // If the next way point is 0(circular), get the last way point in the vector
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size() - 1;
	}

  // Distance in x, y
	double n_x = maps_x[next_wp] - maps_x[prev_wp];
	double n_y = maps_y[next_wp] - maps_y[prev_wp];
  // Difference between current x,y and previous way point
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// Find the projection of x onto n
	double proj_norm = (x_x * n_x + x_y * n_y)/(n_x * n_x + n_y * n_y);
	double proj_x = proj_norm * n_x;
	double proj_y = proj_norm * n_y;

  // Get frenet d
	double frenet_d = distance(x_x, x_y, proj_x, proj_y);

	// See if d value is positive or negative by comparing it to a center point
	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x, center_y, x_x, x_y);
	double centerToRef = distance(center_x, center_y, proj_x, proj_y);

  // Adjust 'd' based on center reference
	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// Calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i+1], maps_y[i+1]);
	}
	frenet_s += distance(0, 0, proj_x, proj_y);

  // Return the 's' and 'd' co-ordinates
 	return {frenet_s, frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
  // Initialize previous way point to be -1
	int prev_wp = -1;

	while(s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size()-1)))
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

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
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

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
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
          // j[1] is the data JSON object
        	// Main car's localization Data
        	double car_x = j[1]["x"];
        	double car_y = j[1]["y"];
        	double car_s = j[1]["s"];
        	double car_d = j[1]["d"];
        	double car_yaw = j[1]["yaw"];
        	double car_speed = j[1]["speed"];

        	// Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];

        	// Previous path's end s and d values
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          // NOTE: 1. (x, y) are in global map co-ordinates
          //       2. (s, d) are the Frenet co-ordinates of that car
          //       3. (vx, vy) are the velocity components in reference to the global map
        	auto sensor_fusion = j[1]["sensor_fusion"];

        	json msgJson;

          // Vectory of next_x_vals and next_y_vals
        	vector<double> next_x_vals;
        	vector<double> next_y_vals;

        	// TODO: Define a path made up of (x,y) points that the car will
          //       visit sequentially every .02 seconds
        #if(EXPERIMENT == 1)
          double dist_inc = 0.3;
          // Move 50 times a second a distance of 0.5m per move.
          // Velocity = 15m/s which is around 33.554mph
          for(int i = 0; i < 50; i++)
          {
            // Next position based on constant yaw whatever the angle the car is at
            next_x_vals.push_back(car_x + (dist_inc * i)*cos(deg2rad(car_yaw)));
            next_y_vals.push_back(car_y + (dist_inc * i)*sin(deg2rad(car_yaw)));
          }
        #elif(EXPERIMENT == 2)
          double pos_x;
          double pos_y;
          double angle;

          // Previous path size
          int path_size = previous_path_x.size();

          // Append values to the list till the previous path size
          for(int i = 0; i < path_size; i++)
          {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // If there was no previous path
          if(path_size == 0)
          {
            // Use the car's position as the location
            pos_x = car_x;
            pos_y = car_y;
            // And the car's angle as the yaw
            angle = deg2rad(car_yaw);
          }
          else
          {
            // Get the last position from the previous path
            pos_x = previous_path_x[path_size-1];
            pos_y = previous_path_y[path_size-1];

            // Position before the last position
            double pos_x2 = previous_path_x[path_size-2];
            double pos_y2 = previous_path_y[path_size-2];

            // Calculate the angle betweent the 2 points
            angle = atan2(pos_y-pos_y2, pos_x-pos_x2);
          }

          double dist_inc = 0.3;
          // Populate the remaining points at an angle
          for(int i = 0; i < 50-path_size; i++)
          {
            // Next points with a changing angle between the points
            next_x_vals.push_back(pos_x + (dist_inc)*cos(angle + (i+1)*(pi()/100)));
            next_y_vals.push_back(pos_y + (dist_inc)*sin(angle + (i+1)*(pi()/100)));
            pos_x += (dist_inc)*cos(angle + (i+1)*(pi()/100));
            pos_y += (dist_inc)*sin(angle + (i+1)*(pi()/100));
          }
        #endif

          // The next_x_vals and next_y_vals are in global co-ordinates
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
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t)
  {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
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
  h.run();
}
