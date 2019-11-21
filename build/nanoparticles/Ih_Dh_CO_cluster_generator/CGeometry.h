#ifndef __CGEOMETRY__
#define __CGEOMETRY__
/*
	CGeometry: Abstract base class which various geometric classes implement

*/

#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<math.h>
#include<sstream>

#include "Struct.h";

#define PI 3.141592653589793238462643383279502884197 //BIG PIE :)
#define DEG2RAD PI/180.0

struct vec3d //a 3d vector with double precission
{
	double x;
	double y;
	double z;
};

class CGeometry
{
public:
	CGeometry();
	~CGeometry();

	//virtual functions
	virtual void generateCoordinates(const int &shell_num, std::string element) = 0;
	virtual void generateSingleShell(const int &shell_num, std::string element) = 0;
	virtual void setRadius(float radius) = 0;
	virtual void setShells(int shell) = 0;

	//other methods
	//every time we add a coordinate we will convert add to m_Cluster;
	void add_coordinate(vec3d temp, std::string element); 
	void clearShell();

	//accessors
	//GetCluster(){ return m_Cluster;}
        clusterData GetShell() { return m_Cluster; }

	//helper public methods
	void point_on_3D_line(const vec3d& b_vec, const vec3d& p2, vec3d& new_point, double nu);
private:
	clusterData m_Cluster; //we hold our cluster data structure here
};

/*
//#define RADIUS_MULTIPLIER 2.7

struct vec3d //a 3d vector with double precission
{
	double x;
	double y;
	double z;
};
class CGeometry
{
public:
	CGeometry();
	~CGeometry();
	void write_coordinaates(string filename);
	void print_coordinates();
	void set_element(string &s){element = s;}
	void set_element(string &s, string &s2, int &coreShell){element = s, element2 = s2, coreShellSize = coreShell, isBimettalic = true; }
	virtual void setRadius(float radius) = 0;
	virtual void generateCoordinates(const int &shell_num) = 0;
	virtual void generateSingleShell(const int &shell_num) = 0;
	void point_on_3D_line(const vec3d& b_vec, const vec3d& p2, vec3d& new_point, double nu);
	void add_coordinate(vec3d temp){coordinates.push_back(temp);}
	vec3d get_coordinate(int i){return coordinates[i];}
	friend ostream& operator<<(ostream& os, const CGeometry& c);
private:
	bool isBimettalic;
	int coreShellSize;
	string element, element2;
};
*/
#endif
