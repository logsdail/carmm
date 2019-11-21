/***************************************************************************
 *   Copyright (C) 2006 by Benjamin Curley                                 *
 *   curley@tc.bham.ac.uk						 						   *	
 *   rhodan@blueyonder.co.uk		                       	               *
 * 																		   *
 *   GeometryMaker: Creates simple 12 vectrex geometries for metallic      *
 *		    clusters						  							   *
 ***************************************************************************/
#ifndef __CICO__
#define __CICO__

#include "CGeometry.h"
#include <cmath>

using namespace std;

class CIcoGen : public CGeometry
{
public:
	CIcoGen();
	~CIcoGen();
	void generate_vertices(int shell_num);
	void generate_edges(int shell_num);
	void generate_faces(int shell_num);
	void clear_arrays(){vertices.clear(); face_atoms.clear(); edge_atoms.clear();}
	virtual void generateCoordinates(const int &shell_num, std::string element);
	virtual void generateSingleShell(const int &shell_num, std::string element);
	void combine_arrays(std::string element);
	void setRadius(float radius){m_fRadius = radius;}
	void setShells(int shells){m_iShells = shells;}
private:
	float m_fRadius;
	int m_iShells;
	vector<vec3d> vertices;
	vector<vec3d> face_atoms;
	vector<vec3d> edge_atoms;

};

class CCubGen : public CGeometry
{
public:
	CCubGen(){};
	~CCubGen(){};
	void generate_vertices(int shell_num); 
	void generate_edges(int shell_num);
	void generate_face111(int shell_num);
	void generate_face100(int shell_num);   
	virtual void generateCoordinates(const int &shell_num, std::string element);
	virtual void generateSingleShell(const int &shell_num, std::string element);
	void combine_arrays(std::string element);
	void setRadius(float radius){m_fRadius = radius;}
	void setShells(int shells){m_iShells = shells;}
private:
	float m_fRadius;
	int m_iShells;
	vector<vec3d> vertices;
	vector<vec3d> face_atoms;
	vector<vec3d> edge_atoms;
	vector<vec3d> single_face;
};

class COctGen : public CGeometry
{
public:
        COctGen(){};
        ~COctGen(){};
        void generate_vertices(int shell_num);
        void generate_edges(int shell_num);
        void generate_face100(int shell_num);
	void generateOrigin(int shell_num, std::string element);
        virtual void generateCoordinates(const int &shell_num, std::string element);
        virtual void generateSingleShell(const int &shell_num, std::string element);
        void combine_arrays(int shell_num, std::string element);
        void setRadius(float radius){m_fRadius = radius;}
	void setShells(int shells){m_iShells = shells;}
private:
        float m_fRadius;
	int m_iShells;
        vector<vec3d> vertices;
        vector<vec3d> face_atoms;
        vector<vec3d> edge_atoms;
        vector<vec3d> single_face;
};

class CDecGen : public CGeometry
{
public:
	CDecGen(){};
	~CDecGen(){};
	void generate_vertices(int shell_num);
	void generate_edges(int shell_num);
	void generate_face111(int shell_num);
	void generate_face100(int shell_num);
	virtual void generateCoordinates(const int &shell_num, std::string element);
	virtual void generateSingleShell(const int &shell_num, std::string element);
	void combine_arrays(std::string element);
	void setRadius(float radius){m_fRadius = radius;}	
	void setShells(int shells){m_iShells = shells;}
private:
	float m_fRadius;
	int m_iShells;
	vector<vec3d> vertices;
	vector<vec3d> face_atoms;
	vector<vec3d> edge_atoms;
	vector<vec3d> single_face;
};


class CPentBi : public CGeometry
{
public:
	CPentBi(){};
	~CPentBi(){};

	void generate_vertices(int shell_num);
	void generate_edges(int shell_num);
	void generate_face111(int shell_num);
	void generate_face100(int shell_num);
	virtual void generateCoordinates(const int& shell_num, std::string element);
	virtual void generateSingleShell(const int& shell_num, std::string element);
	void combine_arrays(std::string element);
	void setRadius(float radius){m_fRadius = radius;}
	void setShells(int shells){m_iShells = shells;}
private:
	float m_fRadius;
	int m_iShells;
	vector<vec3d> vertices;
	vector<vec3d> face_atoms;
	vector<vec3d> edge_atoms;
	vector<vec3d> single_face;
};

#endif
